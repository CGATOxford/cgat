
##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
"""
===========================
Pipeline IDR
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
========

This pipeline is a wrapper for the Irreproducibility Discovery Pipeline
outlined by Anshul Kundaje.
(see https://sites.google.com/site/anshulkundaje/projects/idr)

It is split into two stages and performs several basic steps:

Fist stage
----------

1) Pre-processing bamfiles
Bamfiles are filtered prior to peak-calling. With the option to remove
non-uniquely mapping reads, remove duplicates, and mask genomic regions in a
specified bedfile.
Bamfiles are variously pooled by EXPERIMENT (tissue-condition-agg )and split,
to provide pseudoreplicates of individual samples and pooled-pseudoreplicates.

2) Peak calling is carried out using SPP on three sets of samples:
i) individual samples, ii) pseudoreplicates, and iii) pooled-pseudoreplicates.
Three options exist for using control files in peak-calling:
i) All peak calling may be carried out against a single control per experiment,
in which case the control must be labelled <tissue>-<control>-R1.
ii) Peaks may be called against pooled control files for each EXPERIMENT (as
recommended for IDR) in which case all three peak-calling steps will be carried
out using pooled controls labelled <tissue>-<condition>-R0.
iii) Peak calling for individual samples may be carried out using controls
matched by replicate. In which case peak calling for pseudoreplicates will be
carried out against the same control (i.e. pseudoreplication of the control
file does not take place), and peak calling for the pooled pseudoreplicates is
carried out against a pooled control.
Peakcalling for IDR is necessarily lax and the option exists to specify the
number of peaks called using either -npeaks or -fdr. (It may be necessary to
use a customized version of run_spp.R for one or other of these options,
depending on which version of SPP is being run, this option is also provided.)

3) IDR analysis is carried out using WrapperIDR (which us a python wrapper for
batch-consistency-analysis.r) for i) all pairwise comparisons of individual
replicates within an EXPERIMENT, ii) pairs of psuedoreplicates for each sample,
and iii) paired pooled pseudoreplicates for each EXPERIMENT. These steps
provide IDR analysis on origanl replicates, IDR analysis on self-
pseudoreplicates, and IDR analysis on pooled-pseudoreplicates, respectively.

Second stage
------------

4) Following IDR analysis, a second stage of the pipeline (generatePeakSets)
may be run. Libraries that fail in terms of 'max_numPeaks_Rep' or
'self-consistency' may be specified for exclusion in the config file.
Subsequently, all non-excluded bamfiles are pooled per EXPERIMENT so that both
a conservative peak set and optimal peak set may be derived.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information
how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

Sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini`
file (see :ref:`PipelineDocumenation`). To start with, use the files supplied
with the :ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the
configuration variable :py:data:`annotations_database`
and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software
to be in the path:

+--------------------+-------------------+------------------------------------+
|*Program*           |*Version*          |*Purpose*                           |
+--------------------+-------------------+------------------------------------+
|                    |                   |                                    |
+--------------------+-------------------+------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data\
/pipeline_IDR.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_IDR.tgz
   tar -xvzf pipeline_IDR.tgz
   cd pipeline_IDR
   python <srcdir>/pipeline_IDR.py make full

.. note::
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""

from ruffus import *

import sys
import os
import itertools
import re
import sqlite3
import glob
import shutil

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineUtilities as PU
import CGATPipelines.PipelineIDR as IDR

##########################################################################
##########################################################################
##########################################################################
# Pipeline configuration
##########################################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__")

##########################################################################
##########################################################################
# Helper functions mapping tracks to conditions, etc
##########################################################################

import CGATPipelines.PipelineTracks as PipelineTracks

Sample = PipelineTracks.AutoSample

# define tracks based on all samples in .bamfile that are not input or index
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob(os.path.join(PARAMS["location_bamfiles"], "*.bam")),
    "(\S+).bam",
    exclude=[".+input.+"])


@files(None, None)
def printTracks(infile, outfile):
    P.warn("\n\n\n\nprinting tracks:")
    for track in EXPERIMENTS:
        print "\t"
        print track


def get_peak_caller_parameters(peak_caller_id):
    """
    Returns a dictionary of config file parameters for the chosen peak caller
    (an attempt to keep access to PARAMS out of associated pipeline script).
    """
    caller_parameters = {}
    caller_prefix = peak_caller_id + "_options"
    for key, value in PARAMS.iteritems():
        if re.match(caller_prefix, key):
            caller_parameters[key] = value

    return caller_parameters

##########################################################################
##########################################################################
##########################################################################


def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

##########################################################################
##########################################################################
##########################################################################
# Process Bamfiles
##########################################################################


@follows(mkdir("bamfiles_filtered"))
@transform(os.path.join(PARAMS["location_bamfiles"], "*.bam"),
           regex("(.+)/(.+).bam"),
           r"./bamfiles_filtered/\2.sentinel")
def filterBamfiles(infile, sentinel):
    """
    Pre-process bamfiles prior to peak calling.
    i) sort bamfiles
    ii) remove unmapped readswith bam2bam.py
    iii) remove non-uniquely mapping reads with bam2bam.py (optional)
    iv) remove duplicates with Picards MarkDuplicates (optional)
    v) remove reads from masked regions with bedtools intersect (optional)
    vi) index
    """

    # create tempfile for Picard's MarkDuplicates
    picard_tmp = picard_tmp = P.getTempDir(PARAMS["scratchdir"])

    outfile = P.snip(sentinel, ".sentinel") + ".bam"

    # ensure bamfile is sorted,
    statement = ["samtools sort @IN@ @OUT@", ]

    # remove unmapped reads
    statement.append("python %(scriptsdir)s/bam2bam.py"
                     " --filter=mapped"
                     " --log=%(outfile)s.log"
                     " < @IN@.bam"
                     " > @OUT@")

    # remove non-uniquely mapping reads, if requested
    if PARAMS["filter_remove_non_unique"]:
        statement.append("python %(scriptsdir)s/bam2bam.py"
                         " --filter=unique"
                         " --log=%(outfile)s.log"
                         " < @IN@"
                         " > @OUT@")

    # remove duplicates, if requested
    if PARAMS["filter_remove_duplicates"]:
        statement.append("MarkDuplicates"
                         " INPUT=@IN@"
                         " ASSUME_SORTED=true"
                         " REMOVE_DUPLICATES=true"
                         " QUIET=false"
                         " OUTPUT=@OUT@"
                         " METRICS_FILE=/dev/null"
                         " VALIDATION_STRINGENCY=SILENT"
                         " TMP_DIR=%(picard_tmp)s"
                         " 2> %(outfile)s.log")

    # mask regions, if intervals supplied
    if PARAMS["filter_mask_intervals"]:
        mask = PARAMS["filter_mask_intervals"]
        statement.append("bedtools intersect"
                         " -abam @IN@"
                         " -b %(mask)s"
                         " -wa"
                         " -v"
                         " > @OUT@")

    statement.append("mv @IN@ %(outfile)s")
    statement.append("samtools index %(outfile)s")

    job_options = "-l mem_free=10G"
    statement = P.joinStatements(statement, infile)

    P.run()
    P.touch(sentinel)
    shutil.rmtree(picard_tmp)


@follows(filterBamfiles, mkdir("bamfiles_pseudoreplicates"))
@transform([os.path.join("./bamfiles_filtered", x.asFile() + ".sentinel")
            for x in TRACKS],
           regex("./bamfiles_filtered/(.+).sentinel"),
           r"./bamfiles_pseudoreplicates/\1.sentinel")
def splitBamfiles(infile, sentinel):
    """
    For all tracks, split the filtered bamfile in two using pysam
    """
    infile = P.snip(infile, ".sentinel") + ".bam"
    outfile = P.snip(sentinel, ".sentinel")
    params = '2'
    try:
        module = P.snip(IDR.__file__, ".py")
    except ValueError:
        module = P.snip(IDR.__file__, ".pyc")

    P.submit(module,
             "splitBam",
             params,
             infile,
             outfile)

    P.touch(sentinel)


# AH: remove replicate requirement for input tracks
# AH: avoid merging if only one repliacte
@follows(mkdir("bamfiles_pooled"))
@collate(filterBamfiles,
         regex(r"(.+)/(.+)-input-.*.sentinel"),
         r"./bamfiles_pooled/\2-input-R0.sentinel")
def poolInputBamfiles(infiles, sentinel):
    """
    Merge filtered input files for each tissue, with the option of excluding
    undesirable libraries.
    """
    infiles = [P.snip(x, ".sentinel") + ".bam" for x in infiles]
    outfile = P.snip(sentinel, ".sentinel") + ".bam"
    bad_samples = PARAMS["filter_remove_inputs"].split(",")

    if len(infiles) > 1:
        to_merge = IDR.filterBadLibraries(infiles, bad_samples)
        IDR.mergeBams(to_merge, outfile)
    else:
        os.symlink(infiles[0], outfile)
        os.symlink(infiles[0] + ".bai", outfile + ".bai")

    P.touch(sentinel)


@follows(mkdir("bamfiles_pooled"))
@collate(filterBamfiles,
         regex(r"(.+)/(.+)-((?!input).*)-R[0-9]+.sentinel"),
         r"./bamfiles_pooled/\2-\3-R0.sentinel")
def poolSampleBamfiles(infiles, sentinel):
    """
    Merge filtered sample files for each tissue
    """
    infiles = [P.snip(x, ".sentinel") + ".bam" for x in infiles]
    outfile = P.snip(sentinel, ".sentinel") + ".bam"

    IDR.mergeBams(infiles, outfile)

    P.touch(sentinel)


@follows(mkdir("./bamfiles_pooled_pseudoreplicates"))
@transform(poolSampleBamfiles,
           regex("(.+)/(.+).sentinel"),
           r"./bamfiles_pooled_pseudoreplicates/\2.sentinel")
def splitPooledBamfiles(infile, sentinel):
    infile = P.snip(infile, ".sentinel") + ".bam"
    outfile = P.snip(sentinel, ".sentinel")
    params = '2'
    try:
        module = P.snip(IDR.__file__, ".py")
    except ValueError:
        module = P.snip(IDR.__file__, ".pyc")

    P.submit(module,
             "splitBam",
             params,
             infile,
             outfile)

    P.touch(sentinel)

##########################################################################
##########################################################################
##########################################################################
# Run Peak calling for IDR
##########################################################################


@follows(filterBamfiles,
         poolInputBamfiles,
         mkdir("peakfiles_individual_replicates"))
@transform([os.path.join("./bamfiles_filtered", x.asFile() + ".sentinel")
            for x in TRACKS],
           regex("(.+)/(.+).sentinel"),
           r"./peakfiles_individual_replicates/\2_regionPeak.sentinel")
def callPeaksOnIndividualReplicates(infile, outfile):
    infile = P.snip(infile, ".sentinel") + ".bam"
    # fetch peak calling parameters
    PARAMS_PEAKCALLER = get_peak_caller_parameters(
        PARAMS["options_peak_caller"])

    # call peaks
    IDR.callIDRPeaks(infile,
                     outfile,
                     PARAMS["options_peak_caller"],
                     PARAMS["options_control_type"],
                     PARAMS_PEAKCALLER)

    P.touch(outfile)


@follows(callPeaksOnIndividualReplicates)
@merge("./peakfiles_individual_replicates/*.narrowPeak.gz",
       "./peakfiles_individual_replicates/"
       "peakcalling_summary_individual_replicates.tsv")
def summarizePeaksForIndividualReplicates(infiles, outfile):
    outf = IOTools.openFile(outfile, "w")
    outf.write("Sample_id\t"
               "Experiment\t"
               "Tissue\t"
               "Condition\t"
               "Replicate\t"
               "n_peaks\n")
    IDR.countPeaks(infiles, outf)


@transform(summarizePeaksForIndividualReplicates,
           regex("(.+)/(.+).tsv"),
           r"\2.load")
def loadPeakSummaryForIndividualReplicates(infile, outfile):
    P.load(infile, outfile)


@follows(splitBamfiles,
         poolInputBamfiles,
         mkdir("peakfiles_pseudoreplicates"))
@transform("./bamfiles_pseudoreplicates/*.bam",
           regex("(.+)/(.+).bam"),
           r"./peakfiles_pseudoreplicates/\2_regionPeak.sentinel")
def callPeaksOnPseudoreplicates(infile, outfile):
    # fetch peak calling parameters
    PARAMS_PEAKCALLER = get_peak_caller_parameters(
        PARAMS["options_peak_caller"])

    # call peaks on pseudoreplicates
    IDR.callIDRPeaks(infile,
                     outfile,
                     PARAMS["options_peak_caller"],
                     PARAMS["options_control_type"],
                     PARAMS_PEAKCALLER,
                     pseudoreplicate=True)


@follows(callPeaksOnPseudoreplicates)
@merge("./peakfiles_pseudoreplicates/*.narrowPeak.gz",
       "./peakfiles_pseudoreplicates/"
       "peakcalling_summary_pseudoreplicates.tsv")
def summarizePeaksForPseudoreplicates(infiles, outfile):
    outf = IOTools.openFile(outfile, "w")
    outf.write("Sample_id\t"
               "Experiment\t"
               "Tissue\t"
               "Condition\t"
               "Pseudoreplicate\t"
               "n_peaks\n")
    IDR.countPeaks(infiles, outf)


@transform(summarizePeaksForPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"\2.load")
def loadPeakSummaryForPseudoreplicates(infile, outfile):
    P.load(infile, outfile)


@follows(splitPooledBamfiles,
         poolInputBamfiles,
         mkdir("peakfiles_pooled_pseudoreplicates"))
@transform("./bamfiles_pooled_pseudoreplicates/*.bam",
           regex("(.+)/(.+).bam"),
           r"./peakfiles_pooled_pseudoreplicates/\2.regionPeak.sentinel")
def callPeaksOnPooledPseudoreplicates(infile, outfile):
    # fetch peak calling parameters
    PARAMS_PEAKCALLER = get_peak_caller_parameters(
        PARAMS["options_peak_caller"])

    # call peaks on pseudoreplicates
    IDR.callIDRPeaks(infile,
                     outfile,
                     PARAMS["options_peak_caller"],
                     PARAMS["options_control_type"],
                     PARAMS_PEAKCALLER,
                     pseudoreplicate=True)

    P.touch(outfile)


@follows(callPeaksOnPooledPseudoreplicates)
@merge("./peakfiles_pooled_pseudoreplicates/*.narrowPeak.gz",
       "./peakfiles_pooled_pseudoreplicates/"
       "peakcalling_summary_pooled_pseudoreplicates.tsv")
def summarizePeaksForPooledPseudoreplicates(infiles, outfile):
    outf = IOTools.openFile(outfile, "w")
    outf.write("Sample_id\t"
               "Experiment\t"
               "Tissue\t"
               "Condition\t"
               "Pseudoreplicate\t"
               "n_peaks\n")
    IDR.countPeaks(infiles, outf)


@transform(summarizePeaksForPooledPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"\2.load")
def loadPeakSummaryForPooledPseudoreplicates(infile, outfile):
    P.load(infile, outfile)


##########################################################################
##########################################################################
##########################################################################
# Run IDR
##########################################################################


@follows(callPeaksOnIndividualReplicates,
         mkdir("./idr_individual_replicates"))
@collate("./peakfiles_individual_replicates/*.regionPeak.gz",
         regex(r"(.+)/(.+)-(.+)-(R[0-9]+)_VS_(.+).regionPeak.gz"),
         r"./idr_individual_replicates/\2-\3.idr")
def runIDROnIndividualReplicates(infiles, outfile):
    """
    Run IDR consecutively for each pairwise combination of a particular
    EXPERIMENT
    """
    # set IDR parameters (HACK!) WrapperIDR is in /ifs/devel/CGAT
    chr_table = os.path.join(PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_contigs"])
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")

    # iterate through pairwise combinations of infiles
    for infile1, infile2 in itertools.combinations(infiles, 2):
        # get statement
        statement = IDR.getIDRStatement(infile1,
                                        infile2,
                                        outfile,
                                        PARAMS["idr_options_overlap_ratio"],
                                        PARAMS["idr_options_ranking_measure"],
                                        chr_table,
                                        idr_script)

        # run
        E.info("applyIDR: processing %s and %s" % (infile1, infile2))
        job_options = "-l mem_free=5G"
        P.run()
#        print "\n" + statement + "\n"


@transform(runIDROnIndividualReplicates,
           regex("(.+).idr"),
           add_inputs(r"\1-*-uri.sav"),
           r"\1_batch-consistency.pdf")
def plotBatchConsistencyForIndividualReplicates(infiles, outfile):
    # HACK!
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")
    statement = IDR.getIDRPlotStatement(infiles, outfile, idr_script)
    P.run()
#    print statement


@follows(callPeaksOnPseudoreplicates,
         mkdir("./idr_pseudoreplicates"))
@collate("./peakfiles_pseudoreplicates/*.regionPeak.gz",
         regex("./peakfiles_pseudoreplicates/"
               "(.+)(_00|_01)_VS_(.+).regionPeak.gz"),
         r"idr_pseudoreplicates/\1_pseudoreplicates.idr")
def runIDROnPseudoreplicates(infiles, outfile):
    """
    Run IDR analysis on pseudoreplicates for each TRACK
    """
    # set IDR parameters
    chr_table = os.path.join(PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_contigs"])
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")

    # get statement
    statement = IDR.getIDRStatement(infiles[0],
                                    infiles[1],
                                    outfile,
                                    PARAMS["idr_options_overlap_ratio"],
                                    PARAMS["idr_options_ranking_measure"],
                                    chr_table,
                                    idr_script)

    # run
    E.info("applyIDR: processing %s and %s" % (infiles[0], infiles[1]))
    job_options = "-l mem_free=5G"
    P.run()
#    print  "\n" + statement + "\n"


@collate(runIDROnPseudoreplicates,
         regex("(.+)_pseudoreplicates.idr"),
         add_inputs(r"\1*-uri.sav"),
         r"\1_batch-consistency.pdf")
def plotBatchConsistencyForPseudoreplicates(infiles, outfile):
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")
    statement = IDR.getIDRPlotStatement(infiles[0], outfile, idr_script)
    P.run()


@follows(callPeaksOnPooledPseudoreplicates,
         mkdir("./idr_pooled_pseudoreplicates"))
@collate("./peakfiles_pooled_pseudoreplicates/*.regionPeak.gz",
         regex(r"./peakfiles_pooled_pseudoreplicates/"
               "(.+)(_00|_01)_VS_(.+).regionPeak.gz"),
         r"./idr_pooled_pseudoreplicates/\1_pooled_pseudoreplicate.idr")
def runIDROnPooledPseudoreplicates(infiles, outfile):
    """
    Run IDR analysis on pooled pseudoreplicates for each EXPERIMENT
    """
    # set IDR parameters
    chr_table = os.path.join(PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_contigs"])
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")

    # get statement
    statement = IDR.getIDRStatement(infiles[0],
                                    infiles[1],
                                    outfile,
                                    PARAMS["idr_options_overlap_ratio"],
                                    PARAMS["idr_options_ranking_measure"],
                                    chr_table,
                                    idr_script)

    # run
    E.info("applyIDR: processing %s and %s" % (infiles[0], infiles[1]))
    job_options = "-l mem_free=5G"
    P.run()
#    print "\n" + statement + "\n"


@collate(runIDROnPooledPseudoreplicates,
         regex("(.+)_pooled_pseudoreplicate.idr"),
         add_inputs(r"\1*-uri.sav"),
         r"\1_batch-consistency.pdf")
def plotBatchConsistencyForPooledPseudoreplicates(infiles, outfile):
    idr_script = os.path.join(os.path.dirname(P.__file__), "WrapperIDR.py")
    statement = IDR.getIDRPlotStatement(infiles[0], outfile, idr_script)
    P.run()
#    print "\n" + statement + "\n"

##########################################################################
##########################################################################
##########################################################################
# Post Process IDR
##########################################################################


@follows(plotBatchConsistencyForIndividualReplicates)
@transform(runIDROnIndividualReplicates,
           regex("(.+).idr"),
           add_inputs(r"\1-*-npeaks-aboveIDR.txt"),
           r"\1_npeaks_aboveIDR.tsv")
def combineIDROnIndividualReplicates(infiles, outfile):
    tables = infiles[1:]
    headers = [os.path.basename(x) for x in tables]
    headers = ",".join([P.snip(x, "-npeaks-aboveIDR.txt") for x in headers])
    tables = " ".join(tables)

    to_cluster = False
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --skip-titles"
                 " --headers=%(headers)s"
                 " --log=%(outfile)s.log"
                 " --stdout=%(outfile)s"
                 " %(tables)s")
    P.run()


@transform(combineIDROnIndividualReplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadIDROnIndividualReplicates(infile, outfile):
    P.load(infile, outfile)


@follows(plotBatchConsistencyForPseudoreplicates)
@collate("./idr_pseudoreplicates/*-npeaks-aboveIDR.txt",
         regex("(.+)-R(.+)_vs_(.+)-npeaks-aboveIDR.txt"),
         r"\1_pseudoreplicates_npeaks_aboveIDR.tsv")
def combineIDROnPseudoreplicates(infiles, outfile):
    headers = [os.path.basename(x) for x in infiles]
    headers = ",".join([P.snip(x, "-npeaks-aboveIDR.txt") for x in headers])
    tables = " ".join(infiles)

    to_cluster = False
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --skip-titles"
                 " --headers=%(headers)s"
                 " --log=%(outfile)s.log"
                 " --stdout=%(outfile)s"
                 " %(tables)s")
    P.run()


@transform(combineIDROnPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadIDROnPseudoreplicates(infile, outfile):
    P.load(infile, outfile)


@follows(plotBatchConsistencyForPooledPseudoreplicates)
@collate("./idr_pooled_pseudoreplicates/*-npeaks-aboveIDR.txt",
         regex("(.+)-R(.+)_vs_(.+)-npeaks-aboveIDR.txt"),
         r"\1_pooled_pseudoreplicates_npeaks_aboveIDR.tsv")
def combineIDROnPooledPseudoreplicates(infiles, outfile):
    headers = [os.path.basename(x) for x in infiles]
    headers = ",".join([P.snip(x, "-npeaks-aboveIDR.txt") for x in headers])
    tables = " ".join(infiles)

    to_cluster = False
    statement = ("python %(scriptsdir)s/combine_tables.py"
                 " --columns=1"
                 " --skip-titles"
                 " --headers=%(headers)s"
                 " --log=%(outfile)s.log"
                 " --stdout=%(outfile)s"
                 " %(tables)s")
    P.run()


@transform(combineIDROnPooledPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadIDROnPooledPseudoreplicates(infile, outfile):
    P.load(infile, outfile)


@follows(runIDROnIndividualReplicates)
@merge("./idr_individual_replicates/*-overlapped-peaks.txt",
       "./idr_individual_replicates/individual_replicates.nPeaks.tsv")
def findNPeaksForIndividualReplicates(infiles, outfile):
    idr_thresh = PARAMS["idr_options_inter_replicate_threshold"]
    try:
        module = P.snip(IDR.__file__, ".py")
    except ValueError:
        module = P.snip(IDR.__file__, ".pyc")

    P.submit(module,
             "findNPeaks",
             params=[str(idr_thresh), ],
             infiles=infiles,
             outfiles=outfile)


@transform(findNPeaksForIndividualReplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadNPeaksForIndividualReplicates(infile, outfile):
    P.load(infile, outfile)


@follows(runIDROnPseudoreplicates)
@merge("./idr_pseudoreplicates/*-overlapped-peaks.txt",
       "./idr_pseudoreplicates/pseudoreplicates.nPeaks.tsv")
def findNPeaksForPseudoreplicates(infiles, outfile):
    idr_thresh = PARAMS["idr_options_self_consistency_threshold"]
    try:
        module = P.snip(IDR.__file__, ".py")
    except ValueError:
        module = P.snip(IDR.__file__, ".pyc")

    P.submit(module,
             "findNPeaks",
             params=[str(idr_thresh), ],
             infiles=infiles,
             outfiles=outfile)


@transform(findNPeaksForPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadNPeaksForPseudoreplicates(infile, outfile):
    P.load(infile, outfile)


@follows(runIDROnPooledPseudoreplicates)
@merge("./idr_pooled_pseudoreplicates/*-overlapped-peaks.txt",
       "./idr_pooled_pseudoreplicates/pooled_pseudoreplicates.nPeaks.tsv")
def findNPeaksForPooledPseudoreplicates(infiles, outfile):
    idr_thresh = PARAMS["idr_options_pooled_consistency_threshold"]
    try:
        module = P.snip(IDR.__file__, ".py")
    except ValueError:
        module = P.snip(IDR.__file__, ".pyc")

    P.submit(module,
             "findNPeaks",
             params=[str(idr_thresh), ],
             infiles=infiles,
             outfiles=outfile)


@transform(findNPeaksForPooledPseudoreplicates,
           regex("(.+)/(.+).tsv"),
           r"./\2.load")
def loadNPeaksForPooledPseudoreplicates(infile, outfile):
    P.load(infile, outfile)

##########################################################################
##########################################################################
##########################################################################
# Generate consensus peak set
##########################################################################


@follows(mkdir("bamfiles_final"))
@collate(filterBamfiles,
         regex(r"(.+)/(.+)-((?!input).*)-R[0-9]+.sentinel"),
         r"./bamfiles_final/\2-\3-R0.sentinel")
def reMergeBamfiles(infiles, sentinel):
    infiles = [P.snip(x, ".sentinel") + ".bam" for x in infiles]
    outfile = P.snip(sentinel, ".sentinel") + ".bam"
    bad_samples = PARAMS["options_to_remove"].split(",")

    to_merge = IDR.filterBadLibraries(infiles, bad_samples)

    IDR.mergeBams(to_merge, outfile)
    P.touch(sentinel)


@follows(reMergeBamfiles, mkdir("peakfiles_final"))
@transform("./bamfiles_final/*.bam",
           regex("(.+)/(.+).bam"),
           r"peakfiles_final/\2.narrowPeak.sentinel")
def callPeaksOnPooledReplicates(infile, outfile):
    # fetch peak calling parameters
    PARAMS_PEAKCALLER = get_peak_caller_parameters(
        PARAMS["options_peak_caller"])

    # call peaks on pseudoreplicates
    IDR.callIDRPeaks(infile,
                     outfile,
                     PARAMS["options_peak_caller"],
                     PARAMS["options_control_type"],
                     PARAMS_PEAKCALLER,
                     pseudoreplicate=False)

    P.touch(outfile)


@follows(callPeaksOnPooledReplicates,
         loadNPeaksForIndividualReplicates,
         loadNPeaksForPseudoreplicates,
         loadNPeaksForPooledPseudoreplicates,
         mkdir("peakfiles_final_conservative"),
         mkdir("peakfiles_final_optimum"))
@split("./peakfiles_final/*.narrowPeak.gz",
       regex("(.+)/(.+)-R0_VS_(.+)-R0_peaks.narrowPeak.gz"),
       [r"peakfiles_final_conservative/\2.narrowPeak.gz",
        r"peakfiles_final_optimum/\2.narrowPeak.gz"])
def generatePeakSets(infile, outfiles):
    outf_con, outf_opt = outfiles

    # retrieve maximum number of peaks obtained from inter-replicate IDR
    # (table created by loadNPeaksForIndividualReplicates)
    statement = ("SELECT"
                 " Experiment,"
                 " max(n_peaks) AS nPeaks"
                 " FROM individual_replicates_nPeaks"
                 " GROUP BY experiment")
    df = PU.fetch_DataFrame(statement)
    # reassign experiment as index
    df = df.set_index("Experiment")

    # retrieve number of peaks obtained from pooled_pseudoreplicate IDR
    # (table created by loadNPeaksForPooledPseudoreplicates)
    statement = ("SELECT"
                 " Experiment,"
                 " n_peaks AS nPeaks"
                 " FROM pooled_pseudoreplicates_nPeaks")
    df2 = PU.fetch_DataFrame(statement)

    # reassign experiment as index
    df2 = df2.set_index("Experiment")

    # split the infile name to obtain experiment
    sample_id = os.path.basename(infile).split("_VS_")[0]
    sample = sample_id.split("-")
    experiment = "_".join([sample[0], sample[1]])

    # retrieve max_numPeaks for experiment
    nPeaks = int(df.loc[experiment])
    # retrieve numPeaks_Rep0 for experiment
    nPeaks_rep0 = int(df2.loc[experiment])
    # retrieve maximumn of the two
    nPeaks_max = max(nPeaks, nPeaks_rep0)

    # establish which column to sort by
    if PARAMS["idr_options_ranking_measure"] == "signal.value":
        sort_statement = "sort -k7nr,7nr"
    elif PARAMS["idr_options_ranking_measure"] == "p.value":
        sort_statement = "sort -k8nr,8nr"
    elif PARAMS["idr_options_ranking_measure"] == "q.value":
        sort_statement = "sort -k9nr,9nr"
    else:
        raise ValueError("Unrecognised ranking_measure"
                         " %s don't know which column"
                         " to sort on" % PARAMS["idr_options_ranking_measure"])

    # sort infile by column and write top nPeaks to outfile (conservative)
    ignore_pipe_errors = True
    statement = ("zcat %(infile)s |"
                 " %(sort_statement)s |"
                 " head -%(nPeaks)s |"
                 " gzip > %(outf_con)s")
    P.run()

    # sort infile by column and write top nPeaks_max to outfile (optimum)
    ignore_pipe_errors = True
    statement = ("zcat %(infile)s |"
                 " %(sort_statement)s |"
                 " head -%(nPeaks_max)s |"
                 " gzip > %(outf_opt)s")
    P.run()

##########################################################################
##########################################################################
##########################################################################
# Clean Up!
##########################################################################


@merge(("./bamfiles_filtered/*.bam",
        "./bamfiles_final/*.bam",
        "./bamfiles_pooled/*.bam",
        "./bamfiles_pseudoreplicates/*.bam",
        "./bamfiles_pooled_pseudoreplicates/*.bam"),
       "./bamfiles_removed.sentinel")
def removeBamfiles(infiles, outfile):
    for bamfile in infiles:
        bam_index = bamfile + ".bai"
        os.unlink(bamfile)
        if os.path.exists(bam_index):
            os.unlink(bam_index)
    P.touch(outfile)


@transform(("./bamfiles_filtered/*.bdg",
            "./bamfiles_pooled/*.bdg",
            "./bamfiles_pseudoreplicates/*.bdg",
            "./bamfiles_pooled_pseudoreplicates/*.bdg"),
           suffix(".bdg"),
           ".bw")
def convertBedGraph(infile, outfile):
    contig_file = os.path.join(PARAMS["annotations_dir"], "contigs.tsv")
    statement = ("bedGraphToBigWig %(infile)s %(contig_file)s %(outfile)s")
    P.run()

##########################################################################
##########################################################################
##########################################################################
# primary targets
##########################################################################


@follows(splitBamfiles,
         poolInputBamfiles,
         splitPooledBamfiles)
def preProcessBamfiles():
    pass


@follows(callPeaksOnIndividualReplicates,
         callPeaksOnPseudoreplicates,
         callPeaksOnPooledPseudoreplicates,
         loadPeakSummaryForIndividualReplicates,
         loadPeakSummaryForPseudoreplicates,
         loadPeakSummaryForPooledPseudoreplicates)
def callPeaks():
    pass


@follows(plotBatchConsistencyForIndividualReplicates,
         plotBatchConsistencyForPseudoreplicates,
         plotBatchConsistencyForPooledPseudoreplicates)
def runIDR():
    pass


@follows(loadIDROnIndividualReplicates,
         loadIDROnPseudoreplicates,
         loadIDROnPooledPseudoreplicates,
         loadNPeaksForIndividualReplicates,
         loadNPeaksForPseudoreplicates,
         loadNPeaksForPooledPseudoreplicates)
def summarizeIDR():
    pass


@follows(preProcessBamfiles,
         callPeaks,
         runIDR,
         summarizeIDR,
         generatePeakSets)
def full():
    pass


@follows(removeBamfiles,
         convertBedGraph)
def cleanUp():
    pass

##########################################################################
##########################################################################
##########################################################################
# primary targets
##########################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish():
    '''publish report and data.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit(P.main(sys.argv))
