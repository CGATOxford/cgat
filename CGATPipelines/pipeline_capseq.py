"""
===================
CAPseq pipeline
===================

:Author: David Sims 
:Release: $Id: pipeline_capseq.py 2900 2011-05-24 14:38:00Z david $
:Date: |today|
:Tags: Python

The CAPseq pipeline imports reads from one or more CAPseq experiments and
performs the following tasks:

   * Align reads to the genome using Bowtie
   * Call peaks using MACS
   * Compare intervals across tracks

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_capseq.ini` file. The pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>` documentation).
The configuration file is organized by section and the variables are documented within 
the file. In order to get a local configuration file in the current directory, type::

    python <codedir>/pipeline_cpg.py config

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.


Input
-----

Reads
++++++

Input are :file:`.fastq.gz`-formatted files. The files should be
labeled in the following way::

   sample-condition-replicate.fastq.gz

Note that neither ``sample``, ``condition`` or ``replicate`` should contain 
``_`` (underscore) and ``.`` (dot) characters as these are used by the pipeline
to delineate tasks.

Requirements
------------

The pipeline requires the information from the following pipelines:

:doc:`pipeline_annotations`

set the configuration variables:
   :py:data:`annotations_database` 
   :py:data:`annotations_dir`

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie              |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|MACS                |14                 |peak finding                                    |
+--------------------+-------------------+------------------------------------------------+
|Picard              |>=1.4              |Dupliate removal, mapping stats                 |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |interval comparison                             |
+--------------------+-------------------+------------------------------------------------+


Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_capseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_capseq.tgz
   tar -xvzf pipeline_capseq.tgz
   cd pipeline_capseq.dir
   python <srcdir>/pipeline_capseq.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   bowtie
      bowtie_ - a read mapper

.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml


Code
====

"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import gzip
import sqlite3
import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import cStringIO
import numpy
import CGAT.Masker as Masker
import fileinput
import CGAT.Experiment as E
import logging as L
import PipelineChipseq as PIntervals
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMapping as PipelineMapping
from ruffus import *
from rpy2.robjects import r as R
import rpy2.robjects as ro

###################################################################
###################################################################
###################################################################
# Pipeline configuration
import CGAT.Pipeline as P
P.getParameters("pipeline_capseq.ini")
PARAMS = P.PARAMS
USECLUSTER = True

###################################################################
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.Sample3

TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    [x for x in glob.glob(
        "*.export.txt.gz") if PARAMS["tracks_control"] not in x],
    "(\S+).export.txt.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.sra") if PARAMS["tracks_control"] not in x],
        "(\S+).sra" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob(
            "*.fastq.gz") if PARAMS["tracks_control"] not in x],
        "(\S+).fastq.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob(
            "*.fastq.1.gz") if PARAMS["tracks_control"] not in x],
        "(\S+).fastq.1.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob(
            "*.csfasta.gz") if PARAMS["track_control"] not in x],
        "(\S+).csfasta.gz")
for X in TRACKS:
    print "TRACK=", X, "\n"

TRACKS_CONTROL = PipelineTracks.Tracks(Sample).loadFromDirectory(
    [x for x in glob.glob("*.export.txt.gz") if PARAMS["tracks_control"] in x],
    "(\S+).export.txt.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.sra") if PARAMS["tracks_control"] in x],
        "(\S+).sra" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.fastq.gz") if PARAMS["tracks_control"] in x],
        "(\S+).fastq.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob(
            "*.fastq.1.gz") if PARAMS["tracks_control"] in x],
        "(\S+).fastq.1.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.csfasta.gz") if PARAMS["track_control"] in x],
        "(\S+).csfasta.gz")
for X in TRACKS_CONTROL:
    print "TRACK_CONTROL=", X, "\n"


def getControl(track):
    '''return appropriate control for a track'''
    n = track.clone()
    n.condition = PARAMS["tracks_control"]
    return n

###################################################################
###################################################################
###################################################################
# aggregate per experiment
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
# aggregate per condition
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition",))
# aggregate over tissues and replicates
SAMPLES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", "replicate"))
# aggregate per tissue
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue",))
# compound targets : all experiments
TRACKS_MASTER = EXPERIMENTS.keys() + CONDITIONS.keys()
# compound targets : correlation between tracks
TRACKS_CORRELATION = TRACKS_MASTER + list(TRACKS)

###################################################################
###################################################################
###################################################################
# Section 1: MAP READS


@follows(mkdir("bam"))
@transform(("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz"),
           regex(r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"),
           r"bam/\1.bam")
def buildBAM(infile, outfile):
    '''map reads with bowtie'''
    to_cluster = True
    track = P.snip(os.path.basename(outfile), ".bam")
    job_options = "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    m = PipelineMapping.Bowtie()
    reffile = PARAMS["samtools_genome"]
    statement = m.build((infile,), outfile)
    P.run()

#########################################################################


@transform(buildBAM, suffix(".bam"), ".alignstats")
def buildPicardAlignStats(infile, outfile):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = True
    track = P.snip(os.path.basename(infile), ".bam")
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals(
    )
    P.run()

############################################################


@merge(buildPicardAlignStats, "bam/picard_align_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".alignstats")
        if not os.path.exists(f):
            E.warn("File %s missing" % f)
            continue
        lines = [
            x for x in open(f, "r").readlines() if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        for i in range(1, len(lines)):
            outf.write("%s\t%s" % (track, lines[i]))
    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable(outfile)
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s'''
    P.run()
    os.unlink(tmpfilename)

#########################################################################


@transform(buildBAM, suffix(".bam"), ".gcstats")
def buildPicardGCStats(infile, outfile):
    '''Gather BAM file GC bias stats using Picard '''
    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".bam")
    statement = '''CollectGcBiasMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s 
                   OUTPUT=%(outfile)s CHART_OUTPUT=%(outfile)s.pdf SUMMARY_OUTPUT=%(outfile)s.summary VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################


@merge(buildPicardGCStats, "bam/picard_gcbias_stats.load")
def loadPicardGCStats(infiles, outfile):
    '''Merge Picard insert size stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".dedup.gcstats")
        if not os.path.exists(f):
            E.warn("File %s missing" % f)
            continue
        lines = [
            x for x in open(f, "r").readlines() if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable(outfile)
    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()
    os.unlink(tmpfilename)

#########################################################################


@transform(buildBAM, suffix(".bam"), ".readstats")
def buildBAMStats(infile, outfile):
    '''Count number of reads mapped, duplicates, etc. using bam2stats.py'''
    to_cluster = USECLUSTER
    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''python %(scriptsdir)s/bam2stats.py --force 
                   --output-filename-pattern=%(outfile)s.%%s < %(infile)s > %(outfile)s'''
    P.run()

#########################################################################


@merge(buildBAMStats, "bam/bam_stats.load")
def loadBAMStats(infiles, outfile):
    '''Import bam statistics into SQLite'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join([P.snip(os.path.basename(x), ".readstats")
                      for x in infiles])
    filenames = " ".join(["<( cut -f 1,2 < %s)" % x for x in infiles])
    tablename = P.toTable(outfile)
    E.info("loading bam stats - summary")
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s"""
    P.run()

    for suffix in ("nm",):
        E.info("loading bam stats - %s" % suffix)
        filenames = " ".join(["%s.%s" % (x, suffix) for x in infiles])
        tname = "%s_%s" % (tablename, suffix)

        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()


#########################################################################
@transform(buildBAM, suffix(".bam"), ".dedup.bam")
def dedup(infiles, outfile):
    '''Remove duplicate alignments from BAM files.'''
    to_cluster = USECLUSTER
    track = P.snip(outfile, ".bam")
    dedup_method = PARAMS["dedup_method"]
    if dedup_method == 'samtools':
        statement = '''samtools rmdup %(infiles)s %(outfile)s; ''' % locals()
    elif dedup_method == 'picard':
        statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s 
                           METRICS_FILE=%(track)s.dupstats REMOVE_DUPLICATES=true 
                           VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################


@merge(dedup, "bam/picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = open('dupstats.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".dedup.bam")
        statfile = P.snip(f, ".bam") + ".dupstats"
        if not os.path.exists(statfile):
            E.warn("File %s missing" % statfile)
            continue
        lines = [x for x in open(
            statfile, "r").readlines() if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))

    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable(outfile)
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s '''
    P.run()

############################################################


@transform(dedup, suffix(".dedup.bam"), ".dedup.mapped.bam")
def removeUnmapped(infile, outfile):
    ''' Remove unmapped reads from BAM file'''
    to_cluster = True
    statement = '''samtools view -F 4 -bh %(infile)s > %(outfile)s; 
                   samtools index %(outfile)s; ''' % locals()
    P.run()

############################################################


@follows(removeUnmapped)
@files([(["bam/%s.dedup.mapped.bam" % x for x in SAMPLES[y]] + ["bam/%s.dedup.mapped.bam" % [getControl(x) for x in SAMPLES[y]][0]],
         ["bam/%s.norm.bam" % x for x in SAMPLES[y]] + ["bam/%s.norm.bam" %
                                                        [getControl(x) for x in SAMPLES[y]][0]]
         ) for y in SAMPLES])
def normaliseBAMs(infiles, outfiles):
    '''Within each "sample" where a sample represents a biological object in real life and a combination of
    tissue and replicate labels in the pipeline, Sample reads from the BAM for each track such that all have the same
    number of reads.'''
    to_cluster = True

    # Count reads in chip file
    countfile = {}
    reads = {}

    for infile in infiles:

        countfile[infile] = infile.replace(".bam", ".count")
        statement = ''' samtools idxstats %s | awk '{s+=$3} END {print s}' > %s ''' % (
            infile, countfile[infile])
        print statement
        P.run()
        fh = open(countfile[infile], "r")
        reads[infile] = int(fh.read())
        fh.close()

    smallest_bam = [x for x in reads if reads[x] == min(reads.values())][0]

    for infile in infiles:
        if infile == smallest_bam:
            smallest_bam_out = "%s.norm.bam" % P.snip(
                smallest_bam, ".dedup.mapped.bam")
            statement = ''' cp %(smallest_bam)s %(smallest_bam_out)s; cp %(smallest_bam)s.bai %(smallest_bam_out)s.bai ''' % locals(
            )
            P.run()
        else:
            PIntervals.buildSimpleNormalizedBAM((infile, countfile[infile]), "%s.norm.bam" % P.snip(
                infile, ".dedup.mapped.bam"), reads[smallest_bam])

############################################################


@follows(normaliseBAMs, mkdir("merged_bams"))
@files([(["bam/%s.norm.bam" % y for y in EXPERIMENTS[x]],
         "merged_bams/%s.merge.bam" % str(x).replace("-agg", "")) for x in EXPERIMENTS if len(EXPERIMENTS[x]) > 1])
def mergeReplicateBAMs(infiles, outfile):
    '''Merge normalised BAM files for all replicates, then sort and index. '''
    track = P.snip(outfile, ".merge.bam")
    in_list = " ".join(infiles)
    statement = '''samtools merge %(track)s.rep.bam %(in_list)s;
                   samtools sort  %(track)s.rep.bam %(track)s.merge; 
                   samtools index %(outfile)s;
                   rm %(track)s.rep.bam; ''' % locals()
    P.run()

############################################################


@transform(mergeReplicateBAMs, suffix(".bam"), ".bw")
def getMergedBigWig(infile, outfile):
    '''Merge multiple BAM files per replicate to produce a single non peak-shifted bigwig file'''

    to_cluster = True
    job_options = " -l mem_free=2G"
    statement = '''python %(scriptsdir)s/bam2wiggle.py --output-format=bigwig %(infile)s %(outfile)s '''
    P.run()

############################################################


@follows(normaliseBAMs)
@files([(["bam/%s.norm.bam" % y for y in EXPERIMENTS[x]],
         "merged_bams/%s.shift.merge.bw" % str(x).replace("-agg", "")) for x in EXPERIMENTS])
def getMergedBigWigPeakShift(infiles, outfile):
    '''Merge multiple BAM files per replicate to produce a single peak-shifted bigwig file'''
    expt = P.snip(os.path.basename(outfile), ".merge.bw").replace("-agg", "")
    in_list = " --bamfile=".join(infiles)

    offsets = []
    for t in infiles:
        track = P.snip(os.path.basename(t), ".norm.bam")
        fn = "macs/with_input/%s.macs" % track
        if os.path.exists(fn):
            offsets.append(str(PIntervals.getPeakShiftFromMacs(fn)))

    shifts = " --shift=".join(offsets)
    statement = '''python %(scriptsdir)s/bam2wiggle.py 
                      --output-format=bigwig
                      %(in_list)s
                      %(shifts)s > %(outfile)s'''
    P.run()

############################################################


@follows(mergeReplicateBAMs, mkdir("macs"), mkdir("macs/merged"))
@files([("merged_bams/%s.merge.bam" % str(x).replace("-agg", ""),
         "macs/merged/%s.merged.macs" % str(x).replace("-agg", "")) for x in EXPERIMENTS if len(EXPERIMENTS[x]) > 1])
def runMacsMerged(infile, outfile):
    '''Run MACS for peakshifted wig generation'''
    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".merge.bam")
    statement = '''cd macs/merged/; 
                   macs14 -t ../../%(infile)s 
                          --name=%(track)s.merged
                          --format=BAM
                          --wig -S
                          %(macs_options)s 
                   >& ../../%(outfile)s;
                   cd ../..;'''
    P.run()

############################################################


@follows(runMacsMerged)
@transform("macs/merged/*/treat/*.merged_treat_afterfiting_all.wig.gz", suffix(".wig.gz"), ".macs.norm.wig.gz")
def normaliseWig(infile, outfile):
    '''Run MACS for peakshifted wig generation'''
    track = P.snip(
        os.path.basename(infile), ".merged_treat_afterfiting_all.wig.gz")
    to_cluster = USECLUSTER
    dbhandle = sqlite3.connect(PARAMS["database"])
    cc = dbhandle.cursor()
    query = '''SELECT sum(PF_HQ_ALIGNED_READS)
               FROM picard_align_stats
               WHERE track LIKE "%(track)s%%"''' % locals()
    cc.execute(query)
    total_aligned_reads = cc.fetchone()[0]
    norm = total_aligned_reads / 1e7
    cc.close()

    # Read in wig and divide each value by norm
    outs = gzip.open(outfile, "wb")
    ins = gzip.open(infile, "rb")
    for line in ins:
        if re.match("^\d+", line):
            p, s = line.split()
            pos = int(p)
            score = float(s)
            norm_score = score / norm
            outs.write("%i\t%0.2f\n" % (pos, norm_score))
        else:
            outs.write(line)
    ins.close()
    outs.close()

############################################################
############################################################
############################################################
# BUILD INTERVALS USING CONTROL SAMPLE


@follows(normaliseBAMs, mkdir("macs"), mkdir("macs/with_input"))
@files([(("bam/%s.norm.bam" % x,
          "bam/%s.norm.bam" % getControl(x)),
         "macs/with_input/%s.macs" % x) for x in TRACKS])
def runMACS(infiles, outfile):
    '''Run MACS for peak detection using control sample.'''
    infile, controlfile = infiles
    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".norm.bam")
    dir = os.path.dirname(infile)

    # change to macs directory and run MACS from there so that wig files end
    # up in right place
    statement = '''cd macs/with_input; 
                   macs14 -t ../../%(infile)s 
                          -c ../../%(controlfile)s
                          --name=%(track)s
                          --format=BAM
                          --diag
                          --wig -S
                          %(macs_options)s 
                   >& ../../%(outfile)s;
                   cd ../..;'''
    P.run()

############################################################


@transform(runMACS, regex(r"macs/with_input/(\S+).macs"),
           inputs((r"macs/with_input/\1.macs", r"bam/\1.norm.bam")),
           r"macs/with_input/\1_macs.load")
def loadMACS(infiles, outfile):
    '''Load MACS intervals into database and filter on fold change and p-value'''
    infile, bamfile = infiles
    PIntervals.loadMACS(infile, outfile, bamfile)

############################################################


@merge(runMACS, "macs/with_input/macs.summary")
def summarizeMACS(infiles, outfile):
    '''Parse MACS summary statistics from log file'''
    PIntervals.summarizeMACS(infiles, outfile)

############################################################


@transform(summarizeMACS, suffix(".summary"), "_summary.load")
def loadMACSSummary(infile, outfile):
    '''load macs summary into database'''
    P.load(infile, outfile, "--index=track")

############################################################


@follows(mkdir("intervals"))
@transform(loadMACS, regex(r"macs/with_input/(\S+)_macs.load"), r"intervals/\1.bed")
def exportIntervalsAsBed(infile, outfile):
    '''Export MACS intervals from database as BED file and filter on fold change'''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed(infile, outfile, fc)

############################################################


@follows(runMACS)
@transform("macs/with_input/*_MACS_wiggle/treat/*.wig.gz", suffix(".wig.gz"), ".bw")
def convertWigToBigWig(infile, outfile):
    '''Convert MACS wig file to bigwig file'''
    statement = '''wigToBigWig %(infile)s %(samtools_genome)s.fa.fai %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
# BUILD INTERVALS WITHOUT CONTROL SAMPLE


@follows(normaliseBAMs, mkdir("macs"), mkdir("macs/no_input"))
@files([("bam/%s.norm.bam" % x,
         "macs/no_input/%s.solo.macs" % x) for x in TRACKS])
def runMACSsolo(infile, outfile):
    '''Run MACS for peak detection.'''
    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".norm.bam")
    statement = '''cd macs/no_input/; 
                   macs14 -t ../../%(infile)s 
                          --name=%(track)s.solo
                          --format=BAM
                          --diag
                          --wig -S
                          %(macs_options)s 
                   >& ../../%(outfile)s;
                   cd ../..;'''
    P.run()

############################################################


@transform(runMACSsolo, regex(r"macs/no_input/(\S+).solo.macs"),
           inputs((r"macs/no_input/\1.solo.macs", r"bam/\1.norm.bam")),
           r"macs/no_input/\1.solo_macs.load")
def loadMACSsolo(infiles, outfile):
    infile, bamfile = infiles
    PIntervals.loadMACS(infile, outfile, bamfile)

############################################################


@merge(runMACSsolo, "macs/no_input/macs_solo.summary")
def summarizeMACSsolo(infiles, outfile):
    '''run MACS for peak detection.'''
    PIntervals.summarizeMACSsolo(infiles, outfile)

############################################################


@transform(summarizeMACSsolo, suffix(".summary"), "_summary.load")
def loadMACSsoloSummary(infile, outfile):
    '''load macs summary.'''
    P.load(infile, outfile, "--index=track")

############################################################


@transform(loadMACSsolo, regex(r"macs/no_input/(\S+).macs.load"), r"intervals/\1.bed")
def exportIntervalsAsBedsolo(infile, outfile):
    '''Export list of intervals passing fold change threshold to file '''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed(infile, outfile, fc)

############################################################


@follows(runMACSsolo)
@transform("macs/no_input/*_MACS_wiggle/treat/*.wig.gz", suffix(".wig.gz"), ".bw")
def convertWigToBigWigSolo(infile, outfile):
    '''Convert MACS wig file to bigwig file'''
    statement = '''wigToBigWig %(infile)s %(samtools_genome)s.fa.fai %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
# Merge nearby intervals


@transform((exportIntervalsAsBed), suffix(".bed"), ".merged.bed")
def mergeIntervals(infile, outfile):
    '''Merge intervals less than n bases apart in each dataset and update foldchange scores'''
    d = PARAMS["intervals_merge_dist"]
    method = PARAMS["intervals_foldchange_merge_method"]
    PIntervals.mergeIntervalsWithScores(infile, outfile, d, method)

############################################################


@transform(mergeIntervals, suffix(".merged.bed"), ".merged.cleaned.bed")
def sanitiseIntervals(infile, outfile):
    '''sanatise so that intervals do not exceed contig length'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s -L %(outfile)s.log > %(outfile)s;'''
    P.run()

############################################################


@transform((sanitiseIntervals), suffix(".merged.cleaned.bed"), ".merged.cleaned.bed.load")
def loadMergedIntervals(infile, outfile):
    '''load combined intervals.

    Also, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval

    If *replicates* is true, only replicates will be considered
    for the counting. Otherwise the counts aggregate both replicates
    and conditions.
    '''

    # Write header to output file
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    headers = ("contig", "start", "end", "interval_id", "nPeaks",
               "PeakCenter", "Length", "AvgVal", "PeakVal", "nProbes", "Fold")
    tmpfile.write("\t".join(headers) + "\n")
    contig, start, end, interval_id, npeaks, peakcenter, length, avgval, peakval, nprobes = "", 0, 0, 0, 0, 0, 0, 0, 0, 0

    # Get SAM file and Macs offset
    samfiles, offsets = [], []
    track = P.snip(os.path.basename(infile), ".merged.cleaned.bed")
    base_track = track.replace(".solo", "")

    fn = "bam/%s.norm.bam" % track
    assert os.path.exists(
        fn), "could not find bamfile %s for track %s" % (fn, track)
    samfiles.append(pysam.Samfile(fn,  "rb"))
    if track.find("solo") > -1:
        fn = "macs/no_input/%s.macs" % track
    else:
        fn = "macs/with_input/%s.macs" % track
    if os.path.exists(fn):
        offsets.append(PIntervals.getPeakShiftFromMacs(fn))

    # Loop over input Bed file and calculate stats for merged intervals
    c = E.Counter()
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, int_id, fc = line[:-1].split()[:5]
        start, end = int(start), int(end)
        interval_id = c.input

        npeaks, peakcenter, length, avgval, peakval, nprobes = PIntervals.countPeaks(
            contig, start, end, samfiles, offsets)

        # nreads can be 0 if the intervals overlap only slightly
        # and due to the binning, no reads are actually in the overlap region.
        # However, most of these intervals should be small and have already be deleted via
        # the merge_min_interval_length cutoff.
        # do not output intervals without reads.
        if nprobes == 0:
            c.skipped_reads += 1

        c.output += 1
        tmpfile.write("\t".join(
            map(str, (contig, start, end, int_id, npeaks, peakcenter, length, avgval, peakval, nprobes, fc))) + "\n")

    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_macs_merged_intervals" % track

    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --index=interval_id
                       --index=contig,start 
                       --table=%(tablename)s
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink(tmpfile.name)
    L.info("%s\n" % str(c))

############################################################
############################################################
############################################################
# Assess background (non-peak) binding


@follows(sanitiseIntervals, mkdir("background"))
@files([(("bam/%s.norm.bam" % x,
          "intervals/%s.merged.cleaned.bed" % x),
         "background/%s.bg" % x) for x in TRACKS])
def getBackground(infiles, outfile):
    '''Count the number of reads in the bamfile used for MACS that do not overlap an interval'''
    bam, bed = infiles
    to_cluster = True
    statement = '''intersectBed -abam %(bam)s -b %(bed)s -v -bed | wc -l > %(outfile)s; '''
    statement += '''intersectBed -abam %(bam)s -b %(bed)s -u -bed | wc -l >> %(outfile)s;'''
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################


@transform(getBackground, suffix(".bg"), ".bg.load")
def loadBackground(infile, outfile):
    '''load background into database'''
    track = P.snip(os.path.basename(infile), ".bg")
    header = "out_peaks,in_peaks"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=%(track)s_background
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
# Assess effect of altering fold change threshold


@follows(exportIntervalsAsBed, mkdir("foldchange"))
@files([("macs/with_input/%s.macs" % x,
         "foldchange/%s.foldchange" % x) for x in TRACKS])
def thresholdFoldChange(infile, outfile):
    '''Assess interval overlap between conditions at different fold change thresholds '''

    # Connect to DB
    dbhandle = sqlite3.connect(PARAMS["database"])

    # Make bed files for different fold change thresholds
    track = P.snip(os.path.basename(infile), ".macs").replace("-", "_")
    macsdir = os.path.dirname(infile)
    foldchange = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    for fc in foldchange:
        cc = dbhandle.cursor()
        query = '''SELECT contig, start, end, interval_id, fold FROM %(track)s_macs_intervals WHERE fold > %(fc)i ORDER by contig, start;''' % locals(
        )
        cc.execute(query)
        outbed = "foldchange/" + track + ".fc" + str(fc) + ".bed"
        outs = open(outbed, "w")

        for result in cc:
            contig, start, end, interval_id, fold = result
            outs.write("%s\t%i\t%i\t%s\t%i\n" %
                       (contig, start, end, str(interval_id), fold))
        cc.close()
        outs.close()

        statement = '''echo %(fc)i >> %(outfile)s; cat %(outbed)s | wc -l >> %(outfile)s; '''
        P.run()

    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################


@transform(thresholdFoldChange, suffix(".foldchange"), ".foldchange.load")
def loadFoldChangeThreshold(infile, outfile):
    '''Load intervals overlapping chipseq into database '''
    track = P.snip(os.path.basename(infile), ".foldchange").replace(
        ".", "_").replace("-", "_")
    header = "threshold,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_foldchange
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################


@transform(thresholdFoldChange, suffix(".foldchange"), ".shared.foldchange")
def sharedIntervalsFoldChangeThreshold(infile, outfile):
    '''identify shared intervals between datasets at different foldchange thresholds'''

    # Open foldchange file and read
    fc_file = open(infile, "r")
    foldchange = []
    for line in fc_file:
        threshold, interval_count = line.split()
        foldchange.append(threshold)

    in_track = P.snip(os.path.basename(infile), ".foldchange")
    in_dir = os.path.dirname(infile)
    out_dir = in_dir.replace("macs", "fc")
    try:
        os.mkdir(out_dir)
    except OSError:
        pass

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    # for each foldchange
    for fc in foldchange:

        in_bed = "foldchange/" + \
            in_track.replace("-", "_") + ".fc" + str(fc) + ".bed"

        # For each track
        for track in TRACKS:
            if (str(track) != in_track):
                compare_bed = "foldchange/" + \
                    str(track).replace("-", "_") + ".fc" + str(fc) + ".bed"
                statement = '''echo %(track)s %(fc)s >> %(outfile)s; intersectBed -a %(in_bed)s -b %(compare_bed)s -u | wc -l >> %(outfile)s; '''
                P.run()

    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; sed -i '{s/ /\\t/g}' %(outfile)s; '''
    P.run()

############################################################


@transform(sharedIntervalsFoldChangeThreshold, suffix(".shared.foldchange"), ".shared.foldchange.load")
def loadSharedIntervalsFoldChangeThreshold(infile, outfile):
    '''Load intervals overlapping other tracks into database '''
    track = P.snip(os.path.basename(infile), ".shared.foldchange").replace(
        ".", "_").replace("-", "_")
    header = "track,threshold,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_foldchange_shared
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()


############################################################
############################################################
############################################################
# Calculate replicated intervals
@follows(sanitiseIntervals, mkdir("replicated_intervals"))
@files([(["intervals/%s.merged.cleaned.bed" % y for y in EXPERIMENTS[x]],
         "replicated_intervals/%s.rep.bed" % str(x).replace("-agg", "")) for x in EXPERIMENTS if len(EXPERIMENTS[x]) > 1])
def replicatedIntervals(infiles, outfile):
    '''Combine replicates between experiments.
       First all intervals are merged across replicates 
       and then merged intervals that do not overlap an interval in all replicates are removed. '''

    expt = P.snip(os.path.basename(outfile), ".rep.bed").replace("-agg", "")
    outdir = os.path.dirname(outfile)
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    # merge all files to get total set of intervals
    in_list = " ".join(infiles)
    statement = '''cat %(in_list)s | mergeBed -i stdin -scores max > %(outdir)s/%(expt)s.merged.bed; 
                   cat %(outdir)s/%(expt)s.merged.bed > %(outfile)s;''' % locals()
    P.run()

    # Remove intervals that do not overlap intervals in all reps
    for i in infiles:
        statement = '''intersectBed -a %(outfile)s -b %(i)s -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; '''
        P.run()


############################################################
@transform(replicatedIntervals, suffix(".rep.bed"), ".rep.bed.load")
def loadReplicatedIntervals(infile, outfile):
    '''load replicated intervals.

    Also, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval    '''

    track = P.snip(os.path.basename(infile), ".rep.bed")
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]

    # Write header to output file
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    headers = ("contig", "start", "end", "interval_id", "nPeaks",
               "PeakCenter", "Length", "AvgVal", "PeakVal", "nProbes", "Fold")
    tmpfile.write("\t".join(headers) + "\n")
    contig, start, end, interval_id, npeaks, peakcenter, length, avgval, peakval, nprobes = "", 0, 0, 0, 0, 0, 0, 0, 0, 0

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists(
            fn), "could not find bamfile %s for track %s" % (fn, str(t))
        samfiles.append(pysam.Samfile(fn,  "rb"))
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists(fn):
            offsets.append(PIntervals.getPeakShiftFromMacs(fn))

    # Loop over input Bed file and calculate stats for replicated intervals
    c = E.Counter()
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, fc = line[:-1].split()[:4]
        start, end = int(start), int(end)
        interval_id = c.input

        npeaks, peakcenter, length, avgval, peakval, nprobes = PIntervals.countPeaks(
            contig, start, end, samfiles, offsets)

        # nreads can be 0 if the intervals overlap only slightly
        # and due to the binning, no reads are actually in the overlap region.
        # However, most of these intervals should be small and have already be deleted via
        # the merge_min_interval_length cutoff.
        # do not output intervals without reads.
        if nprobes == 0:
            c.skipped_reads += 1

        c.output += 1
        tmpfile.write("\t".join(map(str, (contig, start, end, interval_id,
                      npeaks, peakcenter, length, avgval, peakval, nprobes, fc))) + "\n")

    tmpfile.close()
    tmpfilename = tmpfile.name

    # Load into database
    tablename = "%s_replicated_intervals" % track
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --index=interval_id
                       --index=contig,start 
                       --table=%(tablename)s
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink(tmpfile.name)
    L.info("%s\n" % str(c))

############################################################


@transform(loadReplicatedIntervals, suffix(".rep.bed.load"), ".replicated.bed")
def exportReplicatedIntervalsAsBed(infile, outfile):
    '''export locations for all intervals.'''
    dbhandle = sqlite3.connect(PARAMS["database"])
    track = P.snip(os.path.basename(infile), ".rep.bed.load").replace("-", "_")
    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id FROM %s_replicated_intervals ORDER by contig, start" % track
    cc.execute(statement)

    # Write to bed file
    outs = open(outfile, "w")
    for result in cc:
        contig, start, stop, interval_id = result
        outs.write("%s\t%i\t%i\t%i\n" % (contig, start, stop, interval_id))
    cc.close()
    outs.close()

############################################################
# Assess effect of altering fold change threshold


@follows(mkdir("foldchange"), mkdir("foldchange/replicated"))
@transform(loadReplicatedIntervals, regex(r"replicated_intervals/(\S+).rep.bed.load"), r"\1.replicated.foldchange")
def replicatedFoldChange(infile, outfile):
    '''Assess interval overlap between conditions at different fold change thresholds '''

    # Connect to DB
    dbhandle = sqlite3.connect(PARAMS["database"])

    # Make bed files for different fold change thresholds
    track = P.snip(os.path.basename(infile), ".rep.bed.load").replace("-", "_")
    macsdir = os.path.dirname(infile)
    foldchange = [6, 8, 10, 12]

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    for fc in foldchange:
        cc = dbhandle.cursor()
        query = '''SELECT contig, start, end, interval_id, fold FROM %(track)s_replicated_intervals WHERE fold > %(fc)i ORDER by contig, start;''' % locals(
        )
        cc.execute(query)
        outbed = "foldchange/replicated/" + track + ".fc" + str(fc) + ".bed"
        outs = open(outbed, "w")

        for result in cc:
            contig, start, end, interval_id, fold = result
            outs.write("%s\t%i\t%i\t%s\t%i\n" %
                       (contig, start, end, str(interval_id), fold))
        cc.close()
        outs.close()

        statement = '''echo %(fc)i >> %(outfile)s; cat %(outbed)s | wc -l >> %(outfile)s; '''
        P.run()

    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
############################################################
############################################################
# Compare intervals


@transform(sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.overlap")
def pairwiseIntervals(infile, outfile):
    '''identify overlapping intervals for each pair of datasets'''
    to_cluster = True
    in_track = P.snip(os.path.basename(infile), ".merged.cleaned.bed")
    statement = '''echo "track" > %(outfile)s; echo "overlap" >> %(outfile)s;'''
    P.run()
    for track in TRACKS:
        if (str(track) != in_track):
            statement = '''echo %(track)s >> %(outfile)s; intersectBed -a %(infile)s -b intervals/%(track)s.merged.cleaned.bed -u | wc -l >> %(outfile)s; '''
            P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################


@transform(pairwiseIntervals, suffix(".merged.overlap"), ".merged.overlap.load")
def loadPairwiseIntervals(infile, outfile):
    '''Load overlapping intervals into database '''
    track = P.snip(os.path.basename(infile), ".overlap").replace(
        ".", "_").replace("-", "_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_overlap
                   > %(outfile)s '''
    P.run()

############################################################


@transform(sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.unique.bed")
def uniqueIntervals(infile, outfile):
    '''identify unique intervals for each dataset'''
    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip(os.path.basename(infile), ".merged.cleaned.bed")
    for track in TRACKS:
        if str(track) != in_track:
            statement = '''intersectBed -a %(outfile)s -b intervals/%(track)s.merged.cleaned.bed -v > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s '''
            P.run()

############################################################


@transform(uniqueIntervals, suffix(".merged.unique.bed"), ".merged.unique.bed.load")
def loadUniqueIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip(os.path.basename(infile), ".unique.bed").replace(
        ".", "_").replace("-", "_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################


@transform(sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.shared.bed")
def sharedIntervals(infile, outfile):
    '''identify shared intervals between datasets'''
    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip(os.path.basename(infile), ".merged.cleaned.bed")
    for track in TRACKS:
        if str(track) != in_track:
            statement = '''intersectBed -a %(outfile)s -b intervals/%(track)s.merged.cleaned.bed -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; '''
            P.run()

############################################################


@transform(sharedIntervals, suffix(".merged.shared.bed"), ".merged.shared.bed.load")
def loadSharedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip(os.path.basename(infile), ".shared.bed").replace(
        ".", "_").replace("-", "_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################


@transform(exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.unique.bed")
def uniqueReplicatedIntervals(infile, outfile):
    '''identify unique intervals for each dataset'''
    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip(os.path.basename(infile), ".replicated.bed")
    for track in EXPERIMENTS:
        track = str(track).replace("-agg", "")
        if track != in_track:
            statement = '''intersectBed -a %(outfile)s -b replicated_intervals/%(track)s.replicated.bed -v > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s '''
            P.run()


###########################################################
@transform(uniqueReplicatedIntervals, suffix(".replicated.unique.bed"), ".replicated.unique.bed.load")
def loadUniqueReplicatedIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip(os.path.basename(infile), ".unique.bed").replace(
        ".", "_").replace("-", "_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################


@transform(exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.shared.bed")
def sharedReplicatedIntervals(infile, outfile):
    '''identify shared intervals between datasets'''
    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip(os.path.basename(infile), ".replicated.bed")
    for track in EXPERIMENTS:
        track = str(track).replace("-agg", "")
        if track != in_track:
            statement = '''intersectBed -a %(outfile)s -b replicated_intervals/%(track)s.replicated.bed -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; '''
            P.run()

############################################################


@transform(sharedReplicatedIntervals, suffix(".replicated.shared.bed"), ".replicated.shared.bed.load")
def loadSharedReplicatedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip(os.path.basename(infile), ".shared.bed").replace(
        ".", "_").replace("-", "_")
    header = "contig,start,stop,interval_id"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################


@transform((exportReplicatedIntervalsAsBed, sharedReplicatedIntervals, uniqueReplicatedIntervals), suffix(".bed"), ".length")
def exportCAPseqReplicatedLength(infile, outfile):
    '''Export length of CAPseq intervals'''
    statement = '''cat %(infile)s | awk '{print $3-$2}' > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
# REPORTS


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''
    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''
    E.info("updating documentation")
    P.run_report(clean=False)


@files("report.log", "publish.log")
def publish_report(infile, outfile):
    '''Link bed, bam, wig and report files to web '''
    publish_dir = PARAMS["publish_dir"]
    species = PARAMS["genome"]
    report_dir = os.path.join(publish_dir, species)
    bam_dir = os.path.join(publish_dir, "bam")
    bed_dir = os.path.join(publish_dir, "bed")
    wig_dir = os.path.join(publish_dir, "wig")
    length_dir = os.path.join(publish_dir, "length")
    working_dir = os.getcwd()
    statement = '''cp -rf report/html/* %(report_dir)s > %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/bam/*.norm.bam %(bam_dir)s >> %(outfile)s;'''
    statement += '''cp -sn %(working_dir)s/macs/with_input/*/treat/*.bw %(wig_dir)s >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/merged_bams/*.bw %(wig_dir)s >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.bed %(bed_dir)s/replicated >> %(outfile)s;'''
    statement += '''cp -sn %(working_dir)s/intervals/*solo*.bed %(bed_dir)s/no_input >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/intervals/*.merged.cleaned.bed %(bed_dir)s/replicates >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.length %(length_dir)s >> %(outfile)s; '''
    P.run()

############################################################
############################################################
############################################################
# Replication


@files("pipeline_capseq.ini", "conf.py")
def replicateMappingAndCalling(infile, outfile):
    '''Link to BAM files and MACS results from previous run'''
    # must copy and modify pipeline.ini first
    prev_run = PARAMS["previous_run"]
    statement = '''cp -ps %(prev_run)s/*.fastq.gz .;'''
    statement += '''cp %(prev_run)s/conf.py %(prev_run)s/sphinxreport.ini .; '''
    statement += '''mkdir bam; cp -psr %(prev_run)s/bam/* bam/; '''
    statement += '''mkdir macs; cp -psr %(prev_run)s/macs/* macs/; '''
    statement += '''mkdir external_bed; cp %(prev_run)s/external_bed/*.bed external_bed/;'''
    statement += '''rm *.load; rm bam/*.load; rm macs/*/*.load;'''
    P.run()

############################################################
############################################################
############################################################
# Pipeline organisation


@follows(buildBAM,
         buildPicardAlignStats, loadPicardAlignStats,
         buildBAMStats, loadBAMStats,
         dedup, loadPicardDuplicateStats)
def mapReads():
    '''Align reads to target genome.'''
    pass


@follows(removeUnmapped, normaliseBAMs)
def NormaliseBAMFiles():
    '''Normalise BAMS'''
    pass


@follows(mergeReplicateBAMs,
         getMergedBigWig,
         runMacsMerged,
         normaliseWig)
def exportBigWigs():
    '''Export merged BAM files per replictae in bigwig format'''
    pass


@follows(runMACS, loadMACS,
         summarizeMACS, loadMACSSummary,
         exportIntervalsAsBed)
def buildIntervalsMacs():
    '''Find peaks using MACS using input as control'''
    pass


@follows(runMACSsolo, loadMACSsolo,
         summarizeMACSsolo, loadMACSsoloSummary,
         exportIntervalsAsBedsolo)
def buildIntervalsMacsNoControl():
    '''Find peaks using MACS without control sample'''
    pass


@follows(mergeIntervals, loadMergedIntervals)
def mergePeaks():
    '''Merge nearby intervals within a single track'''
    pass


@follows(getBackground, loadBackground)
def background():
    '''Assess level of background (non-peak) binding'''
    pass


@follows(thresholdFoldChange, loadFoldChangeThreshold,
         sharedIntervalsFoldChangeThreshold,
         loadSharedIntervalsFoldChangeThreshold)
def FoldChangeThreshold():
    '''Assess number of intervals above different fold change thresholds'''
    pass


@follows(replicatedIntervals, loadReplicatedIntervals)
def findReplicatedIntervals():
    '''Intersect replicates to identify replicated intervals '''
    pass


@follows(pairwiseIntervals, loadPairwiseIntervals,
         uniqueIntervals, loadUniqueIntervals,
         sharedIntervals, loadSharedIntervals)
def comparePeaks():
    '''Compare intervals across tracks'''
    pass


@follows(uniqueReplicatedIntervals, loadUniqueReplicatedIntervals,
         sharedReplicatedIntervals, loadSharedReplicatedIntervals,
         exportCAPseqReplicatedLength)
def compareReplicatedPeaks():
    '''Compare intervals across tracks'''
    pass


@follows(mapReads,
         NormaliseBAMFiles,
         exportBigWigs,
         buildIntervalsMacs,
         buildIntervalsMacsNoControl,
         mergePeaks,
         background,
         FoldChangeThreshold,
         findReplicatedIntervals,
         compareReplicatedPeaks,
         comparePeaks)
def full():
    '''run the full pipeline.'''
    pass


@follows(mapReads,
         NormaliseBAMFiles,
         buildIntervalsMacs,
         buildIntervalsMacsNoControl,
         mergePeaks,
         background,
         FoldChangeThreshold,
         comparePeaks)
def chipseq():
    '''run the pipeline on chipseq data'''
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
