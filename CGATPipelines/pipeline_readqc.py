2##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Tildon Grant Belgard
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
====================
ReadQc pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

The readqc pipeline imports unmapped reads from one or more
fastq and performs basic quality control steps:

   1. per position quality
   2. per read quality
   3. duplicates

For further details see http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/

Additionaly, optional bias analysis can be included through the configuration
file. This analysis is designed to help identify sequence contexts which bias
gene expression and asssess the consistency in the biases between samples.

Bias analysis utilised Sailfish to estimate transcript abundance. This
requires a multi-fasta transcripts file.

For further details see http://www.cs.cmu.edu/~ckingsf/software/sailfish/


Individual tasks are enabled in the configuration file.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|fastqc              |>=0.9.0            |read quality control                            |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the quality of the sequence archive

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_readqc.tgz
   tar -xvzf pipeline_readqc.tgz
   cd pipeline_readqc
   python <srcdir>/pipeline_readqc.py make full


Code
====

"""

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random
import cStringIO

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import sys
import os
import re
import shutil
import string
import itertools
import math
import glob
import time
import gzip
import collections
import random
import numpy
import pandas
import sqlite3
from scipy.stats import linregress
from pandas import DataFrame
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.FastaIterator as FastaIterator
import CGAT.Tophat as Tophat
import rpy2.robjects as ro
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.Stats as Stats
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.Pipeline as P
import CGAT.Fastq as Fastq
import CGAT.CSV2DB as CSV2DB

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])
PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
# define input files
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

# AH: I would put these into the same configuration section
# Include optional bias analysis for RNA-Seq
BIAS_ANALYSIS = P.isTrue("bias_analysis")
# AH: No need to define a global variable here for "transcripts"
#     Keep globals to a minimum.
transcripts = PARAMS["sailfish_transcripts"]

#########################################################################
#########################################################################
#########################################################################


@follows(mkdir(PARAMS["exportdir"]), mkdir(os.path.join(PARAMS["exportdir"], "fastqc")))
@transform(INPUT_FORMATS,
           REGEX_FORMATS,
           r"\1.fastqc")
def runFastqc(infiles, outfile):
    '''convert sra files to fastq and check mapping qualities are in solexa format. 
    Perform quality control checks on reads from .fastq files.'''
    to_cluster = True
    m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"])
    statement = m.build((infiles,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################

def FastqcSectionIterator(infile):
    data = []
    for line in infile:
        if line.startswith(">>END_MODULE"):
            yield name, status, header, data
        elif line.startswith(">>"):
            name, status = line[2:-1].split("\t")
            data = []
        elif line.startswith("#"):
            header = "\t".join([x for x in line[1:-1].split("\t") if x != ""])
        else:
            data.append(
                "\t".join([x for x in line[:-1].split("\t") if x != ""]))


@jobs_limit(1, "db")
@transform(runFastqc, suffix(".fastqc"), "_fastqc.load")
def loadFastqc(infile, outfile):
    '''load FASTQC stats.'''

    track = P.snip(infile, ".fastqc")

    filename = os.path.join(
        PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")

    for fn in glob.glob(filename):
        prefix = os.path.basename(os.path.dirname(fn))
        results = []

        for name, status, header, data in FastqcSectionIterator(IOTools.openFile(fn)):
            # do not collect basic stats, see loadFastQCSummary
            if name == "Basic Statistics":
                continue

            parser = CSV2DB.buildParser()
            (options, args) = parser.parse_args([])
            options.tablename = prefix + "_" + re.sub(" ", "_", name)
            options.allow_empty = True

            inf = cStringIO.StringIO("\n".join([header] + data) + "\n")
            CSV2DB.run(inf, options)
            results.append((name, status))

        # load status table
        parser = CSV2DB.buildParser()
        (options, args) = parser.parse_args([])
        options.tablename = prefix + "_status"
        options.allow_empty = True

        inf = cStringIO.StringIO(
            "\n".join(["name\tstatus"] + ["\t".join(x) for x in results]) + "\n")
        CSV2DB.run(inf, options)

    P.touch(outfile)


def collectFastQCSections(infiles, section):
    '''iterate over all fastqc files and extract a particular section.'''
    results = []

    for infile in infiles:

        track = P.snip(infile, ".fastqc")

        filename = os.path.join(
            PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")

        for fn in glob.glob(filename):
            prefix = os.path.basename(os.path.dirname(fn))
            for name, status, header, data in FastqcSectionIterator(IOTools.openFile(fn)):
                if name == section:
                    results.append((track, status, header, data))

    return results


@merge(runFastqc, "status_summary.tsv.gz")
def buildFastQCSummaryStatus(infiles, outfile):
    '''load fastqc status summaries into a single table.'''

    outf = IOTools.openFile(outfile, "w")
    first = True
    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        filename = os.path.join(
            PARAMS["exportdir"], "fastqc", track + "*_fastqc", "fastqc_data.txt")

        for fn in glob.glob(filename):
            prefix = os.path.basename(os.path.dirname(fn))
            results = []

            names, stats = [], []
            for name, status, header, data in FastqcSectionIterator(IOTools.openFile(fn)):
                stats.append(status)
                names.append(name)

            if first:
                outf.write("track\tfilename\t%s\n" % "\t".join(names))
                first = False

            outf.write("%s\t%s\t%s\n" %
                       (track, os.path.dirname(fn), "\t".join(stats)))
    outf.close()


@merge(runFastqc, "basic_statistics_summary.tsv.gz")
def buildFastQCSummaryBasicStatistics(infiles, outfile):
    '''load fastqc summaries into a single table.'''

    data = collectFastQCSections(infiles, "Basic Statistics")

    outf = IOTools.openFile(outfile, "w")
    first = True
    for track, status, header, rows in data:
        rows = [x.split("\t") for x in rows]
        if first:
            headers = [row[0] for row in rows]
            outf.write("track\t%s\n" % "\t".join(headers))
            first = False
        outf.write("%s\t%s\n" % (track, "\t".join([row[1] for row in rows])))
    outf.close()


@transform((buildFastQCSummaryStatus, buildFastQCSummaryBasicStatistics),
           suffix(".tsv.gz"), ".load")
def loadFastqcSummary(infile, outfile):
    P.load(infile, outfile, options="--index=track")

####################################################
# bias analysis
####################################################
# AH: sections not Pep8 conformant
# AH: put all the conditional into one if statement


@active_if(BIAS_ANALYSIS)
@transform(PARAMS["sailfish_transcripts"],
           regex("(\S+)"),
           "index/transcriptome.sfi")
def indexForSailfish(infile, outfile):
    '''create a sailfish index'''

    outdir = P.snip(outfile, "/transcriptome.sfi")
    kmer = int(PARAMS["sailfish_kmer_size"])
    tmp = P.getTempFilename()

    statement = '''gunzip -c %(infile)s > %(tmp)s;
                   module load bio/sailfish;
                   sailfish index -t %(tmp)s 
                   -k %(kmer)i -o %(outdir)s;
                   rm -f %(tmp)s'''

    P.run()

@follows(indexForSailfish, mkdir("quantification"))
@transform(INPUT_FORMATS,
           REGEX_FORMATS,
           add_inputs(indexForSailfish),
           r"quantification/\1/\1_quant.sf")
def runSailfish(infiles, outfile):
    '''quantify abundance'''

    to_cluster = True
    job_options = "-pe dedicated %i -R y" % PARAMS["sailfish_threads"]

    infile, index = infiles
    index = P.snip(index, "/transcriptome.sfi")

    sample = P.snip(os.path.basename(outfile), "_quant.sf")
    outdir = "quantification/%(sample)s" % locals()
    
    m = PipelineMapping.Sailfish(strand=PARAMS["sailfish_strandedness"],
                                 orient=PARAMS["sailfish_orientation"],
                                 threads=PARAMS["sailfish_threads"])

    statement = m.build((infile,), outfile)

    P.run()

@follows(runSailfish)
@merge(runSailfish,
       "quantification/summary.tsv.gz")
def mergeSailfishResults(infiles,outfile):

    statement='''python %(scriptsdir)s/combine_tables.py
              --glob quantification/*/*quant.sf --columns 1 --take 7 
              --use-file-prefix -v 0| gzip > %(outfile)s'''
    P.run()


# AH: output a compressed file (.tsv.gz)
@active_if(BIAS_ANALYSIS)
@transform(PARAMS["sailfish_transcripts"],
           regex("(\S+)"),
           "transcripts_attributes.tsv.gz")
# take multifasta transcripts file and output file of attributes
def characteriseTranscripts(infile,outfile):

    statement = '''zcat %(infile)s | 
                python %(scriptsdir)s/fasta2table.py
                --split-fasta-identifier --section=dn -v 0
                | gzip > %(outfile)s'''
    P.run()


# where should this code be moved to?
@active_if(BIAS_ANALYSIS)
@follows(characteriseTranscripts)
@transform(characteriseTranscripts,
           regex("transcripts_attributes.tsv.gz"),
           add_inputs(mergeSailfishResults),
           ["quantification/binned_means_correlation.tsv",
            "quantification/binned_means_gradients.tsv"])
def summariseBias(infiles, outfiles):

    transcripts, expression = infiles
    out_correlation, out_gradient = outfiles

    atr = pandas.read_csv(transcripts, sep='\t')
    exp = pandas.read_csv(expression, sep='\t', compression="gzip")
    atr["length"] = numpy.log2(atr["length"])
 
    log_exp = numpy.log2(exp.ix[:,1:]+0.1)
    log_exp["id"] = exp[["Transcript"]]

    bias_factors = list(atr.columns[1:])
    samples = list(exp.columns[1:])

    merged = atr.merge(log_exp, left_index="id", right_index="id")

    def lin_reg_grad(x, y):
        slope, intercept, r, p, stderr = linregress(x, y)
        return slope

    def aggregate_by_factor(df, attribute, sample_names, bins, function):

        temp_dict = dict.fromkeys(sample_names, function)
        temp_dict[attribute] = function

        means_df = merged.groupby(pandas.qcut(df.ix[:, attribute], bins))
        means_df = means_df.agg(temp_dict).sort(axis=1)
        
        corr_matrix = means_df.corr(method='pearson')
        corr_matrix = corr_matrix[corr_matrix.index != attribute]

        factor_gradients = []
        for sample in samples:
            factor_gradients.append(lin_reg_grad(y=means_df[sample],
                                                 x=means_df[factor]))

        return means_df, corr_matrix, factor_gradients

    corr_matrices = {}
    gradient_lists = {}

    for factor in bias_factors:
        means_binned, corr_matrix, gradients = aggregate_by_factor(
            merged, factor, samples, PARAMS["bias_bin"], numpy.mean)
        outfile_means = "%s%s%s" % ("quantification/means_binned_",
                                    factor, ".tsv")
        means_binned.to_csv(outfile_means, sep="\t",
                            index=False, float_format='%.4f')

        corr_matrices[factor] = list(corr_matrix[factor])
        gradient_lists[factor] = gradients

    corr_matrix_df = DataFrame.from_dict(corr_matrices,
                                         orient='columns', dtype=None)
    corr_matrix_df["sample"] = sorted(samples)

    gradient_df = DataFrame.from_dict(gradient_lists,
                                      orient='columns', dtype=None)
    gradient_df["sample"] = sorted(samples)

    corr_matrix_df.to_csv(out_correlation, sep="\t",
                          index=False, float_format='%.6f')

    gradient_df.to_csv(out_gradient, sep="\t",
                       index=False, float_format='%.6f')

@follows(summariseBias)
@transform(summariseBias,
           regex("quantification/(\S+).tsv"),
           r"quantification/\1.load")
def loadBiasSummary(infiles, outfiles):
    for file in glob.glob("quantification/*.tsv"):
        P.load(file, P.snip(file, ".tsv")+".load")

# else:
#    @follows(loadFastqc)
#    def plotBias():
#        pass

#########################################################################


@follows(loadFastqc, loadFastqcSummary, loadBiasSummary)
def full():
    pass

@follows(loadBiasSummary)
def bias():
    pass


#########################################################################


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


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

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
