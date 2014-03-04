"""
================================
Optimizing cufflinks parameters
================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

The cufflinks optimization pipeline attempts to assess the quality of a variety of transcript
assemblies using various sets of parameters.


Overview
==========

Building transcripts is an important prerequisite for a number of projects that we undertake,
especially where the emphasis is on transcript structure or the identification of non-coding
transcripts (often lowly expressed).

In addition to the general issues  associated with transcript building, we are often working on
a variety of data from different library types and as such are not aware of the optimal set of
parameters for a transcript assembly from the outset of an analysis. This is especially important when 
we are faced with such things as high numbers of intronic reads from a non-polyA selected library
or low depth. Ideally we want to optimize the number of spliced reads that are incorporated into
transcript models, thus using the maximal amount of information that is present in any given datset.


Another issue in transcript assembly is the time that is required to run an assembly on an entire
transcriptome, creating a bottleneck in parameter optimization.

Given the above, the cufflinks optimization pipeline uses a reduced set of alignments (from a tophat run)
to optimize a variety of user specified parameters. It performs the following tasks

    * Reduction of input bam files into chr19 only bamfiles

    * Runs the rnaseq transcript building pipeline for all combinations
      of user specified parameters

    * Collects metrics in a single database csvdb from each run (metrics provided by the transcript
      building pipeline) and assesses how many reads contribute to transcripts (inc. spliced reads)

    * It also assesses the ratio of single vs. multi exon transfrags as a measure of overall transcriptome
      quality.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.

Input
-----

Mapped reads
++++++++++++

The principal input of this pipeline is a collection of reads mapped to a reference genome.
Mapped reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.bam

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. 

  
Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|cufflinks_          |>=1.3.0            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_cufflinks_optimization.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_cufflinks_optimization.tgz
   tar -xvzf pipeline_cufflinks_optimization.tgz
   cd pipeline_cufflinks_optimization.dir
   python <srcdir>/pipeline_cufflinks_optimization.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   cufflinks
      cufflinks_ - transcriptome analysis

.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html

Code
====

"""

# load modules
from ruffus import *

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

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
import operator

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Tophat as Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import pysam

import CGAT.Expression as Expression

import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineRnaseq as PipelineRnaseq
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Stats as Stats

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'annotations_dir': "",
        'paired_end': False})

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters(PARAMS["annotations_dir"],
                                      "pipeline_annotations.py")

###################################################################
###################################################################
###################################################################
# get options that are to be tested
cufflinks_options = {}
if "cufflinks_test_options" in PARAMS:
    options = P.asList(PARAMS["cufflinks_test_options"])
    for option in options:
        if option == "--pre-mrna-fraction" \
                or option == "--small-anchor-fraction" \
                or option == "--max-multiread-fraction":
            cufflinks_options[option] = [0, 0.5, 0.75, 1]
        elif option == "--min-isoform-fraction":
            cufflinks_options[option] = [0.05, 0.1, 0.5, 1]
        elif option == "--junc-alpha":
            cufflinks_options[option] = [0.001, 0.01, 0.1]
        elif option == "--min-frags-per-transfrag":
            cufflinks_options[option] = [1, 5, 10]
        elif option == "--overhang-tolerance":
            cufflinks_options[option] = [0, 2, 5, 8]
        elif option == "--overlap-radius":
            cufflinks_options[option] = [50, 100, 200]
        else:
            raise ValueError(
                "pipeline_cufflinks_optimization does not support parameter %s" % option)

if len(cufflinks_options) == 0:
    raise ValueError("no options to optimize specified")

#####################
# get input tracks
#####################

TRACKS = glob.glob("*.accepted.bam")

#############################################################
# ask the user to decide whether to continue with the number
# of assemblies
#############################################################

c = 0
for x in itertools.product(*cufflinks_options.values()):
    c += 1
#raw_input("pipeline_cufflinks_optimization will run %i trancript assemblies: hit enter to continue\n""" % c)

###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################


def updateFile(filename):
    '''
    create empty file for updating purposes
    '''

    outf = open(filename, "w")
    outf.write("file created for ruffus update")
    outf.close()

###################################################################
###################################################################
###################################################################


def options_generator(cufflinks_options):
    '''
    returns a generator for dealing with the cufflinks parameters
    for directory names and .ini file creation
    '''
    for option_values in itertools.product(*cufflinks_options.values()):
        yield " ".join(map(str, reduce(operator.add, zip(cufflinks_options.keys(), list(option_values)))))

###################################################################
###################################################################
###################################################################


def getDirectoryNames(options_generator):
    '''
    get the names for the directories
    '''
    for options in options_generator:
        yield "_".join(options.split(" ")).replace("--", "")

##############################
# produce list of directories
##############################

DIRECTORIES = [x for x in getDirectoryNames(
    options_generator(cufflinks_options))]

###################################################################
###################################################################
###################################################################


def getLogFileNames(options_generator):
    '''
    get filenames for log files
    '''
    for options in options_generator:
        yield "_".join(options.split(" ")).replace("--", "") + ".log"

###################################################################
###################################################################
###################################################################


@transform(TRACKS, suffix(".bam"), ".%s.bam" % PARAMS["chromosome"])
def reduceBamToChr19(infile, outfile):
    '''
    reduce the dataset for parameter testing.
    '''
    bam = pysam.Samfile(infile, "rb")
    outbam = pysam.Samfile(outfile, "wb", template=bam)
    for alignment in bam.fetch(PARAMS["chromosome"]):
        outbam.write(alignment)
    bam.close()
    outbam.close()

###################################################################
###################################################################
###################################################################


@transform(reduceBamToChr19, suffix(".bam"), ".bam.bai")
def indexBam(infile, outfile):
    '''
    index the reduced bam file
    '''
    pysam.index(infile)

###################################################################
###################################################################
###################################################################


@follows(indexBam, mkdir([directory for directory in getDirectoryNames(options_generator(cufflinks_options))]))
@split(indexBam, [logfile for logfile in getLogFileNames(options_generator(cufflinks_options))])
def createLogFiles(infile, outfiles):
    '''
    creates sentinel files
    '''
    for outfile in outfiles:
        updateFile(outfile)

###################################################################
###################################################################
###################################################################


@follows(createLogFiles)
@transform(indexBam, suffix(".bai"), add_inputs([logfile for logfile in getLogFileNames(options_generator(cufflinks_options))]), ".log")
def linkBamToWorkingDirs(infiles, outfile):
    '''
    symlink the bam file and index to the working directories
    for execution of the transcript building pipeline
    '''

    bamfile = P.snip(infiles[0], ".bai")
    indexfile = infiles[0]
    directories = [P.snip(logfile, ".log") for logfile in infiles[1]]

    for directory in directories:
        os.symlink(os.path.abspath(bamfile), os.path.join(directory, bamfile))
        os.symlink(
            os.path.abspath(indexfile), os.path.join(directory, indexfile))
    updateFile(outfile)

###################################################################
###################################################################
###################################################################


@transform(createLogFiles, regex(r"(\S+).log"), r"\1/pipeline.ini")
def createConfigFiles(infile, outfile):
    '''
    create all of the relevant .ini files in each working
    directory in order to execute the transcript building
    '''
    # test options for cufflinks
    cuff_opts = P.snip(infile, ".log").split("_")
    cuff_options = []
    for opt in cuff_opts:
        # not ideal to do my length but all I can think of at the moment
        if len(opt) > 6:
            cuff_options.append("--" + opt)
        else:
            cuff_options.append(opt)
    cuff_options = " ".join(cuff_options)

    options = PARAMS["cufflinks_options"]
    # directory for output config
    outdir = P.snip(infile, ".log")

    outf = open(os.path.join(outdir, "pipeline.ini"), "w")
    config_headers = []
    lines = []
    for line in open("pipeline.ini").readlines():
        lines.append(line)
        if line.find("[cufflinks]") != -1:
            outf.write("[cufflinks]\n\n# general cufflinks options\n\noptions=%s %s   \n" % (
                options, cuff_options))
        elif "[cufflinks]\n" in lines and "[cuffdiff\n]" not in lines:
            if line.find("options=") != -1:
                continue
            else:
                outf.write(line)
        else:
            outf.write(line)
    outf.close()

###################################################################
###################################################################
###################################################################


@follows(createLogFiles)
@files([os.path.join(PARAMS["pipeline_rnaseqtranscripts_dir"], x) for x in ["conf.py", "sphinxreport.ini"]], [os.path.join(directoryname, y) for directoryname, y in itertools.product(getDirectoryNames(options_generator(cufflinks_options)), ["conf.py", "sphinxreport.ini"])])
def linkToPipelineRnaseqTranscriptsConfigFiles(infiles, outfiles):
    '''
    produces links to the configuration files for report building
    '''
    for outfile in outfiles:
        if outfile.find("conf") != -1:
            os.symlink(infiles[0], outfile)
        elif outfile.find("sphinxreport") != -1:
            os.symlink(infiles[1], outfile)
        else:
            raise ValueError("cannot find outfile %s" % outfile)

###################################################################
###################################################################
###################################################################


@follows(linkBamToWorkingDirs, linkToPipelineRnaseqTranscriptsConfigFiles)
@transform(createConfigFiles, regex(r"(\S+)pipeline.ini"), r"\1report.log")
def executePipelineRnaseqTranscripts(infile, outfile):
    '''
    executes the transcripts building pipeline in each directory
    '''

    directory = os.path.dirname(infile)
    statement = '''cd ./%(directory)s; 
                   python %(scriptsdir)s/pipeline_rnaseqtranscripts.py -v5 -p10 make full'''
    P.run()

    statement = '''cd ./%(directory)s; 
                   python %(scriptsdir)s/pipeline_rnaseqtranscripts.py -v5 -p10 make build_report'''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(executePipelineRnaseqTranscripts)
@transform(os.path.join(DIRECTORIES[0], "reference.gtf.gz"), regex(r"(\S+).gtf.gz"), "reference.%s.gtf.gz" % PARAMS["chromosome"])
def filterReferenceGtfForChr19(infile, outfile):
    '''
    for reporting purposes get the chr19 filtered reference gtf file;
    gives an idea of the sort of numbers that are expected from
    the transcript assembly
    '''
    chro = PARAMS["chromosome_chr"]
    statement = '''zcat %(infile)s | grep %(chro)s | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(executePipelineRnaseqTranscripts)
@split([os.path.join(directory, "abinitio_lincrna.gtf.gz") for directory in getDirectoryNames(options_generator(cufflinks_options))], regex(r"(\S+).gtf.gz"), r"\1.*_exon.gtf.gz")
def splitMultiAndSingleExonLincRna(infile, outfiles):
    '''
    pulls out the multi-exonic and the single exonic lincRNA transcripts
    from the lincrna.gtf.gz
    '''

    inf = gzip.open(infile)
    multi = gzip.open(P.snip(infile, ".gtf.gz") + ".multi_exon.gtf.gz", "w")
    single = gzip.open(P.snip(infile, ".gtf.gz") + ".single_exon.gtf.gz", "w")

    for entry in GTF.transcript_iterator(GTF.iterator(inf)):
        if len(entry) > 1:
            for exon in entry:
                multi.write("\t".join(map(str, [exon.contig, exon.source, exon.feature, exon.start, exon.end, ".", exon.strand, "."]))
                            + "\t" + exon.attributes + "\n")
        elif len(entry) == 1:
            for exon in entry:
                single.write("\t".join(map(str, [exon.contig, exon.source, exon.feature, exon.start, exon.end, ".", exon.strand, "."]))
                             + "\t" + exon.attributes + "\n")

    for outfile in outfiles:
        outf = P.snip(outfile, ".gz")
        if not os.path.exists(outfile):
            statement = '''gzip %(outf)s'''
            P.run()

###################################################################
###################################################################
###################################################################


@follows(executePipelineRnaseqTranscripts)
@transform([os.path.join(directory, "abinitio_lincrna.gtf.gz") for directory in getDirectoryNames(options_generator(cufflinks_options))], regex(r"(\S+).gtf.gz"), r"\1.count")
def countMultiAndSingleExonLincRna(infile, outfile):
    '''
    outputs the transcript and gene counts for lincRNA transcripts
    '''
    outf = open(outfile, "w")
    outf.write(
        "no_multi_exon_transcripts\tno_single_exon_transcripts\tproportion_single\n")
    inf = GTF.iterator(IOTools.openFile(infile))
    c_multi = 0
    c_single = 0
    for gtfs in GTF.transcript_iterator(inf):
        if len(gtfs) > 1:
            c_multi += 1
        elif len(gtfs) == 1:
            c_single += 1
    outf.write(
        "\t".join(map(str, [c_multi, c_single, float(c_single) / (c_multi + c_single)])))

###################################################################
###################################################################
###################################################################


@follows(executePipelineRnaseqTranscripts)
@transform([os.path.join(directory, P.snip(bam, ".bam") + ".%s.bam" % PARAMS["chromosome"]) for directory, bam in itertools.product(getDirectoryNames(options_generator(cufflinks_options)), TRACKS)], suffix(".accepted.%s.bam" % PARAMS["chromosome"]), ".summary")
def summariseReadsContributingToTranscripts(infile, outfile):
    '''
    for each run of the pipeline and for each track, count the proportion
    of reads that contribute to the resulting transcript models
    '''
    gtf = P.snip(infile, ".bam") + ".gtf.gz"
    statement = '''python %(scriptsdir)s/bam2transcriptContribution.py -b %(infile)s -g %(gtf)s -o %(outfile)s --log=%(outfile)s.log'''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(executePipelineRnaseqTranscripts, splitMultiAndSingleExonLincRna)
@transform([os.path.join(directory, "abinitio_lincrna.multi_exon.gtf.gz") for directory in getDirectoryNames(options_generator(cufflinks_options))], suffix(".gtf.gz"), ".stats")
def summariseExonCountsAndLengthOfMultiExonicLincRNA(infile, outfile):
    '''
    summarizes some basic statistics on the length and number of exons 
    for each set of parameter values
    '''
    outf = open(outfile, "w")
    outf.write("transcript_id\tno_exons\ttranscriptlength\n")
    inf = GTF.iterator(IOTools.openFile(infile))
    for gtfs in GTF.transcript_iterator(inf):
        outf.write("\t".join((gtfs[0].transcript_id, str(len(gtfs)), str(
            sum([x.end - x.start for x in gtfs])))) + "\n")
    outf.close()

###################################################################
# report building
# the data that is generated in each subdirectory is collated into
# a central database in the working directory
###################################################################


@transform(countMultiAndSingleExonLincRna, regex(r"(\S+).count"), r"\1.load")
def loadCountSingleAndMultiExonLincRNA(infile, outfile):
    '''
    load the counts for the multi and single exon lincRNA
    '''
    tablename = P.toTable(outfile.replace("/", "_")) + ".count"
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log < %(infile)s > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(summariseReadsContributingToTranscripts, regex(r"(\S+).summary"), r"\1.summary.load")
def loadSummariseReadsContributingToTranscripts(infile, outfile):
    '''
    loads the summary of reads contributing to transcripts
    '''
    tablename = P.toTable(outfile.replace("/", "_"))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log < %(infile)s > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(summariseExonCountsAndLengthOfMultiExonicLincRNA, regex(r"(\S+).stats"), r"\1.load")
def loadNumberExonsLengthSummaryStats(infile, outfile):
    '''
    load the table of exon counts and transcript lengths
    '''
    tablename = P.toTable(outfile.replace("/", "_")) + "_stats"
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log < %(infile)s > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
# @follows(executePipelineRnaseqTranscripts)
# @transform([os.path.join(directory, P.snip(bam, ".bam") + ".chr19.class.tsv.gz") for directory, bam in itertools.product(getDirectoryNames(options_generator(cufflinks_options)), TRACKS)]
#            , regex(r"(\S+).class.tsv.gz"), r"\1.class.load" )
# def loadTranscriptClasses(infile, outfile):
#     '''
#     load the transcript class info for each assembly
#     '''
#     tablename = P.toTable(outfile.replace("/", "_"))
#     statement = '''zcat %(infile)s | python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log > %(outfile)s'''
#     P.run()

###################################################################
###################################################################
###################################################################


@files(loadNumberExonsLengthSummaryStats, "NumberExonsLength.tsv")
def outputAllNumberExonsLengthSummaryStats(infiles, outfile):
    '''
    outputs a flat file containing the exon lengths and number of 
    exons for multi-exon transcripts for each track
    '''
    dbh = connect()
    cc = dbh.cursor()
    tablenames = [P.snip(x, ".load").replace(".", "_").replace(
        "/", "_").replace("-", "_") + "_stats" for x in infiles]
    outf = open(outfile, "w")
    outf.write("\t".join(["track", "no_exons", "transcript_length"]) + "\n")
    for table in tablenames:
        length = cc.execute(
            "SELECT avg(transcriptlength) FROM %s" % table).fetchone()
        no_exons = cc.execute(
            "SELECT avg(no_exons) FROM %s" % table).fetchone()
        outf.write(
            table + "\t" + "\t".join(map(str, (no_exons[0], length[0]))) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@files(loadCountSingleAndMultiExonLincRNA, "MultiSingleLincRNACounts.tsv")
def outputAllMultiAndSingleExonLincRNACounts(infiles, outfile):
    '''
    output the counts for multi and single exon counts
    '''
    dbh = connect()
    cc = dbh.cursor()
    tablenames = [P.snip(x, ".load").replace(".", "_").replace(
        "/", "_").replace("-", "_") + "_count" for x in infiles]
    outf = open(outfile, "w")
    outf.write("multi_exon_count\tsingle_exon_count\tproportion_single\n")
    for table in tablenames:
        multi = cc.execute(
            "SELECT no_multi_exon_transcripts FROM %s" % table).fetchone()
        single = cc.execute(
            "SELECT no_single_exon_transcripts FROM %s" % table).fetchone()
        proportion = cc.execute(
            "SELECT proportion_single FROM %s" % table).fetchone()
        outf.write(
            "\t".join(map(str, (multi[0], single[0], proportion[0]))) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@files(loadSummariseReadsContributingToTranscripts, "ReadsToTranscripts.tsv")
def outputReadsContributingToTranscripts(infiles, outfile):
    '''
    output the proportion of reads contributing to transcripts
    '''
    dbh = connect()
    cc = dbh.cursor()
    tablenames = [P.snip(x, ".load").replace(".", "_").replace(
        "/", "_").replace("-", "_") for x in infiles]
    outf = open(outfile, "w")
    outf.write("proportion_reads_to_transcripts\n")
    for table in tablenames:
        proportion = cc.execute(
            "SELECT percent_spliced_alignments_in_transcripts FROM %s" % table).fetchone()
        outf.write(str(proportion[0]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@files([outputAllNumberExonsLengthSummaryStats, outputAllMultiAndSingleExonLincRNACounts, outputReadsContributingToTranscripts], "All_stats_combined.tsv")
def buildAllStats(infiles, outfile):
    '''
    paste stats together
    '''
    statement = '''paste %s > %s''' % (
        " ".join([infile for infile in infiles]), outfile)
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildAllStats, suffix(".tsv"), ".ranked")
def addRanksToAllStats(infile, outfile):
    '''
    rank each track by each variable
    '''
    # NB reverse the ranks as they go from lowest to highest
    R('''data <- read.table("%s", header = T, stringsAsFactors = F) ''' %
      infile)
    R('''data$rank_no_exons <- rev(rank(data$no_exons, ties.method = "min"))''')
    R('''data$rank_transcript_length <- rev(rank(data$transcript_length, ties.method = "min"))''')
    R('''data$rank_multi_exon_count <- rev(rank(data$multi_exon_count, ties.method = "min"))''')
    R('''data$rank_single_exon_count <- rank(data$single_exon_count, ties.method = "min")''')
    R('''data$rank_proportion_reads_to_transcripts <- rev(rank(data$proportion_reads_to_transcripts, ties.method = "min"))''')
    R('''data$rank_sum <- apply(cbind(data$rank_no_exons, data$rank_transcript_length, data$rank_multi_exon_count, data$rank_proportion_single, data$rank_proportion_reads_to_transcripts), 1, sum)''')
    R('''data$rank_all <- rank(data$rank_sum)''')
    R('''write.table(data[order(data$rank_sum),], file = "%s", row.names = F, sep = "\t")''' %
      outfile)

###################################################################
###################################################################
###################################################################


@transform(addRanksToAllStats, suffix(".ranked"), ".load")
def loadRankedAllStats(infile, outfile):
    '''
    load data with ranks in place
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################
## TARGETS ##
#############


@follows(loadNumberExonsLengthSummaryStats, loadSummariseReadsContributingToTranscripts, loadCountSingleAndMultiExonLincRNA)
def summary():
    pass

#------------------------------------------------------------------


@follows(loadRankedAllStats)
def rankData():
    pass

###################################################################
###################################################################
###################################################################


@follows(summary, rankData)
def full():
    pass

###################################################################
###################################################################
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)

###################################################################
###################################################################
###################################################################


@follows(build_report)
def publish():
    '''publish files.'''
    # publish web pages
    patterns = [(re.compile("/ifs/projects/proj012/report/html"),
                 "http://www.cgat.org/downloads/KW0ok5WWly/report"),
                ]

    P.publish_report(patterns=patterns)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
