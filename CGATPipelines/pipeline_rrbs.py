##########################################################################
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

"""====================
rrbs pipeline
====================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


######## UPDATE THIS SECTION ############
(See the readqc pipeline for details of the readqc functions
(target ``readqc``).)

The purpose of this pipeline is to pre-process reads (target ``full``).

Implemented tasks are:

   * :meth:`removeContaminants` - remove contaminants from read sets
   * :meth:`trim` - trim reads by a certain amount
   * :meth:`filter` - filter reads by quality score
   * :meth:`sample` - sample a certain proportion of reads

Individual tasks are enabled in the configuration file.
Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

Removing contaminants
---------------------

Use the task :meth:`removeContaminants` to remove contaminants from read
sets.

Contaminant sequences are listed in the file
:file:`contaminants.fasta`.  If not given, a file with standard
Illumina adapators will be created to remove adaptor contamination.

The task will create output files called :file:`nocontaminants-<infile>`.

The pipeline can then be re-run in order to add stats on the
contaminant-removed files.

.. note::

   Colour space filtering has not been implemented yet.

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while
``replicate`` denotes the :term:`replicate` within an
:term:`experiment`. The ``suffix`` determines the file type.  The
following suffixes/file types are possible:


fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+---------------+----------+------------------------------------------------+
|*Program*      |*Version* |*Purpose*                                       |
+---------------+----------+------------------------------------------------+
|fastqc         |>=0.9.0   |read quality control                            |
+---------------+----------+------------------------------------------------+
|sra-tools      |          |extracting reads from .sra files                |
+---------------+----------+------------------------------------------------+
|picard         |>=1.38    |bam/sam files. The .jar files need to be in your|
|               |          | CLASSPATH environment variable.                |
+---------------+----------+------------------------------------------------+

Pipeline output
===============

The major output is a set of HTML pages and plots reporting on the
quality of the sequence archive

Example
=======

Example data is available at ?
To run the example, simply unpack and untar::

TODO
====


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
import itertools
import glob
import cStringIO
import sqlite3
import CGAT.Experiment as E
import string
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator
import CGATPipelines.PipelineMapping as PipelineMapping
import CGAT.Pipeline as P
import CGAT.Fastq as Fastq
import CGAT.CSV2DB as CSV2DB
import CGAT.Fasta as Fa
import CGATPipelines.PipelineRrbs as RRBS
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import warnings

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

scriptsdir = PARAMS["general_scriptsdir"]


#########################################################################
#########################################################################
#########################################################################
# define input files
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

#########################################################################
#########################################################################
#########################################################################
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise


def connect():
    '''connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''
    dbh = sqlite3.connect(PARAMS["database"])

    return dbh


#########################################################################
# Read mapping
#########################################################################


SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(?P<suffix>fastq.1.gz|fastq.gz)")

#########################################################################
# summarise read 3'
#########################################################################


@follows(mkdir("sequence_characteristics.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"sequence_characteristics.dir/\1.\g<suffix>_start.tsv")
def summariseReadStart(infile, outfile):
    statement = '''zcat %(infile)s |
                paste - - - - | cut -f2 | cut -c1-3 | sort | uniq -c |
                sort -nk1 | awk -F' ' 'BEGIN{total=0; sum=0}
                {total+=$1; OFS"\\t";
                if($2=="CGG"||$2=="TGG"||$2=="CGA"||$2=="TGA")
                {sum+=$1; print $1, $2}}
                END {print total-sum,"others"}' > %(outfile)s ''' % locals()
    P.run()


@merge(summariseReadStart,
       "sequence_characteristics.dir/read_start_summary.tsv")
def combineReadStartSummaries(infiles, outfile):

    infile_list = " ".join(infiles)

    statement = '''echo -e
            "file\\treads\\tsequence\\tsample\\tcondition\\trep"
            > %(outfile)s;
            python %%(scriptsdir)s/combine_tables.py -v0  -a CAT -t
            %(infile_list)s|
            sed -e 's/sequence_characteristics.dir\///g'
            -e 's/.fastq.*start.tsv//g'|
            awk '{OFS="\\t"; split ($1, a, "-");
            print $1,$2,$3,a[1],a[2],a[3]}'
            >> %(outfile)s;''' % locals()

    print statement
    P.run()


@transform(combineReadStartSummaries,
           suffix(".tsv"),
           ".load")
def loadStartSummary(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()

#########################################################################
# Read mapping
#########################################################################


@follows(mkdir("bismark.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"bismark.dir/\1.\g<suffix>_bismark_bt2.bam")
def mapReadsWithBismark(infile, outfile):
    '''map reads with bismark'''

    # can this handle paired end?
    # is appears bismark uses twice as many CPUs as expeceted!
    job_options = "-l mem_free=%s " % PARAMS["bismark_memory"]
    job_threads = PARAMS["bismark_threads"] * 2
    outdir = "bismark.dir"
    bismark_options = PARAMS["bismark_options"]
    m = PipelineMapping.Bismark()
    statement = m.build((infile,), outfile)
    P.run()

#########################################################################
# Call Methylation
#########################################################################


@follows(mkdir("methylation.dir"))
@transform(mapReadsWithBismark,
           regex("bismark.dir/(\S+).bam"),
           r"methylation.dir/\1.bismark.cov")
def callMethylationStatus(infile, outfile):

    if infile.endswith((".fastq.gz_bismark_bt2.bam",
                        ".fastq.gz_bismark_bt.bam")):
        options = " --single-end "
    else:
        options = " --paired-end "

    if PARAMS["bismark_extraction_options"]:
        options += PARAMS["bismark_extraction_options"]

    # not using P.snip as two extensions in outfiles
    CG = ("methylation.dir/CpG_context_" +
          P.snip(os.path.basename(outfile), ".bismark.cov") + ".txt")
    CHG = re.sub("CpG", "CHG", CG)
    CHH = re.sub("CpG", "CHH", CG)

    outdir = "methylation.dir"
    index_dir = PARAMS["bismark_index_dir"]
    genome = PARAMS["bismark_genome"]

    statement = '''bismark_methylation_extractor %(options)s
                --comprehensive --output %(outdir)s --counts
                --cytosine_report --bedGraph
                --genome_folder %(index_dir)s/%(genome)s/
                %(infile)s; gzip -f %(CG)s; gzip -f  %(CHG)s; gzip -f %(CHH)s
                ''' % locals()
    P.run()


@follows(mkdir("plots.dir"))
@transform(callMethylationStatus,
           regex("methylation.dir/(\S+).bismark.cov"),
           r"plots.dir/\1.read_position.tsv")
def plotReadBias(infile, outfile):
    job_options = "-l mem_free=1G"

    m_bias_infile = P.snip(infile, ".bismark.cov") + ".M-bias.txt"

    print m_bias_infile

    RRBS.plotReadBias(m_bias_infile, outfile,
                      submit=True, jobOptions=job_options)

#########################################################################
# Sort and index bams
#########################################################################
# this is done here rather than during mapping as bismark requires read sorted
# bam files, not coordinate sorted
# can this be removed entirely? Nope, unfortunately not. Filtering is
# now done as a postprocessing step in Bismark though


@follows(callMethylationStatus)
@transform(mapReadsWithBismark,
           suffix(".bam"),
           ".sorted.bam")
def sortAndIndexBams(infile, outfile):
    sort_out = P.snip(outfile, ".bam")
    statement = '''samtools sort %(infile)s %(sort_out)s;
                   samtools index %(outfile)s;''' % locals()
    print statement
    P.run()


########################################################################
########################################################################
########################################################################
# RRBS CpGI coverage summary
########################################################################
# this section currently takes an input bed file with the locations of
# CpG Islands and computes coverage across these regions
# an alternative approach would be to load the CpG bed into database and query
# (see below)

@follows(mkdir("coverage.dir"))
@originate("coverage.dir/cpgIslands.bed")
def makeCpgIslandsBed(outfile):
    infile = PARAMS["methylation_summary_cpgislands"]
    out = open(outfile, "w")
    with open(infile, "r") as f:
        for line in f.readlines():
            # this assumes location of req. values
            contig, start, end = line.split()[1:4]
            if not contig == "chrom":
                out.write("%s\t%s\t%s\n" % (contig, start, end))
    out.close()


# @follows(makeCpgIslandsBed)
# @transform(sortFilterAndIndexBams,
#            regex("bismark.dir/(\S+).sorted.bam"),
#            add_inputs(makeCpgIslandsBed),
#            r"coverage.dir/\1.cpg.coverage")
# def overlapCpgIslands(infiles, outfile):
#     infile, bed = infiles
#
#     job_options = "-l mem_free=1G"
#     statement = '''samtools view -b %(infile)s|
#                    coverageBed -hist -d -abam stdin -b %(bed)s
#                    > %(outfile)s''' % locals()
#     print statement
#     P.run()

########################################################################
# RRBS alternative CpGI coverage summary
########################################################################


@subdivide(makeCpgIslandsBed,
           regex("coverage.dir/(\S+).bed"),
           [r"coverage.dir/\1_1based.tsv",
            r"coverage.dir/\1_1based.load"])
def make1basedCpgIslands(infile, outfiles):

    outfile, loadfile = outfiles

    out = open(outfile, "w")
    out.write("%s\t%s\t%s\n" % ("contig", "position", "cpgi"))

    with open(infile, "r") as f:
        lines = f.readlines()
        for line in lines:
            contig, start, stop = line.split()
            for position in [x for x in range(int(start), int(stop)+2)]:
                out.write("%s\t%s\t%s\n" % (contig, position, "CpGIsland"))
    out.close()

    dbh = connect()
    tablename = P.toTable(loadfile)
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(outfile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(loadfile)s''' % locals()
    P.run()


@transform(make1basedCpgIslands,
           regex("coverage.dir/(\S+).load"),
           r"coverage.dir/cpg_table.load")
# not currently implemented
def joinAllCpGsAndCpGIslands(infile, outfile):

    tablename = P.toTable(outfile)
    statement = '''sqlite3 csvdb
                  "CREATE %(tablename)s FROM
                   SELECT *"'''


#########################################################################
# load bismark coverage and join tables in sql
#########################################################################


@transform(callMethylationStatus,
           suffix(".bismark.cov"),
           r".bismark.subset10.cov")
# this is just for testing BiSeq etc, delete in future
def subsetCoverage(infile, outfile):
    statement = '''awk '($5+$6)>=10' %(infile)s > %(outfile)s''' % locals()
    P.run()
##########################################################################


# note the output of this function is 1-based coordinates
# need to remember this!
@originate("methylation.dir/cpg-locations-1.cov")
def findCpGs(outfile):
    genome_infile = PARAMS["methylation_summary_genome_fasta"]
    job_options = "-l mem_free=2G"

    RRBS.fasta2CpG(genome_infile, outfile, submit=True, jobOptions=job_options)


@follows(findCpGs)
@merge([callMethylationStatus,
        findCpGs],
       "methylation.dir/cpgs_meth.tsv")
def mergeCoverage(infiles, outfile):
    cpgs_infile = infiles[-1]
    coverage_infiles = infiles[:-1]

    job_options = "-l mem_free=48G"
    job_threads = 2

    RRBS.mergeAndDrop(cpgs_infile, coverage_infiles, outfile,
                      submit=True, jobOptions=job_options)


@transform(mergeCoverage,
           suffix("_meth.tsv"),
           add_inputs(make1basedCpgIslands),
           "_meth_cpgi.tsv")
def addCpGIs(infiles, outfile):
    infile, CpGI_load, CpGI = infiles
    # very memory intensive!
    job_options = "-l mem_free=20G"
    job_threads = 2

    RRBS.pandasMerge(infile, CpGI, outfile, merge_type="left",
                     left=['contig', 'position'],
                     right=['contig', 'position'],
                     submit=True, jobOptions=job_options)


@transform(addCpGIs,
           suffix(".tsv"),
           ".load")
def loadMergeCoverage(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    scriptsdir = PARAMS["general_scriptsdir"]
    job_options = "-l mem_free=23G"
    job_threads = 2

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


#########################################################################
#########################################################################
#########################################################################
# this is an entirely different approach to generating the coverage table
# using sqlite3. Currently not correctly implemented (requires multiple
# inner joins and union in sqlite3 command)
# @jobs_limit(1)
# @follows(makeCoverageHeader)
# @transform(callMethylationStatus,
#           suffix(".cov"),
#           add_inputs(makeCoverageHeader),
#           ".cov.load")
# def loadAndIndexBismarkCoverage(infiles, outfile):
#    infile, header = infiles
#    dbh = connect()
#    tablename = P.toTable(outfile)
#    scriptsdir = PARAMS["general_scriptsdir"]
#    index = "idx" + re.sub('_fastq_.+', "", tablename)
#
#
#    statement = '''cat %(header)s %(infile)s |
#                python %(scriptsdir)s/csv2db.py
#                --table %(tablename)s --retry --ignore-empty
#                 > %(outfile)s;
#                sqlite3 csvdb "CREATE INDEX %(index)s ON
#                %(tablename)s(contig, start)" ''' % locals()
#    P.run()
#
#
# @merge(loadAndIndexBismarkCoverage,
#       "methylation.dir/csvdb_join_bismark_coverage.sentinel")
# def joinBismarkCoverage(infiles, outfile):
#    print infiles
#    print outfile
#
#    new_columns = ["a.contig", "a.start"]
#    tables = [" FROM "]
#    letters = string.ascii_lowercase
#    n = 0
#    for cov_file in infiles:
#        track = re.sub('[-\.]', "_", os.path.basename(cov_file))
#        column_track = re.sub('\_fastq\_.+', "", track)
#        table_track = P.snip(track, "_load")
#        letter = letters[n]
#
#        new_columns.append(letter +
#                           ".count_meth AS %(column_track)s_meth" % locals())
#        new_columns.append(letter +
#                           ".count_unmeth AS %(column_track)s_unmeth" %
#                           locals())
#        if n == 0:
#            tables.append("%(table_track)s AS %(letter)s" % locals())
#        else:
#            tables.append('''FULL OUTER JOIN
#                          %(table_track)s AS %(letter)s
#                          ON a.contig = %(letter)s.contig
#                          AND a.start = %(letter)s.start''' % locals())
#        n += 1
#
#    new_columns = ", ".join(new_columns)
#    tables = " ".join(tables)
#
#    statement = '''sqlite3 csvdb "CREATE TABLE merged_bismark_cov AS SELECT'''
#    statement += new_columns
#    statement += tables
#    statement += ''' " '''
#
#    print statement
#    P.run()
#########################################################################
#########################################################################
#########################################################################


#########################################################################
# Summarise methylation
#########################################################################

# all these functions should be rewitten to take the output
# from mergeCoverage as input to clean up pipeline


@transform(callMethylationStatus,
           regex("methylation.dir/(\S+).bismark.cov"),
           r"methylation.dir/\1.coverage.tsv")
def summariseCoverage(infile, outfile):

    CG = ("methylation.dir/CpG_context_" +
          P.snip(os.path.basename(infile), ".bismark.cov") + ".txt.gz")

    statement = '''zcat %(CG)s | awk -F'\t' '{print $3+$4;}'| sort -n |
                   uniq -c | cut -f1 |awk '{print $1}' | sort -n | uniq -c|
                   awk -F' ' '{OFS="\\t";} {print $1,$2;}'
                   > %(outfile)s''' % locals()
    print statement
    P.run()


@follows(summariseCoverage)
@merge(summariseCoverage,
       "methylation.dir/coverage.tsv")
def concatenateCoverage(infiles, outfile):
    '''concatenate coverage output'''

    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''echo -e "file\\tfreq\\tcov\\tsample\\tcondition\\trep"
                > %(outfile)s;
                python %(scriptsdir)s/combine_tables.py -v0
                --glob="methylation.dir/*.coverage.tsv" -a CAT -t|
                sed -e 's/methylation.dir\///g'
                -e 's/.fastq.*coverage.tsv//g'|
                awk '{OFS="\\t"; split ($1, a, "-");
                print $1,$2,$3,a[1],a[2],a[3]}'
                >> %(outfile)s;''' % locals()
    P.run()


@transform(concatenateCoverage,
           suffix(".tsv"),
           ".load")
def loadCoverage(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


@follows(summariseCoverage)
@merge(callMethylationStatus,
       "methylation.dir/coverage_overlap.tsv")
def summariseCpGOverlap(infiles, outfile):
    infile_list = []
    for x in infiles:
        if x.endswith(".cov"):
            infile_list.append(x)
    print(len(infile_list))
    infile_list = " ".join(infile_list)
    coverage_range = [1, 2, 5, 10, 20, 30]
    coverage_range = " ".join(map(str, coverage_range))

    statement = '''echo -e "CpGs\\toverlaps\\tthreshold" > %(outfile)s;
                for x in %(coverage_range)s;
                do cat %(infile_list)s | awk -v threshold=$x -F'\\t'
                '{OFS="\\t"; if ($5+$6>threshold) print $1,$2}'| sort| uniq -c|
                awk  -F' '  '{OFS="\\t"; print $1}'| sort| uniq -c|
                awk -v threshold=$x -F' ' '{OFS="\\t";print $1,$2,threshold}'
                >> %(outfile)s; done''' % locals()
    print infile_list
    P.run()


@transform(summariseCpGOverlap,
           suffix(".tsv"),
           ".load")
def loadCpGOverlap(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


@transform(summariseCoverage,
           regex("methylation.dir/(\S+).coverage.tsv"),
           r"methylation.dir/\1.reads_by_threshold.tsv")
def summariseRemainingReadsbyThreshold(infile, outfile):

    statement = '''sum_reads=$(awk -F'\\t' 'BEGIN {total=0} {total+=($1*$2)}
                END {print total}' %(infile)s);
                echo -e "1\\t$sum_reads\\t1" >> %(outfile)s;
                awk -v old_total=$sum_reads 'BEGIN {total=old_total}
                {total -= ($1*$2)}
                {OFS="\\t"; print $2+1,total,total/old_total}'
                %(infile)s >> %(outfile)s''' % locals()
    P.run()


@merge(summariseRemainingReadsbyThreshold,
       "methylation.dir/reads_remaining_by_threshold.tsv")
def concatenateRemainingReads(infiles, outfile):
    '''concatenate coverage output'''

    # change statement to use infiles rather than glob
    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''echo -e
                "file\\tthreshold\\treads\\tpercentage\\tsample\\tcondition\\trep"
                > %(outfile)s;
                python %(scriptsdir)s/combine_tables.py -v0
                --glob="methylation.dir/*.reads_by_threshold.tsv" -a CAT -t|
                sed -e 's/methylation.dir\///g'
                -e 's/.fastq.*.reads_by_threshold.tsv//g'|
                awk '{OFS="\\t"; split ($1, a, "-");
                print $1,$2,$3,$4,a[1],a[2],a[3]}'
                >> %(outfile)s;''' % locals()
    P.run()


@transform(concatenateRemainingReads,
           suffix(".tsv"),
           ".load")
def loadRemainingReads(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


########################################################################
# RRBS Gene meta-profile
########################################################################
# this section makes plots of coverage across TSS
# should other meta-profiles be included?


@follows(makeCpgIslandsBed)
@transform(sortAndIndexBams,
           regex("bismark.dir/(\S+).sorted.bam"),
           add_inputs(PARAMS["annotation_genes"]),
           r"coverage.dir/\1.profile.png")
def makeGeneProfiles(infiles, outfile):
    outname = P.snip(os.path.basename(outfile), ".png")
    infile, genes_gtf = infiles
    scriptsdir = PARAMS["general_scriptsdir"]

    # ensures only large RAM nodes used
    job_options = "-l mem_free=24G"

    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                  --bamfile=%(infile)s --gtffile=%(genes_gtf)s
                  --method=tssprofile --reporter=gene
                  -P %(outname)s
                  --normalization=total-sum
                  --normalize-profile=area;
                  checkpoint;
                  for file in %(outname)s*;
                  do mv $file coverage.dir/.; done'''
    print statement
    P.run()

########################################################################
########################################################################
########################################################################
# Produce summary plots of methylation between samples, CpGI vs. non-CpGI etc


@transform(addCpGIs,
           regex("methylation.dir/(\S+)_meth_cpgi.tsv"),
           r"plots.dir/\1_covered_meth_cpgi.tsv")
def subsetCpGsToCovered(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.subsetToCovered(infile, outfile,
                         submit=True, jobOptions=job_options)


@transform(subsetCpGsToCovered,
           regex("plots.dir/(\S+)_covered_meth_cpgi.tsv"),
           r"plots.dir/\1_covered_with_means.tsv")
def addTreatmentMeans(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.addTreatmentMean(infile, outfile,
                          submit=True, jobOptions=job_options)


@transform(addTreatmentMeans,
           regex("plots.dir/(\S+)_covered_with_means.tsv"),
           r"plots.dir/\1_summary_plots.log")
def makeSummaryPlots(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.summaryPlots(infile, outfile,
                      submit=True, jobOptions=job_options)

########################################################################
########################################################################
########################################################################
# DMR analysis
# modularise BiSeq. For now, code is contained here
# add functionality to deal with no replicate experimental designs


@merge(subsetCoverage,
       "clusters.tsv")
def runBiSeq(infiles, outfile):
    # use this when following callMethylationStatus(i.e full run)
    # cov_infiles = filter(lambda x: 'Liver' in x, infiles)
    cov_infiles = infiles
    # implement properly
    if len(cov_infiles) >= 4:
        pass
    else:
        warnings.warn("nope, not gonna work")
        out = open(outfile, "w")
        out.write("You don't have replicates! This wont work")
        out.close()
        return 0
    print ("how did I get here given that I only have %s samples"
           % len(cov_infiles))

    basenames = [os.path.basename(x) for x in cov_infiles]
    samples = [re.sub(r".fastq\S+", "", x) for x in basenames]
    groups = [x.split("-")[1] for x in samples]
    c = '","'.join(cov_infiles)
    s = '","'.join(samples)
    g = '","'.join(groups)
    print ("basenames: %(basenames)s\ngroups: %(groups)s\nsamples: %(samples)s"
           % locals())
    print "c: %(c)s\ns: %(s)s\ng: %(g)s " % locals()
    # load BiSeq package
    r.library("BiSeq")
    r.library("ggplot2")
    grdevices = importr('grDevices')

    rcode = ('raw = readBismark(files = c("%(c)s"),DataFrame(group = c("%(g)s"),\
    row.names = c("%(s)s")))' % locals())

    r(rcode)

    r('colData(raw)$group <- factor(colData(raw)$group)')
    r('rrbs_rel <- rawToRel(raw)')
    # parameterise this step
    print "working fine to here"
    r('print(head(methReads(raw)))')
    r('print(head(totalReads(raw)))')
    r('print(rowData(raw))')
    r('print(colData(raw))')
    r('print(colData(raw)$group)')

    r('rrbs.clust.unlim <- BiSeq::clusterSites(\
    object = raw,groups = colData(raw)$group,perc.samples = 1,\
    min.sites = 5,max.dist = 50,mc.cores=4)')
    print "made it here"
    r('ind.cov <- totalReads(rrbs.clust.unlim) > 0')
    r('quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)')
    r('rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)')
    # emperically derive value for h (this is the bandwidth)
    r('predictedMeth <- predictMeth(object = rrbs.clust.lim,h=25,mc.cores=4)')
    r('temp_df=data.frame(y=methLevel(predictedMeth)\
    [order(rowData(predictedMeth)),1])')
    r('temp_df["x"]=\
    (methReads(rrbs.clust.lim)/totalReads(rrbs.clust.lim))[,1]')
    r('temp_df["cov"] = totalReads(rrbs.clust.lim)[,1]')
    r('p=ggplot(temp_df,aes(x, y, size=cov))+geom_point(aes(colour=cov))')

    grdevices.png(file="/ifs/projects/proj034/mapping_bismark_0.12.5/test.png")
    r('print(p)')
    grdevices.dev_off()

    r('betaResults <- betaRegression(formula = ~group, link = "probit",\
    object = predictedMeth, type = "BR",mc.cores=8)')
    r('print(head(betaResults))')
    r('predictedMethNull <- predictedMeth')
    r('print((predictedMethNull))')
    r('print(head(methLevel(predictedMethNull)))')
    r('colData(predictedMethNull)$group.null <- rep(c(1,2), 3)')
    r('print((predictedMethNull))')
    r('betaResultsNull <- betaRegression(formula = ~group.null,\
    link = "probit",object = predictedMethNull, type = "BR",mc.cores=4)')
    r('vario <- makeVariogram(betaResultsNull)')
    r('length_vario = length(vario$variogram[[1]][,2])')
    r('estimated_sill =\
    median(vario$variogram[[1]][-(1:(length_vario-50)),2])')
    r('vario.sm <- smoothVariogram(vario, sill = estimated_sill)')

    grdevices.png(
        file="/ifs/projects/proj034/mapping_bismark_0.12.5/test_vario.png")
    r('plot(vario$variogram[[1]],ylim=c(0,5))')
    r('lines(vario.sm$variogram[,c("h", "v.sm")],col = "red", lwd = 1.5)')
    grdevices.dev_off()

    r('print(names(vario.sm))')
    r('print(head(vario.sm$variogram))')
    r('print(head(vario.sm$pValsList[[1]]))')
    r('vario.aux <- makeVariogram(betaResults, make.variogram=FALSE)')
    r('vario.sm$pValsList <- vario.aux$pValsList')
    r('print(head(vario.sm$pValsList[[1]]))')
    r('locCor <- estLocCor(vario.sm)')
    r('clusters.rej <- testClusters(locCor,FDR.cluster = 0.5)')
    r('print(clusters.rej$clusters.reject)')
    r('clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = 0.5)')
    r('print(clusters.trimmed)')
    r('DMRs <- findDMRs(clusters.trimmed,max.dist = 100,diff.dir = TRUE)')
    r('print(DMRs)')
    # concatenate pvalue lists into one table to write out
    r('concat_vario = do.call(rbind, vario.sm$pValsList)')
    r('row_p = rowData(predictedMeth)')
    r('df_ranges<- data.frame(seqnames=seqnames(row_p), starts=start(row_p),\
    ends=end(row_p))')
    r('ranges_pred_meth_df=cbind(df_ranges,methLevel(predictedMeth))')
    r('write.table(concat_vario, "clusters_vario_table_test.tsv",\
    sep="\t",quote = F,row.names = F)')
    r('write.table(ranges_pred_meth_df, "clusters.tsv", sep="\t",\
    quote = F,row.names = F)')

########################################################################
#########################################################################
# spike-in analysis
########################################################################


@follows(mkdir("power.dir"))
@merge(callMethylationStatus,
       "power.dir/power.out")
def generateClusterSpikeIns(infiles, outfile):

    # include some way of restricting power analysis to subset of samples
    job_options = "-l mem_free=23G"

    RRBS.spikeInClusters(infiles, outfile,
                         submit=True, jobOptions=job_options)


@follows(generateClusterSpikeIns)
@transform(generateClusterSpikeIns,
           suffix(".out"),
           ".analysis.out")
def clusterSpikeInsPowerAnalysis(infiles, outfile):

    job_options = "-l mem_free=23G"

    RRBS.spikeInClustersAnalysis(infiles, outfile,
                                 submit=True, jobOptions=job_options)


@transform(clusterSpikeInsPowerAnalysis,
           suffix(".analysis.out"),
           ".plot.out")
def clusterSpikeInsPowerPlot(infiles, outfile):

    job_options = "-l mem_free=23G"

    RRBS.spikeInClustersPlot(infiles, outfile, groups=["Saline", "Dex"],
                             submit=True, jobOptions=job_options)

########################################################################
#########################################################################


@follows(generateClusterSpikeIns)
def test():
    pass


@follows(runBiSeq)
def biseq():
    pass


@follows(loadCoverage,
         loadCpGOverlap,
         loadRemainingReads,
         loadStartSummary,
         findCpGs,
         subsetCoverage,
         sortAndIndexBams,
         makeCpgIslandsBed,
         makeGeneProfiles,
         make1basedCpgIslands,
         loadMergeCoverage,
         makeSummaryPlots,
         mergeCoverage,
         plotReadBias)
def full():
    pass


@follows(callMethylationStatus)
def callMeth():
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
