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

.sra
   Single-end reads
   (Paired-end reads not yet checked)

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
document design file
document power analysis options
parameterise power analysis?
clean up bash code for coverage analysis, i.e replace with python

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
import CGATPipelines.PipelineTracks as PipelineTracks
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
    r"(\S+)-(\S+)-(\S+).(?P<suffix>fastq.1.gz|fastq.gz|sra)")

Sample = PipelineTracks.AutoSample
Sample.attributes = ('tissue', 'condition', 'replicate')
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    [y for x in SEQUENCESUFFIXES for y in glob.glob(x)],
    "(\S+).(fastq.1.gz|fastq.gz|sra)")

EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("tissue", "condition"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
REPLICATES = PipelineTracks.Aggregate(TRACKS, labels=("replicate", ))


print [x for x in EXPERIMENTS]
print [x for x in CONDITIONS]
print [x for x in REPLICATES]


#########################################################################
# summarise read 3'
#########################################################################


@follows(mkdir("sequence_characteristics.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"sequence_characteristics.dir/\1-\2-\3.\g<suffix>_start.tsv")
def summariseReadStart(infile, outfile):
    # this only works for fastq files. Fails with .sra files
    # this function and the next section should be replaced with a call to
    # fastq-dump if the file ends with .sra and then use the functions of
    # the fastq module to count the first bases in the fastq records.
    # for now, create empty outfile
    if infile.endswith(".sra"):
        P.touch(outfile)
    else:
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

    P.run()


@transform(combineReadStartSummaries,
           suffix(".tsv"),
           ".load")
def loadStartSummary(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()

#########################################################################
# Read mapping
#########################################################################


@follows(mkdir("bismark.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"bismark.dir/\1-\2-\3_bismark_bt2.bam")
def mapReadsWithBismark(infile, outfile):
    '''map reads with bismark'''

    # can this handle paired end?
    # it appears bismark uses twice as many CPUs as expeceted!
    job_options = "-l mem_free=%s " % PARAMS["bismark_memory"]
    job_threads = (PARAMS["bismark_threads"] * 2) + 1
    outdir = "bismark.dir"
    bismark_options = PARAMS["bismark_options"]
    m = PipelineMapping.Bismark()
    statement = m.build((infile,), outfile)
    # print statement
    P.run()

#########################################################################
# Call Methylation
#########################################################################


@follows(mkdir("methylation.dir"))
@transform(mapReadsWithBismark,
           regex("bismark.dir/(\S+).bam"),
           r"methylation.dir/\1.bismark.cov")
def callMethylationStatus(infile, outfile):

    if infile.endswith(("bismark_bt2.bam",
                        "bismark_bt.bam")):
        options = " --single-end "
    else:
        options = " --paired-end "

    if PARAMS["bismark_extraction_options"]:
        options += PARAMS["bismark_extraction_options"]

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
                      submit=True, job_options=job_options)

#########################################################################
# Sort and index bams
#########################################################################
# this is done here rather than during mapping as bismark requires read sorted
# bam files, not coordinate sorted
# the infiles from bismark
# both the original Bismark bams and the sorted bams are used downstream
# so currently both bams are kept (any way around this?)


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


########################################################################
# RRBS alternative CpGI coverage summary
########################################################################

@subdivide(makeCpgIslandsBed,
           regex("coverage.dir/(\S+).bed"),
           r"coverage.dir/\1_1based.tsv")
# r"coverage.dir/\1_1based.load"
def make1basedCpgIslands(infile, outfile):

    # outfile, loadfile = outfiles

    out = open(outfile, "w")
    out.write("%s\t%s\t%s\n" % ("contig", "position", "cpgi"))

    with open(infile, "r") as f:
        lines = f.readlines()
        for line in lines:
            contig, start, stop = line.split()
            for position in [x for x in range(int(start), int(stop)+2)]:
                out.write("%s\t%s\t%s\n" % (contig, position, "CpGIsland"))
    out.close()
    # this file takes hours(!) to load and it's not being used from csvdb
    # dbh = connect()
    # tablename = P.toTable(loadfile)

    # statement = '''cat %(outfile)s |
    #            python %%(scriptsdir)s/csv2db.py
    #            --table %(tablename)s --retry --ignore-empty
    #             > %(loadfile)s''' % locals()
    # P.run()


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

# change to subset10 after testing biseq
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

    RRBS.fasta2CpG(genome_infile, outfile,
                   submit=True, job_options=job_options)


@follows(findCpGs)
@merge([callMethylationStatus,
        findCpGs],
       "methylation.dir/cpgs_meth.tsv")
def mergeCoverage(infiles, outfile):
    cpgs_infile = infiles[-1]
    coverage_infiles = infiles[:-1]
    # very memory intensive! - find out why and re-code
    job_options = "-l mem_free=48G"
    job_threads = 2

    RRBS.mergeAndDrop(cpgs_infile, coverage_infiles, outfile,
                      submit=True, job_options=job_options)


@transform(mergeCoverage,
           suffix("_meth.tsv"),
           add_inputs(make1basedCpgIslands),
           "_meth_cpgi.tsv")
def addCpGIs(infiles, outfile):
    infile, CpGI = infiles
    # memory intensive!
    job_options = "-l mem_free=20G"
    job_threads = 2

    RRBS.pandasMerge(infile, CpGI, outfile, merge_type="left",
                     left=['contig', 'position'],
                     right=['contig', 'position'],
                     submit=True, job_options=job_options)


@transform(addCpGIs,
           suffix(".tsv"),
           ".load")
def loadMergeCoverage(infile, outfile):

    dbh = connect()
    tablename = P.toTable(outfile)
    job_options = "-l mem_free=23G"
    job_threads = 2

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


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

    statement = '''echo -e "file\\tfreq\\tcov\\tsample\\tcondition\\trep"
                > %(outfile)s;
                python %%(scriptsdir)s/combine_tables.py -v0
                --glob="methylation.dir/*.coverage.tsv" -a CAT -t|
                sed -e 's/methylation.dir\///g'
                -e 's/_bismark.*.coverage.tsv//g'|
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

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
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

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
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

    statement = '''echo -e
                "file\\tthreshold\\treads\\tpercentage\\tsample\\tcondition\\trep"
                > %(outfile)s;
                python %%(scriptsdir)s/combine_tables.py -v0
                --glob="methylation.dir/*.reads_by_threshold.tsv" -a CAT -t|
                sed -e 's/methylation.dir\///g'
                -e 's/_bismark.*.reads_by_threshold.tsv//g'|
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

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
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

    # ensures only large RAM nodes used
    job_options = "-l mem_free=24G"

    statement = '''python %%(scriptsdir)s/bam2geneprofile.py
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
           r"methylation.dir/\1_covered_meth_cpgi.tsv")
def subsetCpGsToCovered(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.subsetToCovered(infile, outfile,
                         submit=True, jobOptions=job_options)


@transform(subsetCpGsToCovered,
           regex("methylation.dir/(\S+)_covered_meth_cpgi.tsv"),
           r"methylation.dir/\1_covered_with_means.tsv")
def addTreatmentMeans(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.addTreatmentMean(infile, outfile,
                          submit=True, job_options=job_options)


@transform(addTreatmentMeans,
           regex("methylation.dir/(\S+)_covered_with_means.tsv"),
           r"plots.dir/\1_summary_plots")
def makeSummaryPlots(infile, outfile):

    job_options = "-l mem_free=48G"

    RRBS.summaryPlots(infile, outfile,
                      submit=True, job_options=job_options)
    P.touch(outfile)

########################################################################
########################################################################
########################################################################
# DMR analysis
# modularise BiSeq. For now, code is contained here
# add functionality to deal with no replicate experimental designs


@follows(mkdir("biseq.dir"))
@merge(subsetCoverage,
       "biseq.dir/clusters.tsv")
def runBiSeq(infiles, outfile):
    job_options = "-l mem_free=10G -pe dedicated 10"
    RRBS.runBiSeq(infiles, outfile, submit=True, job_options=job_options)

########################################################################
#########################################################################
# spike-in analysis
########################################################################


@follows(mkdir("power.dir"))
@transform(subsetCpGsToCovered,
           regex("methylation.dir/(\S+)_covered_meth_cpgi.tsv"),
           "power.dir/spike_ins.tsv")
def generateClusterSpikeIns(infile, outfile):
    # parametrise binning in pipeline.ini
    job_options = "-l mem_free=4G"

    statement = '''cat %(infile)s |
    python %%(scriptsdir)s/data2spike.py --method=spike
    --design-tsv-file=design.tsv --difference-method=relative
    --spike-shuffle-column-suffix=-perc
    --spike-keep-column-suffix=-meth,-unmeth
    --spike-minimum=100 --spike-maximum=100
    --spike-output-method=seperate
    --spike-cluster-maximum-distance=150
    --spike-cluster-minimum-size=10 --spike-iterations=50
    --spike-type=cluster --spike-change-bin-min=-100
    --spike-change-bin-max=100 --spike-change-bin-width=10
    --spike-initial-bin-min=0 --spike-initial-bin-max=100
    --spike-initial-bin-width=100 --spike-subcluster-min-size=1
    --spike-subcluster-max-size=9 --spike-subcluster-bin-width=1
    > %(outfile)s_tmp; mv %(outfile)s_tmp %(outfile)s''' % locals()
    P.run()


@split(generateClusterSpikeIns,
       "power.dir/spike_ins_subframe.tsv_*")
def splitSpikeClustersDataframe(infile, outfiles):
    outprefix = "power.dir/spike_ins_subframe.tsv"
    df = pd.read_csv(infile, comment='#', sep="\t")
    df['cluster_id'] = df['contig'].apply(
        lambda x: x.split("_")[-1]).astype(float)
    split_at = 1000
    for sub in range(0, int(max(df['cluster_id'])), split_at):
        temp_df = df[df['cluster_id'] >= sub]
        temp_df = temp_df[temp_df['cluster_id'] < sub+split_at]
        temp_df.to_csv("_".join((outprefix, str(sub))),
                       header=True, index=False, sep="\t")


@transform(splitSpikeClustersDataframe,
           regex("power.dir/spike_ins_subframe.tsv_(\d+)$"),
           add_inputs("design.tsv"),
           r"power.dir/spike_ins_subframe.tsv_\1_M3D.stats")
def runM3DSpikeClusters(infiles, outfile):
    job_options = "-l mem_free=4G -pe dedicated 1"
    infile, design = infiles
    RRBS.calculateM3DStat(infile, outfile, design,
                          submit=True, job_options=job_options)


@merge([runM3DSpikeClusters, "design.tsv"],
       "power.dir/spike_in_M3D_stat_merged_between.tsv")
def calculateM3DSpikeClustersPvalue(infiles, outfile):
    job_options = "-l mem_free=4G -pe dedicated 1"
    design = infiles[-1]
    infiles = infiles[:-1]
    RRBS.calculateM3DSpikepvalue(infiles, outfile, design,
                                 submit=True, job_options=job_options)
    P.touch(outfile)


###########################################################################

@transform(splitSpikeClustersDataframe,
           regex("power.dir/spike_ins_subframe.tsv_(\d+)"),
           r"power.dir/spike_ins_subframe.tsv_\1_BiSeq.stats")
def runBiSeqSpikeClusters(infile, outfile):
    pass
    job_options = "-l mem_free=10G -pe dedicated 6"
    # RRBS.calculateBiSeqStat(infile, outfile,
    # submit=True, job_options=job_options)


@follows(generateClusterSpikeIns)
@transform(generateClusterSpikeIns,
           suffix(".tsv"),
           ".analysis.out")
def clusterSpikeInsPowerAnalysis(infiles, outfile):

    job_options = "-l mem_free=23G"

    RRBS.spikeInClustersAnalysis(infiles, outfile,
                                 submit=True, job_options=job_options)


# @transform(clusterSpikeInsPowerAnalysis,
#           suffix(".analysis.out"),
#           ".M3D_plot_list.out")
# def clusterSpikeInsPowerPlotM3D(infiles, outfile):
#
#    job_options = "-l mem_free=23G"
#
#    RRBS.spikeInClustersPlotM3D(infiles, outfile, groups=["Saline", "Dex"],
#                                submit=True, job_options=job_options)
# P.touch(outfile)
#
# @transform(generateClusterSpikeIns,
#           suffix(".out"),
#           ".Biseq_plot_list.out")
# def clusterSpikeInsPowerPlotBiSeq(infiles, outfile):
#
#    job_options = "-l mem_free=48G -pe dedicated 8"
#
#    RRBS.spikeInClustersAnalysisBiSeq(infiles, outfile,
#                                      submit=True, job_options=job_options)


########################################################################
# to do: utilise farm.py to split task up
# will need to write a new iterator function for farm.py
# to keep clusters together

@follows(mkdir("subframes.dir"))
@split(subsetCpGsToCovered,
       "subframes.dir/cluster_subframe_*.tsv")
def splitClustersDataframe(infile, outfiles):
    outprefix = "subframes.dir/cluster_subframe_"
    suffix = ".tsv"

    job_options = "-l mem_free=8G -pe dedicated 1"
    RRBS.splitDataframeClusters(infile, outprefix, suffix,
                                submit=True, job_options=job_options)


####################################################################
# this allows all against all for pairs in EXPERIMENTS
# to do:
# re-code so that user can supply a list of pairs in the pipeline.ini
# e.g PARAMS['pair_x]=y,z where multiple pairs can be given
# these pairs would then overide the itertools combinations


@follows(mkdir("M3D"))
@subdivide(splitClustersDataframe,
           formatter(),
           "M3D/{basename[0]}_M3D_stats_*.tsv",
           "M3D/{basename[0]}_M3D_stats_",
           "design.tsv")
def runM3D(infile, outfile, root, design):
    job_options = "-l mem_free=4G -pe dedicated 1"
    groups = [x for x in itertools.combinations(EXPERIMENTS, 2)]

    # **code repeated - refactor**
    for pair in groups:
        pair = [re.sub("-agg", "", str(x)) for x in pair]
        pair1, pair2 = pair
        pair1_split = pair1.split("-")
        pair2_split = pair2.split("-")
        # only want pairs with one difference
        # e.g treatment or tissue but not both
        if not (pair1_split[0] != pair2_split[0] and
                pair1_split[1] != pair2_split[1]):
            outfile = ("%(root)s%(pair1)s_vs_%(pair2)s.tsv"
                       % locals())
            if pair1_split[0] != pair2_split[0]:
                groups = [pair1_split[0], pair2_split[0]]
            elif pair1_split[1] != pair2_split[1]:
                groups = [pair1_split[1], pair2_split[1]]
            else:
                E.error("This pair does not contain any comparisons: %(pair)s"
                        % locals())

            RRBS.calculateM3DStat(infile, outfile, design, pair=pair,
                                  groups=groups, submit=True,
                                  job_options=job_options)


# @follows(mkdir("M3D_plots.dir"))
# @merge([runM3D, ["design.tsv", "M3D_plots.dir"]],
#       "M3D_plots.dir/cluster_stat_merged.tsv")

@follows(mkdir("M3D_plots.dir"))
@collate(runM3D,
         regex(r"M3D/cluster_subframe_(\d+)_M3D_stats_(\S+)_vs_(\S+).tsv"),
         r"M3D_plots.dir/\2_vs_\3_cluster_stat_merged_between.tsv",
         r"\2", r"\3")
def calculateM3DClustersPvalue(infiles, outfile, pair1, pair2):
    job_options = "-l mem_free=4G -pe dedicated 1"
    infiles = infiles[:-1]
    print "pair1: %s" % pair1
    print "pair2: %s" % pair2
    pair = [pair1, pair2]

    print infiles, outfile, pair
    RRBS.calculateM3Dpvalue(infiles, outfile, pair,
                            submit=True, job_options=job_options)


@transform(calculateM3DClustersPvalue,
           suffix(".tsv"),
           ".load")
def loadM3DClusters(infile, outfile):
    dbh = connect()
    tablename = P.toTable(outfile)

    statement = '''cat %(infile)s |
                python %%(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s''' % locals()
    P.run()


@transform(calculateM3DClustersPvalue,
           suffix(".tsv"),
           "_summary.tsv")
def summariseM3D(infile, outfile):
    ''' summarise the number of cluster passing threshold'''
    # adjusted p-value threshold
    threshold = 0.05
    print infile, outfile, threshold
    RRBS.summariseM3D(infile, outfile, threshold, submit=True)


@merge(summariseM3D,
       ["M3D_plots.dir/summary_table.tsv", "M3D_plots.dir/summary_table.load"])
def combineM3Dsummaries(infiles, outfiles):
    ''' combine M3D summary tables'''
    outfile1, outfile2 = outfiles
    print outfile1, outfile2

    tablename = P.toTable(outfile2)

    statement = ''' python %%(scriptsdir)s/combine_tables.py -v0  -a file
                    --glob=M3D_plots.dir/*between_summary.tsv > %(outfile1)s;
                    cat %(outfile1)s |
                    python %%(scriptsdir)s/csv2db.py
                    --table %(tablename)s --retry --ignore-empty
                    > %(outfile2)s''' % locals()
    print statement
    P.run()


#########################################################################


@follows(calculateM3DSpikeClustersPvalue)
def power():
    pass


@follows(loadM3DClusters,
         combineM3Dsummaries)
def M3D():
    pass


# add this to full once start summarising functions have been refactored
@follows(loadStartSummary)
def startSummary():
    pass


@follows(loadStartSummary,
         loadCoverage,
         loadCpGOverlap,
         loadRemainingReads,
         findCpGs,
         subsetCoverage,
         sortAndIndexBams,
         makeCpgIslandsBed,
         make1basedCpgIslands,
         loadMergeCoverage,
         makeSummaryPlots,
         mergeCoverage,
         plotReadBias,
         power,
         M3D)
def full():
    pass


@follows(mapReadsWithBismark)
def mapReads():
    pass


@follows(mapReadsWithBismark)
def test():
    pass


@follows(runBiSeq)
def biseq():
    pass


@follows(callMethylationStatus,
         plotReadBias)
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
