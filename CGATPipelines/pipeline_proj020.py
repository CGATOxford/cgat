"""
===================
Project 20 pipeline
===================

:Author: David Sims 
:Release: $Id: pipeline_proj020.py 2900 2013-06-03 14:38:00Z david $
:Date: |today|
:Tags: Python

The project 020 pipeline processes reads from a gene trap screen experiment to identify trapped genes:

   * Align reads to the gene trap vector using BLAST


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The pipeline looks for a configuration file in several places:

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
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.


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

Glossary
========

.. glossary::

   bowtie
      bowtie_ - a read mapper

.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml

Code
====

"""
import string
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import re
import glob
import os
import shutil
import collections
import gzip
import sqlite3
import cStringIO
import fileinput
import CGAT.Fastq as fq
import logging as L
import CGAT.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping
from ruffus import *


###################################################################
###################################################################
###################################################################
# Pipeline configuration
import CGAT.Pipeline as P
P.getParameters("pipeline.ini")
PARAMS = P.PARAMS
USECLUSTER = True

###################################################################
###################################################################
###################################################################
# Count raw reads


@transform("*.fastq.1.gz", regex(r"(\S+).fastq.1.gz"), r"\1.nreads")
def countReads(infile, outfile):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()

###################################################################


@merge(countReads, "raw_read_counts.txt")
def mergeReadCounts(infiles, outfile):
    '''Merge counts from all file and load into database'''
    outf = open(outfile, "w")
    outf.write("track\tnreads\n")
    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".nreads")
        rc = open(infile, "r")
        line1 = rc.readline()
        rc.close()
        header, nreads = line1.split()
        outf.write("%(track)s\t%(nreads)s\n" % locals())
    outf.close()

###################################################################


@transform(mergeReadCounts, regex(r"raw_read_counts.txt"), "raw_read_counts.load")
def loadReadCounts(infile, outfile):
    '''load read counts into database'''
    statement = """cat %(infile)s | python %(scriptsdir)s/csv2db.py 
                         --database=%(database)s
                         --table=raw_read_counts
                         --index=track
                 > %(outfile)s; """
    P.run()

###################################################################
###################################################################
###################################################################
# Section 1: MAP READS TO GENE TRAP VECTOR


@follows(mkdir("blast"))
@transform("*.fastq.*.gz", regex(r"(\S+).fastq.(\S+).gz"), r"blast/\1.\2.fasta")
def convertToFasta(infile, outfile):
    '''Convert fastq to fasta using fastx toolkit'''
    to_cluster = True
    qual_score_offset = PARAMS["blast_qual_score_offset"]
    statement = """zcat %(infile)s | fastq_to_fasta -Q %(qual_score_offset)s -n > %(outfile)s """
    P.run()

###################################################################


@transform(convertToFasta, regex(r"(\S+).fasta"), r"\1.blast")
def runBlast(infile, outfile):
    '''map reads to gene trap vector with blast'''
    to_cluster = True
    blast_db = PARAMS["blast_db"]
    statement = """blastn -query %(infile)s -db %(blast_db)s -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" 
                   | python %(scriptsdir)s/blast2table.py --alignment-format=emissions > %(outfile)s"""
    P.run()

###################################################################


@transform(runBlast, regex(r"(\S+).blast"), r"\1.blast.load")
def loadBlast(infile, outfile):
    '''load blast results into database'''
    P.load(infile, outfile, "--index=query_nid --index=pid --index=query_ali")


###################################################################
###################################################################
###################################################################
# Section2: Identify reads containing expected mm Pax2 sequence at the
# start of read1
@follows(mkdir("grep"))
@transform("*.fastq.1.gz", regex(r"(\S+).fastq.1.gz"), r"grep/\1.match")
def grepPrimers(infile, outfile):
    '''count occurences of decreasing primer substrings at start of reads '''
    to_cluster = False
    primer = "a"
    if infile.find("_b.") > 0:
        primer = "b"
    if primer == "a":
        primer_seq = PARAMS["grep_primer_a"]
    else:
        primer_seq = PARAMS["grep_primer_b"]

    for i in range(len(primer_seq), 5, -1):
        primer_subseq = primer_seq[:i]
        statement = '''echo "%(primer_subseq)s" >> %(outfile)s; zcat %(infile)s | grep ^%(primer_subseq)s | wc -l >> %(outfile)s;'''
        P.run()

    # reformat out file
    statement = '''echo "Total reads" >> %(outfile)s; echo `zcat %(infile)s |  wc -l` / 4 | bc >> %(outfile)s;
                   sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

###################################################################


@follows(mkdir("filtered"))
@transform("*.fastq.1.gz", regex(r"(\S+).fastq.1.gz"), (r"filtered/\1.reconciled.fastq.1.gz", r"filtered/\1.reconciled.fastq.2.gz"))
def filterReadsByPrimerMatch(infile, outfiles):
    '''Filter out reads where the start of read 1 does not match primer sequence (14bp)'''
    to_cluster = True
    primer = "a"
    if infile.find("_b.") > 0:
        primer = "b"
    if primer == "a":
        primer_seq = PARAMS["grep_primer_a"]
    else:
        primer_seq = PARAMS["grep_primer_b"]
    grep_filter_length = PARAMS["grep_filter_length"]
    primer_subseq = primer_seq[:grep_filter_length]

    track = P.snip(os.path.basename(infile), ".fastq.1.gz")
    infile2 = track + ".fastq.2.gz"
    outfile1, outfile2 = outfiles
    tempfile = "filtered/" + track + ".filtered.fastq.1.gz"

    # filter by primer match
    fastq_in = open(infile, "r")
    fastq_out = open(tempfile, "wb")
    for read in fq.iterate(fastq_in):
        if read.seq[:grep_filter_length] == primer_subseq:
            fastq_out.writeln("@" + read.id)
            fastq_out.writeln(read.seq)
            fastq_out.writeln("+")
            fastq_out.writeln(read.qual)
    fastq_in.close()
    fastq_out.close()

    # reconcile read pairs
    statement = '''python %(scriptsdir)s/fastqs2fastq.py --method=reconcile %(tempfile)s %(infile2)s --output-pattern=filtered/%(track)s.reconciled.fastq.%%i.gz'''
    P.run()

###################################################################


@follows(mkdir("filtered"))
@transform("*.fastq.1.gz", regex(r"(\S+).fastq.1.gz"), r"filtered/\1.cutadapt.read1.fastq.gz")
def trimGeneTrapVectorRead1(infile, outfile):
    '''trim primer sequence from 5` end of read1 to identify enriched fragments containing gene traps '''
    to_cluster = True
    track = P.snip(os.path.basename(infile), ".fastq.1.gz")
    primer = "a"
    if infile.find("_b.") > 0:
        primer = "b"
    if primer == "a":
        primer_seq = PARAMS["trim_primer_seq_a"]
    else:
        primer_seq = PARAMS["trim_primer_seq_b"]
    trim_options = PARAMS["trim_options"]

    #contaminant_file = PARAMS["trim_contaminants"]
    #adaptors = []
    # for entry in FastaIterator.FastaIterator( IOTools.openFile( contaminant_file ) ):
    #    adaptors.append( "-a %s" % entry.sequence )
    #adaptors= " ".join(adaptors)

    statement = '''cutadapt -b %(primer_seq)s  
                       %(trim_options)s
                       --untrimmed-output=filtered/%(track)s.untrimmed.read1.fastq.gz
                       --too-short-output=filtered/%(track)s.tooshort.read1.fastq.gz
                       -o %(outfile)s
                       %(infile)s
                       &> filtered/%(track)s.cutadapt.read1.log;'''
    P.run()

###################################################################


@follows(mkdir("filtered"))
@transform("*.fastq.2.gz", regex(r"(\S+).fastq.2.gz"), r"filtered/\1.cutadapt.read2.fastq.gz")
def trimGeneTrapVectorRead2(infile, outfile):
    '''trim primer sequence from 3` end read2 to enable alignment '''
    to_cluster = True

    track = P.snip(os.path.basename(infile), ".fastq.2.gz")
    adaptor_seq = PARAMS["trim_pax2_exon3"]
    trim_options = PARAMS["trim_options"]
    adaptor_rc = adaptor_seq[::-1]
    adaptor_seq = adaptor_seq.translate(string.maketrans("ACTGN", "TGACN"))

    statement = '''cutadapt -a %(adaptor_seq)s -a %(adaptor_rc)s
                       %(trim_options)s
                       --too-short-output=filtered/%(track)s.tooshort.read2.fastq.gz
                       -o %(outfile)s
                       %(infile)s
                       &> filtered/%(track)s.cutadapt.read2.log;'''
    P.run()

###################################################################


@collate((trimGeneTrapVectorRead1, trimGeneTrapVectorRead2), regex(r"(\S+).cutadapt.read(\S+).fastq.gz"), (r"\1.cutadapt.reconciled.read1.fastq.gz", r"\1.cutadapt.reconciled.read2.fastq.gz"))
def reconcileReadPairs(infiles, outfiles):
    '''Remove read2 if read1 did not contain match to gene trap'''
    to_cluster = True
    read1, read2 = infiles
    out1, out2 = outfiles
    track = P.snip(os.path.basename(read1), ".cutadapt.read1.fastq.gz")
    # reconcile read pairs
    statement = '''python %(scriptsdir)s/fastqs2fastq.py --method=reconcile %(read1)s %(read2)s --output-pattern=filtered/%(track)s.cutadapt.reconciled.read%%i.fastq.gz'''
    P.run()

###################################################################


@transform(reconcileReadPairs, regex(r"filtered/(\S+).cutadapt.reconciled.read1.fastq.gz"), r"filtered/\1.cutadapt.reconciled.read1.nreads")
def countTaggedReads(infiles, outfile):
    '''count number of reads in input files.'''
    to_cluster = True
    read1, read2 = infiles
    m = PipelineMapping.Counter()
    statement = m.build((read1,), outfile)
    P.run()

###################################################################


@merge(countTaggedReads, "tagged_read_counts.txt")
def mergeTaggedReadCounts(infiles, outfile):
    '''Merge counts from all file and load into database'''
    outf = open(outfile, "w")
    outf.write("track\tnreads\n")
    for infile in infiles:
        track = P.snip(
            os.path.basename(infile), ".cutadapt.reconciled.read1.nreads")
        rc = open(infile, "r")
        line1 = rc.readline()
        rc.close()
        header, nreads = line1.split()
        outf.write("%(track)s\t%(nreads)s\n" % locals())
    outf.close()

###################################################################


@transform(mergeTaggedReadCounts, regex(r"tagged_read_counts.txt"), "tagged_read_counts.load")
def loadTaggedReadCounts(infile, outfile):
    '''load read counts into database'''
    statement = """cat %(infile)s | python %(scriptsdir)s/csv2db.py 
                         --database=%(database)s
                         --table=tagged_read_counts
                         --index=track
                 > %(outfile)s; """
    P.run()

###################################################################
# align read 2 against transcriptome


@follows(reconcileReadPairs, mkdir("transcriptome"))
@transform("filtered/*.cutadapt.reconciled.read2.fastq.gz", regex(r"filtered/(\S+).cutadapt.reconciled.read2.fastq.gz"), r"transcriptome/\1.read2.bam")
def alignReadsToTranscriptome(infile, outfile):
    '''map reads to transcriptome with bowtie'''
    to_cluster = True
    track = P.snip(os.path.basename(outfile), ".bam")
    job_options = "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    m = PipelineMapping.Bowtie()
    reffile = PARAMS["bowtie_transcriptome"]
    bowtie_options = PARAMS["bowtie_options"]
    statement = m.build((infile,), outfile)
    P.run()

###################################################################
###################################################################
# align read2 against genome using tophat using pipeline_mapping.py
# use tophat directory
# python %(scriptsdir)s/pipeline_mapping.py -v 5 -p 10 make full
###################################################################
###################################################################

###################################################################
# Alignment stats


@transform("tophat/tophat.dir/*.bam", suffix(".bam"), ".alignstats")
def buildPicardAlignStats(infile, outfile):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = True
    track = P.snip(os.path.basename(infile), ".bam")
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(picard_genome)s ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals(
    )
    P.run()

############################################################


@merge(buildPicardAlignStats, "tophat/tophat.dir/picard_align_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = open("tophat/tophat.dir/picard_align_stats.tsv", "w")
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

#########################################################################


@transform("tophat/tophat.dir/*.bam", suffix(".bam"), ".dedup.stats")
def dupStats(infile, outfile):
    '''Remove duplicate alignments from BAM files.'''
    to_cluster = USECLUSTER
    track = P.snip(infile, ".bam")
    statement = '''MarkDuplicates INPUT=%(infile)s ASSUME_SORTED=true OUTPUT=/dev/null 
                   METRICS_FILE=%(outfile)s VALIDATION_STRINGENCY=SILENT; ''' % locals()
    P.run()

#########################################################################


@merge(dupStats, "tophat/tophat.dir/picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = open('tophat/tophat.dir/dupstats.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".dedup.stats")
        statfile = f
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

###################################################################
# quantify reads over transcripts


@follows(mkdir("htseq-counts"))
@transform("tophat/tophat.dir/*.bam", regex(r"tophat/tophat.dir/(\S+).bam"), r"htseq-counts/\1.htseq-counts.gz")
def quantitateWithHTSeq(infile, outfile):
    '''Use htseq script to count reads overlapping genes'''
    gtf = PARAMS["htseq_gtf"]
    options = PARAMS["htseq_options"]
    bamfile = infile
    to_cluster = True
    temp_out = "%s_%s" % (os.path.basename(bamfile), os.path.basename(gtf))
    statement = ''' samtools sort -on %(bamfile)s %(temp_out)s.sorted | samtools view -
                   | grep -v chrM
                   | htseq-count - %(gtf)s %(options)s
                   | gzip -c > %(outfile)s 2> %(outfile)s.log '''
    P.run()

###################################################################


@merge(quantitateWithHTSeq, r"htseq-counts/agg-agg-agg.tsv.gz")
def combineSampleCounts(infiles, outfile):
    '''Merge the genecounts for all samples'''
    to_cluster = True
    inlist = " ".join(infiles)
    statement = ''' python %(scriptsdir)s/combine_tables.py
                        -c 1
                        -k 2
                        -t
                        --use-file-prefix
                        --regex-filename='(\S+).cutadapt.reconciled.read2.tophat.htseq-counts.gz'
                        -L %(outfile)s.log
                         %(inlist)s > %(outfile)s '''
    P.run()

###################################################################


@follows(convertToFasta,
         runBlast,
         loadBlast)
def blast():
    '''run the pipeline'''
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
