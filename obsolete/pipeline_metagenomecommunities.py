"""
=====================================================
Community analysis of metgenomic shotgun sequencing
=====================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Pipeline_metagenomecommunities.py takes as input a set of fastq
files from a shotgun sequencing experiment of environmental samples
and assesses community structure and function.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions
(:term:`experiment`) with one or more biological and/or technical
replicates (:term:`replicate`). A :term:`replicate` within each
:term:`experiment` is a :term:`track`.

Community profiling
--------------------

The pipeline uses various tools for assessing the abundance of taxa 
in an environmental sample. To assess relative abundance of taxa in
a community, we use metaphlan as it is easy to use and because it performs
alignments against a reduced set of clade-specific marker genes it is relatively
fast. Nevertheless, we also use lcamapper (from MEGAN), which in contrast to metaphlan
attempts to assign every read to a taxonomic group using the lowest common ancestor
approach.


Functional profiling
---------------------

Whether DNA-seq or RNA-seq is used, functional profiling can be performed. Pathway-centric
and gene-centric approaches are used and implemented using lcamapper.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineDocumenation`). To
start with, use the files supplied with the :ref:`Example` data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the
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
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|blastn              |>=2.2.25           |Simlilarity searching algorithm (nucleotides)   |
+--------------------+-------------------+------------------------------------------------+
|blastp              |>=2.2.25           |Simlilarity searching algorithm (proteins)      |
+--------------------+-------------------+------------------------------------------------+
|metaphlan           |                   |Community profiling                             |
+--------------------+-------------------+------------------------------------------------+
|lcamapper           |                   |Community profiling + fuctional profiling       |
+--------------------+-------------------+------------------------------------------------+



Pipeline output
===============

TODO::

Additional outputs are stored in the database file :file:`csvdb`.

Glossary
========

.. glossary::

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

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMetagenomeAssembly as PipelineMetagenomeAssembly
import CGAT.FastaIterator as FastaIterator
import CGAT.Metaphlan as Metaphlan
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import pysam
import CGAT.Fastq as Fastq

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["pipeline.ini"])


PARAMS = P.PARAMS

###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

# collect fastq.gz tracks
TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
    glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        glob.glob("*.fastq.1.gz"), "(\S+).fastq.1.gz")

ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

###################################################################
# sequence files as input
###################################################################
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz",
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz")
SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fastq$|fastq.gz|fastq.1.gz)")

###################################################################
# connecting to database
###################################################################


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''
    dbh = sqlite3.connect(PARAMS["database"])
    return dbh

###################################################################
###################################################################
###################################################################
# load number of reads
###################################################################
###################################################################
###################################################################


@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"\1.nreads")
def countReads(infile, outfile):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build((infile,), outfile)
    P.run()


@merge(countReads, "reads_summary.load")
def loadReadCounts(infiles, outfile):
    '''load read counts into database.'''

    to_cluster = False
    outf = P.getTempFile()
    outf.write("track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile(infile).readlines()
        nreads = int(lines[0][:-1].split("\t")[1])
        outf.write("%s\t%i\n" % (track, nreads))
    outf.close()
    inname = outf.name

    tablename = P.toTable(outfile)
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log 
                  < %(inname)s > %(outfile)s'''
    P.run()
    os.unlink(outf.name)

###################################################################
###################################################################
###################################################################
# Preprocessing reads for community analysis
###################################################################
###################################################################
###################################################################


@follows(mkdir("fasta.dir"))
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"fasta.dir/\1.fa")
def preprocessReads(infile, outfile):
    '''
    create merged fasta file for use with metaphlan 
    '''
    # check for second read in the pair
    if infile.endswith(".fastq.gz"):
        E.info("converting fastq file to fasta file")
        outf = open(outfile, "w")
        for fastq in Fastq.iterate(IOTools.openFile(infile)):
            outf.write("%s\n%s\n" % (">" + fastq.identifier, fastq.seq))
        outf.close()
    elif infile.endswith(".1.gz"):
        read2 = P.snip(infile, ".1.gz") + ".2.gz"
        assert os.path.exists(read2), "file does not exist %s" % read2

        log = infile.replace("fastq.", "")
        statement = '''python %(scriptsdir)s/fastqs2fasta.py 
                   -a %(infile)s 
                   -b %(read2)s 
                   --log=%(log)s.log 
                   > %(outfile)s'''
        P.run()

###################################################################
###################################################################
###################################################################
# Estimate taxonomic relative abundances using metaphlan
###################################################################
###################################################################
###################################################################


@follows(mkdir("metaphlan.dir"))
@transform(preprocessReads, regex("(\S+)/(\S+).fa"), r"metaphlan.dir/\2.readmap")
def buildMetaphlanReadmap(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via 
    blastn searching
    '''
    to_cluster = True

    # at present the pipeline will take a set of files
    # and compute the abundances of different taxonomic groups
    # based on ALL reads i.e. paired data are combined into
    # a single file for analysis
    if PARAMS["metaphlan_executable"] == "bowtie2":
        assert os.path.exists(
            PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(
            PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for blast?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"
    statement = PipelineMetagenomeAssembly.Metaphlan().build(
        infile, method="read_map")
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildMetaphlanReadmap, suffix(".readmap"), ".readmap.load")
def loadMetaphlanReadmaps(infile, outfile):
    '''
    load the metaphlan read maps
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@merge(loadMetaphlanReadmaps, "metaphlan.dir/taxonomic.counts")
def countMetaphlanTaxonomicGroups(infiles, outfile):
    '''
    count the total number of species that
    were found by metaphlan
    '''
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\tcount\n")
    taxons = ["_order", "class", "family",
              "genus", "kingdom", "phylum", "species"]
    dbh = connect()
    cc = dbh.cursor()
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_readmap")
        for taxon in taxons:
            count = cc.execute(
                """SELECT COUNT(DISTINCT %s) FROM %s""" % (taxon, table)).fetchone()[0]
            outf.write("\t".join([track, taxon, str(count)]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@follows(mkdir("metaphlan.dir"))
@transform(preprocessReads, regex("(\S+)/(\S+).fa"), r"metaphlan.dir/\2.relab")
def buildMetaphlanRelativeAbundance(infile, outfile):
    '''
    metaphlan is a program used in metagenomics. It assigns
    reads to clades based on specific genetic markers via 
    blastn searching
    '''
    to_cluster = True
    # at present the pipeline will take a set of files
    # and compute the abundances of different taxonomic groups
    # based on ALL reads i.e. paired data are combined into
    # a single file for analysis
    if PARAMS["metaphlan_executable"] == "bowtie2":
        assert os.path.exists(
            PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(
            PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"

    statement = PipelineMetagenomeAssembly.Metaphlan().build(
        infile, method="rel_ab")
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildMetaphlanRelativeAbundance, suffix(".relab"), ".relab.load")
def loadMetaphlanRelativeAbundances(infile, outfile):
    '''
    load the metaphlan relative abundances
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@merge(loadMetaphlanRelativeAbundances, "metaphlan.dir/taxonomic.abundances")
def buildMetaphlanTaxonomicAbundances(infiles, outfile):
    '''
    build a file that combines taxonomic abundances 
    from each sample
    '''
    dbh = connect()
    cc = dbh.cursor()
    outf = open(outfile, "w")
    outf.write("track\ttaxon_level\ttaxon\tabundance\tidx\n")
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_relab")
        for data in cc.execute("""SELECT taxon_level, taxon, rel_abundance FROM %s""" % table).fetchall():
            idx = track.split("_")[1]
            outf.write(
                "\t".join([track, data[0], data[1], str(data[2]), idx]) + "\n")
    outf.close()

#########################################
# metaphlan target
#########################################


@follows(loadMetaphlanRelativeAbundances,
         buildMetaphlanTaxonomicAbundances,
         countMetaphlanTaxonomicGroups,
         loadMetaphlanReadmaps)
def metaphlan():
    pass

###################################################################
###################################################################
###################################################################
# Assign reads to taxa using Lowest Common Ancestor (LCA)
###################################################################
###################################################################
###################################################################


@follows(mkdir("lca.dir"))
@transform(preprocessReads, regex("(\S+)/(\S+).fa"), r"lca.dir/\2.blast.gz")
def runBlastOnRawSequences(infile, outfile):
    '''
    A blast alignment is the first step in analysis of 
    community structure and function using the LCA approach
    implemented in MEGAN. For large numbers of reads this
    may become infeasible in terms of running time and the
    shear size of the blast output. TODO: add blast options
    to parameters available
    '''
    db = PARAMS["megan_db"]
    evalue = PARAMS["megan_evalue"]
    blast_options = PARAMS["megan_blast_options"]
    statement = '''cat %(infile)s | python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" 
                    --log=%(outfile)s.log 
                    --chunksize=100 "blastx -db %(db)s %(blast_options)s -evalue %(evalue)s" 
                  | gzip > %(outfile)s; checkpoint
                  '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(runBlastOnRawSequences, suffix(".blast.gz"), ".lca.gz")
def runLCA(infile, outfile):
    '''
    run the lowest common ancestor algorithm
    on the blast output to assign reads to
    taxa - from mtools. Runs with defaults at
    the moment. At this point we also build
    the functional profile using KEGG - this 
    comes with mtools package
    '''
    track = P.snip(outfile, ".lca.gz")
    outf_tax = P.snip(outfile, ".gz")
    outf_kegg = P.snip(outfile, ".lca.gz") + ".kegg.gz"
    statement = '''lcamapper.sh 
                   -k
                   -i %(infile)s
                   -o %(outf_tax)s > %(outfile)s.log
                   ; gzip %(track)s.blast-kegg.txt
                   ; gzip %(outf_tax)s
                   ; mv %(track)s.blast-kegg.txt.gz %(outf_kegg)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".level.count")
def countLcaPerLevelTaxa(infile, outfile):
    '''
    count the number of taxa as found using LCA algorithm
    '''
    statement = '''zcat %(infile)s |
                   python %(scriptsdir)s/lca2table.py
                   --summarise=level-counts
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(countLcaPerLevelTaxa, suffix(".count"), ".count.load")
def loadCountLcaPerLevelTaxa(infile, outfile):
    '''
    load taxa level counts
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".taxa.count")
def countLcaTaxa(infile, outfile):
    '''
    count the number of taxa as found using LCA algorithm
    '''
    statement = '''zcat %(infile)s |
                   python %(scriptsdir)s/lca2table.py
                   --summarise=taxa-counts
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(countLcaTaxa, suffix(".count"), ".count.load")
def loadCountLcaTaxa(infile, outfile):
    '''
    load taxa level counts
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@transform(runLCA, suffix(".lca.gz"), ".taxa.gz")
def buildLCA(infile, outfile):
    '''
    tabulate LCA output into nice format. Per read assignment
    '''
    statement = '''zcat %(infile)s 
                   | python %(scriptsdir)s/lca2table.py
                     --summarise=individual
                     --log=%(outfile)s.log
                   | sed -e 's/order/_order/g'
                   | gzip > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildLCA, suffix(".taxa.gz"), ".taxa.readcounts")
def countContributingReads(infile, outfile):
    '''
    count number of reads with a taxnomic assignment
    '''
    levels = [
        "phylum", "class", "order", "family", "genus", "species"]
    result = collections.OrderedDict()
    for level in levels:
        result[level] = 0
    inf = IOTools.openFile(infile)
    header = inf.readline().split("\t")

    # column indices
    indices = [3, 5, 7, 9, 11, 13]
    total = 0
    for line in inf.readlines():
        total += 1
        data = line[:-1].split("\t")
        phylum, _class, order, family, genus, species = [
            data[i] for i in indices]
        if phylum != "NA":
            result["phylum"] += 1
        if _class != "NA":
            result["class"] += 1
        if order != "NA":
            result["order"] += 1
        if family != "NA":
            result["family"] += 1
        if genus != "NA":
            result["genus"] += 1
        if species != "NA":
            result["species"] += 1
    outf = open(outfile, "w")
    outf.write("level\tn_reads\tpct_reads\n")
    for level, count in result.iteritems():
        nreads, prop = count, float(count) / total
        outf.write("\t".join([level, str(nreads), str(prop * 100)]) + "\n")
    outf.close()

###################################################################
###################################################################
###################################################################


@transform(countContributingReads, suffix(".readcounts"), ".readcounts.load")
def loadCountContributingReads(infile, outfile):
    '''
    load contributing read counts
    '''
    P.load(infile, outfile)

#########################################
# LCA target
#########################################


@follows(buildLCA,
         loadCountLcaPerLevelTaxa,
         loadCountLcaTaxa,
         loadCountContributingReads)
def lca():
    pass

###################################################################
###################################################################
###################################################################
# Counting KEGG associations
###################################################################
###################################################################
###################################################################


@follows(mkdir("kegg.dir"))
@transform(PARAMS.get("kegg_tre"),
           regex("(\S+)/(\S+).tre"),
           add_inputs(PARAMS.get("kegg_map")),
           r"kegg.dir/\2.tsv")
def buildKeggTable(infiles, outfile):
    '''
    build kegg table mapping KO identifiers to pathways. This is
    based on the file that was downloaded with mtools (D.Huson)
    '''
    keggtre, keggmap = infiles
    statement = '''python %(scriptsdir)s/keggtre2table.py 
                   --kegg-tre=%(keggtre)s
                   --map=%(keggmap)s 
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(buildKeggTable, suffix(".tsv"), ".load")
def loadKeggTable(infile, outfile):
    '''
    load KEGG table
    '''
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################


@follows(mkdir("kegg.dir"))
@transform(runLCA,
           regex("(\S+)/(\S+).lca.gz"),
           add_inputs(buildKeggTable),
           r"kegg.dir/\2.kegg.counts")
def countKeggAssociations(infiles, outfile):
    '''
    count the number of reads associted with Kegg pathways
    '''
    infile = infiles[0].replace(".lca", ".kegg")
    kegg_table = infiles[1]
    level = PARAMS.get("kegg_level")
    statement = '''zcat %(infile)s |
                   python %(scriptsdir)s/lcakegg2counts.py
                   --kegg-table=%(kegg_table)s
                   --level=%(level)s
                   --method=proportion
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(countKeggAssociations, suffix(".counts"), ".counts.load")
def loadCountKeggAssociations(infile, outfile):
    '''
    load counts of KO associations
    '''
    P.load(infile, outfile, "--header=pathway,p_annotated_reads")

###################################################################
###################################################################
###################################################################


@follows(mkdir("kegg.dir"))
@transform(runLCA,
           regex("(\S+)/(\S+).lca.gz"),
           add_inputs(buildKeggTable),
           r"kegg.dir/\2.kegg.ko.counts")
def countKeggGenes(infiles, outfile):
    '''
    count the number of reads associted with KO identifiers
    '''
    infile = infiles[0].replace(".lca", ".kegg")
    kegg_table = infiles[1]
    level = PARAMS.get("kegg_level")
    statement = '''zcat %(infile)s |
                   python %(scriptsdir)s/lcakegg2counts.py
                   --kegg-table=%(kegg_table)s
                   --level=D
                   --method=proportion
                   --log=%(outfile)s.log
                   > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(countKeggGenes, suffix(".ko.counts"), ".ko.counts.load")
def loadCountKeggGenes(infile, outfile):
    '''
    load counts of KO associations
    '''
    P.load(infile, outfile, "--header=KO,p_annotated_reads")

#########################################
# kegg target
#########################################


@follows(loadCountKeggAssociations)
def kegg():
    pass

#########################################
# full target
#########################################


@follows(loadReadCounts, metaphlan, lca, kegg)
def full():
    pass

#########################################
#########################################
#########################################
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
