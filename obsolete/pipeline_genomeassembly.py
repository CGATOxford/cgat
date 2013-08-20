###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 
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
#################################################################################
"""
==========================
Genome assembly pipeline
==========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

The genome assembly pipeline takes reads from one or more NGS experiment and
assembles into contigs / scaffolds

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

Assembly stategy
----------------

The strategy used for assembling short reads derived from an NGS experiment into contigs and 
scaffolds depends largely on the origin of the data. For single genomes with an assumed
uniform read distribution it is possible to use single genome assembly algorithms. However,
as the complexity of of a dataset increases e.g. for matagenomic data where genomes from multiple
species need to be assembled, there are further considerations and more specialist tools
are required.

The ability of any genome assembler to accurately reconstruct long contigs will depend on multiple factors.

*The length of sequenced reads.
*The depth of sequencing.
*The k-mer lengths used to contruct de bruijn graphs (for de brijn-based assemblers)

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

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

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|SOAPdenovo2         |>=2.1.2            |Single genome assembler                         |
+--------------------+-------------------+------------------------------------------------+
|meta-velvetg        |>=1.2.02           |Meta-genome assembler                           |
+--------------------+-------------------+------------------------------------------------+
|velvet              |                   |Single genome assembler                         |
+--------------------+-------------------+------------------------------------------------+
|idba                |>=1.1.0            |Meta-genome assembler                           |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main output is the genome assembly - output as a .fasta formatted file.
Additional outputs are stored in the database file :file:`csvdb`.

Glossary
========

.. glossary::

   SOAPdenovo2
     SOAPdenovo2 _ - a single genome assembly tool

   meta-velvetg
     meta-velvetg_ - a meta-genome assembly tool

   velvet
     velvet_ - a single genome assembly tool. Functions from velvet are also
                used in the first set of operations in meta-velvet

   idba
     idba_ - a meta-genome assembly tool

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
import CGAT.GFF as GFF
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineGenomeAssembly as PipelineGenomeAssembly
import CGAT.FastaIterator as FastaIterator
import metaphlan_utils
import CGATPipelines.PipelineMapping as PipelineMapping

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    "pipeline.ini" )


PARAMS = P.PARAMS

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# collect fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" )
            
ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
## Global flags
###################################################################
ASSEMBLERS = P.asList( PARAMS["general_assemblers"] )
METAGENOME = "meta-velvet" in ASSEMBLERS or "ibda" in ASSEMBLERS or "cortex_var" in ASSEMBLERS

ASSEMBLERS = P.asList( PARAMS["assemblers"] )
MAPPER = PARAMS["coverage_mapper"]
BOWTIE = MAPPER == "bowtie"
BWA = MAPPER == "bwa"

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''
    dbh = sqlite3.connect( PARAMS["database"] )
    return dbh

###################################################################
###################################################################
###################################################################
# genome assembly
###################################################################
###################################################################
###################################################################
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz"
                 , "*.fastq","*.fastq.gz", "*.fastq.1.gz")
SEQUENCEFILES_REGEX = regex( r"(\S+).(fastq.1.gz|fasta|fasta.gz|fasta.1.gz|fastq|fastq.gz|fastq.1.gz)")
 
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## load number of reads                                                                                                                                                                                                                      
###################################################################                                                                                                                                                                          
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"\1.nreads" )
def countReads( infile, outfile ):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build( (infile,), outfile )
    P.run()

@merge(countReads, "reads_summary.load" )
def loadReadCounts( infiles, outfile ):
    '''load read counts into database.'''

    outf = P.getTempFile()
    outf.write( "track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile( infile ).readlines()
        nreads = int( lines[0][:-1].split("\t")[1])
        outf.write( "%s\t%i\n" % (track,nreads))
    outf.close()

    P.load( outf.name, outfile )

    os.unlink(outf.name)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## assemble reads with meta-velvet
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@follows(mkdir("metavelvet.dir"))
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX
            , r"metavelvet.dir/\1.contigs.fa")
def runMetavelvet(infile, outfile):
    '''
    run meta-velvet on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G"
    statement = PipelineGenomeAssembly.Metavelvet().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit(1, "R")
@transform(runMetavelvet
           , suffix(".contigs.fa")
           , ".stats.pdf")
def plotCoverageHistogram(infile, outfile):
    '''
    plot the coverage over kmers
    '''
    inf = P.snip(infile, ".contigs.fa") + ".stats.txt"
    outf = P.snip(inf, ".txt") + ".pdf"
    R('''library(plotrix)''')
    R('''data = read.table("%s", header=TRUE)''' % inf)
    R('''pdf("%s", height = 7, width = 7 )''' % outf)
    R('''weighted.hist(data$short1_cov, data$lgth, breaks=seq(0, 200, by=1))''')
    R["dev.off"]()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runMetavelvet
           , suffix(".contigs.fa")
           , ".stats.load")
def loadMetavelvetRawStats(infile, outfile):
    '''
    load the assembly stats for meta-velvet
    '''
    inf = P.snip(infile, ".contigs.fa") + ".stats.txt"
    P.load(inf, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runMetavelvet, suffix(".contigs.fa"), ".summary.tsv")
def buildMetavelvetStats(infile, outfile):
    '''
    build metavelvet stats:
    N50
    Number of scaffolds 
    Total scaffold length
    '''
    PipelineGenomeAssembly.contig_to_stats(infile, outfile, PARAMS)
    
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetavelvetStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadMetavelvetStats(infile, outfile):
    '''
    load the metavelvet stats
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## assemble reads with idba
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS)
@follows(mkdir("idba.dir"))
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX
            , r"idba.dir/\1.scaffold.fa")
def runIdba(infile, outfile):
    '''
    run idba on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G"
    statement = PipelineGenomeAssembly.Idba().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runIdba, suffix(".scaffold.fa"), ".summary.tsv")
def buildIdbaStats(infile, outfile):
    '''
    build idba stats:
    N50
    Number of scaffolds 
    Total scaffold length
    '''
    PipelineGenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildIdbaStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadIdbaStats(infile, outfile):
    '''
    load the idba stats
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@follows(mkdir("ray.dir"))
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX
            , r"ray.dir/\1.contigs.fa")
def runRay(infile, outfile):
    '''
    run Ray on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G"
    statement = PipelineGenomeAssembly.Ray().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runRay, suffix(".contigs.fa"), ".summary.tsv")
def buildRayStats(infile, outfile):
    '''
    build ray stats:
    N50
    Number of scaffolds 
    Total scaffold length
    max length
    '''
    PipelineGenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildRayStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadRayStats(infile, outfile):
    '''
    load the ray stats
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
ASSEMBLY_TARGETS = []
assembly_targets = {"metavelvet": runMetavelvet
                    , "idba": runIdba
                    , "ray": runRay
                    }
for x in ASSEMBLERS:
    ASSEMBLY_TARGETS.append(assembly_targets[x])

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
STATS_TARGETS = []
stats_targets = {"metavelvet": loadMetavelvetStats
                    , "idba": loadIdbaStats
                    , "ray": loadRayStats
                    }
for x in ASSEMBLERS:
    STATS_TARGETS.append(stats_targets[x])

@follows(*STATS_TARGETS)
def contig_stats():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@split(STATS_TARGETS, "*/contig.summary.tsv")
def buildContigSummary(infiles, outfile):
    '''
    merge the contig summary statistics
    '''
    stats = collections.defaultdict(list)
    for filepath in infiles:
        dirname = os.path.dirname(filepath)
        stats[dirname].append(os.path.basename(filepath))

    N = PARAMS["scaffold_n"]

    # connect to database
    dbh = connect()
    cc = dbh.cursor()
    for dirname in stats.keys():
        outfname = os.path.join(dirname, "contig.summary.tsv")
        outf = open(outfname, "w")
        outf.write("track\tnscaffolds\tscaffold_length\tN%i\tmean_length\tmedian_length\tmax_length\n" % N)
        for infile in stats[dirname]:
            track = P.snip(infile.split(dirname.split(".dir")[0])[1][1:], ".summary.load")
            table = P.toTable(infile)
            data = cc.execute("""SELECT nscaffolds
                                 , scaffold_length
                                 , N50
                                 , mean_length
                                 , median_length
                                 , max_length FROM %s""" % table).fetchone()
            outf.write("\t".join(map(str, [track, data[0], data[1], data[2], data[3], data[4], data[5]])) + "\n")
        outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigSummary, suffix(".tsv"), ".load")
def loadContigSummary(infile, outfile):
    '''
    load contig summary stats for each assembler
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + os.path.basename(infile) + ".load"
    P.load(infile, outname)
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(ASSEMBLY_TARGETS, suffix(".fa"), ".lengths.tsv")
def buildContigLengths(infile, outfile):
    '''
    output lengths for each contig in each of the assemblies
    '''
    PipelineGenomeAssembly.build_scaffold_lengths(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigLengths, suffix(".lengths.tsv"), ".load")
def loadContigLengths(infile, outfile):
    '''
    load contig lengths
    '''
    outname = P.snip(os.path.dirname(infile), ".dir") + "_" + P.snip(os.path.basename(infile), ".tsv") + ".load"
    P.load(infile, outname)
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## build indices for mapping - this is for coverage analysis
###################################################################                                                                                                                                                                          
@active_if(BOWTIE)
@transform(ASSEMBLY_TARGETS, suffix(".fa"), ".ebwt")
def buildAssemblyBowtieIndices(infile, outfile):
    '''
    build bowtie indices
    '''
    outbase = TRACKS.getTracks()[0]
    directory = os.path.dirname(infile)
    statement = '''bowtie-build -f %(infile)s %(directory)s/%(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## map with bowtie
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@follows(buildAssemblyBowtieIndices, runMetavelvet)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"metavelvet.dir/\1.bam")
def mapReadsWithBowtieAgainstMetavelvetContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    PARAMS["bowtie_index_dir"] = "metavelvet.dir"
    PARAMS["genome"] = TRACKS.getTracks(infile)[0].split(".")[0]
        
    infile, reffile = infile,  os.path.join("metavelvet.dir", TRACKS.getTracks(infile)[0])
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS)
@follows(buildAssemblyBowtieIndices, runIdba)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"idba.dir/\1.bam")
def mapReadsWithBowtieAgainstIdbaContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    PARAMS["bowtie_index_dir"] = "idba.dir"
    PARAMS["genome"] = TRACKS.getTracks(infile)[0].split(".")[0]
        
    infile, reffile = infile,  os.path.join("idba.dir", TRACKS.getTracks(infile)[0])
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@follows(buildAssemblyBowtieIndices, runRay)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"ray.dir/\1.bam")
def mapReadsWithBowtieAgainstRayContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    PARAMS["bowtie_index_dir"] = "ray.dir"
    PARAMS["genome"] = TRACKS.getTracks(infile)[0].split(".")[0]
        
    infile, reffile = infile,  os.path.join("ray.dir", TRACKS.getTracks(infile)[0])
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@transform(buildContigLengths, suffix(".lengths.tsv")
           , add_inputs(mapReadsWithBowtieAgainstMetavelvetContigs)
           , ".coverage")
def buildCoverageOverMetavelvetContigs(infiles, outfile):
    '''
    build histograms of the coverage over each of the contigs
    '''
    inf = infiles[0].find("metavelvet") != -1 
    if inf:
        size = infiles[0]
        bam = infiles[1]
        statement = '''bedtools genomecov -split -ibam %(bam)s -g %(size)s -d > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS)
@transform(buildContigLengths, suffix(".lengths.tsv")
           , add_inputs(mapReadsWithBowtieAgainstIdbaContigs)
           , ".coverage")
def buildCoverageOverIdbaContigs(infiles, outfile):
    '''
    build histograms of the coverage over each of the contigs
    '''
    inf = infiles[0].find("idba") != -1 
    if inf:
        size = infiles[0]
        bam = infiles[1]
        statement = '''bedtools genomecov -split -ibam %(bam)s -g %(size)s -d > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@transform(buildContigLengths, suffix(".lengths.tsv")
           , add_inputs(mapReadsWithBowtieAgainstRayContigs)
           , ".coverage")
def buildCoverageOverRayContigs(infiles, outfile):
    '''
    build histograms of the coverage over each of the contigs
    '''
    inf = infiles[0].find("ray") != -1 
    if inf:
        size = infiles[0]
        bam = infiles[1]
        statement = '''bedtools genomecov -split -ibam %(bam)s -g %(size)s -d > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
COVERAGE_TARGETS = []
coverage_targets = {"metavelvet": buildCoverageOverMetavelvetContigs
                    , "idba": buildCoverageOverIdbaContigs
                    , "ray": buildCoverageOverRayContigs}

for x in ASSEMBLERS:
    COVERAGE_TARGETS.append(coverage_targets[x])

@follows(*COVERAGE_TARGETS)
def coverage():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## annotate metagenomic reads with metaphlan
###################################################################                                                                                                                                                                          
@follows(mkdir("metaphlan.dir"))
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX
            , r"metaphlan.dir/\1.readmap")
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
        assert os.path.exists(PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for blast?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"
    statement = PipelineGenomeAssembly.Metaphlan().build(infile, method="read_map")
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetaphlanReadmap, suffix(".readmap"), ".readmap.load")
def loadMetaphlanReadmaps(infile, outfile):
    '''
    load the metaphlan read maps
    '''
    P.load(infile,outfile)

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
    taxons = ["_order", "class", "family", "genus", "kingdom", "phylum", "species"]
    dbh = connect()
    cc = dbh.cursor()
    for infile in infiles:
        table = P.toTable(infile)
        track = P.snip(table, "_readmap")
        for taxon in taxons:
            count = cc.execute("""SELECT COUNT(DISTINCT %s) FROM %s""" % (taxon, table)).fetchone()[0]
            outf.write("\t".join([track, taxon, str(count)]) + "\n")
    outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(mkdir("metaphlan.dir"))
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX
            , r"metaphlan.dir/\1.relab")
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
        assert os.path.exists(PARAMS["metaphlan_db"] + ".1.bt2"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + ".1.bt2"
        method = "--bowtie2db"
    elif PARAMS["metaphlan_executable"] == "blast":
        assert os.path.exists(PARAMS["metaphlan_db"] + "nin"), "missing file %s: Are you sure you have the correct database for bowtie2?" % PARAMS["metaphlan_db"] + "nin"
        method = "--blastdb"

    statement = PipelineGenomeAssembly.Metaphlan().build(infile, method="rel_ab")
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildMetaphlanRelativeAbundance, suffix(".relab"), ".relab.load")
def loadMetaphlanRelativeAbundances(infile, outfile):
    '''
    load the metaphlan relative abundances
    '''
    P.load(infile,outfile)

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
            # remove unassigned species at this point
            idx = track.split("_")[1]
            if data[1].find("unclassified") != -1: continue
            outf.write("\t".join([track, data[0], data[1], str(data[2]), idx]) + "\n")
    outf.close()

@follows( coverage
          , contig_stats
          , loadMetaphlanRelativeAbundances)
def full():
    pass


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )




