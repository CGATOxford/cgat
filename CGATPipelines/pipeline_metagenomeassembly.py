/###############################################################################
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
=============================
Metagenome assembly pipeline
=============================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

The metagenome assembly pipeline takes reads from one or more NGS experiments and
assembles into contigs / scaffolds. Genes present on contigs are predicted using ORF
perdiction software.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

Assembly stategy
----------------

While there exist many tools for assembling reads from single genomes, only recently has 
software been specifically developed (or extended) to allow for the assembly of metagenomes.
The major factor affecting the ability to assemble a metagenome is the diversity of the sample.
Accurate assembly of long contigs is difficult in the presence of many species at differential
abundances. This is in constrast to single genome assembly where reads are (or should be) uniformly
sampled across the genome.

This pieline therefore uses a range of available software for the assembly of metagenomes.

Considerations
--------------

Metagenomics is a young and rapidly developing field. There is, as yet, no gold standard for 
assembly. It is likely that the presence of multiple, highly related species will lead to 
the assembly of schimeric contigs i.e. contigs derived from more than one species. It is 
generally considered that longer K-mer lengths used in the construction of the de-bruijn
graph (for de-bruijn graph assemblers) will result in fewer chimeras. Nevertheless, longer
k-mers may also result in more, short contigs being produced as a result of a neccessity for
a greater overlap between nodes in the graph. The length of k-mer chosen is also dependent on
the length of reads that you are trying to assemble - longer reads means you can use longer
k-mers. Which k-mer to use in the assembly process is therefore dependent on the data used
and the expected complexity of the sample. We make no effort here to advise on k-mer length.

  
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
|ray meta            |>=2.2.0            |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|meta-velvet         |>=1.2.02           |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|idba_ud             |>=1.1.0            |Metagenome assembler                            |
+--------------------+-------------------+------------------------------------------------+
|metaphlan           |>=1.7.7            |Taxon abundance estimator                       |
+--------------------+-------------------+------------------------------------------------+
|MetaGeneMark        |>=1.0.0            |ORF prediction                                  |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |>=2.17.0           |BED interval analysis suite                     |
+--------------------+-------------------+------------------------------------------------+
|bowtie              |>=0.12.7           |Short read alignment algorithm                  |
+--------------------+-------------------+------------------------------------------------+
|bowtie2             |>=2.0.0            |Short read alignment algorithm                  |
+--------------------+-------------------+------------------------------------------------+
|blastn              |>=2.2.25           |Simlilarity searching algorithm (nucleotides)   |
+--------------------+-------------------+------------------------------------------------+
|blastp              |>=2.2.25           |Simlilarity searching algorithm (proteins)      |
+--------------------+-------------------+------------------------------------------------+


Pipeline output
===============

The main output is the genome assembly - output as a fasta formatted file.
Additional outputs include taxon abundance estimation (metaphlan) and ORF
predictions (MetaGeneMark).

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
import CGAT.GFF as GFF
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
ASSEMBLERS = P.asList( PARAMS["assemblers"] )
MAPPER = PARAMS["coverage_mapper"]
BOWTIE = MAPPER == "bowtie"
BOWTIE2 = MAPPER == "bowtie2"

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
SEQUENCEFILES_REGEX = regex(r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fastq$|fastq.gz|fastq.1.gz)")
 
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
## preprocess reads for metaphlan and IDBA
###################################################################                                                                                                                                                                          
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"\1.fa")
def preprocessReads(infile, outfile):
    '''
    create merged fasta file for use with metaphlan 
    and IDBA
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
    
        statement = '''python %(scriptsdir)s/fastqs2fasta.py 
                   -a %(infile)s 
                   -b %(read2)s 
                   --log=%(infile)s.log 
                   > %(outfile)s'''
        P.run()

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
    statement = PipelineMetagenomeAssembly.Metavelvet().build(infile)
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
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)
    
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
@transform(preprocessReads, regex("(\S+).fa"), r"idba.dir/\1..contigs.fa")
def runIdba(infile, outfile):
    '''
    run idba on each track
    '''
    to_cluster = True
    job_options = " -l mem_free=30G"
    statement = PipelineMetagenomeAssembly.Idba().build(infile)
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runIdba, suffix(".contigs.fa"), ".summary.tsv")
def buildIdbaStats(infile, outfile):
    '''
    build idba stats:
    N50
    Number of scaffolds 
    Total scaffold length
    '''
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

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
    job_options=" -l mem_free=30G -pe dedicated 5"
    statement = PipelineMetagenomeAssembly.Ray().build(infile)
    P.run()

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
@transform(ASSEMBLY_TARGETS, suffix(".contigs.fa"), ".filtered.contigs.fa")
def filterContigs(infile, outfile):
    '''
    filter contigs if specified in .ini file. If not specified
    then the pipeline will not remove any but will produce a new 
    outfile - this is not space efficient and SHOULD BE CHANGED
    '''
    if not PARAMS["filter"]:
        length = 0
    else:
        length = PARAMS["filter"]
    
    PipelineMetagenomeAssembly.filterContigs(infile, outfile, length)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".summary.tsv")
def buildContigStats(infile, outfile):
    '''
    build contig stats:
    N50
    Number of scaffolds 
    Total scaffold length
    max length
    '''
    PipelineMetagenomeAssembly.contig_to_stats(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigStats, regex("(\S+).dir/(\S+).tsv"), r"\1.dir/\1-\2.load")
def loadContigStats(infile, outfile):
    '''
    load the contig stats
    '''
    P.load(infile, outfile)

@follows(loadContigStats)
def contig_stats():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@split(loadContigStats, "*/contig.summary.tsv")
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
@transform(filterContigs, suffix(".fa"), ".lengths.tsv")
def buildContigLengths(infile, outfile):
    '''
    output lengths for each contig in each of the assemblies
    '''
    PipelineMetagenomeAssembly.build_scaffold_lengths(infile, outfile, PARAMS)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildContigLengths, suffix(".lengths.tsv"), ".lengths.load")
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
## gene finding using MetaGeneMark
###################################################################                                                                                                                                                                          
@transform(filterContigs, suffix(".fa"), ".genes.tsv")
def findGenesUsingMetaGeneMark(infile, outfile):
    '''
    Use the MetaGeneMark hmm software to predict genes.
    Output is tsv - similar to gff format but with
    sequences
    '''
    to_cluster=True
    mparams = PARAMS["metagenemark_model_params"]
    statement = '''gmhmmp -a -d -f G -m %(mparams)s -o %(outfile)s %(infile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.gff.gz")
def parseGenesGff(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format gff 
                                    --log=%(outfile)s.log 
                                     | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.fasta.gz")
def parseGenesFasta(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format fasta 
                                    --log=%(outfile)s.log 
                                    | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(findGenesUsingMetaGeneMark, regex("(\S+).tsv"), r"\1.aa.gz")
def parseGenesAa(infile, outfile):
    '''
    parse the genes file
    '''
    to_cluster = True
    statement = '''cat %(infile)s | python %(scriptsdir)s/formatMetagenemark.py 
                                    --format aa 
                                    --log=%(outfile)s.log 
                                    | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(mkdir("genes.dir"))
@transform(parseGenesAa, regex("(\S+).dir/(\S+).aa.gz"), r"genes.dir/\1_\2.blast.result")
def runBlastOnAminoAcidSequences(infile, outfile):
    '''
    look for homology with known genes
    '''
    to_cluster = True
    db = PARAMS["blastp_db"]
    evalue = PARAMS["blastp_evalue"]
    if PARAMS["blastp_ungapped"]:
        ungapped = "-ungapped"
    else:
        ungapped = ""

    statement = '''zcat %(infile)s 
                  | python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" --chunksize=10000 "blastp -db %(db)s -evalue %(evalue)s
                  -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq'"
                  | python %(scriptsdir)s/blast2table.py 
                  --alignment-format=emissions > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(runBlastOnAminoAcidSequences, suffix(".result"), r".result.load")  
def loadBlastOnAminoAcidSequences(infile, outfile):
    '''
    load blastp results
    '''
    P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform([parseGenesGff,parseGenesFasta,parseGenesAa]
           ,regex("(\S+).dir/(\S+).genes.(\S+).gz")
           , r"\1.dir/\1_\2.genes.\3.tsv.gz")
def buildGeneTables(infile, outfile):
    '''
    build gene tables
    '''
    to_cluster = True
    if infile.endswith(".gff.gz"): 
        outf = gzip.open(outfile, "w")
        outf.write("chr\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
        for line in gzip.open(infile).readlines():
            outf.write(line)
        outf.close()
    else:
        statement = '''zcat %(infile)s | python %(scriptsdir)s/fasta2table.py 
                                         -s sequence 
                                         --log=%(outfile)s.log | gzip > %(outfile)s'''
        P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
jobs_limit(1, "db")
@transform(buildGeneTables, regex("(\S+)/(\S+).genes.(\S+).tsv.gz"), r"\1/\2.genes.\3.load")
def loadGeneTables(infile, outfile):
     '''
     load genes from metagenemaek analysis
     '''
     if infile.find("gff") != -1:
         P.load(infile, outfile)
     else:
         P.load(infile, outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(loadGeneTables)
def metagenemark():
    pass

@follows(loadBlastOnAminoAcidSequences)
def geneSimilarity():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## build indices for mapping - this is for coverage analysis
###################################################################                                                                                                                                                                          
@active_if(BOWTIE)
@transform(filterContigs, suffix(".fa"), ".ebwt")
def buildAssemblyBowtieIndices(infile, outfile):
    '''
    build bowtie indices
    '''
    to_cluster = True
    outbase = P.snip(infile, ".fa")
    directory = os.path.dirname(infile)
    statement = '''bowtie-build -f %(infile)s %(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if(BOWTIE2)
@transform(filterContigs, suffix(".fa"), ".bt2")
def buildAssemblyBowtie2Indices(infile, outfile):
    '''
    build bowtie indices
    '''
    to_cluster = True
    outbase = P.snip(infile, ".fa")
    statement = '''bowtie2-build -f %(infile)s %(outbase)s'''
    P.run()
    P.touch(outfile)

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## map with bowtie
###################################################################                                                                                                                                                                          
bowtie_index = {"bowtie":buildAssemblyBowtieIndices
                , "bowtie2":buildAssemblyBowtie2Indices}
BOWTIE_INDEX = bowtie_index[MAPPER]

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("metavelvet" in ASSEMBLERS)
@follows(BOWTIE_INDEX)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"metavelvet.dir/\1.filtered.bam")
def mapReadsWithBowtieAgainstMetavelvetContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    to_cluster = True

    bowtie_index_dir = "metavelvet.dir"
    genome = re.search(".*R[0-9]*", infile).group(0) + ".filtered.contigs"
    infile, reffile = infile,  os.path.join(bowtie_index_dir, genome) + ".fa"
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("idba" in ASSEMBLERS)
@follows(BOWTIE_INDEX)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"idba.dir/\1.filtered.bam")
def mapReadsWithBowtieAgainstIdbaContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    to_cluster = True
    bowtie_index_dir = "idba.dir"
    genome = re.search(".*R[0-9]*", infile).group(0) + ".filtered.contigs"
    infile, reffile = infile,  os.path.join(bowtie_index_dir, genome) + ".fa"
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@active_if("ray" in ASSEMBLERS)
@follows(BOWTIE_INDEX)
@transform(SEQUENCEFILES, SEQUENCEFILES_REGEX, r"ray.dir/\1.filtered.bam")
def mapReadsWithBowtieAgainstRayContigs(infile, outfile):
    '''
    map reads against contigs with bowtie
    '''
    to_cluster = True
    bowtie_index_dir = "ray.dir"
    genome = re.search(".*R[0-9]*", infile).group(0) + ".filtered.contigs"
    infile, reffile = infile,  os.path.join(bowtie_index_dir, genome) + ".fa"
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
ALIGNMENT_TARGETS = []
alignment_targets = {"metavelvet":mapReadsWithBowtieAgainstMetavelvetContigs
                     , "idba":mapReadsWithBowtieAgainstIdbaContigs
                     , "ray":mapReadsWithBowtieAgainstRayContigs}
for x in ASSEMBLERS:
    ALIGNMENT_TARGETS.append(alignment_targets[x])

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform("ray.dir/sim-15M250BP-R1.bam"
           , regex("(\S+).dir/(\S+).bam")
           , r"\1.dir/\1_\2.unique")
def countUniqueAlignments(infile, outfile):
    '''
    count the number of uniquely mapping reads
    '''
    samfile = pysam.Samfile(infile)
    c_unique = 0
    c_total = 0
    E.info("iterating over bam file")
    for alignment in samfile.fetch():
        c_total += 1
        nh = [x[1] for x in alignment.tags if x[0] == "NH"]
        print nh
        # hack - need to discover why some reads do not contain
        # nh tag - maybe if not mapped
        if len(nh) == 0: continue
        if nh[0] == 1:
            c_unique += 1
            print c_unique
    E.info("writing results for %s" % infile)
    outf = open(outfile, "w")
    outf.write("total_reads\tpunique\n%i\t%f\n" % (c_total, float(c_unique/c_total)))
    outf.close()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(ALIGNMENT_TARGETS
           , regex("(\S+).dir/(\S+).bam")
           , r"\1.dir/\1_\2.picard_stats")
def buildPicardStats(infile, outfile):
    '''build alignment stats using picard.
    Note that picards counts reads but they are in fact alignments.
    '''
    reffile = P.snip(infile, ".bam") + ".contigs.fa"
    PipelineMappingQC.buildPicardAlignmentStats( infile, 
                                                 outfile,
                                                 reffile )

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit( 1, "db" )
@merge( buildPicardStats, "picard_stats.load" )
def loadPicardStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(*ALIGNMENT_TARGETS)
@transform(buildContigLengths, suffix(".lengths.tsv"), ".coverage.gz")
def buildCoverageOverContigs(infile, outfile):
    '''
    build histograms of the coverage over each of the contigs
    '''
    to_cluster = True
    size = infile
    bam = P.snip(infile, ".contigs.lengths.tsv") + ".bam"
    statement = '''bedtools genomecov -ibam %(bam)s -g %(size)s -d | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@transform(buildCoverageOverContigs, suffix(".gz"), ".stats.gz")
def buildCoverageStats(infile, outfile):
    '''
    build coverage statistics - mean and standard deviation
    '''
    job_options = " -l mem_free=30G"
    to_cluster = True
    statement = '''zcat %(infile)s | python %(scriptsdir)s/coverage2stats.py --log=%(outfile)s.log | gzip > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@jobs_limit(1, "db")
@transform(buildCoverageStats, suffix(".gz"), ".load")
def loadCoverageStats(infile, outfile):
    '''
    load coverage stats
    '''
    tablename = P.toTable(P.snip(os.path.dirname(infile), ".dir") + "_%s" % os.path.basename(outfile))
    statement = '''zcat %(infile)s | python %(scriptsdir)s/csv2db.py 
                -t %(tablename)s 
                --index=contig 
                --log=%(outfile)s.log > %(outfile)s'''
    P.run()

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
@follows(loadCoverageStats)
def coverage():
    pass

###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
###################################################################                                                                                                                                                                          
## annotate metagenomic reads with metaphlan
###################################################################                                                                                                                                                                          
@follows(mkdir("metaphlan.dir"))
@transform(preprocessReads
           , regex("(\S+).fa")
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
    statement = PipelineMetagenomeAssembly.Metaphlan().build(infile, method="read_map")
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
@transform(preprocessReads
           , regex("(\S+).fa")
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

    statement = PipelineMetagenomeAssembly.Metaphlan().build(infile, method="rel_ab")
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

#########################################
# taxonomic classification targets
#########################################
@follows(loadMetaphlanRelativeAbundances
         , buildMetaphlanTaxonomicAbundances
         , countMetaphlanTaxonomicGroups
         , loadMetaphlanReadmaps)
def metaphlan():
    pass


####################
# full targets
####################
@follows( loadReadCounts
          , loadContigLengths
          , loadPicardStats
          , coverage
          , contig_stats
          , metaphlan
          , metagenemark)
def full():
    pass

####################
# report building
####################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''
    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''
    E.info( "updating documentation" )
    P.run_report( clean = False )


if __name__== "__main__":
    
    if sys.argv [3] == "plotToFile":
        pipeline_printout_graph ("flowchart.jpg", "jpg", [full],
                                 user_colour_scheme = {"colour_scheme_index" :6})
    else:
        sys.exit( P.main(sys.argv) )
    



