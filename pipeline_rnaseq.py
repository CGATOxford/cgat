################################################################################
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
#################################################################################
"""
====================
RNA-Seq pipeline
====================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The RNA-Seq pipeline imports unmapped reads from one or more
RNA-Seq experiments and performs the basic RNA-Seq analysis steps:

   1. Map reads to genome
   2. Build transcript models
   3. Estimate differential expression

This pipeline works on a single genome.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

The pipeline performs the following tasks:

   * analyse each experiment:
      * for each replicate
          * map reads using tophat for each term:`replicate` separately. 
          * predict splice isoforms and expression levels with :term:`cufflinks`.
          * estimate expression levels of reference gene set with :term:`cufflinks`.
          * annotate isoforms in replicates with genomic annotations
      * compare isoforms in replicates within each :term:`experiment` (:term:`cuffcompare`) 
          and to reference gene set.
      * summary statistics on reproducibility within each experiment
   * build a combined gene set including the reference gene set and isoforms predicted by :term:`cufflinks`.
      * compare all isoforms in all experiments+isoforms (:term:`cuffcompare`) to each other 
         and the reference gene set
        * summary statistics on isoforms with respect to gene set
   * estimate differential expression levels of transcripts
      * different gene sets
         * reference gene set
         * combined gene set
         * novel gene set
      * different methods
         * :term:`DESeg` (tag counting)
         * :term:`cuffdiff``
      * summary statistics on differential expression

Mapping strategy
----------------

The best strategy for mapping and transcriptome assembly depends on the length of your reads.
With short reads, detecting novel splice-junctions is a difficult task. In this case it will be
best to rely on a set of known splice-junctions. Longer reads map more easily across splice-junctions.

From the tophat manual::

   TopHat finds splice junctions without a reference annotation. By first mapping RNA-Seq reads 
   to the genome, TopHat identifies potential exons, since many RNA-Seq reads will contiguously 
   align to the genome. Using this initial mapping, TopHat builds a database of possible splice 
   junctions, and then maps the reads against this junction to confirm them.

   Short read sequencing machines can currently produce reads 100bp or longer, but many exons are 
   shorter than this, and so would be missed in the initial mapping. TopHat solves this problem 
   by splitting all input reads into smaller segments, and then mapping them independently. The segment 
   alignments are "glued" back together in a final step of the program to produce the end-to-end read alignments.

   TopHat generates its database of possible splice junctions from three sources of evidence. The first 
   source is pairings of "coverage islands", which are distinct regions of piled up reads in the 
   initial mapping. Neighboring islands are often spliced together in the transcriptome, so 
   TopHat looks for ways to join these with an intron. The second source is only used when 
   TopHat is run with paired end reads. When reads in a pair come from different exons of a 
   transcript, they will generally be mapped far apart in the genome coordinate space. When 
   this happens, TopHat tries to "close" the gap between them by looking for subsequences of 
   the genomic interval between mates with a total length about equal to the expected distance 
   between mates. The "introns" in this subsequence are added to the database. The third, and 
   strongest, source of evidence for a splice junction is when two segments from the same read #
   are mapped far apart, or when an internal segment fails to map. With long (>=75bp) reads, 
   "GT-AG", "GC-AG" and "AT-AC" introns be found ab initio. With shorter reads, TopHat only 
   reports alignments across "GT-AG" introns

Thus, in order to increase the sensitivity of splice-site detection, it might be best to derive a set of 
splice-junctions using all reads. This is not done automatically, but can be done manually by 
adding a file with junctions to the ``tophat_options`` entry in the configuration file.

The pipeline supplies tophat with a list of all coding exons to facilitate mapping across known
splice-junctions. If they are prioritized, I do not know.

Transcripts are built individually for each :term:`track`. This seems to be the most rigorous way
as there might be conflicting transcripts between replicates and merging the sets might confuse transcript
reconstruction. Also, conflicting transcripts between replicates give an idea of the variability of the data. 
However, if there are only few reads,  there might be a case for building transcript models using reads 
from all replicates of an experiment. However, there is no reason to merge reads between experiments.

LncRNA
--------

One of the main benefits of RNASeq over microarrays is that novel transcripts can be detected. A particular
interest are currently novel long non-coding RNA. Unfortunately, it seems that these transcripts
are often expressed at very low levels and possibly in a highly regulated manner, for example only in
certain tissues. On top of their low abundance, they frequently seem to be co-localized with protein
coding genes, making it hard to distinguish them from transcription artifacts.

Success in identifying lincRNA will depend a lot on your input data. Long, paired-end reads are likely 
to lead to success. Unfortunately, many exploratory studies go for single-ended, short read data.
With such data, identification of novel spliced transcripts will be rare and the set of novel transcripts
is likely to contain many false-positives.

The pipeline constructs a set of novel lncRNA in the following manner:
   1. All transcript models overlapping protein coding transcripts are removed.
   2. Overlapping lncRNA on the same strand are merged.

Artifacts in lncRNA analysis
++++++++++++++++++++++++++++

There are several sources of artifacts in lncRNA analysis

Read mapping errors 
~~~~~~~~~~~~~~~~~~~

Mapping errors are identifyable as sharp peaks in the
coverage profile. Mapping errors occur if the true location of a read has more mismatches
than the original location or it maps across an undetected splice-site. Most of the 
highly-expressed lncRNA are due to mapping errors. Secondary locations very often overlap
highly-expressed protein-coding genes. These errors are annoying for two reasons: they
provide false positives, but at the same time prevent the reads to be counted towards
the expression of the true gene.

They can be detected in two ways:

1. via a peak-like distribution of reads which should result in a low entropy of start position 
density. Note that this possibly can remove transcripts that are close to the length of a single
read.
          
2. via mapping against known protein coding transcripts. However, getting this mapping right
is hard for two reasons. Firstly, mapping errors usually involve reads aligned with mismatches.
Thus, the mapping has to be done either on the read-level (computationally expensive), or 
on the transcript level after variant calling. (tricky, and also computationally expensive). 
Secondly, as cufflinks extends transcripts generously, only a part of a transcript might actually
be a mismapped part. Distinguishing partial true matches from random matches will be tricky.

Read mapping errors can also be avoided by

1. Longer read lengths
2. Strict alignment criteria
3. A two-stage mapping process.

Fragments
~~~~~~~~~

As lncRNA are expressed at low levels, it is likely that only a partial transcript
can be observed.


Differential expression
-----------------------

The quality of the rnaseq data (read-length, paired-end) determines the
quality of transcript models. For instance, if reads are short (35bp) and/or 
reads are not paired-ended, transcript models will be short and truncated.
In these cases it might be better to concentrate the analysis on only previously
known transcript models. 

The pipeline offers various sets for downstream analysis of differential expression.

1. A set of previously known transcripts (:file:`reference.gtf.gz`). Use this set if 
   only interested in the transcription of previously known transcripts or read length
   does not permit transcript assembly. This does include all transcripts within the
   ENSEMBL gene set, including processed but untranscribed transcripts, transcripts with
   retained introns, pseudogenes, etc.
   
2. A set of previously known protein coding transcripts (:file:`refcoding.gtf.gz`).
   This set is derived from (:file:`reference.gtf.gz`) but only includes exons of 
   transcripts that are protein coding.

3. An ab-initio gene set (:file:`abinitio.gtf.gz`). The ab-initio set is built by running :term:`cuffcompare` on
   the combined individual :term:`cufflinks` results. Transcripts that have been observed in only
   one :term:`track` are removed (removed transcripts end up in :file:`removed.gtf.gz`) in order
   to exclude partial transcripts. Use this set if reads are of good length and/or are paired-ended.

4. A set of novel transcribed loci (:file:`novel.gtf.gz`). This gene set is derived from the
   set of ab-initio transcripts. All ab-initio transcripts overlapping protein coding transcripts 
   in :file:`refcoding.gtf.gz` are removed. Overlapping transcripts are merged into a single 
   transcript/gene. This removes individual transcript structure, but retains constitutive 
   introns. This set retains transcripts that are only observed in a single experiment. It also
   includes known non-coding transcripts, so a locus might not necessarily be "novel".

Transcripts are the natural choice to measure expression of. However other quantities
might be of interest. Some quantities are biological meaningful, for example differential 
expression from a promotor shared by several trancripts. Other quantities might no biologically
meaningful but are necessary as a technical comprise.
For example, the overlapping transcripts might be hard to resolve and thus might need to be
aggregated per gene. Furthermore, functional annotation is primarily associated with genes
and not individual transcripts. The pipeline attempts to measure transcription and differential
expression for a variety of entities following the classification laid down by :term:`cuffdiff`:

isoform
   Transcript level
gene
   Gene level, aggregates several isoform/transcripts
tss
   Transcription start site. Aggregate all isoforms starting from the same :term:`tss`.
cds
   Coding sequence expression. Ignore reads overlapping non-coding parts of transcripts (UTRs, etc.). Requires
   annotation of the cds and thus only available for :file:`reference.gtf.gz`.

Methods differ in their ability to measure transcription on all levels. 
 
.. todo::
   add promoters and splicing output

Overprediction of differential expression for low-level expressed transcripts with :term:`cuffdiff`
is a `known problem <http://seqanswers.com/forums/showthread.php?t=6283&highlight=fpkm>`_.

Estimating coverage
-------------------

An important question in RNASeq analysis is if the sequencing has been done to sufficient depth.
The questions split into two parts:

   * What is the minimum abundant transcript that should be detectable with the number of reads
     mapped? See for example `PMID: 20565853 <http://www.ncbi.nlm.nih.gov/pubmed/20565853>`

   * What is the minimum expression change between two conditions that can be reliably inferred?
     See for examples `PMID: 21498551 <http://www.ncbi.nlm.nih.gov/pubmed/21498551?dopt=Abstract>`

These questions are difficult to answer due to the complexity of RNASeq data: Genes have multiple
transcripts, transcript/gene expression varies by orders of magnitude and a large fraction of
reads might stem from repetetive RNA.

See `figure 4 <http://www.nature.com/nbt/journal/v28/n5/full/nbt.1621.html>`_ from the cufflinks paper
to get an idea about the reliability of transcript construction with varying sequencing depth.

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

sra
   Short-Read Archive format. Reads will be extracted using the :file:`fastq-dump` tool.

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

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie_             |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|tophat_             |>=1.3.1            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|cufflinks_          |>=1.0.3            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+
|R/DESeq             |                   |differential expression                         |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.42             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|bamstats_           |>=1.22             |from CGR, Liverpool                             |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

For each :term:`experiment` there will be the following tables:

<track>_cuffcompare_benchmark 
   results from comparing gene models against reference gene set
   primary key: track

<track>_cuffcompare_transcripts
   transcript expression values (FPKMs)
   primary key: track+transfrag_id
   foreign key: transfrag_id

<track>_cuffcompare_tracking
   tracking information linking loci against transcripts.
   primary key: transfrag_id, 
   foreign key: locus_id

<track>_cuffcompare_tracking
   locus information (number of transcripts within locus per track)
      primary key: locus_id

Differential gene expression results
-------------------------------------

Differential expression is estimated for different genesets
with a variety of methods. Differential expression can be defined
for various levels.

<geneset>_<method>_<level>_diff
    Results of the pairwise tests for differential expression
    primary keys: track1, track2

<geneset>_<method>_<level>_levels
    Expression levels

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_rnaseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseq.tgz
   tar -xvzf pipeline_rnaseq.tgz
   cd pipeline_rnaseq
   python <srcdir>/pipeline_rnaseq.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   cufflinks
      cufflinks_ - transcriptome analysis

   tophat
      tophat_ - a read mapper to detect splice-junctions

   bowtie
      bowtie_ - a read mapper


   
.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011

Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import GFF, GTF, IOTools, IndexedFasta
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

import PipelineGeneset
import PipelineMapping
import PipelineRnaseq
import PipelineMappingQC
import Stats

# levels of cuffdiff analysis
# (no promotor and splice -> no lfold column)
CUFFDIFF_LEVELS = ("gene", "cds", "isoform", "tss" )

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ],
    defaults = {
        'annotations_dir' : "",
        'paired_end' : False } )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "*.sra" ), "(\S+).sra" ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" ) +\
            PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" )

ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
##################################################################
##################################################################
## genesets build - defined statically here, but could be parsed
## from configuration options
## Needs to done in turn to be able to turn off potentially empty gene sets
## such as refnoncoding
GENESETS = ("novel", "abinitio", "reference", "refcoding", "refnoncoding", "lincrna" )

# reference gene set for QC purposes
REFERENCE= "refcoding"

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

USECLUSTER=True

#########################################################################
#########################################################################
#########################################################################
def writePrunedGTF( infile, outfile ):
    '''remove various gene models from a gtf file.
    '''
    to_cluster = USECLUSTER

    cmds = []

    rna_file = os.path.join( PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_rna_gff"] )

    if "geneset_remove_repetetive_rna" in PARAMS:
        
        cmds.append( '''python %s/gtf2gtf.py
        --remove-overlapping=%s
        --log=%s.log''' % (PARAMS["scriptsdir"], 
                           rna_file, outfile ) )
                     
    if "geneset_remove_contigs" in PARAMS:
        cmds.append( '''awk '$1 !~ /%s/' ''' % PARAMS["geneset_remove_contigs"] )

    cmds = " | ".join( cmds )

    if infile.endswith(".gz"):
        uncompress = "zcat"
    else:
        # wastefull
        uncompress = "cat"
        
    if outfile.endswith(".gz"):
        compress = "gzip"
    else:
        compress = "cat"

    # remove \0 bytes within gtf file
    statement = '''%(uncompress)s %(infile)s 
    | %(cmds)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=contig+gene --log=%(outfile)s.log
    | %(compress)s > %(outfile)s'''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def mergeAndFilterGTF( infile, outfile, logfile ):
    '''sanitize transcripts file for cufflinks analysis.

    Merge exons separated by small introns (< 5bp).

    Transcript will be ignored that
       * have very long introns (max_intron_size) (otherwise, cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)

    This method preserves all features in a gtf file (exon, CDS, ...).

    returns a dictionary of all gene_ids that have been kept.
    '''

    max_intron_size =  PARAMS["max_intron_size"]

    c = E.Counter()

    outf = gzip.open( outfile, "w" )

    E.info( "filtering by contig and removing long introns" )    
    contigs = set(IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) ).getContigs())

    rx_contigs = None
    if "geneset_remove_contigs" in PARAMS:
        rx_contigs = re.compile( PARAMS["geneset_remove_contigs"] )
        E.info( "removing contigs %s" % PARAMS["geneset_remove_contigs"] )

    rna_index = None
    if "geneset_remove_repetetive_rna" in PARAMS:
        rna_file = os.path.join( PARAMS["annotations_dir"],
                                 PARAMS_ANNOTATIONS["interface_rna_gff"] )
        if not os.path.exists( rna_file ):
            E.warn( "file '%s' to remove repetetive rna does not exist" % rna_file )
        else:
            rna_index = GFF.readAndIndex( GFF.iterator( IOTools.openFile( rna_file, "r" ) ) )
            E.info( "removing ribosomal RNA in %s" % rna_file )
    
    gene_ids = {}

    logf = IOTools.openFile( logfile, "w" )
    logf.write( "gene_id\ttranscript_id\treason\n" )

    for all_exons in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile )) ):

        c.input += 1
        
        e = all_exons[0]
        # filtering 
        if e.contig not in contigs:
            c.missing_contig += 1
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "missing_contig" )) + "\n" )
            continue

        if rx_contigs and rx_contigs.search(e.contig):
            c.remove_contig += 1
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "remove_contig" )) + "\n" )
            continue

        if rna_index and all_exons[0].source != 'protein_coding':
            found = False
            for exon in all_exons:
                if rna_index.contains( e.contig, e.start, e.end ):
                    found = True
                    break
            if found:
                logf.write( "\t".join( (e.gene_id, e.transcript_id, "overlap_rna" )) + "\n" )
                c.overlap_rna += 1
                continue
                
        is_ok = True

        # keep exons and cds separate by grouping by feature
        all_exons.sort( key = lambda x: x.feature )
        new_exons = []

        for feature, exons in itertools.groupby( all_exons, lambda x: x.feature ):

            tmp = sorted( list(exons) , key = lambda x: x.start )

            gene_ids[tmp[0].transcript_id] = tmp[0].gene_id
            
            l, n = tmp[0], []

            for e in tmp[1:]:
                d = e.start - l.end
                if d > max_intron_size:
                    is_ok = False
                    break
                elif d < 5:
                    l.end = max(e.end, l.end)
                    c.merged += 1
                    continue

                n.append( l )
                l = e

            n.append( l )
            new_exons.extend( n )

            if not is_ok: break

        if not is_ok: 
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "bad_transcript" )) + "\n" )
            c.skipped += 1
            continue

        new_exons.sort( key = lambda x: x.start )

        for e in new_exons:
            outf.write( "%s\n" % str(e) )
            c.exons += 1

        c.output += 1


    outf.close()

    L.info( "%s" % str(c) )
    
    return gene_ids

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "reference.gtf.gz" )
def buildReferenceGeneSet( infile, outfile ):
    '''sanitize ENSEMBL transcripts file for cufflinks analysis.

    Merge exons separated by small introns (< 5bp).

    Removes unwanted contigs according to configuration
    value ``geneset_remove_contigs``.

    Removes transcripts overlapping ribosomal genes if 
    ``geneset_remove_repetitive_rna`` is set. Protein
    coding transcripts are not removed.

    Transcript will be ignored that
       * have very long introns (max_intron_size) (otherwise, cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)
       
    The result is run through cuffdiff in order to add the p_id and tss_id tags
    required by cuffdiff. 

    This will only keep sources of the type 'exon'. It will also remove
    any transcripts not in the reference genome.

    Cuffdiff requires overlapping genes to have different tss_id tags.

    This gene is the source for most other gene sets in the pipeline.
    '''

    tmpfilename = P.getTempFilename( "." )
    tmpfilename2 = P.getTempFilename( "." )
    tmpfilename3 = P.getTempFilename( "." )

    gene_ids = mergeAndFilterGTF( infile, tmpfilename, "%s.removed.gz" % outfile )
    
    #################################################
    E.info( "adding tss_id and p_id" )

    # The p_id attribute is set if the fasta sequence is given.
    # However, there might be some errors in cuffdiff downstream:
    #
    # cuffdiff: bundles.cpp:479: static void HitBundle::combine(const std::vector<HitBundle*, std::allocator<HitBundle*> >&, HitBundle&): Assertion `in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id()' failed.
    #
    # I was not able to resolve this, it was a complex
    # bug dependent on both the read libraries and the input reference gtf files
    statement = '''
    cuffcompare -r <( gunzip < %(tmpfilename)s )
         -T 
         -s %(cufflinks_genome_dir)s/%(genome)s.fa
         -o %(tmpfilename2)s
         <( gunzip < %(tmpfilename)s )
         <( gunzip < %(tmpfilename)s )
    > %(outfile)s.log
    '''
    P.run()

    #################################################
    E.info( "resetting gene_id and transcript_id" )

    # reset gene_id and transcript_id to ENSEMBL ids
    # cufflinks patch: 
    # make tss_id and p_id unique for each gene id
    outf = IOTools.openFile( tmpfilename3, "w")
    map_tss2gene, map_pid2gene = {}, {}
    inf = IOTools.openFile( tmpfilename2 + ".combined.gtf" )

    def _map( gtf, key, val, m ):
        if val in m:
            while gene_id != m[val]:
                val += "a" 
                if val not in m: break
        m[val] = gene_id

        gtf.setAttribute( key, val )
        
    for gtf in GTF.iterator( inf ):
        transcript_id = gtf.oId
        gene_id = gene_ids[transcript_id] 
        gtf.setAttribute( "transcript_id", transcript_id )
        gtf.setAttribute( "gene_id", gene_id )
        
        # set tss_id
        try: tss_id = gtf.tss_id
        except AttributeError: tss_id = None
        try: p_id = gtf.p_id
        except AttributeError: p_id = None
        
        if tss_id: _map( gtf, "tss_id", tss_id, map_tss2gene)
        if p_id: _map( gtf, "p_id", p_id, map_pid2gene)

        outf.write( str(gtf) + "\n" )
        
    outf.close()

    # sort gtf file
    PipelineGeneset.sortGTF( tmpfilename3, outfile )

    os.unlink( tmpfilename )
    # make sure tmpfilename2 is NEVER empty
    assert tmpfilename2
    for x in glob.glob( tmpfilename2 + "*" ): os.unlink( x )
    os.unlink( tmpfilename3 )

#########################################################################
#########################################################################
#########################################################################
@transform( buildReferenceGeneSet, 
            suffix("reference.gtf.gz"),
            "refnoncoding.gtf.gz" )
def buildNoncodingGeneSet( infile, outfile ):
    '''build a new gene set without protein coding 
    transcripts.'''
    
    to_cluster = True
    statement = '''
    zcat %(infile)s | awk '$2 == "lincRNA" || $2 == "processed_transcript"' | gzip > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "reference_with_cds.gtf.gz" )
def buildReferenceGeneSetWithCDS( infile , outfile ):
    '''build a new gene set without protein coding 
    transcripts.'''
    
    mergeAndFilterGTF( infile, outfile, "%s.removed.gz" % outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildReferenceGeneSet, 
            suffix("reference.gtf.gz"),
            "%s.gtf.gz" % REFERENCE )
def buildCodingGeneSet( infile, outfile ):
    '''build a gene set with only protein coding 
    transcripts.

    Genes are selected via their gene biotype in the GTF file.
    Note that this set will contain all transcripts of protein
    coding genes, including processed transcripts.

    This set includes UTR and CDS.
    '''
    
    to_cluster = True
    statement = '''
    zcat %(infile)s | awk '$2 == "protein_coding"' | gzip > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildReferenceGeneSet, 
            suffix("reference.gtf.gz"),
            "refcodingtranscripts.gtf.gz" )
def buildCodingTranscriptSet( infile, outfile ):
    '''build a gene set with only protein coding transcripts.

    Protein coding transcripts are selected via the ensembl
    transcript biotype
    '''

    dbh = connect()

    statement = '''SELECT DISTINCT transcript_id FROM transcript_info WHERE transcript_biotype = 'protein_coding' '''
    cc = dbh.cursor()
    transcript_ids = set( [x[0] for x in cc.execute(statement)] )
    
    inf = IOTools.openFile( infile )
    outf = IOTools.openFile( outfile, 'w')
    
    for g in GTF.iterator( inf ):
        if g.transcript_id in transcript_ids:
            outf.write( str(g) + "\n" )
    
    outf.close()
    inf.close()


#########################################################################
#########################################################################
#########################################################################
@transform( buildCodingGeneSet, 
            suffix( "%s.gtf.gz" % REFERENCE ),
            "%s.gff.gz" % REFERENCE )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.
    '''
    PipelineGeneset.buildGeneRegions( infile, outfile )
    

@transform( buildGeneRegions,
            suffix( "%s.gff.gz" % REFERENCE ),
            "%s.terminal_exons.bed.gz" % REFERENCE )
def buildTerminalExons( infile, outfile ):
    '''output terminal truncated exons.'''
    
    size = 50 

    outf1 = IOTools.openFile( outfile, "w" )
        
    for gff in GTF.flat_gene_iterator( GTF.iterator_filtered( GTF.iterator( IOTools.openFile( infile )),
                                                              feature = "exon" )):
        gene_id, contig, strand = gff[0].gene_id, gff[0].contig, gff[0].strand

        gff.sort( key = lambda x: x.start )

        if strand == "-":
            exon = gff[0].start, gff[0].end
            cut_exon = gff[0].start, gff[0].start + size
        elif strand == "+":
            exon = gff[-1].start, gff[-1].end
            cut_exon = gff[-1].end - size, gff[-1].end
        else:
            continue

        outf1.write( "%s\t%i\t%i\t%s\t%i\t%s\n" % (contig, cut_exon[0], cut_exon[1], gene_id, 0, strand ) )
        
    outf1.close()

    '''
    zcat refcoding_terminal_exons.bed.gz | python ~/cgat/bed2gff.py --as-gtf | python ~/cgat/gtf2table.py --counter=read-coverage --bam-file=heart-library-R1.accepted.bam >& exon.coverage
    '''

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_geneset_flat_gtf"]),
        "introns.gtf.gz" )
def buildIntronGeneModels(infile, outfile ):
    '''build protein-coding intron-transcipts.

    Intron-transcripts are the reverse complement of transcripts.

    Only protein coding genes are taken.

    10 bp are truncated on either end of an intron and need
    to have a minimum length of 100.

    Introns from nested genes might overlap, but all exons
    are removed.
    '''

    to_cluster = True
    
    filename_exons = os.path.join( PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_exons_gtf"] )

    statement = '''gunzip
        < %(infile)s 
        | awk '$2 == "protein_coding"'
	| python %(scriptsdir)s/gtf2gtf.py --sort=gene 
	| python %(scriptsdir)s/gtf2gtf.py 
               --exons2introns 
               --intron-min-length=100 
               --intron-border=10 
               --log=%(outfile)s.log
        | python %(scriptsdir)s/gff2gff.py
               --crop=%(filename_exons)s
               --log=%(outfile)s.log
	| python %(scriptsdir)s/gtf2gtf.py 
              --set-transcript-to-gene 
              --log=%(outfile)s.log 
	| perl -p -e 's/intron/exon/'
        | gzip
        > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".junctions")
def buildJunctions( infile, outfile ):
    '''build file with splice junctions from gtf file.
    
    A junctions file is a better option than supplying a GTF
    file, as parsing the latter often fails. See:

    http://seqanswers.com/forums/showthread.php?t=7563

    '''
    
    outf = IOTools.openFile( outfile, "w" )
    njunctions = 0
    for gffs in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile, "r" ) )):
        
        gffs.sort( key = lambda x: x.start )
        end = gffs[0].end
        for gff in gffs[1:]:
            # subtract one: these are not open/closed coordinates but the 0-based coordinates
            # of first and last residue that are to be kept (i.e., within the exon).
            outf.write( "%s\t%i\t%i\t%s\n" % (gff.contig, end - 1, gff.start, gff.strand ) )
            end = gff.end
            njunctions += 1
            
    outf.close()

    if njunctions == 0:
        E.warn( 'no junctions found in gene set' )
        return
    else:
        E.info( 'found %i junctions before removing duplicates' % njunctions )


    # make unique
    statement = '''mv %(outfile)s %(outfile)s.tmp; 
                   cat < %(outfile)s.tmp | sort | uniq > %(outfile)s;
                   rm -f %(outfile)s.tmp; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".fa")
def buildReferenceTranscriptome( infile, outfile ):
    '''build reference transcriptome. 

    The reference transcriptome contains all known 
    protein coding transcripts.

    The sequences include both UTR and CDS.

    '''

    to_cluster = USECLUSTER

    statement = '''
    zcat %(infile)s
    | awk '$3 == "exon"'
    | python %(scriptsdir)s/gff2fasta.py
        --is-gtf
        --genome=%(genome_dir)s/%(genome)s
        --log=%(outfile)s.log
    | perl -p -e "if (/^>/) { s/ .*$// }"
    | python %(scriptsdir)s/sequence2sequence.py -v 0
    | fold 
    > %(outfile)s;
    checkpoint; 
    samtools faidx %(outfile)s
    ''' 

    P.run()
    
    prefix = P.snip( outfile, ".fa" )

    # build raw index
    statement = '''
    bowtie-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
    '''

    P.run()

    # build color space index
    statement = '''
    bowtie-build -C -f %(outfile)s %(prefix)s_cs >> %(outfile)s.log 2>&1
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra",
             "*.csfasta.gz",
             "*.csfasta.F3.gz" ),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)"), 
            add_inputs( buildReferenceTranscriptome ), 
            r"\1.trans.bam" )
def mapReadsWithBowtieAgainstTranscriptome( infiles, outfile ):
    '''map reads from short read archive sequence using bowtie against
    transcriptome data.
    '''

    # Mapping will permit up to one mismatches. This is sufficient
    # as the downstream filter in rnaseq_bams2bam requires the
    # number of mismatches less than the genomic number of mismatches.
    # Change this, if the number of permitted mismatches for the genome
    # increases.

    # Output all valid matches in the best stratum. This will 
    # inflate the file sizes due to matches to alternative transcripts
    # but otherwise matches to paralogs will be missed (and such
    # reads would be filtered out).
    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    to_cluster = USECLUSTER
    m = PipelineMapping.BowtieTranscripts()
    infile, reffile = infiles
    prefix = P.snip( reffile, ".fa" )
    bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_options"] 
    statement = m.build( (infile,), outfile ) 
    P.run()


#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra",
             "*.csfasta.gz",
             "*.csfasta.F3.gz",
             ),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)"), 
            add_inputs( buildJunctions), 
            r"\1.genome.bam" )
def mapReadsWithTophat( infiles, outfile ):
    '''map reads from .fastq or .sra files.

    A list with known splice junctions is supplied.
    '''
    job_options= "-pe dedicated %i -R y" % PARAMS["tophat_threads"]
    
    if "--butterfly-search" in PARAMS["tophat_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_options += " -l mem_free=50G"

    to_cluster = USECLUSTER
    m = PipelineMapping.Tophat()
    infile, reffile = infiles
    tophat_options = PARAMS["tophat_options"] + " --raw-juncs %(reffile)s " % locals()
    statement = m.build( (infile,), outfile ) 
    P.run()


#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@merge( (mapReadsWithTophat, buildJunctions), "junctions.fa" )
def buildJunctionsDB( infiles, outfile ):
    '''build a database of all junctions.'''

    to_cluster = USECLUSTER
    outfile_junctions = outfile + ".junctions.bed.gz"
    min_anchor_length = 3
    read_length = 50

    tmpfile = P.getTempFile( "." )

    for infile in infiles:
        if infile.endswith(".bam"):
            junctions_file = P.snip( infile, ".bam" ) + ".junctions.bed.gz"
            columns = (0,1,2,5)
        else:
            junctions_file = infile
            columns = (0,1,2,3)

        if not os.path.exists( junctions_file ):
            E.warn( "can't find junctions file '%s'" % junctions_file )
            continue

        inf = IOTools.openFile( junctions_file )
        for line in inf:
            if line.startswith("#"): continue
            if line.startswith("track"): continue
            data = line[:-1].split("\t")
            try:
                tmpfile.write( "\t".join( [data[x] for x in columns] ) + "\n" )
            except IndexError:
                raise IndexError( "parsing error in line %s" % line)

    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
    sort %(tmpfilename)s | gzip > %(outfile_junctions)s
    '''

    P.run()

    os.unlink( tmpfilename )

    E.info( "building junctions database" )
    statement = '''
    juncs_db %(min_anchor_length)i %(read_length)i 
              <( zcat %(outfile_junctions)s )
              /dev/null /dev/null 
              %(bowtie_index_dir)s/%(genome)s.fa
              > %(outfile)s
              2> %(outfile)s.log
    '''
    P.run()

    E.info( "indexing junctions database" )

    prefix = P.snip( outfile, ".fa" )

    # build raw index
    statement = '''
    bowtie-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
    '''

    P.run()

    # build color space index
    statement = '''
    bowtie-build -C -f %(outfile)s %(prefix)s_cs >> %(outfile)s.log 2>&1
    '''

    P.run()

if "tophat_add_separate_junctions" in PARAMS and PARAMS["tophat_add_separate_junctions"]:
#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
    @transform( ("*.fastq.1.gz", 
                 "*.fastq.gz",
                 "*.sra",
                 "*.csfasta.gz",
                 "*.csfasta.F3.gz" ),
                regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)"), 
                add_inputs( buildJunctionsDB, 
                            os.path.join(PARAMS["annotations_dir"],
                                         PARAMS_ANNOTATIONS["interface_contigs"])), 
                r"\1.junc.bam" )
    def mapReadsWithBowtieAgainstJunctions( infiles, outfile ):
        '''map reads from short read archive sequence using bowtie against
        junctions.
        '''

        # Mapping will permit up to one mismatches. This is sufficient
        # as the downstream filter in rnaseq_bams2bam requires the
        # number of mismatches less than the genomic number of mismatches.
        # Change this, if the number of permitted mismatches for the genome
        # increases.

        # Output all valid matches in the best stratum. This will 
        # inflate the file sizes due to matches to alternative transcripts
        # but otherwise matches to paralogs will be missed (and such
        # reads would be filtered out).
        job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
        to_cluster = USECLUSTER
        m = PipelineMapping.BowtieJunctions()
        infile, reffile, contigsfile = infiles
        prefix = P.snip( reffile, ".fa" )
        bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_options"] 
        statement = m.build( (infile,), outfile ) 
        P.run()
else:
    @transform( ("*.fastq.1.gz", 
                 "*.fastq.gz",
                 "*.sra",
                 "*.csfasta.gz",
                 "*.csfasta.F3.gz" ),
                regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)"), 
                r"\1.junc.bam" )
    def mapReadsWithBowtieAgainstJunctions( infiles, outfile ):
        P.touch(outfile)
        
############################################################
############################################################
############################################################

@follows( mkdir( os.path.join( PARAMS["exportdir"], "fastqc" ) ) )
@transform( mapReadsWithTophat, suffix( ".bam"), ".fastqc" )
def buildFastQCReport( infile, outfile ):
    '''run fastqc on aligned reads.'''
    
    to_cluster = USECLUSTER
    statement = '''fastqc --outdir=%(exportdir)s/fastqc %(infile)s >& %(outfile)s''' 
    P.run()

############################################################
############################################################
############################################################
@collate( (mapReadsWithTophat, 
           mapReadsWithBowtieAgainstTranscriptome, 
           mapReadsWithBowtieAgainstJunctions),
          regex(r"(.+)\..*.bam"),  
          add_inputs( buildCodingGeneSet ), 
          r"\1.accepted.bam" )
def buildBAMs( infiles, outfile):
    '''reconcile genomic and transcriptome matches.
    '''
    genome, transcriptome, junctions, reffile = infiles[0][0], infiles[2][0], infiles[1][0], infiles[0][1]

    outfile_mismapped = P.snip(outfile, ".accepted.bam") + ".mismapped.bam"

    assert genome.endswith( ".genome.bam" )

    to_cluster = USECLUSTER

    options = []
    if "tophat_unique" in PARAMS and PARAMS["tophat_unique"]:
        options.append( "--unique" )

    if "tophat_remove_contigs" in PARAMS and PARAMS["tophat_remove_contigs"]:
        options.append( "--remove-contigs=%s" % PARAMS["tophat_remove_contigs"] )

    if "tophat_remove_rna" in PARAMS and PARAMS["tophat_remove_rna"]:
        options.append( "--filename-regions=<( zcat %s | grep 'repetetive_rna' )" %\
                            (os.path.join( 
                    PARAMS["annotations_dir"], 
                    PARAMS_ANNOTATIONS["interface_genomic_context_bed"]) ) )

    if "tophat_remove_mismapped" in PARAMS and PARAMS["tophat_remove_mismapped"]:
        options.append( "--filename-transcriptome=%(transcriptome)s" % locals() )

    if "tophat_add_separate_junctions" in PARAMS and PARAMS["tophat_add_separate_junctions"]:
        options.append( "--filename-junctions=%(junctions)s" % locals() )

    options = " ".join(options)

    tmpfile = P.getTempFilename()

    prefix = P.snip( outfile, ".bam")

    if os.path.exists( "%(outfile)s.log" % locals() ):
        os.remove( "%(outfile)s.log" % locals() )
        
    statement = '''
    python %(scriptsdir)s/rnaseq_bams2bam.py 
       --force
       --filename-gtf=%(reffile)s
       --filename-mismapped=%(outfile_mismapped)s
       --log=%(outfile)s.log
       --filename-stats=%(outfile)s.tsv
       %(options)s
       %(genome)s
    | samtools sort - %(prefix)s 2>&1 >> %(outfile)s.log;
    checkpoint;
    samtools index %(outfile_mismapped)s 2>&1 >> %(outfile)s.log;
    checkpoint;
    samtools index %(outfile)s 2>&1 >> %(outfile)s.log;
    '''

    P.run()

############################################################
############################################################
############################################################
@transform( buildBAMs, suffix(".accepted.bam"), ".mismapped.bam" )
def buildMismappedBAMs( infile, outfile ):
    '''pseudo target - update the mismapped bam files.'''
    P.touch( outfile )

############################################################
############################################################
############################################################
@transform( (mapReadsWithTophat, buildBAMs, buildMismappedBAMs), 
            suffix(".bam" ), ".picard_stats")
def buildPicardStats( infile, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    PipelineMappingQC.buildPicardAlignmentStats( infile, outfile,
                                                 os.path.join( PARAMS["cufflinks_genome_dir"],
                                                               PARAMS["genome"] + ".fa" ) )

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "bamstats" ) ) )
@transform( (buildMismappedBAMs, mapReadsWithTophat, buildBAMs ), 
            suffix(".bam" ), ".bam.report")
def buildBAMReports( infile, outfile ):
    '''build alignment stats using bamstats

    '''
    to_cluster = USECLUSTER

    # requires a large amount of memory to run.
    # only use high-mem machines
    job_options = "-l mem_free=32G"

    # xvfb-run  -f ~/.Xauthority -a 
    track = P.snip( infile, ".bam" )

    # use a fake X display in order to avoid problems with
    # no X connection on the cluster
    xvfb_command = P.which("xvfb-run" )

    # permit multiple servers using -a option
    if xvfb_command: xvfb_command+= " -a "

    # bamstats can not accept a directory as output, hence cd to exportdir
    statement = '''
    cd %(exportdir)s/bamstats;
    %(xvfb_command)s
    %(cmd-run)s bamstats -i ../../%(infile)s -v html -o %(track)s.html 
                         --qualities --mapped --lengths --distances --starts
    >& ../../%(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
@merge( buildPicardStats, "picard_stats.load" )
def loadPicardStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

############################################################
############################################################
############################################################
@merge( mapReadsWithTophat, "tophat_stats.tsv" )
def buildTophatStats( infiles, outfile ):

    def _select( lines, pattern ):
        x = re.compile(pattern)
        for line in lines:
            r = x.search( line )
            if r: 
                g = r.groups()
                if len(g) > 1: return g
                else: return g[0]

        raise ValueError( "pattern '%s' not found %s" % (pattern, lines ))

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("track",
                            "reads_in",
                            "reads_removed", 
                            "reads_out", 
                            "junctions_loaded", 
                            "junctions_found", 
                            "possible_splices" ) ) + "\n" )
    
    for infile in infiles:
        
        track = P.snip( infile, ".bam" )
        indir = infile + ".logs" 

        fn = os.path.join( indir, "prep_reads.log" )
        lines = open( fn ).readlines()
        reads_removed, reads_in = map(int, _select( lines, "(\d+) out of (\d+) reads have been filtered out" ) )
        reads_out = reads_in - reads_removed
        prep_reads_version = _select( lines, "prep_reads (.*)$" )
        
        fn = os.path.join( indir, "reports.log" )
        lines = open( fn ).readlines()
        tophat_reports_version = _select( lines, "tophat_reports (.*)$" )
        junctions_loaded = int( _select( lines, "Loaded (\d+) junctions") )
        junctions_found = int( _select( lines, "Found (\d+) junctions from happy spliced reads") )

        fn = os.path.join( indir, "segment_juncs.log" )
        lines = open( fn ).readlines()
        if len(lines) > 0:
            segment_juncs_version =  _select( lines, "segment_juncs (.*)$" )
            possible_splices = int( _select( lines, "Reported (\d+) total possible splices") )
        else:
            segment_juncs_version = "na"
            possible_splices = ""

        # fix for paired end reads - tophat reports pairs, not reads
        if PARAMS["paired_end"]:
            reads_in *= 2
            reads_out *= 2
            reads_removed *= 2

        outf.write( "\t".join( map(str, (track,
                                         reads_in, reads_removed, reads_out, 
                                         junctions_loaded, junctions_found, possible_splices ) ) ) + "\n" )

    outf.close()

############################################################
############################################################
############################################################
@transform( buildTophatStats, suffix(".tsv"), ".load" )
def loadTophatStats( infile, outfile ):
    P.load( infile, outfile )
        
############################################################
############################################################
############################################################
@merge( buildBAMs, "mapping_stats.load" )
def loadMappingStats( infiles, outfile ):

    header = ",".join( [P.snip( x, ".bam") for x in infiles] )
    filenames = " ".join( [ "%s.tsv" % x for x in infiles ] )
    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    P.run()

############################################################
############################################################
############################################################
@transform( (mapReadsWithTophat, buildBAMs, buildMismappedBAMs), 
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER

    rna_file = os.path.join( PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_rna_gff"] )

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --filename-rna=%(rna_file)s
         --remove-rna
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
@transform( (mapReadsWithBowtieAgainstTranscriptome,),
            suffix(".bam"),
            ".readstats" )
def buildTranscriptBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    
    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''

    header = ",".join( [P.snip( x, ".readstats") for x in infiles] )
    filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    P.run()

    for suffix in ("nm", "nh"):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --table=%(tname)s 
                >> %(outfile)s
                """
    
        P.run()

############################################################
############################################################
############################################################
@transform( (mapReadsWithTophat, buildBAMs, buildMismappedBAMs),
            suffix(".bam"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"] ) ),
            ".contextstats" )
def buildContextStats( infiles, outfile ):
    '''build mapping context stats.

    Examines the genomic context to where reads align.

    A read is assigned to the genomic context that it
    overlaps by at least 50%. Thus some reads mapping
    several contexts might be dropped.
    '''

    infile, reffile = infiles

    min_overlap = 0.5

    to_cluster = USECLUSTER
    statement = '''
       python %(scriptsdir)s/rnaseq_bam_vs_bed.py
              --min-overlap=%(min_overlap)f
              --log=%(outfile)s.log
              %(infile)s %(reffile)s
       > %(outfile)s
       '''

    P.run()

############################################################
############################################################
############################################################
@follows( loadBAMStats )
@merge( buildContextStats, "context_stats.load" )
def loadContextStats( infiles, outfile ):
    """load context mapping statistics."""

    header = ",".join( [P.snip( x, ".contextstats") for x in infiles] )
    filenames = " ".join( infiles  )
    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --skip-titles
                   %(filenames)s
                | perl -p -e "s/bin/track/; s/\?/Q/g"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
                """
    P.run()
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    cc = Database.executewait( dbhandle, '''ALTER TABLE %(tablename)s ADD COLUMN mapped INTEGER''' % locals())
    statement = '''UPDATE %(tablename)s SET mapped = 
                                       (SELECT b.mapped FROM bam_stats AS b 
                                            WHERE %(tablename)s.track = b.track)''' % locals()

    cc = Database.executewait( dbhandle, statement )
    dbhandle.commit()
    
#########################################################################
#########################################################################
#########################################################################
@transform( buildBAMs, suffix(".accepted.bam"), r"\1.gtf.gz")
def buildGeneModels(infile, outfile):
    '''build transcript models for each track separately.
    '''

    to_cluster = USECLUSTER    
    job_options= "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]

    track = os.path.basename( P.snip( outfile, ".gtf.gz" ) )

    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    infile = os.path.abspath( infile )
    outfile = os.path.abspath( outfile )

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s           
              --num-threads %(cufflinks_threads)i
              --library-type %(tophat_library_type)s
              --frag-bias-correct %(cufflinks_genome_dir)s/%(genome)s.fa
              --multi-read-correct
              %(cufflinks_options)s
              %(infile)s 
    >& %(outfile)s.log;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s;
    '''

    P.run()

    # version 0.9.3
    #mv genes.expr %(outfile)s.genes.expr;
    #mv transcripts.expr %(outfile)s.transcripts.expr

    shutil.rmtree( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( buildGeneModels,
            suffix(".gtf.gz"),
            add_inputs( buildCodingGeneSet ),
            ".class.tsv.gz" )
def oldClassifyTranscripts( infiles, outfile ):
    '''classify transcripts against a reference gene set.
    '''
    
    infile, reffile = infiles

    statement = '''gunzip 
    < %(infile)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=transcript
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
		--counter=position 
		--counter=classifier 
		--section=exons 
		--counter=length 
		--counter=splice 
		--counter=splice-comparison 
		--log=%(outfile)s.log 
                --filename-format=gff
		--filename-gff=%(annotation)s 
		--genome-file=%(genome_dir)s/%(genome)s"
    | gzip
    > %(outfile)s
    '''
    P.run()


#########################################################################
#########################################################################
#########################################################################
@transform("*.bam", 
           suffix(".bam"), 
           add_inputs(buildCodingGeneSet),
           ".ref.gtf.gz")
def estimateExpressionLevelsInReference(infiles, outfile):
    '''estimate expression levels against a set of reference gene models.
    '''

    to_cluster = USECLUSTER    
    job_options= "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]

    track = os.path.basename( outfile[:-len(".gtf")] )
    
    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    bamfile, gtffile = infiles

    gtffile = os.path.abspath( gtffile )
    bamfile = os.path.abspath( bamfile )
    outfile = os.path.abspath( outfile )

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    # increase max-bundle-length to 4.5Mb due to Galnt-2 in mm9 with a 4.3Mb intron.
    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s      
              --GTF=<(gunzip < %(gtffile)s)
              --num-threads=%(cufflinks_threads)i
              --frag-bias-correct %(cufflinks_genome_dir)s/%(genome)s.fa
              --library-type %(tophat_library_type)s
              %(cufflinks_options)s
              %(bamfile)s 
    >& %(outfile)s.log;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s;
    '''

    P.run()

    shutil.rmtree( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( (estimateExpressionLevelsInReference, buildGeneModels), 
            suffix(".gtf.gz"), 
            "_gene_expression.load")
def loadExpressionLevels( infile, outfile ):
    '''load expression level measurements.'''

    track = P.snip( outfile, "_gene_expression.load" )
    P.load( infile + ".genes.expr",
            outfile = track + "_gene_expression.load",
            options = "--index=gene_id" )

    tablename = track + "_transcript_expression"
    infile2 = infile + ".transcripts.expr"

    statement = '''cat %(infile2)s
    | perl -p -e "s/trans_id/transcript_id/"
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=transcript_id 
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def runCuffCompare( infiles, outfile, reffile ):
    '''run cuffcompare.

    Will create a .tmap and .refmap input file for each input file.
    '''

    to_cluster = USECLUSTER

    tmpdir = P.getTempDir( "." )
    
    cmd_extract = "; ".join( [ "gunzip < %s > %s/%s" % (x,tmpdir,x) for x in infiles ] )

    # note: cuffcompare adds \0 bytes to gtf file - replace with '.'
    statement = '''
        %(cmd_extract)s;
        cuffcompare -o %(outfile)s
                    -s %(cufflinks_genome_dir)s/%(genome)s.fa
                    -r <( gunzip < %(reffile)s)
                    %(inf)s
        >& %(outfile)s.log;
        checkpoint;
        perl -p -e "s/\\0/./g" < %(outfile)s.combined.gtf | gzip > %(outfile)s.combined.gtf.gz;
        checkpoint;
        rm -f %(outfile)s.combined.gtf;
        checkpoint;
        gzip -f %(outfile)s.{tracking,loci};
        '''

    # the following is a hack. I was running into the same problem as described here:
    # http://seqanswers.com/forums/showthread.php?t=5809
    # the bug depended on various things, including the order of arguments
    # on the command line. Until this is resolved, simply try several times
    # with random order of command line arguments.

    for t in range( PARAMS["cufflinks_ntries"] ):
        random.shuffle( infiles )
        inf = " ".join( ["%s/%s" % (tmpdir,x) for x in infiles] )
        try:
            P.run()
            break
        except P.PipelineError, msg:
            E.warn("caught exception - trying again" )

    shutil.rmtree( tmpdir )

#########################################################################
#########################################################################
#########################################################################
@follows( buildGeneModels )
@files( [( ([ "%s.gtf.gz" % y.asFile() for y in EXPERIMENTS[x]], buildCodingGeneSet), 
           "%s.cuffcompare" % x.asFile()) 
         for x in EXPERIMENTS ] )
def compareTranscriptsPerExperiment( infiles, outfile ):
    '''compare transcript models between replicates within each experiment.'''
    infiles, reffile = infiles
    runCuffCompare( infiles, outfile, reffile )

#########################################################################
#########################################################################
#########################################################################
@merge( buildGeneModels, "%s.cuffcompare" % ALL.asFile() )
def compareTranscriptsBetweenExperiments( infiles, outfile ):
    '''compare transcript models between replicates in all experiments.'''
    # needs to be parameterized, unfortunately @merge has no add_inputs
    reffile = "%s.gtf.gz" % REFERENCE
    runCuffCompare( infiles, outfile, reffile )

#########################################################################
#########################################################################
#########################################################################
@transform( (compareTranscriptsBetweenExperiments, 
             compareTranscriptsPerExperiment),
            suffix(".cuffcompare"), 
            "_cuffcompare.load" )
def loadTranscriptComparison( infile, outfile ):
    '''load data from transcript comparison.

    This task creates two tables:

    <track>_benchmark
    <track>_loci

    The following tables are only present if there are
    multiple replicates in a sample:

    <track>_tracking 
    '''
    tracks, result = Tophat.parseTranscriptComparison( IOTools.openFile( infile ))
    tracks = [ P.snip( os.path.basename(x), ".gtf.gz" ) for x in tracks ]

    tmpfile = P.getTempFilename()
    tmpfile2 = P.getTempFilename()
    tmpfile3 = P.getTempFilename()

    outf = open( tmpfile, "w") 
    outf.write( "track\n" )
    outf.write( "\n".join( tracks ) + "\n" )
    outf.close()

    #########################################################
    ## load tracks
    #########################################################
    tablename = P.toTable(outfile) + "_tracks"

    statement = '''cat %(tmpfile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

    L.info( "loaded %s" % tablename )

    #########################################################
    ## load benchmarking data
    #########################################################
    outf = open( tmpfile, "w")
    outf.write( "track\tcontig\t%s\n" % "\t".join( Tophat.CuffCompareResult.getHeaders() ) )

    for track, vv in result.iteritems():
        track = P.snip( os.path.basename(track), ".gtf.gz" )
        for contig, v in vv.iteritems():
            if v.is_empty: continue
            outf.write( "%s\t%s\t%s\n" % (P.quote( track ), contig, str(v) ) )
    outf.close()

    tablename = P.toTable( outfile ) + "_benchmark"

    statement = '''cat %(tmpfile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track
              --index=contig
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

    L.info( "loaded %s" % tablename )

    #########################################################
    ## load tracking and transcripts information
    #########################################################
    outf = open( tmpfile, "w")
    outf.write( "%s\n" % "\t".join( ( "transfrag_id",
                                      "locus_id",
                                      "ref_gene_id",
                                      "ref_transcript_id",
                                      "code", 
                                      "nexperiments" ) ) )

    outf2 = open( tmpfile2, "w")
    outf2.write( "%s\n" % "\t".join( ( "track",
                                       "transfrag_id",
                                       "gene_id",
                                       "transcript_id",
                                       "fmi", 
                                       "fpkm", 
                                       "conf_lo", 
                                       "conf_hi", 
                                       "cov",
                                       "length" ) ) )
    outf3 = open( tmpfile3, "w" )
    outf3.write( "transfrag_id\t%s\n" % "\t".join( [ P.quote( x ) for x in tracks ] ) )

    fn = "%s.tracking.gz" % infile

    if os.path.exists( fn ):
        for transfrag in Tophat.iterate_tracking( IOTools.openFile( fn, "r") ):

            nexperiments = len( [x for x in transfrag.transcripts if x] )

            outf.write( "%s\n" % \
                            "\t".join( (transfrag.transfrag_id, 
                                        transfrag.locus_id, 
                                        transfrag.ref_gene_id,
                                        transfrag.ref_transcript_id,
                                        transfrag.code,
                                        str(nexperiments))))

            outf3.write( "%s" % transfrag.transfrag_id )

            for track, t in zip(tracks, transfrag.transcripts):
                if t:
                    outf2.write("%s\n" % "\t".join( map(str, (track,
                                                              transfrag.transfrag_id ) + t ) ) )

                    outf3.write( "\t%f" % t.fpkm )
                else:
                    outf3.write( "\t" )

            outf3.write( "\n" )
    else:
        E.warn( "no tracking file %s - skipped " )
            
    outf.close()
    outf2.close()
    outf3.close()

    tablename = P.toTable( outfile ) + "_tracking"
    statement = '''cat %(tmpfile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=locus_id
              --index=transfrag_id
              --index=code
              --table=%(tablename)s 
    >> %(outfile)s
    '''

    P.run()
    L.info( "loaded %s" % tablename )

    tablename = P.toTable( outfile ) + "_transcripts"
    statement = '''cat %(tmpfile2)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=transfrag_id
              --index=ref_gene_id
              --index=ref_transcript_id
              --index=transcript_id
              --index=gene_id
              --index=track
              --table=%(tablename)s 
    >> %(outfile)s
    '''

    P.run()
    L.info( "loaded %s" % tablename )

    tablename = P.toTable( outfile ) + "_fpkm"
    statement = '''cat %(tmpfile3)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=transfrag_id
              --table=%(tablename)s 
    >> %(outfile)s
    '''

    P.run()
    L.info( "loaded %s" % tablename )

    #########################################################
    ## load locus information
    #########################################################
    outf = open( tmpfile, "w")
    outf.write( "%s\n" % "\t".join( ( "locus_id",
                                      "contig",
                                      "strand",
                                      "start",
                                      "end", 
                                      "nexperiments", ) + tuple(tracks) ) )

    for locus in Tophat.iterate_locus( IOTools.openFile( "%s.loci.gz" % infile, "r") ):

        counts = [ len(x) for x in locus.transcripts ] 
        nexperiments = len( [x for x in counts if x > 0] )

        outf.write( "%s\t%s\t%s\t%i\t%i\t%i\t%s\n" % \
                        (locus.locus_id, locus.contig, locus.strand, 
                         locus.start, locus.end,
                         nexperiments,
                         "\t".join( map( str, counts) ) ) )
    outf.close()

    tablename = P.toTable( outfile ) + "_loci"

    statement = '''cat %(tmpfile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=locus_id
              --table=%(tablename)s 
    >> %(outfile)s
    '''

    P.run()
    L.info( "loaded %s" % tablename )

    os.unlink( tmpfile )
    os.unlink( tmpfile2 )
    os.unlink( tmpfile3 )

#########################################################################
#########################################################################
#########################################################################
@transform( compareTranscriptsBetweenExperiments, 
            suffix(".cuffcompare"),
            ".gtf.gz" )
def buildAbinitioGeneSet( infile, outfile ):
    '''builds ab-initio gene set.
    
    The ab-initio gene set is derived from the cuffcompare result.

    The following transfrags are removed at this stage:

        * transfrags overlapping RNA genes
        * transfrags on certain contigs (usually: mitochondrial genes)
        
    '''
    infile += ".combined.gtf.gz"
    writePrunedGTF( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows(loadTranscriptComparison)
@merge( (buildAbinitioGeneSet, buildReferenceGeneSet),
        "abinitio.gtf.gz" )
def buildFullGeneSet( infiles, outfile ):
    '''builds a gene set by merging the ab-initio gene set and
    the reference gene set.
    
    The gene set is cleaned in order to permit differential expression
    analysis.
    
    Only transfrags are kept that are:

    1. observed in at least 2 samples to remove partial transfrags that
        are the result of low coverage observations in one sample
    
    see also: http://seqanswers.com/forums/showthread.php?t=3967
    
    Will also build removed.gtf.gz of removed transcripts.
    '''
    abinitio_gtf, reference_gtf = infiles
    keep_gtf = outfile
    remove_gtf = "removed.gtf.gz"

    tablename = P.quote( P.snip( abinitio_gtf, ".gtf.gz") + "_cuffcompare_tracking" )
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()    

    statement = '''SELECT transfrag_id FROM %(tablename)s WHERE nexperiments > 1''' % locals()

    keep = set( [ x[0] for x in cc.execute(statement).fetchall()] )
    
    E.info( "keeping %i transfrags" % len(keep) )

    inf = GTF.iterator( IOTools.openFile( abinitio_gtf ) )
    outf1 = IOTools.openFile( keep_gtf, "w" )
    outf2 = IOTools.openFile( remove_gtf, "w" )
    
    c = E.Counter()
    for gtf in inf:
        c.input += 1
        if gtf.transcript_id in keep:
            c.kept += 1
            outf1.write( "%s\n" % str(gtf ) )
        else:
            c.removed += 1
            outf2.write( "%s\n" % str(gtf ) )
        
    outf1.close()
    outf2.close()

    E.info("%s" % str(c))


#########################################################################
#########################################################################
#########################################################################
@merge( (buildAbinitioGeneSet, buildReferenceGeneSet,
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_repeats_gff"] )),
        "novel.gtf.gz" )
def buildNovelGeneSet( infiles, outfile ):
    '''build a gene set of novel genes by merging the ab-initio gene set and
    the reference gene set.
    
    Ab-initio transcripts are removed based on features in the reference gene set.
    Removal is aggressive  - as soon as one transcript of a 
    gene/locus overlaps, all transcripts of that gene/locus are gone.

    Transcripts that lie exclusively in repetetive sequence are removed, too.

    The resultant set contains a number of novel transcripts. However, these
    transcripts will still overlap some known genomic features like pseudogenes.

    '''

    abinitio_gtf, reference_gtf, repeats_gff = infiles
    
    E.info( "indexing geneset for filtering" )

    sections = ("protein_coding", "lincRNA", "processed_transcript")

    indices = {}
    for section in sections:
        indices[section] = GTF.readAndIndex( 
            GTF.iterator_filtered( GTF.iterator( IOTools.openFile( reference_gtf) ),
                                   source = section ),
            with_value = False ) 
        
    E.info( "build indices for %i features" % len(indices))

    repeats = GFF.readAndIndex( GFF.iterator( IOTools.openFile( repeats_gff) ),
                                with_value = False )

    E.info( "build index for repeats" )

    total_genes, remove_genes = set(), collections.defaultdict( set )
    inf = GTF.iterator( IOTools.openFile( abinitio_gtf ) )
    for gtf in inf:
        total_genes.add( gtf.gene_id )
        for section in sections:
            if indices[section].contains( gtf.contig, gtf.start, gtf.end):
                remove_genes[gtf.gene_id].add( section )

        try:
            for r in repeats.get( gtf.contig, gtf.start, gtf.end ):
                if r[0] <= gtf.start and r[1] >= gtf.end:
                    remove_genes[gtf.gene_id].add( "repeat" )
                    break
        except KeyError:
            pass

    E.info( "removing %i out of %i genes" % (len(remove_genes), len(total_genes)) )

    PipelineRnaseq.filterAndMergeGTF( abinitio_gtf, outfile, remove_genes, merge = True )

#########################################################################
#########################################################################
#########################################################################
@merge( (buildAbinitioGeneSet, buildReferenceGeneSet,
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_repeats_gff"]),
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_pseudogenes_gtf"] ),
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_numts_gtf"] ),
         ), "lincrna.gtf.gz" )
def buildLincRNAGeneSet( infiles, outfile ):
    '''build lincRNA gene set. 
    
    The lincRNA gene set contains all known lincRNA transcripts from
    the reference gene set plus all transcripts in the novel set that
    do not overlap at any protein coding, processed or pseudogene transcripts 
    (exons+introns) in the reference gene set.

    Transcripts that lie exclusively in repetetive sequence are removed, too.

    lincRNA genes are often expressed at low level and thus the resultant transcript
    models are fragmentory. To avoid some double counting in downstream analyses, 
    transcripts overlapping on the same strand are merged.
    
    Transcripts need to have a length of at least 200 bp.

    '''
    
    infile_abinitio, reference_gtf, repeats_gff, pseudogenes_gtf, numts_gtf = infiles

    E.info( "indexing geneset for filtering" )

    input_sections = ("protein_coding", 
                      "lincRNA", 
                      "processed_transcript" )

    indices = {}
    for section in input_sections:
        indices[section] = GTF.readAndIndex( 
            GTF.iterator_filtered( GTF.merged_gene_iterator( GTF.iterator( IOTools.openFile( reference_gtf) )),
                                   source = section ),
            with_value = False )
        
    E.info( "built indices for %i features" % len(indices))

    indices["repeats"] = GFF.readAndIndex( GFF.iterator( IOTools.openFile( repeats_gff) ), with_value = False )
    
    E.info( "added index for repeats" )

    indices["pseudogenes"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( pseudogenes_gtf) ), with_value = False )

    E.info( "added index for pseudogenes" )

    indices["numts"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( numts_gtf) ), with_value = False )

    E.info( "added index for numts" )

    sections = indices.keys()

    total_genes, remove_genes = set(), collections.defaultdict( set )
    inf = GTF.iterator( IOTools.openFile( infile_abinitio ) )

    E.info( "collecting genes to remove" )

    min_length = int(PARAMS["lincrna_min_length"])

    for gtfs in GTF.transcript_iterator( inf ):
        gene_id = gtfs[0].gene_id 
        total_genes.add( gene_id )

        l = sum( [ x.end - x.start for x in gtfs ] )

        if l < min_length:
            remove_genes[gene_id].add( "length" )
            continue

        for section in sections:
            for gtf in gtfs:
                if indices[section].contains( gtf.contig, gtf.start, gtf.end):
                    remove_genes[gene_id].add( section )
        
    E.info( "removing %i out of %i genes" % (len(remove_genes), len(total_genes)) )

    PipelineRnaseq.filterAndMergeGTF( infile_abinitio, outfile, remove_genes, merge = True )

    E.info( "adding known lincRNA set" )

    # add the known lincRNA gene set.
    statement = '''zcat %(reference_gtf)s
    | awk '$2 == "lincRNA"' 
    | gzip 
    >> %(outfile)s
    '''
    P.run()

    # sort
    statement = '''
    mv %(outfile)s %(outfile)s.tmp;
    checkpoint;
    zcat %(outfile)s.tmp
    | python %(scriptsdir)s/gtf2gtf.py --sort=contig+gene --log=%(outfile)s.log
    | gzip > %(outfile)s;
    checkpoint;
    rm -f %(outfile)s.tmp
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( (buildLincRNAGeneSet, buildReferenceTranscriptome ), "lincrna.pseudos.tsv" )
def annotateLincRNA( infiles, outfile ):
    '''align lincRNA against reference transcriptome
    in order to spot pseudogenes.
    '''
    
    linc_fasta, reference_fasta = infiles

    format = ("qi", "qS", "qab", "qae", 
              "ti", "tS", "tab", "tae", 
              "s",
              "pi", 
              "C")
    
    format = "\\\\t".join( ["%%%s" % x for x in format] )

    statement = '''
    zcat %(linc_fasta)s
    | python %(scriptsdir)s/gff2fasta.py 
              --is-gtf 
              --genome=%(genome_dir)s/%(genome)s 
              --log=%(outfile)s.log
    | %(cmd-farm)s --split-at-regex=\"^>(\S+)\" --chunksize=400 --log=%(outfile)s.log
    "exonerate --target %%STDIN%%
              --query %(reference_fasta)s
              --model affine:local
              --score %(lincrna_min_exonerate_score)i
              --showalignment no --showsugar no --showcigar no 
              --showvulgar no
              --bestn 5
              --ryo \\"%(format)s\\n\\" 
    " 
    | grep -v -e "exonerate" -e "Hostname"
    | gzip > %(outfile)s.links.gz
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (buildLincRNAGeneSet,
             buildNovelGeneSet), 
            suffix(".gtf.gz"), 
            "_build_summary.load" )
def loadGeneSetsBuildInformation( infile, outfile ):
    '''load results from gene set filtering into database.'''
    infile += ".summary.tsv.gz" 
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (buildCodingGeneSet,
             buildNoncodingGeneSet,
             buildGeneModels, 
             buildFullGeneSet, 
             buildLincRNAGeneSet,
             buildNovelGeneSet),
            suffix(".gtf.gz"), 
            add_inputs( buildReferenceGeneSetWithCDS ),
            ".class.tsv.gz" )
def classifyTranscripts( infiles, outfile ):
    '''classify transcripts.
    '''
    to_cluster = USECLUSTER
    
    infile, reference = infiles

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gtf2table.py
           --counter=classifier-rnaseq 
           --reporter=transcripts
           --filename-gff=%(reference)s
           --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


## need to change pipeline logic to avoid this duplication
@transform( (compareTranscriptsPerExperiment,
             compareTranscriptsBetweenExperiments ), 
            suffix(".cuffcompare"), 
            add_inputs( buildReferenceGeneSetWithCDS ),
            ".class.tsv.gz" )
def classifyTranscriptsCuffcompare( infiles, outfile ):
    '''classify transcripts.
    '''
    to_cluster = USECLUSTER
    
    infile, reference = infiles

    statement = '''
    zcat %(infile)s.combined.gtf.gz
    | python %(scriptsdir)s/gtf2table.py
           --counter=classifier-rnaseq 
           --reporter=transcripts
           --filename-gff=%(reference)s
           --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (classifyTranscripts, classifyTranscriptsCuffcompare), suffix(".tsv.gz"), ".load" )
def loadClassification( infile, outfile ):
    P.load( infile, outfile, options = "--index=transcript_id --index=match_gene_id --index=match_transcript_id --index=source" )

#########################################################################
#########################################################################
#########################################################################
@merge( (buildGeneModels, 
         buildAbinitioGeneSet, 
         compareTranscriptsPerExperiment,
         compareTranscriptsBetweenExperiments,
         buildFullGeneSet,
         buildReferenceGeneSet,
         buildCodingGeneSet,
         buildNoncodingGeneSet,
         buildLincRNAGeneSet,
         buildNovelGeneSet),
        "geneset_stats.tsv" )
def buildGeneSetStats( infiles, outfile ):
    '''compile gene set statistics.
    '''

    to_cluster = USECLUSTER

    cuffcompare = [ x + ".combined.gtf.gz" for x in infiles if x.endswith("cuffcompare")]
    other = [ x for x in infiles if x.endswith(".gtf.gz")]

    if os.path.exists("removed.gtf.gz"):
        other.append( "removed.gtf.gz" )

    allfiles = " ".join( other + cuffcompare )

    statement = '''
    python %(scriptsdir)s/gff2stats.py --is-gtf
    %(allfiles)s --log=%(outfile)s.log
    | perl -p -e "s/.gtf.gz//"
    | perl -p -e "s/^agg.*cuffcompare.combined/unfiltered/"
    | perl -p -e "s/.cuffcompare.combined//"
    > %(outfile)s
    '''

    P.run()
    
#########################################################################
#########################################################################
#########################################################################
@transform( buildGeneSetStats, suffix(".tsv"), ".load" )
def loadGeneSetStats( infile, outfile ):
    '''load geneset statisticts.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (
        buildReferenceGeneSet,
        buildCodingGeneSet,
        buildAbinitioGeneSet,
        buildFullGeneSet,
        buildNoncodingGeneSet,
        buildLincRNAGeneSet,
        buildNovelGeneSet),
            suffix(".gtf.gz"),
            ".mappability.gz" )
def annotateTranscriptsMappability( infile, outfile ):
    '''classify transcripts with respect to the gene set.
    '''
    # script will be farmed out
    to_cluster = False

    if "geneset_mappability" not in PARAMS or not PARAMS["geneset_mappability"]:
        P.touch(outfile)
        return

    statement = """
    zcat < %(infile)s 
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
    "python %(scriptsdir)s/gtf2table.py 
                --reporter=transcripts
                --counter=bigwig-counts
                --bigwig-file=%(geneset_mappability)s
		--log=%(outfile)s.log"
    | gzip
    > %(outfile)s"""
    
    P.run()

############################################################
@transform( annotateTranscriptsMappability, suffix(".mappability.gz"), "_mappability.load" )
def loadTranscriptsMappability( infile, outfile ):
    '''load interval annotations: genome architecture
    '''
    if "geneset_mappability" not in PARAMS or not PARAMS["geneset_mappability"]:
        P.touch(outfile)
        return

    P.load( infile, outfile, "--index=transcript_id --allow-empty" )

#########################################################################
#########################################################################
#########################################################################
@transform( (buildFullGeneSet,
             buildNovelGeneSet),
            suffix(".gtf.gz"),
            ".annotations.gz" )
def annotateTranscripts( infile, outfile ):
    '''classify transcripts with respect to the gene set.
    '''
    to_cluster = USECLUSTER

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_annotation"] )

    statement = """
    zcat < %(infile)s 
    | python %(scriptsdir)s/gtf2table.py 
                --reporter=transcripts
		--counter=position 
		--counter=classifier
		--section=exons 
		--counter=length 
		--log=%(outfile)s.log 
		--filename-gff=%(annotation_file)s 
		--genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s"""
    
    P.run()

############################################################
@transform( annotateTranscripts, suffix(".annotations"), "_annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
def hasReplicates( track ):
    '''indicator function - return true if track has replicates'''
    replicates = PipelineTracks.getSamplesInTrack( track, TRACKS )
    return len(replicates) > 1

@follows( loadTranscriptComparison, mkdir( os.path.join( PARAMS["exportdir"], "cuffcompare" ) ) )
@files( [ ("%s.cuffcompare" % x.asFile(), "%s.reproducibility" % x.asFile() )
          for x in EXPERIMENTS if hasReplicates( x )] )
def buildReproducibility( infile, outfile ):
    '''all-vs-all comparison between samples.

    Compute correlation between expressed transfrags. Transfrags missing
    from another set are ignored.
    '''

    track = TRACKS.factory( filename = outfile[:-len(".reproducibility")] )
    
    replicates = PipelineTracks.getSamplesInTrack( track, TRACKS )

    dbhandle = sqlite3.connect( PARAMS["database"] )

    tablename = "%s_cuffcompare_fpkm" % track.asTable()
    tablename2 = "%s_cuffcompare_tracking" % track.asTable()

    ##################################################################
    ##################################################################
    ##################################################################
    ## build table correlating expression values
    ##################################################################
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "track1\ttrack2\tcode\tpairs\tnull1\tnull2\tboth_null\tnot_null\tone_null\t%s\n" % "\t".join( Stats.CorrelationTest.getHeaders() ) )

    for rep1, rep2 in itertools.combinations( replicates, 2 ):
        
        track1, track2 = rep1.asTable(), rep2.asTable()

        def _write( statement, code ):
            data = Database.executewait( dbhandle, statement ).fetchall()
            if len(data) == 0: return
            both_null = len( [ x for x in data if x[0] == 0 and x[1] == 0 ] )
            one_null = len( [ x for x in data if x[0] == 0 or x[1] == 0 ] )
            null1 = len( [ x for x in data if x[0] == 0 ] )
            null2 = len( [ x for x in data if x[1] == 0 ] )
            not_null = [ x for x in data if x[0] != 0 and x[1] != 0 ]
            if len(not_null) > 1:
                x,y = zip( *not_null )
                result = Stats.doCorrelationTest( x, y)
            else:
                result = Stats.CorrelationTest()

            outf.write( "%s\n" % "\t".join( map(str, (track1, track2, code, 
                                                      len(data),
                                                      null1, null2, both_null, 
                                                      len(not_null), 
                                                      one_null,
                                                      str(result) ) ) ) )

        for code in PARAMS["reproducibility_codes"]:
            statement = '''SELECT CASE WHEN %(track1)s THEN %(track1)s ELSE 0 END, 
                                  CASE WHEN %(track2)s THEN %(track2)s ELSE 0 END
                       FROM %(tablename)s AS a,
                            %(tablename2)s AS b
                       WHERE a.transfrag_id = b.transfrag_id AND
                             b.code = '%(code)s'
                    '''
            
            _write( statement % locals(), code )

        statement = '''SELECT CASE WHEN %(track1)s THEN %(track1)s ELSE 0 END, 
                                  CASE WHEN %(track2)s THEN %(track2)s ELSE 0 END
                       FROM %(tablename)s AS a
                    '''
        _write( statement % locals(), "*" )


    ##################################################################
    ##################################################################
    ##################################################################
    ## plot pairwise correlations
    ##################################################################
    # plot limit
    lim = 1000

    outdir = os.path.join( PARAMS["exportdir"], "cuffcompare" )
        
    R('''library(RSQLite)''')
    R('''drv = dbDriver( "SQLite" )''' )
    R('''con <- dbConnect(drv, dbname = 'csvdb')''')
    columns = ",".join( [ x.asTable() for x in replicates ] )
    data = R('''data = dbGetQuery(con, "SELECT %(columns)s FROM %(tablename)s")''' % locals())
    R.png( "%(outdir)s/%(outfile)s.pairs.png" % locals())
    R('''pairs( data, pch = '.', xlim=c(0,%(lim)i), ylim=c(0,%(lim)i) )''' % locals())
    R('''dev.off()''')

    for rep1, rep2 in itertools.combinations( replicates, 2 ):
        a,b = rep1.asTable(), rep2.asTable()
        r = R('''r = lm( %(a)s ~ %(b)s, data)''' % locals() )
        R.png( "%(outdir)s/%(outfile)s.pair.%(rep1)s_vs_%(rep2)s.png" % locals())
        R('''plot(data$%(a)s, data$%(b)s, pch='.', xlim=c(0,%(lim)i), ylim=c(0,%(lim)i),)''' % locals() )

        try: R('''abline(r)''')
        except RRuntimeError: pass

        R('''dev.off()''')

#########################################################################
#########################################################################
#########################################################################
@transform( buildReproducibility, suffix(".reproducibility"), "_reproducibility.load" )
def loadReproducibility( infile, outfile ):
    '''load reproducibility results.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
# @files( [ ( ([ "%s.bam" % xx.asFile() for xx in EXPERIMENTS[x] ], 
#              [ "%s.bam" % yy.asFile() for yy in EXPERIMENTS[y] ]),
#             "%s_vs_%s.cuffdiff" % (x.asFile(),y.asFile()) )
#           for x,y in itertools.combinations( EXPERIMENTS, 2) ] )
# def estimateDifferentialExpressionPairwise( infiles, outfile ):
#     '''estimate differential expression using cuffdiff.

#     Replicates are grouped.
#     '''
    
#     to_cluster = USECLUSTER
#     job_options= "-pe dedicated %i -R y" % PARAMS["cuffdiff_threads"]

#     reffile = "reference.gtf.gz"        

#     outdir = outfile + ".dir" 
#     try: os.mkdir( outdir )
#     except OSError: pass

#     reps = "%s    %s" % (",".join( infiles[0]),
#                          ",".join( infiles[1]) )
    
#     statement = '''
#     cuffdiff -o %(outdir)s
#              --verbose
#              -r %(cufflinks_genome_dir)s/%(genome)s.fa
#              --num-threads %(cuffdiff_threads)i
#              <(gunzip < %(reffile)s)
#              %(reps)s
#     >& %(outfile)s
#     '''
#     P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (buildFullGeneSet, 
             buildReferenceGeneSet,
             buildCodingGeneSet,
             buildLincRNAGeneSet,
             buildNoncodingGeneSet,
             buildNovelGeneSet),
            suffix(".gtf.gz"),
            "_geneinfo.load" )
def loadGeneSetGeneInformation( infile, outfile ):
    PipelineGeneset.loadGeneStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (buildFullGeneSet, 
             buildReferenceGeneSet,
             buildCodingGeneSet,
             buildLincRNAGeneSet,
             buildNoncodingGeneSet,
             buildNovelGeneSet),
            suffix(".gtf.gz"),
            "_transcript2gene.load" )
def loadGeneInformation( infile, outfile ):
    PipelineGeneset.loadTranscript2Gene( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (
        buildGeneModels,
        buildFullGeneSet, 
        buildReferenceGeneSet,
        buildCodingGeneSet,
        buildNoncodingGeneSet,
        buildLincRNAGeneSet,
        buildNovelGeneSet),
            suffix(".gtf.gz"),
            "_transcriptinfo.load" )
def loadGeneSetTranscriptInformation( infile, outfile ):
    PipelineGeneset.loadTranscriptStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (buildFullGeneSet, 
             buildReferenceGeneSet,
             buildCodingGeneSet,
             buildNoncodingGeneSet,
             buildLincRNAGeneSet,
             buildNovelGeneSet),
            suffix(".gtf.gz"),
            ".cuffdiff" )
def runCuffdiff( infile, outfile ):
    '''estimate differential expression using cuffdiff.

    Replicates are grouped.
    '''

    to_cluster = USECLUSTER

    outdir = outfile + ".dir" 
    try: os.mkdir( outdir )
    except OSError: pass

    job_options= "-pe dedicated %i -R y" % PARAMS["cuffdiff_threads"]

    # replicates are separated by ","
    reps, labels = [], []
    for group, replicates in EXPERIMENTS.iteritems():
        reps.append( ",".join( [ "%s.accepted.bam" % r.asFile() for r in replicates] ) )
        labels.append( group.asFile() )

    reps = "   ".join( reps )
    labels = ",".join( labels )

    statement = '''date > %(outfile)s; hostname >> %(outfile)s;
    cuffdiff --output-dir %(outdir)s
             --library-type %(tophat_library_type)s
             %(cuffdiff_options)s
             --verbose
             --num-threads %(cuffdiff_threads)i
             --labels %(labels)s
             --FDR %(cuffdiff_fdr)f
             <(gunzip < %(infile)s )
             %(reps)s
    >> %(outfile)s 2>&1;
    date >> %(outfile)s;
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( runCuffdiff, 
            suffix(".cuffdiff"), 
            "_cuffdiff.load" )
def loadCuffdiff( infile, outfile ):
    '''load results from differential expression analysis and produce
    summary plots.

    Note: converts from ln(fold change) to log2 fold change.
   
    The cuffdiff output is parsed. 

    Pairwise comparisons in which one gene is not expressed (fpkm < fpkm_silent)
    are set to status 'NOCALL'. These transcripts might nevertheless be significant.
    '''

    prefix = P.toTable( outfile )
    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    to_cluster = False
    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    # ignore promoters and splicing - no fold change column, but  sqrt(JS)
    for fn, level in ( ("cds_exp.diff", "cds"),
                       ("gene_exp.diff", "gene"),
                       ("isoform_exp.diff", "isoform"),
                       # ("promoters.diff", "promotor"),
                       # ("splicing.diff", "splice"), 
                       ("tss_group_exp.diff", "tss") ):
        
        tablename = prefix + "_" + level + "_diff"

        # max/minimum fold change seems to be (-)1.79769e+308
        # ln to log2: multiply by log2(e)
        statement = '''cat %(indir)s/%(fn)s
        | perl -p -e "s/sample_/track/g; s/value_/value/g; s/yes$/1/; s/no$/0/; s/ln\\(fold_change\\)/lfold/; s/p_value/pvalue/"
        | awk -v OFS='\\t' '/test_id/ {print;next;} 
                              { $10 = $10 * 1.44269504089; 
                                if( $6 == "OK" && ($7 < %(cuffdiff_fpkm_expressed)f || $8 < %(cuffdiff_fpkm_expressed)f )) { $6 = "NOCALL"; };
                                print; } '
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track1
              --index=track2
              --index=test_id
              --table=%(tablename)s 
         >> %(outfile)s.log
         '''
        
        P.run()

    for fn, level in ( ("cds.fpkm_tracking", "cds" ),
                       ("genes.fpkm_tracking", "gene"),
                       ("isoforms.fpkm_tracking", "isoform"),
                       ("tss_groups.fpkm_tracking", "tss") ):

        tablename = prefix + "_" + level + "_levels" 

        statement = '''cat %(indir)s/%(fn)s
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=tracking_id
              --table=%(tablename)s 
         >> %(outfile)s.log
         '''
        
        P.run()

    ## build convenience table with tracks
    tablename = prefix + "_isoform_levels"
    tracks = Database.getColumnNames( dbhandle, tablename )
    tracks = [ x[:-len("_FPKM")] for x in tracks if x.endswith("_FPKM") ]
    
    tmpfile = P.getTempFile()
    tmpfile.write( "track\n" )
    tmpfile.write("\n".join(tracks) + "\n" )
    tmpfile.close()
    
    statement = P.load( tmpfile.name, outfile )
    os.unlink( tmpfile.name )
    
#########################################################################
#########################################################################
#########################################################################
def buildExpressionStats( tables, method, outfile ):
    '''build expression summary statistics.
    
    Creates some diagnostic plots in

    <exportdir>/<method> directory.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    def togeneset( tablename ):
        return re.match("([^_]+)_", tablename ).groups()[0]

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("geneset", "level", "track1", "track2", "tested",
                            "\t".join( [ "status_%s" % x for x in keys_status ] ),
                            "significant",
                            "twofold" ) ) + "\n" )

    all_tables = set(Database.getTables( dbhandle ))
    outdir = os.path.join( PARAMS["exportdir"], method )

    for level in CUFFDIFF_LEVELS:

        for tablename in tables:

            tablename_diff = "%s_%s_diff" % (tablename, level)
            tablename_levels = "%s_%s_diff" % (tablename, level )
            geneset = togeneset( tablename_diff )
            if tablename_diff not in all_tables: continue

            def toDict( vals, l = 2 ):
                return collections.defaultdict( int, [ (tuple( x[:l]), x[l]) for x in vals ] )
            
            tested = toDict( Database.executewait( dbhandle,
                                               """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            status = toDict( Database.executewait( dbhandle,
                                                   """SELECT track1, track2, status, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2,status""" % locals() ).fetchall(), 3 )
            signif = toDict( Database.executewait( dbhandle,
                                                   """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE significant
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            fold2 = toDict( Database.executewait( dbhandle,
                    """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE (lfold >= 1 or lfold <= -1) AND significant
                                    GROUP BY track1,track2,significant""" % locals() ).fetchall() )
            
            for track1, track2 in itertools.combinations( EXPERIMENTS, 2 ):
                outf.write( "\t".join(map(str, (
                                geneset,
                                level,
                                track1,
                                track2,
                                tested[(track1,track2)],
                                "\t".join( [ str(status[(track1,track2,x)]) for x in keys_status]),
                                signif[(track1,track2)],
                                fold2[(track1,track2)] ) ) ) + "\n" )
                
            ###########################################
            ###########################################
            ###########################################
            # plot length versus P-Value
            data = Database.executewait( dbhandle, 
                                         '''SELECT i.sum, pvalue 
                                 FROM %(tablename_diff)s, 
                                 %(geneset)s_geneinfo as i 
                                 WHERE i.gene_id = test_id AND significant'''% locals() ).fetchall()

            # require at least 10 datapoints - otherwise smooth scatter fails
            if len(data) > 10:
                data = zip(*data)

                pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_pvalue_vs_length.png" % locals()
                R.png( pngfile )
                R.smoothScatter( R.log10( ro.FloatVector(data[0]) ),
                                 R.log10( ro.FloatVector(data[1]) ),
                                 xlab = 'log10( length )',
                                 ylab = 'log10( pvalue )',
                                 log="x", pch=20, cex=.1 )

                R['dev.off']()

    outf.close()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "cuffdiff" ) ) )
@transform( loadCuffdiff,
            suffix(".load"), 
            ".plots" )
def buildCuffdiffPlots( infile, outfile ):
    '''create summaries of cufflinks results (including some diagnostic plots)

    Plots are created in the <exportdir>/cuffdiff directory.

    Plots are:

    <geneset>_<method>_<level>_<track1>_vs_<track2>_significance.png
        fold change against expression level
    '''
    ###########################################
    ###########################################
    ## create diagnostic plots
    ###########################################
    outdir = os.path.join( PARAMS["exportdir"], "cuffdiff" )
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    prefix = P.snip( infile, ".load" )

    geneset, method = prefix.split("_")
    
    for level in CUFFDIFF_LEVELS:
        tablename_diff = prefix + "_%s_diff" % level
        tablename_levels = prefix + "_%s_levels" % level
        
        # note that the ordering of EXPERIMENTS and the _diff table needs to be the same
        # as only one triangle is stored of the pairwise results.
        # do not plot "undefined" lfold values (where value1 or value2 = 0)
        # do not plot lfold values where the confidence bounds contain 0.
        for track1, track2 in itertools.combinations( EXPERIMENTS, 2 ):
            statement = """
                        SELECT CASE WHEN d.value1 < d.value2 THEN d.value1 ELSE d.value2 END, d.lfold, d.significant
                        FROM %(tablename_diff)s AS d
                        WHERE track1 = '%(track1)s' AND 
                              track2 = '%(track2)s' AND 
                              status = 'OK' AND
                              value1 > 0 AND 
                              value2 > 0 
                        """ % locals()
            
            data = zip( *Database.executewait( dbhandle, statement ))
            
            pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_significance.png" % locals()
            #ian: Bug fix: moved R.png to after data check so that no plot is started if there is no data
            #     this was leading to R falling over from too many open devices

            if len(data) == 0:
                E.warn( "no plot for %s - %s -%s vs %s" % ( pngfile, level, track1, track2))
                continue

            R.png( pngfile )
            R.plot( ro.FloatVector(data[0]), 
                    ro.FloatVector(data[1]), 
                    xlab = 'min(FPKM)',
                    ylab = 'log2fold',
                    log="x", pch=20, cex=.1,
                    col = R.ifelse( ro.IntVector(data[2]), "red", "black" ) )

            R['dev.off']()

    P.touch( outfile )
        
#########################################################################
#########################################################################
#########################################################################
@merge( loadCuffdiff,
        "cuffdiff_stats.tsv" )
def buildCuffdiffStats( infiles, outfile ):
    tablenames = [P.toTable(x) for x in infiles ]
    buildExpressionStats( tablenames, "cuffdiff", outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildCuffdiffStats,
            suffix(".tsv"), 
            ".load" )
def loadCuffdiffStats( infile, outfile ):
    '''import cuffdiff results.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
def getLibrarySizes( infiles ):
    
    vals = []
    for infile in infiles:
        assert infile.endswith( ".readstats")
        val, cont = [ x[:-1].split("\t") for x in open(infile).readlines() if re.search( "\tmapped", x ) ][0]
        vals.append(int(val))
        
    return vals

#########################################################################
#########################################################################
#########################################################################
@merge( loadExpressionLevels,
        "genelevel_fpkm_tagcounts.tsv.gz")
def buildFPKMGeneLevelTagCounts( infiles, outfile ):
    '''build tag counts using normalized counts from tophat.

    These are gene-length normalized count levels.

    They are obtained by multiplying the FPKM value
    by the median library size.
    '''
    infiles = [ x for x in infiles if x.endswith( ".ref_gene_expression.load" ) ]

    tracks = [ P.snip( x,".ref_gene_expression.load" ) for x in infiles ]

    # get normalization values
    library_sizes = getLibrarySizes( [ "%s.readstats" % x for x in tracks ] )
    if len(library_sizes) == 0: raise ValueError("could not get library sizes" )

    median_library_size = numpy.median( library_sizes )

    # dbhandle = sqlite3.connect( os.path.join( PARAMS["annotations_dir"],
    #                                           PARAMS_ANNOTATIONS["interface_database"] ) )
    # cc = dbhandle.cursor()    
    # median_gene_length = numpy.median( [ x for x in cc.execute( "SELECT sum FROM gene_stats") ] )

    scale = median_library_size / 1000000.0

    L.info( "normalization: median library size=%i, factor=1.0 / %f" % \
                (median_library_size, scale) )

    # normalize
    results = []
    dbhandle = sqlite3.connect( PARAMS["database"] )

    for track in tracks:
        table = "%s_ref_gene_expression" % P.quote(track)
        statement = "SELECT gene_id, FPKM / %(scale)f FROM %(table)s" % locals()
        results.append( dict( Database.executewait( dbhandle, statement ).fetchall() ) )
    
    outf = IOTools.openFile( outfile, "w" )
    gene_ids = set()
    for x in results: gene_ids.update( x.keys() )

    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    for gene_id in gene_ids:
        outf.write( "%s\t%s\n" % ( gene_id, "\t".join( [str(int(x[gene_id])) for x in results ] ) ) )
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "coding_exons.gtf.gz" )
def buildCodingExons( infile, outfile ):
    '''compile set of protein coding exons.

    This set is used for splice-site validation
    '''

    to_cluster = True
    statement = '''
    zcat %(infile)s 
    | awk '$2 == "protein_coding" && $3 == "CDS"'
    | perl -p -e "s/CDS/exon/" 
    | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log 
    | gzip 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
@transform( buildBAMs,
            suffix(".bam"),
            add_inputs( buildCodingExons ),
            ".exon.validation.tsv.gz" )
def buildExonValidation( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    infile, exons = infiles
    statement = '''cat %(infile)s
    | python %(scriptsdir)s/rnaseq_bam_vs_exons.py
         --filename-exons=%(exons)s
         --force
         --log=%(outfile)s.log
         --output-filename-pattern="%(outfile)s.%%s.gz"
    | gzip
    > %(outfile)s
    '''

    P.run()


############################################################
############################################################
############################################################
@merge( buildExonValidation, "exon_validation.load" )
def loadExonValidation( infiles, outfile ):
    '''merge alignment stats into single tables.'''
    suffix = suffix = ".exon.validation.tsv.gz" 
    P.mergeAndLoad( infiles, outfile, suffix = suffix )
    for infile in infiles:
        track = P.snip( infile, suffix )
        o = "%s_overrun.load" % track 
        P.load( infile + ".overrun.gz", o )

#########################################################################
#########################################################################
#########################################################################
@transform( (buildReferenceGeneSet, 
             buildCodingGeneSet,
             buildNovelGeneSet,
             buildLincRNAGeneSet,
             buildNoncodingGeneSet,
             buildFullGeneSet),
            suffix(".gtf.gz"),
            ".unionintersection.bed.gz" )
def buildUnionIntersectionExons( infile, outfile ):
    '''build union/intersection genes according to Bullard et al. (2010) BMC Bioinformatics.

    Builds a single-segment bed file.
    '''

    to_cluster = USECLUSTER
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --intersect-transcripts --with-utr --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (buildReferenceGeneSet, 
             buildCodingGeneSet,
             buildNovelGeneSet,
             buildLincRNAGeneSet,
             buildNoncodingGeneSet,
             buildFullGeneSet),
            suffix(".gtf.gz"),
            ".union.bed.gz" )
def buildUnionExons( infile, outfile ):
    '''build union genes.

    Exons across all transcripts of a gene are merged.
    They are then intersected between genes to remove any overlap.

    Builds a single-segment bed file.    
    '''

    to_cluster = USECLUSTER
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
# note - needs better implementation, currently no dependency checks.
@follows( buildUnionExons, mkdir( "exon_counts.dir" ) )
@files( [ ( ("%s.accepted.bam" % x.asFile(), "%s.union.bed.gz" % y ),
            ("exon_counts.dir/%s_vs_%s.bed.gz" % (x.asFile(),y ) ) )
          for x,y in itertools.product( TRACKS, GENESETS ) ] )
def buildExonLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = USECLUSTER

    # note: needs to set flags appropriately for
    # single-end/paired-end data sets
    # set filter options
    # for example, only properly paired reads
    paired = False
    if paired:
        flag_filter = "-f 0x2"
    else:
        flag_filter = ""

    # note: the -split option only concerns the stream in A - multiple
    # segments in B are not split. Hence counting has to proceed via
    # single exons - this can lead to double counting if exon counts
    # are later aggregated.

    statement = '''
    samtools view -b %(flag_filter)s -q %(deseq_min_mapping_quality)s %(infile)s
    | coverageBed -abam stdin -b %(exons)s -split
    | sort -k1,1 -k2,2n
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@collate(buildExonLevelReadCounts,
         regex(r"exon_counts.dir/(.+)_vs_(.+)\.bed.gz"),  
         r"\2.exon_counts.load")
def loadExonLevelReadCounts( infiles, outfile ):
    '''load exon level read counts.
    '''
    
    to_cluster = USECLUSTER

    # aggregate not necessary for bed12 files, but kept in
    # ims: edited so that picks up chromosome, start pos and end pos for downstream use.
    src = " ".join( [ "<( zcat %s | cut -f 1,2,3,4,7 )" % x for x in infiles] )

    tmpfile = P.getTempFilename( "." )
    tmpfile2 = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s'''
    
    P.run()

    tracks = [P.snip(x, ".bed.gz" ) for x in infiles ]
    tracks = [re.match( "exon_counts.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

    outf = IOTools.openFile( tmpfile2, "w")
    outf.write( "gene_id\tchromosome\tstart\tend\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        # ims: edit so that now skips five and ens_id is in 3rd index
        genes = list(set([ data[x] for x in range(3,len(data), 5 ) ]))
        # ims: add entries for chromosome, start and ends
        chrom = list(set( [data[x] for x in range(0,len(data), 5 ) ]))
        starts = list(set( [data[x] for x in range(1,len(data), 5 ) ]))
        ends = list(set( [data[x] for x in range(2,len(data), 5 ) ]))
        # ims: edit as value is now in postion 4 and there are 5 columns per line
        values = [ data[x] for x in range(4,len(data), 5 ) ]
        # ims: extra assets for chrom, starts and ends
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        assert len(chrom) == 1, "paste command failed, wrong number of chromosomes per line"
        assert len(starts) == 1, "paste command failed, wrong number of starts per line"
        assert len(ends) == 1, "paste command failed, wrong number of ends per line"
        # ims: add extra coloumns into output
        outf.write( "%s\t%s\t%s\t%s\t%s\n" % (genes[0], chrom[0], starts[0], ends[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    P.load( tmpfile2, outfile )

    os.unlink( tmpfile )
    os.unlink( tmpfile2 )

#########################################################################
#########################################################################
#########################################################################
@follows( buildUnionExons, mkdir( "gene_counts.dir" ) )
@transform( buildBAMs, 
            regex(r"(\S+).accepted.bam"), 
            add_inputs( buildCodingGeneSet ),
            r"gene_counts.dir/\1.gene_counts.tsv.gz" )
def buildGeneLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = USECLUSTER

    statement = '''
    zcat %(exons)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("gene_counts.dir"), buildGeneModels )
@files( [( ([ "%s.accepted.bam" % y.asFile() for y in EXPERIMENTS[x]], buildCodingGeneSet), 
           "gene_counts.dir/%s.gene_counts.tsv.gz" % x.asFile()) 
         for x in EXPERIMENTS ] +\
            [ ( ( ["%s.accepted.bam" % y.asFile() for y in TRACKS], buildCodingGeneSet),
                "gene_counts.dir/%s.gene_counts.tsv.gz" % ALL.asFile()) ] )
def buildAggregateGeneLevelReadCounts( infiles, outfile):
    '''count reads falling into transcripts of protein coding 
       gene models.

    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
       
    '''
    bamfiles, geneset = infiles
    
    to_cluster = USECLUSTER

    bamfiles = ",".join(bamfiles)

    statement = '''
    zcat %(geneset)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(bamfiles)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( ( buildGeneLevelReadCounts,
              buildAggregateGeneLevelReadCounts),
           suffix(".tsv.gz"),
           ".load" )
def loadGeneLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "intron_counts.dir" ) )
@transform( buildBAMs, 
            regex(r"(\S+).accepted.bam"), 
            add_inputs( buildIntronGeneModels ),
            r"intron_counts.dir/\1.intron_counts.tsv.gz" )
def buildIntronLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = USECLUSTER

    statement = '''
    zcat %(exons)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=length
          --prefix="introns_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform(buildIntronLevelReadCounts,
           suffix(".tsv.gz"),
           ".load" )
def loadIntronLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "extension_counts.dir" ) )
@transform( buildBAMs, 
            regex(r"(\S+).accepted.bam"), 
            r"extension_counts.dir/\1.extension_counts.tsv.gz" )
def buildGeneLevelReadExtension( infile, outfile ):
    '''compute extension of cds. 

    Known UTRs are counted as well.
    '''

    to_cluster = USECLUSTER

    cds = os.path.join( PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"] )
    
    territories = os.path.join( PARAMS["annotations_dir"],
                                PARAMS_ANNOTATIONS["interface_territories_gff"] )

    utrs = os.path.join( PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_annotation_gff"] )

    if "geneset_remove_contigs" in PARAMS:
        remove_contigs = '''| awk '$1 !~ /%s/' ''' % PARAMS["geneset_remove_contigs"] 
    else:
        remove_contigs = ""

    statement = '''
    zcat %(cds)s 
    %(remove_contigs)s
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=position
          --counter=read-extension
          --output-filename-pattern=%(outfile)s.%%s.tsv.gz
          --filename-gff=%(territories)s
          --filename-gff=%(utrs)s
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "utr_extension" ) ) )
@transform( buildGeneLevelReadExtension,
            suffix(".tsv.gz"),
            ".plot" )
def plotGeneLevelReadExtension( infile, outfile ):
    '''plot reads extending beyond last exon.'''
    PipelineRnaseq.plotGeneLevelReadExtension( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "utr_extension" ) ) )
@transform( buildGeneLevelReadExtension,
            suffix(".tsv.gz"),
            ".utr.gz" )
def buildUTRExtension( infile, outfile ):
    '''build new utrs.'''
    PipelineRnaseq.buildUTRExtension( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildUTRExtension,
            suffix(".utr.gz"),
            "_utr.load" )
def loadUTRExtension( infile, outfile ):
    P.load( infile, outfile, "--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
@merge( buildUTRExtension, "utrs.bed.gz" )
def buildUTRs( infiles, outfile ):
    '''build new utrs by merging estimated UTR extensions
    from all data sets.
    '''
    
    infiles = " " .join(infiles)

    to_cluster = USECLUSTER

    statement = '''
    zcat %(infiles)s 
    | python %(scriptsdir)s/csv_cut.py contig max_5utr_start max_5utr_end gene_id max_5utr_length strand 
    | awk -v FS='\\t' '$1 != "contig" && $2 != ""' 
    | mergeBed -nms -s 
    > %(outfile)s.5
    '''
    
    P.run()

    statement = '''
    zcat %(infiles)s 
    | python %(scriptsdir)s/csv_cut.py contig max_3utr_start max_3utr_end gene_id max_3utr_length strand 
    | awk -v FS='\\t' '$1 != "contig" && $2 != ""' 
    | mergeBed -nms -s 
    > %(outfile)s.3
    '''
    
    P.run()

    statement = '''
    cat %(outfile)s.5 %(outfile)s.3
    | sort -k 1,1 -k2,2n 
    | gzip 
    > %(outfile)s'''
    
    P.run()

    os.unlink( "%s.5" % outfile )
    os.unlink( "%s.3" % outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("transcript_counts.dir") )
@transform( buildBAMs, 
            regex(r"(\S+).accepted.bam"), 
            add_inputs( buildCodingGeneSet ),
            r"transcript_counts.dir/\1.transcript_counts.tsv.gz" )
def buildTranscriptLevelReadCounts( infiles, outfile):
    '''count reads falling into transcripts of protein coding 
       gene models.

    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
       
    '''
    infile, geneset = infiles
    
    to_cluster = USECLUSTER

    statement = '''
    zcat %(geneset)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=transcripts
          --bam-file=%(infile)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("transcript_counts.dir"), buildGeneModels )
@files( [( ([ "%s.accepted.bam" % y.asFile() for y in EXPERIMENTS[x]], buildCodingGeneSet), 
           "transcript_counts.dir/%s.transcript_counts.tsv.gz" % x.asFile()) 
         for x in EXPERIMENTS ] )
def buildAggregateTranscriptLevelReadCounts( infiles, outfile):
    '''count reads falling into transcripts of protein coding 
       gene models.

    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
       
    .. note:: this step takes very long if multiple bam-files are supplied.
    It has thus been taken out of the pipeline. The aggregate can be derived from summing
    the individual counts anyways.

    '''
    bamfiles, geneset = infiles
    
    to_cluster = USECLUSTER

    bamfiles = ",".join(bamfiles)

    statement = '''
    zcat %(geneset)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=transcripts
          --bam-file=%(bamfiles)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (buildTranscriptLevelReadCounts,
             buildAggregateTranscriptLevelReadCounts),
            suffix(".tsv.gz"),
            ".load" )
def loadTranscriptLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=transcript_id" )

#########################################################################
#########################################################################
#########################################################################
@collate(buildExonLevelReadCounts,
         regex(r"exon_counts.dir/(.+)_vs_(.+)\.bed.gz"),  
         r"\2.exon_counts.tsv.gz")
def aggregateExonLevelReadCounts( infiles, outfile ):
    '''aggregate exon level tag counts for each gene.

    coverageBed adds the following four columns:

    1) The number of features in A that overlapped (by at least one base pair) the B interval.
    2) The number of bases in B that had non-zero coverage from features in A.
    3) The length of the entry in B.
    4) The fraction of bases in B that had non-zero coverage from features in A.

    For bed6: use column 7
    For bed12: use column 13

    This method uses the sum of the number of reads found in all of the exons as the 
    tag count.
    '''
    
    to_cluster = USECLUSTER

    # aggregate not necessary for bed12 files, but kept in
    src = " ".join( [ "<( zcat %s | sort -k4,4 | groupBy -i stdin -g 4 -c 7 -o sum | sort -k1,1)" % x for x in infiles ] )

    tmpfile = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s'''
    
    P.run()

    tracks = [P.snip(x, ".bed.gz" ) for x in infiles ]
    tracks = [re.match( "exon_counts.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ data[x] for x in range(1,len(data), 2 ) ]
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink( tmpfile )

#########################################################################
#########################################################################
#########################################################################
@transform( ( aggregateExonLevelReadCounts),
            suffix(".tsv.gz"),
            ".load" )
def loadAggregateExonLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "deseq" ) ) )
@transform( aggregateExonLevelReadCounts,
            suffix(".exon_counts.tsv.gz"),
            ".deseq")
def runDESeq( infile, outfile ):
    '''estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared 
    cuffdiff.

    Plots are:

    <geneset>_<method>_<level>_<track1>_vs_<track2>_significance.png
        fold change against expression level
    
    '''
    
    to_cluster = USECLUSTER

    outdir = os.path.join( PARAMS["exportdir"], "deseq" )

    geneset, method = outfile.split(".")
    level = "gene"

    # load data 
    R('''suppressMessages(library('DESeq'))''')
    R( '''counts_table <- read.delim( '%s', header = TRUE, row.names = 1, stringsAsFactors = TRUE )''' % infile )

    # get conditions to test
    # note that tracks in R use a '.' as separator
    tracks = R('''colnames(counts_table)''')
    map_track2column = dict( [ (y,x) for x,y in enumerate( tracks ) ] )
    
    sample2condition = [None] * len(tracks)
    conditions = []
    no_replicates = False
    for group, replicates in EXPERIMENTS.iteritems():
        if len(replicates) == 1:
            E.warn( "only one replicate in %s - replicates will be ignored in ALL data sets for variance estimation" % group )
            no_replicates = True

        for r in replicates:
            sample2condition[map_track2column[r.asR()]] = group.asR()
        conditions.append( group )

    ro.globalenv['conds'] = ro.StrVector(sample2condition)

    def build_filename2( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_%(section)s.png" % kwargs
    def build_filename1( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(section)s_%(track)s.png" % kwargs
    def build_filename0( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(section)s.png" % kwargs

    # this analysis follows the 'Analysing RNA-Seq data with the "DESeq" package'
    # tutorial 
    R('''cds <-newCountDataSet( counts_table, conds) ''')
    R('''cds <- estimateSizeFactors( cds )''')

    if no_replicates:
        R('''cds <- estimateVarianceFunctions( cds, method="blind" )''')
    else:
        R('''cds <- estimateVarianceFunctions( cds )''')

    L.info("creating diagnostic plots" ) 
    size_factors = R('''sizeFactors( cds )''')
    R.png( build_filename0( section = "scvplot", **locals() ) )
    R('''scvPlot( cds, ylim = c(0,3))''')
    R['dev.off']()

    R('''vsd <- getVarianceStabilizedData( cds )''' )
    R('''dists <- dist( t( vsd ) )''')
    R.png( build_filename0( section = "heatmap", **locals() ) )
    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )
    R['dev.off']()

    for track in conditions:
        condition = track.asR()
        R.png( build_filename1( section = "fit", **locals() ) )
        R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % condition )
        if not no_replicates:
            R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
            R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
            R['dev.off']()
            R.png( build_filename1( section = "residuals", **locals() ) )
            R('''residualsEcdfPlot( cds, "%s" )''' % condition )
            R['dev.off']()

    L.info("calling differential expression")

    outf = IOTools.openFile( outfile, "w" )
    names = None
    fdr = PARAMS["cuffdiff_fdr"]
    isna = R["is.na"]

    for track1, track2 in itertools.combinations( conditions, 2 ):
        R('''res <- nbinomTest( cds, '%s', '%s' )''' % (track1.asR(),track2.asR()))

        R.png( build_filename2( section = "significance", **locals() ) )
        R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1,
                   col = ifelse( res$padj < %(cuffdiff_fdr)s, "red", "black" ) )''' % PARAMS )
        R['dev.off']()
        if not names:
            names = list(R['res'].names)
            m = dict( [ (x,x) for x in names ])
            m.update( dict(
                    pval = "pvalue", 
                    baseMeanA = "value1", 
                    baseMeanB = "value2",
                    id = "test_id", 
                    log2FoldChange = "lfold") )
            
            header = [ m[x] for x in names ] 
            outf.write( "track1\ttrack2\t%s\tstatus\tsignificant\n" % "\t".join(header))
        else:
            if names != list(R['res'].names):
                raise ValueError( "different column headers in DESeq output: %s vs %s" % (names, list(R['res'].names)))

        rtype = collections.namedtuple( "rtype", names )
        
        for data in zip( *R['res']) :
            d = rtype._make( data )
            outf.write( "%s\t%s\t" % (track1,track2))
            # set significant flag
            if d.padj <= fdr: signif = 1
            else: signif = 0

            # set lfold change to 0 if both are not expressed
            if d.baseMeanA == 0.0 and d.baseMeanB == 0.0:
                d = d._replace( foldChange = 0, log2FoldChange = 0 )

            if isna( d.pval )[0]: status = "FAIL"
            else: status = "OK"

            outf.write( "\t".join( map(str, d) ))
            outf.write("\t%s\t%s\n" % (status, str(signif)))
            
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform( runDESeq,
            suffix(".deseq"), 
            "_deseq.load" )
def loadDESeq( infile, outfile ):
    '''load differential expression results.
    '''
    # add gene level follow convention "<level>_diff"
    
    # if one expression value is 0, the log fc is inf or -inf.
    # substitute with 10

    tablename = P.snip( outfile, ".load") + "_gene_diff" 
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track1
              --index=track2
              --index=test_id
              --table=%(tablename)s 
            > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( loadDESeq, "deseq_stats.tsv" )
def buildDESeqStats( infiles, outfile ):
    tablenames = [P.toTable( x ) for x in infiles ] 
    buildExpressionStats( tablenames, "deseq", outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildDESeqStats,
            suffix(".tsv"), 
            ".load" )
def loadDESeqStats( infile, outfile ):
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
## targets related to exporting results of the pipeline
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "roi" ) ) )
@transform( (loadCuffdiff, loadDESeq),
            regex(r"(.*).load"),
            r"%s/roi/differentially_expressed_\1" % (PARAMS["exportdir"]) )
def buildGeneSetsOfInterest( infile, outfile ):
    '''export gene sets of interest

    * differentially expressed genes

    Regions of interest are exported as :term:`bed` formatted files.
    '''
    
    dbh = connect()
    
    table = P.toTable( infile ) + "_gene_diff"
    track = table[:table.index('_')]

    statement = '''SELECT test_id, track1, track2, info.contig, info.start, info.end, info.strand, lfold
                          FROM %(table)s, 
                               %(track)s_geneinfo AS info
                          WHERE 
                               significant AND
                               info.gene_id = test_id
                 ''' % locals()

    data = Database.executewait( dbh, statement % locals() )

    outfiles = IOTools.FilePool( outfile + "_%s.bed.gz" )

    for test_id, track1, track2, contig, start, end, strand, lfold in data:
        try: 
            lfold = float(lfold)
        except TypeError:
            lfold = 0

        key = "%s_vs_%s" % (track1, track2)
        outfiles.write( key, "%s\t%i\t%i\t%s\t%5.2f\t%s\n" % \
                        (contig, start, end, test_id, lfold, strand))

    outfiles.close()

    P.touch( outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( buildBAMs,
          buildFastQCReport,
          loadTophatStats,
          loadBAMStats,
          loadPicardStats,
          loadContextStats,
          loadMappingStats,
          )
def mapping(): pass

@follows( buildGeneModels,
          loadTranscriptComparison,
          buildAbinitioGeneSet,
          buildReferenceGeneSet,
          buildCodingGeneSet,
          buildFullGeneSet,
          buildLincRNAGeneSet,
          buildNoncodingGeneSet,
          buildNovelGeneSet,
          loadGeneSetsBuildInformation,
          loadClassification,
          loadGeneLevelReadCounts,
          loadIntronLevelReadCounts,
          loadGeneInformation,
          loadGeneSetStats,
          loadGeneSetGeneInformation,
          loadGeneSetTranscriptInformation,
          loadReproducibility,
          loadTranscriptsMappability,
          loadTranscriptLevelReadCounts,
          loadGeneLevelReadCounts,
          loadExonLevelReadCounts,
          )
def genesets(): pass

@follows( buildUTRs, 
          plotGeneLevelReadExtension,
          loadUTRExtension)
def utrs(): pass

@follows( loadCuffdiff,
          loadDESeq,
          buildCuffdiffPlots,
          loadCuffdiffStats,
          loadDESeqStats )
def expression(): pass

@follows( buildGeneSetsOfInterest )
def export(): pass

@follows( loadExonValidation )
def validate(): pass

###################################################################
###################################################################
###################################################################
## export targets
###################################################################
@merge( mapping,  "view_mapping.load" )
def createViewMapping( infile, outfile ):
    '''create view in database for alignment stats.

    This view aggregates all information on a per-track basis.

    The table is built from the following tracks:
    
    tophat_stats: .genome
    mapping_stats: .accepted
    bam_stats: .accepted
    context_stats: .accepted
    picard_stats: .accepted
    '''

    tablename = P.toTable( outfile )
    # can not create views across multiple database, so use table
    view_type = "TABLE"
    
    dbhandle = connect()
    Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

    statement = '''
    CREATE %(view_type)s %(tablename)s AS
    SELECT SUBSTR( b.track, 1, LENGTH(b.track) - LENGTH( '.accepted')) AS track, *
    FROM bam_stats AS b,
          mapping_stats AS m,
          context_stats AS c,          
          picard_stats_alignment_summary_metrics AS a,
          tophat_stats AS t
    WHERE b.track LIKE "%%.accepted" 
      AND b.track = m.track
      AND b.track = c.track
      AND b.track = a.track
      AND SUBSTR( b.track, 1, LENGTH(b.track) - LENGTH( '.accepted')) || '.genome' = t.track
    '''

    Database.executewait( dbhandle, statement % locals() )

    nrows = Database.executewait( dbhandle, "SELECT COUNT(*) FROM view_mapping" ).fetchone()[0]
    
    if nrows == 0:
        raise ValueError( "empty view mapping, check statement = %s" % (statement % locals()) )

    E.info( "created view_mapping with %i rows" % nrows )

    P.touch( outfile )

###################################################################
###################################################################
###################################################################
@follows( createViewMapping )
def views():
    pass

###################################################################
###################################################################
###################################################################
@follows( mapping,
          genesets,
          expression,
          utrs,
          validate,
          export,
          views)
def full(): pass

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "%s/bamfiles" % PARAMS["web_dir"]), 
          mkdir("%s/genesets" % PARAMS["web_dir"]),
          mkdir("%s/classification" % PARAMS["web_dir"]),
          mkdir("%s/differential_expression" % PARAMS["web_dir"]),
          update_report,
          )
def publish():
    '''publish files.'''
    # publish web pages
    P.publish_report()

    # publish additional data
    web_dir = PARAMS["web_dir"]
    project_id = P.getProjectId()

    # directory, files
    exportfiles = {
        "bamfiles" : glob.glob( "*.accepted.bam" ) + glob.glob( "*.accepted.bam.bai" ),
        "genesets": [ "lincrna.gtf.gz", "abinitio.gtf.gz" ],
        "classification": glob.glob("*.class.tsv.gz") ,
        "differential_expression" : glob.glob( "*.cuffdiff.dir" ),
        }
    
    bams = []

    for targetdir, filenames in exportfiles.iteritems():
        for src in filenames:
            dest = "%s/%s/%s" % (web_dir, targetdir, src)
            if dest.endswith( ".bam"): bams.append( dest )
            dest = os.path.abspath( dest )
            if not os.path.exists( dest ):
                os.symlink( os.path.abspath(src), dest )
    
    # output ucsc links
    for bam in bams: 
        filename = os.path.basename( bam )
        track = P.snip( filename, ".bam" )
        print """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/bamfiles/%(filename)s""" % locals()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
