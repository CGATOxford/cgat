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
reconstruction. Also, it is important to detect these conflicting transcripts between replicates
in order to get an idea of the variability of the data. However, if there are only few reads, 
there might be a case for building transcript models using reads from all replicates of an experiment. 
However, there is no reason to merge reads between experiments.

See `figure 4 <http://www.nature.com/nbt/journal/v28/n5/full/nbt.1621.html>`_ from the cufflinks paper
to get an idea about the reliability of transcript construction with varying sequencing depth.

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

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` file (see :ref:`PipelineDocumenation`).

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

The pipeline requires the results from :doc:`pipeline_annotations`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie_             |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|tophat_             |>=1.2.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|cufflinks_          |>=0.9.3            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.12           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+
|R/DESeq             |                   |differential expression                         |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
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

Example data is available at http://www.cgat.org/~andreas/sample_data/rnaseq.tgz.
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

.. todo::

   abstract source of sequences and library type from
   mapping processes.

Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import GTF, IOTools, IndexedFasta
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro

import PipelineGeneset
import PipelineMapping
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
     "pipeline.ini" ] )
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
    glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" )

ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
## genesets build - needs to be defined statically.
GENESETS = ("novel", "abinitio", "reference", "refcoding" )

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
@merge( os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_geneset_gtf"]),
        "reference.gtf.gz" )
def buildReferenceGeneSet( infile, outfile ):
    '''sanitize transcripts file for cufflinks analysis.

    Merge exons separated by small introns.

    Transcript will be ignored that
       * have very long introns (max_intron_size) (otherwise, cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)
       
    The result is run through cuffdiff in order to add the p_id and tss_id tags
    required by cuffdiff. 

    This will only keep sources of the type 'exon'. It will also remove
    any transcripts not in the reference genome.

    Cuffdiff requires overlapping genes to have different tss_id tags.
    '''
    max_intron_size =  PARAMS["max_intron_size"]

    c = E.Counter()

    tmpfilename = P.getTempFilename( "." )
    tmpfilename2 = P.getTempFilename( "." )
    tmpfilename3 = P.getTempFilename( "." )

    tmpf = gzip.open( tmpfilename, "w" )

    E.info( "filtering by contig and removing long introns" )

    contigs = set(IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) ).getContigs())
    
    gene_ids = {}
    for all_exons in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile )) ):

        c.info += 1

        if all_exons[0].contig not in contigs:
            c.missing_contig += 1
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
            L.info( "removing transcript %s" % all_exons[0].transcript_id )
            c.skipped += 1
            continue

        new_exons.sort( key = lambda x: x.start )

        for e in new_exons:
            # add chr prefix 
            tmpf.write( "%s\n" % str(e) )
            c.exons += 1

        c.output += 1

    L.info( "%s" % str(c) )

    tmpf.close()

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

    writePrunedGTF( tmpfilename3, outfile )

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
            "refcoding.gtf.gz" )
def buildCodingGeneSet( infile, outfile ):
    '''build a new gene set with only protein coding 
    transcripts.'''
    
    to_cluster = True
    statement = '''
    zcat %(infile)s | awk '$2 == "protein_coding"' | gzip > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".junctions.gz")
def buildJunctions( infile, outfile ):
    '''build file with splice junctions from gtf file.
    
    A junctions file is a better option than supplying a GTF
    file, as parsing the latter often fails. See:

    http://seqanswers.com/forums/showthread.php?t=7563
    '''
    
    outf = IOTools.openFile( outfile, "w" )
    for gffs in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile, "r" ) )):
        
        end = gffs[0].end
        for gff in gffs[1:]:
            outf.write( "%s\t%i\t%i\t%s\n" % (gff.contig, end, gff.start, gff.strand ) )
            end = gff.end
                        
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".fa")
def buildReferenceTranscriptome( infile, outfile ):
    '''build reference transcriptome
    '''

    to_cluster = USECLUSTER

    statement = '''
    zcat %(infile)s
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

    statement = '''
    bowtie-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra"),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"), 
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
    bowtie_options = "-v 2 --best --strata -a"
    statement = m.build( (infile,), outfile ) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra"),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"), 
            add_inputs( buildJunctions), 
            r"\1.genome.bam" )
def mapReadsWithTophat( infiles, outfile ):
    '''map reads from .fastq or .sra files.

    A list with known splice junctions is supplied.
    '''
    job_options= "-pe dedicated %i -R y" % PARAMS["tophat_threads"]
    to_cluster = USECLUSTER
    m = PipelineMapping.Tophat()
    infile, reffile = infiles
    tophat_options = PARAMS["tophat_options"] + " --raw-juncs <( gunzip < %(reffile)s ) " % locals()
    statement = m.build( (infile,), outfile ) 
    P.run()

############################################################
############################################################
############################################################
@collate( (mapReadsWithTophat, mapReadsWithBowtieAgainstTranscriptome),
          regex(r"(.+)\..*.bam"),  
          add_inputs( buildCodingGeneSet ), 
          r"\1.bam" )
def buildBAMs( infiles, outfile):

    genome, transcriptome, reffile = infiles[0][0], infiles[1][0], infiles[0][1]
    outfile_mismapped = P.snip(outfile, ".bam") + ".mismapped.bam"

    assert genome.endswith( ".genome.bam" )

    to_cluster = USECLUSTER

    statement = '''
    python %(scriptsdir)s/rnaseq_bams2bam.py 
       --force
       --filename-gtf=%(reffile)s
       --filename-mismapped=%(outfile_mismapped)s
       %(transcriptome)s %(genome)s %(outfile)s
    > %(outfile)s.log;
    checkpoint;
    samtools index %(outfile)s 2>&1 >> %(outfile)s.log;
    samtools index %(outfile_mismapped)s 2>&1 >> %(outfile)s.log;
    '''
    P.run()
        
############################################################
############################################################
############################################################
@transform( buildBAMs, suffix(".bam"), ".mismapped.bam" )
def buildMismappedBAMs( infile, outfile ):
    '''pseudo target - update the mismapped bam files.'''
    P.touch( outfile )

############################################################
############################################################
############################################################
@transform( (mapReadsWithTophat, buildBAMs, buildMismappedBAMs), 
            suffix(".bam" ), ".bam.stats")
def buildAlignmentStats( infile, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    to_cluster = USECLUSTER
    
    # replace the SO field from tophat/samtools with coordinate to indicate
    # that the file is sorted by coordinate.
    # naturally - the bam files should be sorted.
    statement = '''
    java -Xmx2g net.sf.picard.analysis.CollectMultipleMetrics
            I=<(samtools view -h %(infile)s | perl -p -e "s/SO:\S+/SO:coordinate/" ) 
            O=%(outfile)s 
            R=%(cufflinks_genome_dir)s/%(genome)s.fa
            ASSUME_SORTED=true
    >& %(outfile)s
    '''
    
    P.run()


############################################################
############################################################
############################################################
@merge( buildAlignmentStats, "alignment_stats.load" )
def loadAlignmentStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    tablename = P.toTable( outfile )


    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( f, ".bam.stats" )
        fn = f + ".alignment_summary_metrics" 
        if not os.path.exists( fn ): 
            E.warn( "file %s missing" % fn )
            continue
        lines = [ x for x in open( fn, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        outf.write( "%s\t%s" % (track,lines[1] ))
        
    outf.close()
    tmpfilename = outf.name

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()

    for suffix, column in ( ("quality_by_cycle_metrics", "cycle"),
                            ("quality_distribution_metrics", "quality") ):

        # some files might be missing - bugs in Picard
        xfiles = [ x for x in infiles if os.path.exists( "%s.%s" % (x, suffix) ) ]

        header = ",".join( [P.snip( x, ".bam.stats") for x in xfiles] )        
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in xfiles ] )

        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --missing=0
                   %(filenames)s
                | python %(scriptsdir)s/csv2db.py
                      --header=%(column)s,%(header)s
                      --replace-header
                      --index=track
                      --table=%(tname)s 
                >> %(outfile)s
                """
    
        P.run()

    os.unlink( tmpfilename )

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
        segment_juncs_version =  _select( lines, "segment_juncs (.*)$" )
        possible_splices = int( _select( lines, "Reported (\d+) total possible splices") )

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
                      --missing=0
                   %(filenames)s
                | python %(scriptsdir)s/csv2db.py
                      --header=%(suffix)s,%(header)s
                      --replace-header
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
@transform( buildBAMs, suffix(".bam"), r"\1.gtf.gz")
def buildGeneModels(infiles, outfile):
    '''build transcript models - run cufflinks on each region seperately
    '''

    infile, infile_mismapped = infiles

    to_cluster = USECLUSTER    
    job_options= "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]

    track = os.path.basename( outfile[:-len(".gtf")] )

    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    infile = os.path.abspath( infile )
    outfile = os.path.abspath( outfile )

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s           
              --reference %(cufflinks_genome_dir)s/%(genome)s.fa
              --num-threads %(cufflinks_threads)i
              %(cufflinks_options)s
              %(infile)s 
    >& %(outfile)s.log;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s;
    mv genes.expr %(outfile)s.genes.expr;
    mv transcripts.expr %(outfile)s.transcripts.expr
    '''

    P.run()

    shutil.rmtree( tmpfilename )


#########################################################################
#########################################################################
#########################################################################
@transform("*.bam", 
           suffix(".bam"), 
           add_inputs(buildReferenceGeneSet),
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
              --reference %(cufflinks_genome_dir)s/%(genome)s.fa
              --num-threads=%(cufflinks_threads)i
              %(cufflinks_options)s
              %(bamfile)s 
    >& %(outfile)s.log;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s;
    mv -f genes.expr %(outfile)s.genes.expr;
    mv -f transcripts.expr %(outfile)s.transcripts.expr
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
        rm -f $(outfile)s.combined.gtf;
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
@files( [( ([ "%s.gtf.gz" % y.asFile() for y in EXPERIMENTS.getTracks(x)], buildCodingGeneSet), 
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
    reffile = "refcoding.gtf.gz"
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

    creates four tables:
    <track>_benchmark
    <track>_tracking
    <track>_transfrags
    <track>_loci
    '''
    tracks, result = Tophat.parseTranscriptComparison( IOTools.openFile( infile ))
    tracks = [ P.snip( os.path.basename(x), ".gtf.gz" ) for x in tracks ]

    tmpfile = P.getTempFilename()
    tmpfile2 = P.getTempFilename()
    tmpfile3 = P.getTempFilename()

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

    for transfrag in Tophat.iterate_tracking( IOTools.openFile( "%s.tracking.gz" % infile, "r") ):

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

    outf.close()
    outf2.close()
    outf3.close()

    tablename = P.toTable( outfile ) + "_tracking"
    statement = '''cat %(tmpfile)s
    | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
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
    
    Transfrags not overlapping previously known annotations are kept
    in order retain novel lincRNA.
    
    Will also build removed.gtf.gz, but not part of split to avoid
    downstream processing.
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
@files( ( ( (buildAbinitioGeneSet, buildCodingGeneSet),
            "novel.gtf.gz" ) ,) )
def buildNovelGeneSet( infiles, outfile ):
    '''builds a gene set by merging the ab-initio gene set and
    the reference gene set.
    
    Ab-initio transcripts overlapping protein coding transcripts in the reference
    gene set are removed (requires the source ``protein_coding`` to be set in the
    reference gtf file). 
    
    The removal is aggressive and works by gene_id - as soon as one transcript of a 
    gene/locus overlaps, all transcripts of that gene/locus are gone.

    Transcripts overlapping on the same strand are merged.

    '''
    abinitio_gtf, reference_gtf = infiles
    
    E.info( "reading index" )

    index = GTF.readAndIndex( 
        GTF.iterator_filtered( GTF.iterator( IOTools.openFile( reference_gtf) ),
                               source = "protein_coding" ) )
    E.info( "indexed genes on %i contigs" % len(index))

    total_genes, remove_genes = set(), set()
    inf = GTF.iterator( IOTools.openFile( abinitio_gtf ) )
    for gtf in inf:
        total_genes.add( gtf.gene_id )
        if index.contains( gtf.contig, gtf.start, gtf.end):
            remove_genes.add( gtf.gene_id )
    
    E.info( "removing %i out of %i genes" % (len(remove_genes), len(total_genes)) )

    tmpfile = P.getTempFile( "." )
    inf = GTF.iterator( IOTools.openFile( abinitio_gtf ) )

    for gtf in inf:
        if gtf.gene_id in remove_genes:
            continue
        
        tmpfile.write( "%s\n" % str(gtf))

    tmpfile.close()
    tmpfilename = tmpfile.name

    # close-by exons need to be merged, otherwise 
    # cuffdiff fails for those on "." strand

    statement = '''
    %(scriptsdir)s/gff_sort pos < %(tmpfilename)s
    | python %(scriptsdir)s/gtf2gtf.py
        --unset-genes="NONC%%06i"
        --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
        --merge-genes
        --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
        --merge-exons
        --merge-exons-distance=5
        --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
        --renumber-genes="NONC%%06i"
        --log=%(outfile)s.log
    | python %(scriptsdir)s/gtf2gtf.py
        --renumber-transcripts="NONC%%06i"
        --log=%(outfile)s.log
    | %(scriptsdir)s/gff_sort genepos 
    | gzip > %(outfile)s
    '''
    P.run()

    os.unlink( tmpfilename )

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
@follows( loadTranscriptComparison, mkdir( os.path.join( PARAMS["exportdir"], "cuffcompare" ) ) )
@transform( compareTranscriptsPerExperiment, 
            suffix(".cuffcompare"), 
            ".reproducibility" )
def buildReproducibility( infile, outfile ):
    '''all-vs-all comparison between samples.

    Compute correlation between expressed transfrags. Transfrags missing
    from another set are ignored.
    '''

    track = TRACKS.factory( filename = outfile[:-len(".reproducibility")] )
    
    replicates = PipelineTracks.getSamplesInTrack( track, TRACKS )

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()    

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
            data = cc.execute( statement ).fetchall()
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
        R('''abline(r)''')
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
        reps.append( ",".join( [ "%s.bam" % r.asFile() for r in replicates] ) )
        labels.append( group.asFile() )

    reps = "   ".join( reps )
    labels = ",".join( labels )

    statement = '''date > %(outfile)s; hostname >> %(outfile)s;
    cuffdiff --output-dir %(outdir)s
             --verbose
             --reference-seq %(cufflinks_genome_dir)s/%(genome)s.fa
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

    Note: converts to log2 fold change.
   
    The cuffdiff output is parsed. Pairwise comparisons
    in which one gene is not expressed (fpkm < fpkm_silent)
    are set to status 'NOCALL'
    '''

    prefix = P.toTable( outfile )
    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    to_cluster = False
    
    # ignore promoters and splicing - no fold change column, but  sqrt(JS)
    for fn, level in ( ("cds_exp.diff", "cds"),
                       ("gene_exp.diff", "gene"),
                       ("isoform_exp.diff", "isoform"),
                       # ("promoters.diff", "promotor"),
                       # ("splicing.diff", "splice"), 
                       ("tss_group_exp.diff", "tss") ):
        
        tablename = prefix + "_" + level + "_diff"

        statement = '''cat %(indir)s/%(fn)s
        | perl -p -e "s/sample_/track/g; s/value_/value/g; s/yes$/1/; s/no$/0/; s/ln\\(fold_change\\)/lfold/; s/p_value/pvalue/"
        | awk -v OFS='\\t' '/test_id/ {print;next;} 
                              { $9 = $9 / log(2); 
                                if( $6 == "OK" && ($7 < %(cuffdiff_fpkm_expressed)f || $8 < %(cuffdiff_fpkm_expressed)f )) { $6 = "NOCALL"; };
                                print; } '
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track1
              --index=track2
              --index=test_id
              --table=%(tablename)s 
         >> %(outfile)s
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
         >> %(outfile)s
         '''
        
        P.run()

#########################################################################
#########################################################################
#########################################################################
def buildExpressionStats( tables, method, outfile ):
    '''build expression summary statistics.
    
    Creates some diagnostic plots in

    <exportdir>/<method> directory.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()    

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
            
            tested = toDict( cc.execute( """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            status = toDict( cc.execute( """SELECT track1, track2, status, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2,status""" % locals() ).fetchall(), 3 )
            signif = toDict( cc.execute( """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE significant
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            fold2 = toDict( cc.execute( """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE (lfold >= 1 or lfold <= -1) AND significant
                                    GROUP BY track1,track2,significant""" % locals() ).fetchall())
            
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
            data = cc.execute('''SELECT i.sum, pvalue 
                                 FROM %(tablename_diff)s, 
                                 %(geneset)s_geneinfo as i 
                                 WHERE i.gene_id = test_id AND significant'''% locals() ).fetchall()
            if len(data):
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
    cc = dbhandle.cursor()    
    
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
            
            data = zip( *cc.execute( statement ))

            pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_significance.png" % locals()
            R.png( pngfile )
            if len(data) == 0:
                E.warn( "no plot for %s - %s -%s vs %s" % ( pngfile, level, track1, track2))
                continue

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
    cc = dbhandle.cursor()    

    for track in tracks:
        table = "%s_ref_gene_expression" % P.quote(track)
        statement = "SELECT gene_id, FPKM / %(scale)f FROM %(table)s" % locals()
        results.append( dict( cc.execute( statement ).fetchall() ) )
    
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
@transform( (buildReferenceGeneSet, 
             buildCodingGeneSet,
             buildNovelGeneSet,
             buildFullGeneSet),
            suffix(".gtf.gz"),
            ".union.bed.gz" )
def buildUnionExons( infile, outfile ):
    '''build union/intersection genes according to Bullard et al. (2010) BMC Bioinformatics.

    Builds a multi-segment bed file
    '''

    to_cluster = USECLUSTER
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --intersect-transcripts --with-utr --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k4,4 -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py --method=block  --log=%(outfile)s.log
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
# note - needs better implementation, currently no dependency checks.
@follows( buildUnionExons, mkdir( "exon_coverage.dir" ) )
@files( [ ( ("%s.bam" % x.asFile(), "%s.union.bed.gz" % y ),
            ("exon_coverage.dir/%s_vs_%s.bed.gz" % (x.asFile(),y ) ) )
          for x,y in itertools.product( TRACKS, GENESETS) ] )
def buildExonCoverage( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = USECLUSTER

    # set filter options
    # for example, only properly paired reads
    flag_filter = "-f 0x2"
    flag_filter = ""

    statement = '''
    samtools view -b %(flag_filter)s %(infile)s
    | bamToBed -i stdin 
    | coverageBed -a stdin -b %(exons)s 
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@collate(buildExonCoverage,
         regex(r"exon_coverage.dir/(.+)_vs_(.+)\.bed.gz"),  
         r"\2.tagcounts.tsv.gz")
def buildRawGeneLevelTagCounts( infiles, outfile ):
    '''aggregate exon level tag counts for each gene.

    coverageBed adds the following four columns:

    1) The number of features in A that overlapped (by at least one base pair) the B interval.
    2) The number of bases in B that had non-zero coverage from features in A.
    3) The length of the entry in B.
    4) The fraction of bases in B that had non-zero coverage from features in A.

    For bed6: use column 7
    For bed12: use column 13
    '''
    
    to_cluster = USECLUSTER

    # aggregate not necessary for bed12 files, but kept in
    src = " ".join( [ "<( zcat %s | groupBy -i stdin -g 4 -c 13 -o sum | sort -k1,1)" % x for x in infiles ] )

    tmpfile = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s'''
    
    P.run()

    tracks = [P.snip(x, ".bed.gz" ) for x in infiles ]
    tracks = [re.match( "exon_coverage.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

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
@follows( mkdir( os.path.join( PARAMS["exportdir"], "deseq" ) ) )
@transform( buildRawGeneLevelTagCounts,
            suffix(".tagcounts.tsv.gz"),
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
    R('''library('DESeq')''')
    R( '''counts_table <- read.delim( '%s', header = TRUE, row.names = 1, stringsAsFactors = TRUE )''' % infile )

    # get conditions to test
    # note that tracks in R use a '.' as separator
    tracks = R('''colnames(counts_table)''')
    map_track2column = dict( [ (y,x) for x,y in enumerate( tracks ) ] )
    
    sample2condition = [None] * len(tracks)
    conditions = []
    for group, replicates in EXPERIMENTS.iteritems():
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

    for track1,track2 in itertools.combinations( conditions, 2 ):
        R('''res <- nbinomTest( cds, '%s', '%s' )''' % (track1.asR(),track2.asR()))

        R.png( build_filename2( section = "significance", **locals() ) )
        R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1,
                   col = ifelse( res$padj < %(cuffdiff_fdr)s, "red", "black" ) )''' % PARAMS )
        R['dev.off']()
        if not names:
            names = list(R['res'].names)
            m = dict( [ (x,x) for x in names ])
            m.update( dict(
                    pval = "pvalue", baseMeanA = "value1", baseMeanB = "value2",
                    id = "test_id", log2FoldChange = "lfold") )
            
            header = [ m[x] for x in names ] 
            outf.write( "track1\ttrack2\t%s\tstatus\tsignificant\n" % "\t".join(header))
        else:
            if names != list(R['res'].names):
                raise ValueError( "different column headers in DESeq output: %s vs %s" % (names, list(R['res'].names)))

        rtype = collections.namedtuple( "rtype", names )
        
        for data in zip( *R['res']) :
            d = rtype._make( data )
            outf.write( "%s\t%s\t" % (track1,track2))
            if d.padj <= fdr: signif = 1
            else: signif = 0
            if isna( d.pval ): status = "OK"
            else: status = "FAIL"

            outf.write( "\t".join( map(str, data) ))
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
@follows( buildBAMs,
          loadTophatStats,
          loadBAMStats,
          loadAlignmentStats,
          loadContextStats,
          )
def mapping(): pass

@follows( buildGeneModels,
          loadTranscriptComparison,
          buildAbinitioGeneSet,
          buildReferenceGeneSet,
          buildCodingGeneSet,
          buildFullGeneSet,
          buildNovelGeneSet,
          loadGeneSetStats,
          loadGeneSetGeneInformation,
          loadGeneSetTranscriptInformation,
          loadReproducibility
          )
def genesets(): pass

@follows( loadCuffdiff,
          loadDESeq,
          buildCuffdiffPlots,
          loadCuffdiffStats,
          loadDESeqStats )
def expression(): pass

@follows( mapping,
          genesets,
          expression)
def full(): pass

def export(): pass

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )

    dirname, basename = os.path.split( os.path.abspath( __file__ ) )
    docdir = os.path.join( dirname, "pipeline_docs", P.snip( basename, ".py" ) )

    # requires libtk, which is not present on the nodes
    to_cluster = False
    
    job_options= "-pe dedicated %i -R y" % PARAMS["report_threads"]

    statement = '''
    rm -rf report _cache _static;
    sphinxreport-build 
           --num-jobs=%(report_threads)s
           sphinx-build 
                    -b html 
                    -d %(report_doctrees)s
                    -c . 
           %(docdir)s %(report_html)s
    > report.log
    '''

    P.run()

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "starting documentation build process from scratch" )

    dirname, basename = os.path.split( os.path.abspath( __file__ ) )
    docdir = os.path.join( dirname, "pipeline_docs", P.snip( basename, ".py" ) )

    # requires libtk, which is not present on the nodes
    to_cluster = False
    
    job_options= "-pe dedicated %i -R y" % PARAMS["report_threads"]

    statement = '''
    sphinxreport-build 
           --num-jobs=%(report_threads)s
           sphinx-build 
                    -b html 
                    -d %(report_doctrees)s
                    -c . 
           %(docdir)s %(report_html)s
    > report.log
    '''

    P.run()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
