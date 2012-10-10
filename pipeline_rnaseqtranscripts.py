###############################################################################
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
=================================
RNA-Seq Transcript Build pipeline
=================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The rnaseq transcript build pipeline attempts build a variety of gene sets
from reads mapped to a reference genome.

This pipeline works on a single genome.

Overview
========

h_The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

The pipeline builds the following genesets for each input track:

   * Direct cufflinks output:

      1. gene set per replicate

      2. gene set per experiment. Data from individual replicates are combined.

      3. a complete gene set (:file:`agg-agg-agg.gtf.gz`). 
         Data from all experiments are combined

   * Derive gene sets:

      1. full geneset (:file:`full.gtf.gz`)
         all transcripts predicted by cuffcompare. Derived from :file:`agg-agg-agg.gtf.gz`

      2. pruned geneset (:file:`pruned.gtf.gz`)
         only transfrags are kept that are present in at least two samples. 
         Derived from :file:`full.gtf.gz`.

      3. novel geneset (:file:`novel.gtf.gz`)
         only transfrags that do not overlap any of the transcripts in the reference
         gene set. This data set is derived from :file:`pruned.gtf.gz`.
 
      4. lincRNA gene set (:file:`lincrna.gtf.gz`)
         

Novel gene set
---------------

The novel gene set is build from the pruned gene set.

Transcripts are removed based on features in the reference gene set. By default, these
are features called ``protein_coding``, ``lincRNA`` or ``processed_transcript``.
Transcripts that lie exclusively in repetetive sequence are removed, too.

Removal is aggressive  - as soon as one transcript of a gene/locus overlaps, all 
transcripts of that gene/locus are gone.

The resultant set contains a number of novel transcripts. However, these
transcripts will still overlap some known genomic features like pseudogenes.

LincRNA
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

Mapped reads
++++++++++++

The principal input of this pipeline is a collection of reads mapped to a reference genome.
Mapped reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.bam

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. 

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
|cufflinks_          |>=1.3.0            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+
|R/DESeq             |                   |differential expression                         |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.42             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|bamstats_           |>=1.22             |from CGR, Liverpool                             |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

The pipeline will build :term:`gtf` formatted files for each of
the geneset builds.

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

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqtranscripts.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqtranscripts.tgz
   tar -xvzf pipeline_rnaseqtranscripts.tgz
   cd pipeline_rnaseq
   python <srcdir>/pipeline_rnaseqtranscripts.py make full

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

import Expression

import PipelineGeneset
import PipelineMapping
import PipelineRnaseq
import PipelineMappingQC
import Stats

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
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
    glob.glob( "*.bam" ), "(\S+).bam" )

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

#########################################################################
#########################################################################
#########################################################################
def writePrunedGTF( infile, outfile ):
    '''remove various gene models from a gtf file.
    '''
#    to_cluster = True

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
         -s %(bowtie_index_dir)s/%(genome)s.fa
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
    transcripts. 
    
    # Nick
    Also remove transcripts that are ambiguous in
    biotype in terms of their non-codingness

    filter the refnoncoding geneset for things that are described in ensembl
    as being:
    Ambiguous_orf
    Retained_intron
    Sense_intronic
    antisense
    Sense_overlapping
    Processed transcript
    '''

    # Nick - changed the filters for refnoncoding set
    to_cluster = True
    statement = '''zcat %(infile)s 
                   | awk '$2 == "lincRNA" || $2 == "non_coding" || $2 == "3prime_overlapping_ncrna" || $2 == "ncRNA_host"' | gzip > %(outfile)s'''                                                                                                                                             
    P.run()

    # statement = '''
    # zcat %(infile)s | awk '$2 == "lincRNA" || $2 == "processed_transcript"' | gzip > %(outfile)s
    # '''
    # P.run()

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
            "refcoding.gtf.gz" )
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
#########################################################################

# Nick - added building of a mask file for omitting certain regions during gene model building

@files(os.path.join(PARAMS["annotations_dir"], "geneset_all.gtf.gz"), "geneset_mask.gtf")
def buildMaskGtf(infile, outfile):
    '''
    This takes ensembl annotations (geneset_all.gtf.gz) and writes out all entries that 
    have a 'source' match to "rRNA" or 'contig' match to "chrM". for use with cufflinks
    '''
    geneset = IOTools.openFile(infile)
    outf = open(outfile, "wb")
    for entry in GTF.iterator(geneset):
        if re.findall("rRNA", entry.source) or re.findall("chrM", entry.contig):
            outf.write("\t".join((map(str,[entry.contig
                              , entry.source
                              , entry.feature
                              , entry.start
                              , entry.end
                              , "."
                              , entry.strand
                              , "."
                              , "transcript_id" + " " + '"' + entry.transcript_id + '"' + ";" + " " + "gene_id" + " " + '"' + entry.gene_id + '"'])))
                               + "\n")

    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform( "*.bam",
            suffix(".bam"), 
            add_inputs( buildMaskGtf, buildReferenceGeneSet ),
            r"\1.gtf.gz")
def buildTranscriptsWithCufflinks(infiles, outfile):
    '''build transcript models for each track separately.
    '''

    infile, mask_file, reference_file = infiles

    to_cluster = True    
    job_options= "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]
    
    track = os.path.basename( P.snip( outfile, ".gtf.gz" ) )

    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    infile = os.path.abspath( infile )
    outfile = os.path.abspath( outfile )

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    genome_file = os.path.abspath( os.path.join( PARAMS["bowtie_index_dir"], PARAMS["genome"] + ".fa" ) )

    options=PARAMS["cufflinks_options"]

    # Nick - added options to mask rRNA and ChrM from gene modle builiding. 
    # Also added options for faux reads. RABT - see cufflinks manual
    if PARAMS["cufflinks_include_mask"]:
        mask_file = os.path.abspath( mask_file )
        options = options + " -M %s" % mask_file # add mask option

    statement = '''mkdir %(tmpfilename)s; '''

    if PARAMS["cufflinks_include_guide"]:
        # add reference for RABT - this is all genes in reference ensembl 
        # geneset so includes known lincRNA and transcribed pseudogenes
        # TODO: remove explicit file reference
        statement += '''zcat %(reference_file)s > %(tmpfilename)s/reference.gtf; ''' % locals()
        options = options + " --GTF-guide %(tmpfilename)s/reference.gtf" % locals()

    statement += '''
        cd %(tmpfilename)s;
                cufflinks 
              --label %(track)s           
              --num-threads %(cufflinks_threads)i
              --library-type %(tophat_library_type)s
              --frag-bias-correct %(genome_file)s
              --multi-read-correct
              %(options)s
              %(infile)s 
        >& %(outfile)s.log;
        perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s;
       '''

    P.run()

    # version 0.9.3
    #mv genes.expr %(outfile)s.genes.expr;
    #mv transcripts.expr %(outfile)s.transcripts.expr

    shutil.rmtree( tmpfilename )

###################################################################
##################################################################
##################################################################
mapToMethodTargets = { 'cufflinks': buildTranscriptsWithCufflinks,
                       # 'scripture': buildTranscriptsWithScripture,
                       }
METHODTARGET = mapToMethodTargets.get( PARAMS["method"], None )

@follows( METHODTARGET )
def transcripts(): pass

#########################################################################
#########################################################################
#########################################################################
def runCuffCompare( infiles, outfile, reffile ):
    '''run cuffcompare.

    Will create a .tmap and .refmap input file for each input file.
    '''

#    to_cluster = True

    tmpdir = P.getTempDir( "." )
    
    cmd_extract = "; ".join( [ "gunzip < %s > %s/%s" % (x,tmpdir,x) for x in infiles ] )

    genome = os.path.join ( PARAMS["bowtie_index_dir"], PARAMS["genome"]) + ".fa"
    genome = os.path.abspath( genome )

    # note: cuffcompare adds \0 bytes to gtf file - replace with '.'
    statement = '''
        %(cmd_extract)s;
        cuffcompare -o %(outfile)s
                    -s %(genome)s
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
    if 0:
        for t in range( PARAMS["cufflinks_ntries"] ):
            random.shuffle( infiles )
            inf = " ".join( ["%s/%s" % (tmpdir,x) for x in infiles] )
            try:
                P.run()
                break
            except P.PipelineError, msg:
                E.warn("caught exception - trying again" )
    else:
        inf = " ".join( ["%s/%s" % (tmpdir,x) for x in infiles] )
        P.run()
                
    shutil.rmtree( tmpdir )

#########################################################################
#########################################################################
#########################################################################
# @files( [( ([ "%s.gtf.gz" % y.asFile() for y in EXPERIMENTS[x]], buildCodingGeneSet), 
#            "%s.cuffcompare" % x.asFile()) 
#          for x in EXPERIMENTS ] )
@collate( METHODTARGET,
          regex( r"(.*)-(.*)-(.*).gtf.gz" ),
          add_inputs( buildCodingGeneSet ),
          r"\1-\2-agg.cuffcompare" )
def compareTranscriptsPerExperiment( infiles, outfile ):
    '''compare transcript models between replicates within each experiment.'''
    reffile = infiles[0][1]
    infiles = [ x[0] for x in infiles ]
    runCuffCompare( infiles, outfile, reffile )

#########################################################################
#########################################################################
#########################################################################
@merge( METHODTARGET, "%s.cuffcompare" % ALL.asFile() )
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
@merge( compareTranscriptsBetweenExperiments, 
        "full.gtf.gz" )
def buildFullGeneSet( infile, outfile ):
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
@merge( (buildFullGeneSet, buildReferenceGeneSet),
        "pruned.gtf.gz" )
def buildPrunedGeneSet( infiles, outfile ):
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
    tables = Database.getTables( dbhandle )
    if tablename in tables:
        cc = dbhandle.cursor()    
        statement = '''SELECT transfrag_id FROM %(tablename)s WHERE nexperiments > 1''' % locals()
        keep = set( [ x[0] for x in cc.execute(statement).fetchall()] )
        E.info( "keeping %i transfrags" % len(keep) )

    else:
        E.warn( "table %s missing - no replicates - keepy all transfrags" % tablename ) 
        keep = None

    inf = GTF.iterator( IOTools.openFile( abinitio_gtf ) )
    outf1 = IOTools.openFile( keep_gtf, "w" )
    outf2 = IOTools.openFile( remove_gtf, "w" )
    
    c = E.Counter()
    for gtf in inf:
        c.input += 1
        if keep == None or gtf.transcript_id in keep:
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
@merge( (buildPrunedGeneSet, buildReferenceGeneSet,
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

    sections = P.asList( PARAMS["novel_features"] )

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
@merge( (buildPrunedGeneSet, buildReferenceGeneSet,
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
@merge( (buildPrunedGeneSet, buildReferenceGeneSet, buildNoncodingGeneSet,
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_pseudogenes_gtf"] ),
         os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_numts_gtf"] ),
         ), "abinitio_lincrna.gtf.gz" )
def buildAbinitioLincRNAGeneSet( infiles, outfile ):
    '''
    build ab initio lincRNA gene set. 
    In contrast to the buildLincRNAGeneSet this lincRNA set does not contain
    the reference noncoding gene set. It is transcripts in the abinitio (pruned) set that
    do not overlap at any protein coding, processed or pseudogene transcripts 
    (exons+introns) in the reference gene set.

    Transcripts overlapping protein coding transcripts on opposite transcripts are retained.
    It also does not filter out repeats as LincRNA transcripts seem to be found in such regions on occasion
    '''
    
    infile_abinitio, reference_gtf, refnoncoding_gtf, pseudogenes_gtf, numts_gtf = infiles

    E.info( "indexing geneset for filtering" )

    input_sections = ("protein_coding", 
                      "processed_pseudogene",
                      "unprocessed_pseudogene",
                      "nonsense_mediated_decay",
                      "retained_intron")
                    
    indices = {}
    for section in input_sections:
        indices[section] = GTF.readAndIndex( 
            GTF.iterator_filtered(  GTF.merged_gene_iterator(GTF.iterator( IOTools.openFile( reference_gtf ) )),
                                   source = section ),
            with_value = True )
        
    E.info( "built indices for %i features" % len(indices))

    indices["numts"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( numts_gtf) ), with_value = True )

    E.info( "added index for numts" )

    indices["pseudogenes"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( pseudogenes_gtf) ), with_value = True )

    E.info( "added index for pseudogenes" )

    noncoding = {}
    noncoding["noncoding"] = GTF.readAndIndex(GTF.iterator(IOTools.openFile(refnoncoding_gtf)), with_value = True)

    E.info("created index for known noncoding exons to avoid filtering")

    sections = indices.keys()

    total_transcripts, remove_transcripts = set(), collections.defaultdict( set )
    transcript_length = collections.defaultdict( int )
    inf = GTF.iterator( IOTools.openFile( infile_abinitio ) )

    E.info( "collecting genes to remove" )

    min_length = int(PARAMS["lincrna_min_length"])

    for gtfs in GTF.transcript_iterator( inf ): 
        l = sum([x.end-x.start for x in gtfs])
        if l < min_length:
            remove_transcripts[gtfs[0].transcript_id].add("length")

        for section in sections:
            for gtf in gtfs:
                transcript_id = gtf.transcript_id   
                total_transcripts.add( transcript_id )
                
                # only filter on things that don't overlap in the known non-coding set
                # this was a problem originally because of the biotype descriptions e.g. MALAT1 was both lincRNA and processed transcript
                if not noncoding["noncoding"].contains(gtf.contig, gtf.start, gtf.end): 
                    if indices[section].contains( gtf.contig, gtf.start, gtf.end):                    

                    # retain antisense transcripts
                        for gtf2 in indices[section].get(gtf.contig, gtf.start, gtf.end):
                            if gtf.strand == gtf2[2].strand:
                                remove_transcripts[transcript_id].add( section )
                    
    E.info( "removing %i out of %i transcripts" % (len(remove_transcripts), len(total_transcripts)) )
    
    outf = open("lincrna_removed.tsv", "w")
    outf.write("transcript_id" + "\t" + "removed" + "\n" )
    for x, y in remove_transcripts.iteritems():
        outf.write("%s\t%s\n" % (x, ",".join(y)))

    # write out transcripts that are not in removed set
    outf = gzip.open(outfile, "w")
    for entry in GTF.iterator(IOTools.openFile(infile_abinitio)):
        if entry.transcript_id in remove_transcripts: continue
        outf.write("%s\n" % str(entry))
    outf.close()

###################################################################
##################################################################
##################################################################
GENESETTARGETS = []
mapToGenesetTargets = { 'full': buildFullGeneSet,
                        'pruned': buildPrunedGeneSet,
                        'novel': buildNovelGeneSet,
                        'lincrna': buildLincRNAGeneSet,
                        'abinitio_lincrna': buildAbinitioLincRNAGeneSet}
                        
for x in P.asList( PARAMS["genesets"]):
    GENESETTARGETS.append( mapToGenesetTargets[x] )

@follows( *GENESETTARGETS )
def genesets(): pass

#########################################################################
#########################################################################
#########################################################################
@merge( GENESETTARGETS + [ buildTranscriptsWithCufflinks, 
          compareTranscriptsPerExperiment,
          compareTranscriptsBetweenExperiments,
          buildReferenceGeneSet,
          buildCodingGeneSet],
        "geneset_stats.tsv" )
def buildGeneSetStats( infiles, outfile ):
    '''compile gene set statistics.
    '''

    to_cluster = True

    infiles = IOTools.flatten( infiles )

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
@transform( GENESETTARGETS + [
        buildReferenceGeneSet,
        buildCodingGeneSet],
            suffix(".gtf.gz"),
            "_geneinfo.load" )
def loadGeneSetGeneInformation( infile, outfile ):
    PipelineGeneset.loadGeneStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( GENESETTARGETS + [
        buildReferenceGeneSet,
        buildCodingGeneSet ],
            suffix(".gtf.gz"),
            "_transcript2gene.load" )
def loadGeneInformation( infile, outfile ):
    PipelineGeneset.loadTranscript2Gene( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( GENESETTARGETS + [
        buildReferenceGeneSet,
        buildCodingGeneSet ],
            suffix(".gtf.gz"),
            "_transcriptinfo.load" )
def loadGeneSetTranscriptInformation( infile, outfile ):
    PipelineGeneset.loadTranscriptStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( GENESETTARGETS + [buildTranscriptsWithCufflinks,],
            suffix(".gtf.gz"), 
            add_inputs( buildReferenceGeneSetWithCDS ),
            ".class.tsv.gz" )
def classifyTranscripts( infiles, outfile ):
    '''classify transcripts.
    '''
    to_cluster = True
    
    infile, reference = infiles


    #IMS: changed to allow different classifiers
    counter = PARAMS['gtf2table_classifier']

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gtf2table.py
           --counter=%(counter)s  
           --reporter=transcripts
           --filename-gff=%(reference)s
           --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()


## need to change pipeline logic to avoid this duplication
@transform( ( compareTranscriptsPerExperiment,
              compareTranscriptsBetweenExperiments),
            suffix(".cuffcompare"), 
            add_inputs( buildReferenceGeneSetWithCDS ),
            ".class.tsv.gz" )
def classifyTranscriptsCuffcompare( infiles, outfile ):
    '''classify transcripts.
    '''
    to_cluster = True
    
    infile, reference = infiles

    #IMS: change to allow different classifiers
    counter = PARAMS['gtf2table_classifier']

    statement = '''
    zcat %(infile)s.combined.gtf.gz
    | python %(scriptsdir)s/gtf2table.py
           --counter=%(counter)s 
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
@transform( (classifyTranscripts, classifyTranscriptsCuffcompare), 
            suffix(".tsv.gz"), 
            ".load" )
def loadClassification( infile, outfile ):
    P.load( infile, outfile, 
            options = "--index=transcript_id --index=match_gene_id --index=match_transcript_id --index=source" )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("transcript_counts.dir") )
@transform( "*.bam",
            regex(r"(.*).bam"),
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
    
    to_cluster = True

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
@transform( buildTranscriptLevelReadCounts,
            suffix(".tsv.gz"),
            ".load" )
def loadTranscriptLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=transcript_id" )

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

    tables = Database.getTables( dbhandle )
    if tablename2 not in tables:
        E.warn( "table %s missing - no replicates" % tablename2 )
        P.touch( outfile )
        return

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
            print statement
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
    tablename = P.toTable(infile)
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log --allow-empty=True < %(infile)s'''

#########################################################################
#########################################################################
#########################################################################
@follows( loadGeneSetStats, 
          loadGeneSetGeneInformation,
          loadGeneInformation,
          loadGeneSetTranscriptInformation,
          loadClassification,
          loadTranscriptLevelReadCounts,
          loadReproducibility,
          )
def qc(): pass

###################################################################
###################################################################
###################################################################
@follows( transcripts, genesets, qc )
def full(): pass

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

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
