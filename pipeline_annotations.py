################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
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
===================
Annotation pipeline
===================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The annotation pipeline imports various annotations and organizes them
for use in other pipelines.

   * a geneset (from ENSEMBL gene sets)
   * repeats (from UCSC repeatmasker tracks)

This pipeline works on a single genome. Annotations are often shared
between versions within the same project or even between projects, hence
this separate pipeline. The output of this pipeline is used by various
other pipelines, for example the the :doc:`pipeline_rnaseq` and the 
:doc:`pipeline_chipseq`.

Overview
========

The pipeline takes as input an ENSEMBL gene set and builds various gene sets
of interest. 

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The :file:`pipeline.ini` needs to be edited so that it points to the
appropriate locations of the auxiliary files. See especially:

1. section ``[ensembl]`` with the location of the ENSEMBL dump
    files (``filename_gtf``, filename_pep``, ``filename_cdna``)

2. section ``[general]`` with the location of the indexed genomic
    fasta files to use and the name of the genome (default=``hg19``),
    see :doc:`IndexedFasta`.

3. section ``[ucsc]`` with the name of the database to use (default=``hg19``).

Input
-----

This script requires no input within the :term:`working directory`. 

Pipeline output
===============

The results of the computation are all stored in an sqlite relational
database file. Additional files are output in the working directory
for use in other pipelines. These are:

csvdb
   An sqlite database with most of the information created by this pipeline

annotation_gff.gz
   A :term:`gff` formatted file annotating the genome with respect to the geneset.
   Annotations are non-overlapping.

geneset_all.gtf.gz
   The full gene set after reconciling with assembly. Chromosomes names are
   renamed to be consistent with the assembly.

geneset_cds.gtf.gz
   A :term:`gtf` formatted file with only the CDS parts of transcripts.
   This set will naturally include only coding transcripts.

geneset_exons.gtf.gz
   A :term:`gtf` formatted file with only the exon parts of transcripts.
   This set includes both coding and non-coding transcripts. Coding 
   transcripts span both the UTR and the CDS.

geneset_flat.gtf.gz
   A :term:`gtf` formatted flattened gene models. All overlapping transcripts
   have been merged. This set includes both coding and non-coding transcripts.

pseudogenes.gtf.gz
   A :term:`gtf` formatted file with pseudogenes. Pseudogenes are either taken from
   the ENSEMBL annotation or processed transcripts with similarity to protein coding
   sequence. As some protein coding genes contain processed transcripts without an ORF,
   Pseudogenes might overlap with protein coding transcripts.

peptides.fasta
   A :term:`fasta` formatted file of peptide sequences of coding transcripts.

cds.fasta
   A :term:`fasta` formatted file of coding sequence of coding transcripts.

cdna.fasta
   A :term:`fasta` formatted file of transcripts including both coding and non-coding parts.

contigs.tsv
   A :term:`tsv` formatted table with contig sizes

tss.bed.gz
   A :term:`bed` formatted file with transcription start sites.

tts.bed.gz
   A :term:`bed` formatted file with transcription termination sites.

promotors.bed.gs
   A :term:`bed` formatted file with promotor regions (fixed witdth segments upstream of 
   transcription start sites).

repeats.bed.gz
   A :term:`bed` formatted file of repetitive sequences (obtained from UCSC repeatmasker tracks).

rna.gff.gz
   A :term:`gff` formatted file of repetitive RNA sequences in the genome
   (obtained from UCSC repeatmasker tracks).

go.tsv.gz
   A list of :term:`GO` assignments for each gene.

goslim.tsv.gz
   A list of :term:`GOSlim` assignments for each gene.

territories.gff.gz
   A :term:`gff` formatted file of non-overlapping gene territories.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz
   tar -xvzf pipeline_annotations.tgz
   cd pipeline_annotations
   python <srcdir>/pipeline_annotations.py make full

Code
====

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *

import sqlite3

# for UCSC import
import MySQLdb

import IndexedFasta, IOTools, GFF, Bed
import PipelineGeneset as PGeneset
import PipelineBiomart as PBiomart
import PipelineDatabase as PDatabase
import PipelineGO

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
PARAMS = P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

if os.path.exists("pipeline_conf.py"):
    E.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

############################################################
############################################################
############################################################
## genome section
############################################################
############################################################
############################################################
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"), 
        PARAMS['interface_contigs'] )
def buildContigSizes( infile, outfile ):
    '''output contig sizes
    '''
    prefix = P.snip( infile, ".fasta" )
    fasta = IndexedFasta.IndexedFasta( prefix )
    outs = IOTools.openFile(outfile, "w" )

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        outs.write( "%s\t%i\n" % ( contig,size) )

    outs.close()

############################################################
############################################################
############################################################
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"), 
        PARAMS["interface_genome_tsv"])
def buildGenomeInformation( infile, outfile ):
    '''compute genome composition information.'''

    to_cluster = True

    statement = '''
    cat %(infile)s
    | python %(scriptsdir)s/fasta2table.py 
        --section=length
        --section=cpg
    | gzip
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( buildGenomeInformation, suffix(".tsv.gz"), ".load" )
def loadGenomeInformation( infile, outfile ):
    '''load genome information.'''
    P.load( infile, outfile )
    
    
############################################################
############################################################
############################################################
## Mappability
@files( os.path.join(PARAMS["gem_dir"],PARAMS["genome"]+".gem"),  PARAMS["genome"]+".mappability" )
def calculateMappability( infile, outfile ):
    '''Calculate mappability using GEM '''
    index = P.snip(infile, ".gem")
    to_cluster = True
    window_size = PARAMS["gem_window_size"]
    threads = PARAMS["gem_threads"]
    mismatches = PARAMS["gem_mismatches"]
    max_indel_length = PARAMS["gem_max_indel_length"]
    statement = '''gem-mappability -t %(threads)s -m %(mismatches)s --max-indel-length %(max_indel_length)s -l %(window_size)s -I %(index)s -o %(outfile)s ''' % locals()
    P.run()

###################################################################
@transform( calculateMappability, suffix(".mappability"), ".mappability.count" )
def countMappableBases( infile, outfile ):
    '''Count mappable bases in genome'''
    to_cluster = True
    statement = '''cat %(infile)s | tr -cd ! | wc -c > %(outfile)s''' % locals()
    P.run()
    
###################################################################
@transform( countMappableBases, suffix(".count"), ".count.load" )
def loadMappableBases( infile, outfile ):
    '''load count of mappable bases in genome'''
    to_cluster = True
    header = "total_mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=total_mappable_bases
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

###################################################################
@transform( calculateMappability, suffix(".mappability"), ".split.log" )
def splitMappabiliyFileByContig( infile, outfile ):
    '''Count mappable bases in genome'''
    to_cluster = True
    track = P.snip( os.path.basename(infile), ".mappability" )
    statement = '''mkdir contigs; 
                   csplit -k -f contigs/contig %(infile)s '/^~[a-zA-Z]/' {100000} > %(outfile)s;
                   rm contigs/contig00;''' % locals()
    P.run()

###################################################################
@follows( splitMappabiliyFileByContig )
@merge( "contigs/contig*", PARAMS["genome"]+"_mappability_per_contig.tsv" )
def countMappableBasesPerContig( infiles, outfile ):
    '''Count mappable bases for each contig'''
    for infile in infiles:
        statement = '''grep '~' %(infile)s | sed s/~//g >> %(outfile)s; cat %(infile)s | tr -cd ! | wc -c >> %(outfile)s'''
        P.run()
    
    statement = '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s;'''
    P.run()

###################################################################
@transform( countMappableBasesPerContig, suffix(".tsv"), ".tsv.load" )
def loadMappableBasesPerContig( infile, outfile ):
    '''load count of mappable bases per contig '''
    to_cluster = True
    header = "contig,mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=mappable_bases_per_contig
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()


############################################################
############################################################
############################################################
## gene set section
############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_geneset_all_gtf'] )
def buildGeneSet( infile, outfile ):
    '''build a gene set - reconciles chromosome names. '''
    to_cluster = True

    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/gff2gff.py 
                  --sanitize=genome 
                  --skip-missing 
                  --genome-file=%(genome_dir)s/%(genome)s 
                  --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_annotation_gff'] )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. In case of overlapping
    genes, only take the longest (in genomic coordinates).
    Genes not on UCSC contigs are removed.
    '''
    PGeneset.buildGeneRegions( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_geneset_flat_gtf'] )
def buildFlatGeneSet( infile, outfile ):
    '''build a flattened gene set.

    All transcripts in a gene are merged into a single transcript. 

    *infile* is an ENSEMBL gtf file.
    '''
    PGeneset.buildFlatGeneSet( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "gene_info.load" )
def loadGeneInformation( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadGeneInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( buildFlatGeneSet, "gene_stats.load" )
def loadGeneStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_cds_gtf"] )
def buildCDSTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS exons are parts of exons are output - UTR's are removed.
    '''
    PGeneset.buildCDS( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_exons_gtf"] )
def buildExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the ENSEMBL gene set.'''
    PGeneset.buildExons( infile, outfile )


############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_coding_exons_gtf"] )
def buildCodingExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PGeneset.buildCodingExons( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_noncoding_exons_gtf"] )
def buildNonCodingExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PGeneset.buildNonCodingExons( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_lincrna_exons_gtf"] )
def buildLincRNAExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PGeneset.buildLincRNAExons( infile, outfile )
    
############################################################
############################################################
############################################################
@transform( buildCDSTranscripts, suffix(".gtf.gz"), "_gtf.load" )
def loadCDSTranscripts( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildCDSTranscripts, "cds_stats.load" )
def loadCDSStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_info.load" )
def loadTranscriptInformation( infile, outfile ):
    '''load transcript information.'''
    
#    PGeneset.loadTranscriptInformation( infile, 
#                                        outfile,
#                                       only_proteincoding = PARAMS["ensembl_only_proteincoding"] )

    tablename = P.toTable( outfile )

    columns = {
        "ensembl_gene_id" : "gene_id",
        "ensembl_transcript_id" : "transcript_id",
        "ensembl_peptide_id" : "protein_id",
        "gene_biotype" : "gene_biotype",
        "transcript_biotype" : "transcript_biotype",
        "source" : "source",
        "status" : "gene_status",
        "transcript_status" : "transcript_status",
        "external_gene_id" : "gene_name",
        "external_transcript_id" : "transcript_name",
        }

    data = PBiomart.biomart_iterator( columns.keys()
                                      , biomart = "ensembl"
                                      , dataset = PARAMS["ensembl_biomart_dataset"]  )

    PDatabase.importFromIterator( outfile
                                  , tablename
                                  , data
                                  , columns = columns 
                                  , indices = ("gene_id", "transcript_id", "protein_id", "gene_name", "transcript_name") )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_pep"], 
           PARAMS["interface_peptides_fasta"] ), ) )
def buildPeptideFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .pep.all.fa.gz file.
    '''
    PGeneset.buildPeptideFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_cdna"], 
           PARAMS["interface_cdna_fasta"] ), ) )
def buildCDNAFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''
    PGeneset.buildCDNAFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@merge( buildCDSTranscripts, PARAMS["interface_cds_fasta"] )
def buildCDSFasta( infile, outfile ):
    '''build cds sequences from peptide and cds file.
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''

    PGeneset.buildCDSFasta( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_pep"], "protein_stats.load" )
def loadProteinStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadProteinStats( infile, outfile )

############################################################
############################################################
############################################################
@files( (buildGeneSet,
         buildPeptideFasta),
        PARAMS["interface_pseudogenes_gtf"] )
def buildPseudogenes( infile, outfile ):
    '''build set of pseudogenes.'''
    PGeneset.buildPseudogenes( infile, outfile )

############################################################
############################################################
############################################################
@files( (None,),
        PARAMS["interface_numts_gtf"] )
def buildNUMTs( infile, outfile ):
    '''build list of NUMTs.'''
    PGeneset.buildNUMTs( infile, outfile )

############################################################
############################################################
############################################################
@merge( (loadProteinStats, loadTranscriptInformation), "seleno.list")
def buildSelenoList( infile, outfile ):
    '''export a list of seleno cysteine transcripts.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    statement = '''
    SELECT DISTINCT transcript_id
    FROM transcript_info as t,
         protein_stats as p
    WHERE p.protein_id = t.protein_id AND
         p.nU > 0
    '''
    outf = open(outfile, "w" )
    outf.write("transcript_id\n")
    outf.write("\n".join( [x[0] for x in cc.execute( statement) ] ) + "\n" )
    outf.close()


############################################################
############################################################
############################################################
@merge( buildFlatGeneSet, PARAMS["interface_territories_gff"] )
def buildGeneTerritories( infile, outfile ):
    '''build gene territories from protein coding genes.'''

    to_cluster=True
    
    statement = '''
    gunzip < %(infile)s
    | awk '$2 == "protein_coding"'
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr
    | python %(scriptsdir)s/gtf2gtf.py --sort=position
    | python %(scriptsdir)s/gtf2gff.py 
          --genome-file=%(genome_dir)s/%(genome)s 
          --log=%(outfile)s.log
          --radius=%(geneset_territories_radius)s
          --method=territories
    | python %(scriptsdir)s/gtf2gtf.py --filter=longest-gene --log=%(outfile)s.log 
    | gzip
    > %(outfile)s '''
    
    P.run()

############################################################
############################################################
############################################################
@merge( buildCodingExonTranscripts, PARAMS["interface_promotors_bed"] )
def buildPromotorRegions( infile, outfile ):
    '''annotate promotor regions from reference gene set.'''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=%(geneset_promotor_size)s --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | gzip 
        > %(outfile)s
    """

    P.run()
    
############################################################
## TRANSCRIPTS
@merge( buildCodingExonTranscripts, PARAMS["interface_transcripts_gtf"] )
def buildTranscripts( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log 
        > %(outfile)s """
    P.run()

############################################################
@follows( buildTranscripts )
@merge( PARAMS["interface_transcripts_gtf"], PARAMS["interface_transcripts_bed"] )
def buildTranscriptsBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        cat %(infile)s
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run() 
        
############################################################
## NON-CODING TRANSCRIPTS
@merge( buildNonCodingExonTranscripts, PARAMS["interface_noncoding_gtf"] )
def buildNoncodingTranscripts( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log 
        > %(outfile)s """
    P.run()

############################################################
@follows(buildNoncodingTranscripts)
@merge( PARAMS["interface_noncoding_gtf"], PARAMS["interface_noncoding_bed"] )
def buildNoncodingTranscriptsBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
############################################################
@follows(buildNoncodingTranscripts)
@merge( PARAMS["interface_noncoding_gtf"], PARAMS["interface_noncoding_genes_bed"] )
def buildNoncodingGenesBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()    

############################################################
@merge( buildLincRNAExonTranscripts, PARAMS["interface_lincrna_gtf"] )
def buildLincRNATranscripts( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log 
        > %(outfile)s """
    P.run()

############################################################
@follows(buildLincRNATranscripts)
@merge( PARAMS["interface_lincrna_gtf"], PARAMS["interface_lincrna_bed"] )
def buildLincRNATranscriptsBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
            
############################################################
############################################################
############################################################
## TRANSCRIPTION START SITES
@merge( buildCodingExonTranscripts, PARAMS["interface_tss_bed"] )
def buildTranscriptTSS( infile, outfile ):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1. '''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()

############################################################
@merge( buildCodingExonTranscripts, PARAMS["interface_tss_gene_bed"] )
def buildGeneTSS( infile, outfile ):
    '''create a single TSS for each gene'''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()

############################################################
@merge( buildCodingExonTranscripts, PARAMS["interface_tss_gene_interval_bed"] )
def buildGeneTSSInterval( infile, outfile ):
    '''create a single interval that encompasses all annotated TSSs for a given gene'''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | sed s/transcript/exon/g | sed s/exon_id/transcript_id/g 
        | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
############################################################
@merge( buildNonCodingExonTranscripts, PARAMS["interface_tss_gene_noncoding_bed"] )
def buildNoncodingGeneTSS( infile, outfile ):
    '''Assign a TSS for each non-coding gene'''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
############################################################
@transform( (buildTranscriptTSS, buildGeneTSS, buildGeneTSSInterval, buildNoncodingGeneTSS), suffix(".bed.gz"), ".extended.bed.gz" )
def ExtendRegion( infile, outfile ):
    '''convert bed to gtf'''
    statement = """gunzip < %(infile)s 
                   | slopBed -i stdin -g %(faidx)s -b 1000  
                   | gzip
                   > %(outfile)s """
    P.run()
    
############################################################
@transform( (buildTranscriptTSS, buildGeneTSS, buildGeneTSSInterval, buildNoncodingGeneTSS, ExtendRegion), suffix(".bed.gz"), ".gtf" )
def convertToGTF( infile, outfile ):
    '''convert bed to gtf'''
    statement = """gunzip < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf  --log=%(outfile)s.log 
                   > %(outfile)s """
    P.run()


############################################################
############################################################
############################################################
## Transcription termination sites
@merge( buildCodingExonTranscripts, PARAMS["interface_tts_bed"] )
def buildTranscriptTTS( infile, outfile ):
    '''annotate transcription termination sites from reference gene set. '''
    statement = """gunzip < %(infile)s 
                   | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gtf2gtf.py --join-exons --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gtf2gff.py --method=tts --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
                   | gzip
                   > %(outfile)s"""
    P.run()

############################################################
@merge( buildCodingExonTranscripts, PARAMS["interface_tts_gene_bed"] )
def buildGeneTTS( infile, outfile ):
    '''annotate transcription termination sites from reference gene set. '''
    statement = """gunzip < %(infile)s 
                   | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gtf2gff.py --method=tts --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
                   | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
                   | gzip
                   > %(outfile)s"""
    P.run()

############################################################
############################################################
############################################################
## UCSC tracks
############################################################
def connectToUCSC():
    dbhandle = MySQLdb.Connect( host = PARAMS["ucsc_host"],
                                user = PARAMS["ucsc_user"] )

    cc = dbhandle.cursor()
    cc.execute( "USE %s " %  PARAMS["ucsc_database"] )

    return dbhandle

def getRepeatsFromUCSC( dbhandle, repclasses, outfile ):
    '''select repeats from UCSC and write to *outfile* in gff format.
    '''

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is 
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [ x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError( "could not find any `rmsk` tables" )

    # now collect repeats
    tmpfile = P.getTempFile(".")
    
    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, strand, '.', '.', 
                      CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\";')
               FROM %(table)s
               WHERE repClass in ('%(repclasses)s') """ % locals() 
        E.debug( "executing sql statement: %s" % sql )
        cc.execute( sql )
        for data in cc.fetchall():
            tmpfile.write( "\t".join(map(str,data)) + "\n" )

    tmpfile.close()

    to_cluster = True

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = '''cat %(tmpfilename)s
        | %(scriptsdir)s/gff_sort pos 
        | python %(scriptsdir)s/gff2gff.py 
            --sanitize=genome 
            --skip-missing 
            --genome-file=%(genome_dir)s/%(genome)s
            --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    '''
    P.run()

    os.unlink( tmpfilename)

############################################################
############################################################
############################################################
@files( ((None, PARAMS["interface_rna_gff"] ), ) )
def importRNAAnnotationFromUCSC( infile, outfile ):
    '''import repetetive RNA from a UCSC formatted file.

    The repetetive RNA are taken from the repeat-masker track.

    The results are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(PARAMS["ucsc_rnatypes"].split(",") )
    dbhandle = connectToUCSC()
    getRepeatsFromUCSC( dbhandle, repclasses, outfile )

############################################################
############################################################
############################################################
@files( ((None, PARAMS["interface_repeats_gff"] ), ) )
def importRepeatsFromUCSC( infile, outfile ):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(PARAMS["ucsc_repeattypes"].split(","))
    dbhandle = connectToUCSC()
    getRepeatsFromUCSC( dbhandle, repclasses, outfile )

#############################################################
@transform( importRepeatsFromUCSC, suffix(".gff.gz"), ".gff.gz.load" )
def loadRepeats( infile, outfile ):
    '''load total repeats length'''

    #statement = """zcat %(infile)s | awk '{print $5-$4}' | awk '{ sum+=$1} END {print sum}' > %(outfile)s; """

    headers = "contig,start,stop,class"
    statement = """zcat %(infile)s | python %(scriptsdir)s/gff2bed.py --name=class | grep -v "#" | cut -f1,2,3,4
                   | python %(scriptsdir)s/csv2db.py 
                         --table=repeats
                         --header=%(headers)s
                         --index=
                 > %(outfile)s; """
    P.run()

#############################################################
@transform( loadRepeats, suffix(".gff.gz.load"), ".counts.load" )
def countTotalRepeatLength( infile, outfile):
    ''' Count total repeat length and add to database '''
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    statement = """create table repeat_length as SELECT sum(stop-start) as total_repeat_length from repeats"""
    cc.execute( statement )
    cc.close()
    
    statement = "touch %(outfile)s"
    P.run()

############################################################
############################################################
############################################################

if 0:
    ############################################################
    ############################################################
    ############################################################
    ## get UCSC tables
    ############################################################
    def getUCSCTracks(infile = PARAMS["filename_ucsc_encode"]):
        '''return a list of UCSC tracks from infile.'''
        tables = []
        with open(infile) as f:
            for line in f:
                if line.startswith("#"): continue
                tablename = line[:-1].strip()
                if tablename == "": continue
                tables.append( tablename )
        return tables

    ############################################################
    ############################################################
    ############################################################
    ## import UCSC encode tracks
    ############################################################
    @posttask( touch_file("ucsc_encode.import") )
    @files( PARAMS["filename_ucsc_encode"], "ucsc_encode.import") 
    def importUCSCEncodeTracks( infile, outfile ):

        statement = '''
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B -e "SELECT * FROM %(tablename)s" %(ucsc_database)s |\
        python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --table=%(tablename)s \
        >> %(outfile)s

        '''

        dbhandle = sqlite3.connect( PARAMS["database"] )

        cc = dbhandle.cursor()
        tables = set( [ x[0] for x in cc.execute( "SELECT name FROM sqlite_master WHERE type='table'") ] )
        cc.close()

        for tablename in getUCSCTracks( infile ):
            if tablename in tables:
                E.info( "skipping %(tablename)s - already exists" % locals())
                continue            

            E.info( "importing %(tablename)s" % locals() )
            P.run()



    ############################################################
    ############################################################
    ############################################################
    ## export UCSC encode tracks as bed
    ############################################################
    @transform( importUCSCEncodeTracks, suffix(".import"), ".bed")
    def exportUCSCEncodeTracks( infile, outfile ):

        dbhandle = sqlite3.connect( PARAMS["database"] )

        outs = open(outfile, "w")
        for tablename in getUCSCTracks():
            outs.write( "track name=%s\n" % tablename )

            cc = dbhandle.cursor()
            statement = "SELECT chrom, chrostart, chroend FROM %s ORDER by chrom, chrostart" % (tablename)
            cc.execute( statement )
            for contig, start, end in cc:
                outs.write("%s\t%i\t%i\n" % (contig, start, end) )
        outs.close()

###################################################################
###################################################################
###################################################################
## gene list analyses
###################################################################
@files( [ (None, PARAMS["interface_go"] ), ] )
def createGO( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGO( infile, outfile )

############################################################
@transform( createGO, 
            regex("(.*)"),
            PARAMS["interface_goslim"])
def createGOSlim( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGOSlim( infile, outfile )

############################################################
@transform( (createGO, createGOSlim), 
            suffix(".tsv.gz"),
            r"\1_assignments.load" )
def loadGOAssignments( infile, outfile ):

    table = P.toTable( outfile )

    statement = '''
    zcat < %(infile)s
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --table=%(table)s
              --index=gene_id
              --index=go_id
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@merge( (importRepeatsFromUCSC, 
         importRNAAnnotationFromUCSC,
         PARAMS["ensembl_filename_gtf"],
         buildFlatGeneSet,
         createGO,
         ),
        PARAMS["interface_genomic_context_bed"] )
def buildGenomicContext( infiles, outfile ):
    '''build a file with genomic context.
    
    The output is a bed formatted file, annotating genomic segments
    according to whether they are any of the ENSEMBL annotations.
    
    It also adds the RNA and repeats annotations from the UCSC.

    The annotations can be partially or fully overlapping.

    Adjacent features (less than 10 bp apart) of the same type are merged.
    '''

    to_cluster = True

    repeats_gff, rna_gff, annotations_gtf, geneset_flat_gff, go_tsv = infiles

    tmpfile = P.getTempFilename( "." )
    tmpfiles = [ "%s_%i" % (tmpfile, x) for x in range( 5 ) ]

    distance=10

    # add ENSEMBL annotations
    statement = """
            gunzip 
            < %(annotations_gtf)s
            | python %(scriptsdir)s/gtf2gtf.py --sort=gene
            | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s 
            | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log 
            | python %(scriptsdir)s/gff2bed.py --name=source --is-gtf --log=%(outfile)s.log
            | sort -k 1,1 -k2,2n
            | python %(scriptsdir)s/bed2bed.py --method=merge --merge-by-name --merge-distance=%(distance)i --log=%(outfile)s.log
            > %(tmpfile)s_0
    """
    P.run()
            
    # rna
    statement = '''
    zcat %(repeats_gff)s %(rna_gff)s 
    | python %(scriptsdir)s/gff2bed.py --name=family --is-gtf -v 0 
    | sort -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py --method=merge --merge-by-name --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_1''' 
    P.run()

    ## add aggregate intervals for repeats
    statement = '''
    zcat %(repeats_gff)s 
    | python %(scriptsdir)s/gff2bed.py --name=family --is-gtf -v 0 
    | awk -v OFS="\\t" '{$4 = "repeats"; print}'
    | sort -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py --method=merge --merge-by-name --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_2''' 
    P.run()

    ## add aggregate intervals for rna
    statement = '''
    zcat %(rna_gff)s 
    | python %(scriptsdir)s/gff2bed.py --name=family --is-gtf -v 0 
    | awk -v OFS="\\t" '{$4 = "repetetive_rna"; print}'
    | sort -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py --method=merge --merge-by-name --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_3 ''' 
    P.run()

    ## add ribosomal protein coding genes
    goids = ("GO:0003735", )
    
    patterns = "-e %s" % ( "-e ".join( goids ) )
    
    statement = ''' 
    zcat %(geneset_flat_gff)s
    | python %(scriptsdir)s/gtf2gtf.py 
        --apply=<(zcat %(go_tsv)s |grep %(patterns)s | cut -f 2 | sort | uniq)
        --filter=gene
        --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py 
        --log=%(outfile)s.log
    | awk -v OFS="\t" '{$4 = "ribosomal_coding"; print}'
    | sort -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py --method=merge --merge-by-name --merge-distance=%(distance)i --log=%(outfile)s.log
    > %(tmpfile)s_4
    '''
    P.run()

    ## sort and merge
    files = " ".join( tmpfiles )
    statement = '''
    sort --merge -k1,1 -k2,2n %(files)s
    | gzip
    > %(outfile)s
    '''
    P.run()

    for x in tmpfiles: os.unlink(x)

@transform( buildGenomicContext, suffix(".bed.gz"), ".tsv")
def buildGenomicContextStats( infile, outfile ):
    '''analysis overlap of genomic contexts.'''
    
    to_cluster= True
    tmpdir = P.getTempDir(".")
    
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/split_file.py
        --pattern-output=%(tmpdir)s/%%s.bed
        --column=4
    > %(outfile)s.log
    '''
    
    P.run()
    
    statement = '''
    python %(scriptsdir)s/diff_bed.py
       %(tmpdir)s/*.bed
    > %(outfile)s
    '''
    P.run()

    
    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
## BED feature files for genes, TSS intervals, 5'/3' flanks and intergenic regions
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_genic_gtf'] )
def buildGeneIntervals( infile, outfile ):
    ''' Merge all transcripts per gene (including utr) to get start and stop 
        coordinates for every protein-coding gene and store in a GTF file'''

    if infile.endswith(".gz"):
        uncompress = "zcat"
    else:
        uncompress = "cat"

    statement = '''%(uncompress)s %(infile)s 
                   | awk '$2 == "protein_coding"' 
                   | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log
                   | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr --log=%(outfile)s.log
                   | grep -v "#" > %(outfile)s;'''
    P.run()

############################################################
@transform(buildGeneIntervals, regex(PARAMS['interface_genic_gtf']), PARAMS['interface_genic_bed'])
def buildGeneBed( infile, outfile ):
    ''' Convert genic GTF to BED'''

    if infile.endswith(".gz"):
        uncompress = "zcat"
    else:
        uncompress = "cat"

    statement = '''%(uncompress)s %(infile)s 
                   | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log
                   | grep -v "#" 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

############################################################
@transform(buildGeneBed, regex(PARAMS['interface_genic_bed']), PARAMS['interface_upstream_flank_bed'] )
def buildUpstreamFlankBed( infile, outfile ):
    ''' build interval upstream of gene start for each entry in bed file'''

    window=PARAMS["geneset_flank"]
    faidx=PARAMS["faidx"]

    statement = '''flankBed -i %(infile)s -g %(faidx)s -l %(window)s -r 0 -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

############################################################
@transform(buildGeneBed, regex(PARAMS['interface_genic_bed']), PARAMS['interface_downstream_flank_bed'] )
def buildDownstreamFlankBed( infile, outfile ):
    ''' build interval downstream of gene start for each entry in bed file'''

    window=PARAMS["geneset_flank"]
    faidx=PARAMS["faidx"]

    statement = '''flankBed -i %(infile)s -g %(faidx)s -l 0 -r %(window)s -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

############################################################
@merge((buildGeneBed, buildUpstreamFlankBed, buildDownstreamFlankBed), PARAMS['interface_intergenic_bed'] )
def buildIntergenicBed( infiles, outfile ):
    ''' Genomic regions not associated with any other features'''
    inlist = " ".join(infiles)

    statement = '''cat %(inlist)s | complementBed -i stdin -g %(faidx)s > %(outfile)s'''
    P.run()

############################################################
@transform((buildUpstreamFlankBed,buildDownstreamFlankBed,buildIntergenicBed), suffix(".bed"), ".gtf" )
def convertBedToGtf( infile, outfile ):
    ''' convert bed files to GTF'''

    statement = '''cat %(infile)s | python %(scriptsdir)s/bed2gff.py --as-gtf | grep -v "#" > %(outfile)s'''
    P.run()

##################################################################
##################################################################
##################################################################
## Primary targets
##################################################################
@follows( buildContigSizes,
          loadGenomeInformation )
def genome():
    '''import information on geneset.'''
    pass

@follows( buildGeneSet,
          buildGeneTerritories,
          loadCDSTranscripts,
          loadTranscriptInformation,
          loadGeneStats,
          loadGeneInformation,
          buildExonTranscripts,
          buildCodingExonTranscripts,
          buildNonCodingExonTranscripts,
          buildPseudogenes,
          buildSelenoList,
          buildTranscripts,
          buildTranscriptsBed,
          buildNoncodingTranscripts,
          buildNoncodingTranscriptsBed)
def geneset():
    '''import information on geneset.'''
    pass

@follows( importRepeatsFromUCSC,
          importRNAAnnotationFromUCSC,
          buildGenomicContext,
          loadRepeats,
          countTotalRepeatLength )
def repeats():
    '''import repeat annotations.'''
    pass

@follows( buildPeptideFasta,
          buildCDSFasta,
          buildCDNAFasta )
def fasta():
    '''build fasta files.'''
    pass

@follows( buildTranscriptTSS, buildGeneTSS,
          buildGeneTSSInterval, buildNoncodingGeneTSS,
          ExtendRegion, convertToGTF,
          buildTranscriptTTS, buildGeneTTS )
def tss():
    '''build promotors.'''
    pass
    
@follows( buildPromotorRegions )
def promotors():
    '''build promotors.'''
    pass

@follows( loadGOAssignments, )
def ontologies():
    '''create and load ontologies'''
    pass

@follows( buildGeneIntervals,
          buildGeneBed,
          buildUpstreamFlankBed,
          buildDownstreamFlankBed,
          buildIntergenicBed,
          convertBedToGtf )
def GenicBedFiles():
    '''build genic bed files.'''
    pass

@follows( calculateMappability, countMappableBases,
          loadMappableBases, splitMappabiliyFileByContig,
          countMappableBasesPerContig, loadMappableBasesPerContig )
def gemMappability():
    '''Count mappable bases in genome'''
    pass

@follows( genome, geneset, repeats, fasta, tss, promotors, ontologies, GenicBedFiles, gemMappability )
def full():
    '''build all targets.'''
    pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish_report():
    '''publish report.'''

    E.info( "publishing report" )
    P.publish_report()


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
