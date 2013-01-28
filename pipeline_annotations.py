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
   Annotations are non-overlapping and are based only on protein coding transcripts.

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

gc_segmentation.bed.gz
   A :term:`bed` formatted file with the genome segmented in regions
   of different G+C content.

cpg.bed.gz
   A list of all CpGs in the genome sequence

gwas_catalog.bed.gz
   Bed intervals with GWAS intervals from the gwas catalog.

gwas_distild.bed.gz
   Bed intervals with GWAS intervals from the DistiLD database.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_annotations.tgz
   tar -xvzf pipeline_annotations.tgz
   cd pipeline_annotations.dir
   python <srcdir>/pipeline_annotations.py make full

Code
====

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import Experiment as E
import Pipeline as P
from ruffus import *
from bx.bbi.bigwig_file import BigWigFile
import sqlite3
# for UCSC import
import MySQLdb
import IndexedFasta, IOTools, GFF, GTF
import PipelineGeneset as PipelineGeneset
import PipelineBiomart as PBiomart
import PipelineDatabase as PDatabase
import PipelineGO
import PipelineUCSC
import PipelineKEGG
import Intervals

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
PARAMS = P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )

if os.path.exists("pipeline_conf.py"):
    E.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    return dbh

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
    '''output contig sizes.
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
        PARAMS['interface_contigs_bed'] )
def buildContigBed( infile, outfile ):
    '''output bed file with contigs
    '''
    prefix = P.snip( infile, ".fasta" )
    fasta = IndexedFasta.IndexedFasta( prefix )
    outs = IOTools.openFile(outfile, "w" )

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        outs.write( "%s\t%i\t%i\n" % ( contig,0,size) )

    outs.close()

############################################################
############################################################
############################################################
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"), 
        (PARAMS['interface_contigs_ungapped_bed'],
         PARAMS['interface_gaps_bed'],
         ))
def buildUngappedContigBed( infile, outfiles ):
    '''output bed file with contigs. 

    Genomic regions with gaps are excluded.
    '''
    prefix = P.snip( infile, ".fasta" )
    fasta = IndexedFasta.IndexedFasta( prefix )
    outs_nogap = IOTools.openFile(outfiles[0], "w" )
    outs_gap = IOTools.openFile(outfiles[1], "w" )
    min_gap_size = PARAMS["gaps_min_size"]

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():

        seq = fasta.getSequence( contig )

        def gapped_regions( seq ):
            is_gap = seq[0] == "N"
            last = 0
            for x,c in enumerate(seq):
                if c == "N":
                    if not is_gap:
                        last = x
                        is_gap = True
                else:
                    if is_gap:
                        yield( last, x )
                        last = x
                        is_gap = False
            if is_gap:
                yield last, size
                
        last_end = 0
        for start, end in gapped_regions(seq):
            if end - start < min_gap_size: continue
            
            if last_end != 0:
                outs_nogap.write( "%s\t%i\t%i\n" % ( contig,last_end,start) )
            outs_gap.write( "%s\t%i\t%i\n" % ( contig,start,end) )
            last_end = end

        if last_end < size:
            outs_nogap.write( "%s\t%i\t%i\n" % ( contig,last_end,size) )

    outs_nogap.close()
    outs_gap.close()

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
    statement = '''gem-mappability -t %(gem_threads)s -m %(gem_mismatches)s 
                                   --max-indel-length %(gem_max_indel_length)s 
                                   -l %(gem_window_size)s 
                                   -I %(index)s -o %(outfile)s ''' 
    P.run()

###################################################################
@transform( calculateMappability, suffix(".mappability"), ".mappability.count" )
def countMappableBases( infile, outfile ):
    '''Count mappable bases in genome'''
    to_cluster = True
    statement = '''cat %(infile)s | tr -cd ! | wc -c > %(outfile)s''' 
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
                   rm contigs/contig00;''' 
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
    statement = '''cat %(infile)s 
                   | python %(scriptsdir)s/csv2db.py
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
def annotateGenome( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. 
    Processed_transcripts tend to cover larger genomic regions
    and often overlap between adjacent protein coding genes.

    In case of overlapping genes, only take the longest 
    (in genomic coordinates).

    Genes not on UCSC contigs are removed.
    '''
    PipelineGeneset.annotateGenome( infile, 
                                    outfile,
                                    only_proteincoding = True )

############################################################
############################################################
############################################################
@follows( annotateGenome )
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_geneset_flat_gtf'] )
def buildFlatGeneSet( infile, outfile ):
    '''build a flattened gene set.

    All transcripts in a gene are merged into a single transcript. 

    *infile* is an ENSEMBL gtf file.
    '''
    PipelineGeneset.buildFlatGeneSet( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "gene_info.load" )
def loadGeneInformation( infile, outfile ):
    '''load the transcript set.'''
    PipelineGeneset.loadGeneInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( buildFlatGeneSet, "gene_stats.load" )
def loadGeneStats( infile, outfile ):
    '''load the transcript set.'''
    PipelineGeneset.loadGeneStats( infile, outfile )

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
    PipelineGeneset.buildCDS( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_exons_gtf"] )
def buildExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildExons( infile, outfile )

############################################################
############################################################
############################################################
@files( buildExonTranscripts, "exon_stats.load" )
def loadExonStats( infile, outfile ):
    '''load the transcript set stats.'''
    PipelineGeneset.loadTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_coding_exons_gtf"] )
def buildCodingExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildCodingExons( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_noncoding_exons_gtf"] )
def buildNonCodingExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildNonCodingExons( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_geneset_lincrna_exons_gtf"] )
def buildLincRNAExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set. '''
    PipelineGeneset.buildLincRNAExons( infile, outfile )
    
############################################################
############################################################
############################################################
@transform( buildCDSTranscripts, suffix(".gtf.gz"), "_gtf.load" )
def loadCDSTranscripts( infile, outfile ):
    '''load the transcript set.'''
    PipelineGeneset.loadTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildCDSTranscripts, "cds_stats.load" )
def loadCDSStats( infile, outfile ):
    '''load the transcript set stats.'''

    PipelineGeneset.loadTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_info.load" )
def loadTranscriptInformation( infile, outfile ):
    '''load transcript information.'''
    
#    PipelineGeneset.loadTranscriptInformation( infile, 
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
        "uniprot_sptrembl": "uniprot_id",
        "uniprot_genename": "uniprot_name",
        }

    data = PBiomart.biomart_iterator( columns.keys()
                                      , biomart = "ensembl"
                                      , dataset = PARAMS["ensembl_biomart_dataset"] )
    
    # The full list of genes from this table is too extensive. The following are removed:
    # 1. Some genes are present as LRGxxx identifiers
    # """LRG stands for Locus Reference Genomic. An LRG is a fixed sequence, 
    # independent of the genome, specifically created for the diagnostic community 
    # to record DNA sequence variation on a fixed framework"""
    # These are removed below:
    data = filter( lambda x: not x['ensembl_gene_id'].startswith( "LRG"), data )

    # 2. Some genes are present on artificial chromosomes such as
    # ENSG00000265928 on HG271_PATCH.
    # To filter these out, the gene ids are cross-checked against those in
    # the ensembl gtf file.
    gene_ids = set()
    with IOTools.openFile( infile ) as inf:
        for gtf in GTF.iterator(inf):
            gene_ids.add( gtf.gene_id )

    data = filter( lambda x: x['ensembl_gene_id'] in gene_ids, data )

    PDatabase.importFromIterator( outfile
                                  , tablename
                                  , data
                                  , columns = columns 
                                  , indices = ("gene_id", "transcript_id", "protein_id", "gene_name", "transcript_name", "uniprot_id") )

    # validate: 1:1 mapping between gene_ids and gene_names
    dbh = connect()
    cc = dbh.cursor()
    data = cc.execute( """SELECT gene_name, count(distinct gene_id) from %(tablename)s 
                        GROUP BY gene_name HAVING count(distinct gene_id) > 1""" % locals())
    l = data.fetchall()
    if len(l) > 0:
        E.warn( "there are %i gene_names mapped to different gene_ids" % len(l))
    for gene_name, counts in l:
        E.info( "ambiguous mapping: %s->%i" % (gene_name, counts))
        

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_synonyms.load" )
def loadTranscriptSynonyms( infile, outfile ):
    '''load table with synonyms for transcript identifiers.'''
    
    tablename = P.toTable( outfile )

    columns = {
        "ensembl_transcript_id" : "transcript_id",
        "external_transcript_id" : "transcript_name",
        "refseq_mrna" : "refseq_id",
        }

    data = PBiomart.biomart_iterator( columns.keys()
                                      , biomart = "ensembl"
                                      , dataset = PARAMS["ensembl_biomart_dataset"] )
    
    PDatabase.importFromIterator( outfile
                                  , tablename
                                  , data
                                  , columns = columns 
                                  , indices = ("transcript_id", "transcript_name", "refseq_id") )


###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_pep"], 
           PARAMS["interface_peptides_fasta"] ), ) )
def buildPeptideFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .pep.all.fa.gz file.
    '''
    PipelineGeneset.buildPeptideFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_cdna"], 
           PARAMS["interface_cdna_fasta"] ), ) )
def buildCDNAFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''
    PipelineGeneset.buildCDNAFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@merge( buildCDSTranscripts, PARAMS["interface_cds_fasta"] )
def buildCDSFasta( infile, outfile ):
    '''build cds sequences from peptide and cds file.
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''

    PipelineGeneset.buildCDSFasta( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_pep"], "protein_stats.load" )
def loadProteinStats( infile, outfile ):
    '''load the transcript set.'''

    PipelineGeneset.loadProteinStats( infile, outfile )

############################################################
############################################################
############################################################
@files( (buildGeneSet,
         buildPeptideFasta),
        PARAMS["interface_pseudogenes_gtf"] )
def buildPseudogenes( infile, outfile ):
    '''build set of pseudogenes.'''
    PipelineGeneset.buildPseudogenes( infile, outfile )

############################################################
############################################################
############################################################
@files( (None,),
        PARAMS["interface_numts_gtf"] )
def buildNUMTs( infile, outfile ):
    '''build list of NUMTs.'''
    PipelineGeneset.buildNUMTs( infile, outfile )

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
        | gzip
        > %(outfile)s """
    P.run()

############################################################
@follows( buildTranscripts )
@merge( PARAMS["interface_transcripts_gtf"], PARAMS["interface_transcripts_bed"] )
def buildTranscriptsBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        zcat %(infile)s
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
        | gzip
        > %(outfile)s """
    P.run()

############################################################
@follows(buildNoncodingTranscripts)
@merge( PARAMS["interface_noncoding_gtf"], PARAMS["interface_noncoding_bed"] )
def buildNoncodingTranscriptsBed( infile, outfile ):
    '''annotate transcripts from reference gene set. '''
    statement = """
        zcat < %(infile)s 
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
        | gzip
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
@follows( buildContigSizes )
@transform( (buildTranscriptTSS, buildGeneTSS, buildGeneTSSInterval, buildNoncodingGeneTSS), 
            suffix(".bed.gz"), 
            ".extended.bed.gz" )
def extendRegion( infile, outfile ):
    '''convert bed to gtf'''
    statement = """gunzip < %(infile)s 
                   | slopBed -i stdin -g %(interface_contigs)s -b 1000  
                   | gzip
                   > %(outfile)s """
    P.run()
    
############################################################
@transform( (buildTranscriptTSS, buildGeneTSS, buildGeneTSSInterval, buildNoncodingGeneTSS, extendRegion), suffix(".bed.gz"), ".gtf" )
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
@files( ((None, PARAMS["interface_rna_gff"] ), ) )
def importRNAAnnotationFromUCSC( infile, outfile ):
    '''import repetetive RNA from a UCSC formatted file.

    The repetetive RNA are taken from the repeat-masker track.

    The results are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(PARAMS["ucsc_rnatypes"].split(",") )
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC( dbhandle, repclasses, outfile )

############################################################
############################################################
############################################################
@files( ((None, PARAMS["interface_repeats_gff"] ), ) )
def importRepeatsFromUCSC( infile, outfile ):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(PARAMS["ucsc_repeattypes"].split(","))
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC( dbhandle, repclasses, outfile )

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
@files( ((None, PARAMS["interface_allrepeats_gff"] ), ) )
def importAllRepeatsFromUCSC( infile, outfile ):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses = None
    dbhandle = PipelineUCSC.connectToUCSC()
    PipelineUCSC.getRepeatsFromUCSC( dbhandle, repclasses, outfile )

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################
## build bed files with mappable regions
############################################################
@transform( os.path.join( PARAMS["ucsc_dir"], 
                      "gbdb",
                      PARAMS["ucsc_database"],
                      "bbi",
                      "*rgMapability*.bw"),
            regex( ".*rgMapabilityAlign(\d+)mer.bw" ),
            add_inputs( os.path.join( PARAMS["genome_dir"],
                                      PARAMS["genome"] + ".fasta" ) ),
            r"mapability_\1.bed.gz" )
def buildMapableRegions( infiles, outfile ): 
    '''build bed files with mappable regions.

    Convert bigwig tracks with mappability information to a
    bed-formatted file that contains only mappable regions of the
    genome.

    A mapable region is more permissive than a mapable position.

    This method assumes that files use the ``CRG Alignability
    tracks``.

    UCSC says:
 
      The CRG Alignability tracks display how uniquely k-mer sequences
      align to a region of the genome. To generate the data, the
      GEM-mappability program has been employed. The method is
      equivalent to mapping sliding windows of k-mers (where k has been
      set to 36, 40, 50, 75 or 100 nts to produce these tracks) back to
      the genome using the GEM mapper aligner (up to 2 mismatches were
      allowed in this case). For each window, a mapability score was
      computed (S = 1/(number of matches found in the genome): S=1 means
      one match in the genome, S=0.5 is two matches in the genome, and
      so on). The CRG Alignability tracks were generated independently
      of the ENCODE project, in the framework of the GEM (GEnome
      Multitool) project.

    For the purpose of these tracks, a region is defined to be un-mapable
    if its maximum mapability score is less than 0.5. 
    
    Unmapable regions that are less than half kmer size are mapable, as
    reads from the left/right mapable positions will extend into the region.
    '''

    infile, fastafile = infiles
    fasta = IndexedFasta.IndexedFasta( P.snip( fastafile, ".fasta" ) )
    contigs = fasta.getContigSizes( with_synonyms = False )
    
    kmersize = int(re.search( ".*Align(\d+)mer.bw", infile ).groups()[0])
    
    E.info( "creating mapable regions bed files for kmer size of %i" % kmersize )

    max_distance = kmersize // 2

    f = open(infile)
    bw = BigWigFile(file=f)
    
    def _iter_mapable_regions( bw, contig, size ):
        
        min_score = PARAMS["ucsc_min_mappability"]

        # there is no iterator access, results are returned as list
        # thus proceed window-wise in 10Mb windows
        window_size = 10000000
        last_start, start = None, None

        for window_start in xrange(0, size, window_size):
            values = bw.get( contig, window_start, window_start + window_size )
            if values == None: continue

            for this_start, this_end, value in values:
                if value < min_score:
                    if start: yield start, this_start
                    start = None
                else:
                    if start == None: 
                        start = this_start

        if start != None:
            yield start, this_end
                
    outf = IOTools.openFile( outfile, "w" )

    for contig, size in contigs.iteritems():

        last_start, last_end = None, None
        for start, end in _iter_mapable_regions( bw, contig, size ):
            if last_start == None: last_start, last_end = start, end
            if start - last_end >= max_distance:
                outf.write( "%s\t%i\t%i\n" % (contig, last_start, last_end ) )
                last_start = start

            last_end = end

        if last_start != None:
            outf.write( "%s\t%i\t%i\n" % (contig, last_start, last_end ) )

    outf.close()

############################################################
############################################################
############################################################
## 
############################################################
@transform( buildMapableRegions, suffix( ".bed.gz"), ".filtered.bed.gz" )
def filterMapableRegions( infile, outfile ):
    '''remove small windows from a mapability track.

    Too many fragmented regions will cause gat to fail as it
    fragments the workspace into too many individual segments.

    The filtering works by merging all segments that are
    within mapability_merge_distance and removing all those
    that are larger than mapabpility_min_segment_size
    '''
    
    to_cluster = True

    statement = '''
    mergeBed -i %(infile)s -d %(mapability_merge_distance)i
    | awk '$3 - $2 >= %(mapability_min_segment_size)i'
    | gzip 
    > %(outfile)s
    '''

    P.run()

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
    PipelineGO.createGOFromENSEMBL( infile, outfile )

############################################################
@transform( createGO, 
            regex("(.*)"),
            PARAMS["interface_goslim"])
def createGOSlim( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGOSlimFromENSEMBL( infile, outfile )

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
@transform( createGO, suffix( ".tsv.gz"), ".paths" )
def buildGOPaths( infile, outfile ):
    '''compute a file with paths of each GO term to the ancestral node.'''
    infile = P.snip( infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOPaths( infile, outfile )

############################################################
############################################################
############################################################
@transform( createGO, suffix( ".tsv.gz"), ".desc.tsv" )
def buildGOTable( infile, outfile ):
    '''build a simple table with GO descriptions in obo.'''
    infile = P.snip( infile, ".tsv.gz") + "_ontology.obo"
    PipelineGO.buildGOTable( infile, outfile )

@files(None,PARAMS['interface_kegg'])
def importKEGGAssignments(infile,outfile):
    ''' import the KEGG annotations from the R KEGG.db 
    annotations package. Note that since KEGG is no longer
    publically availible, this is not up-to-date and maybe removed
    from bioconductor in future releases '''

    biomart_dataset = PARAMS["KEGG_dataset"]
    mart = PARAMS["KEGG_mart"]

    PipelineKEGG.importKEGGAssignments(outfile, mart, biomart_dataset)

@transform(importKEGGAssignments, suffix(".tsv.gz"),"_assignments.load")
def loadKEGGAssignments(infile, outfile):

    P.load(infile,outfile, options = "-i gene_id -i kegg_id")


############################################################
############################################################
############################################################
@follows( loadTranscriptInformation, loadGOAssignments )
@files( [ (None, PARAMS["interface_go_geneontology"]), ] )
def createGOFromGeneOntology( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PipelineGO.createGOFromGeneOntology( infile, outfile )

############################################################
@transform( createGOFromGeneOntology, 
            suffix( ".tsv.gz"), 
            add_inputs(buildGOPaths),
            ".imputed.tsv.gz")
def imputeGO( infiles, outfile ):
    PipelineGO.imputeGO( infiles[0], infiles[1], outfile )

############################################################
############################################################
############################################################
@merge( (buildGeneTerritories, loadGOAssignments),
        ( PARAMS["interface_genomic_function_bed"],
          PARAMS["interface_genomic_function_tsv"],
          ) )
def buildGenomicFunctionalAnnotation( infiles, outfiles ):
    '''output a bed file with genomic regions with functional annotations.

    Each bed entry is a gene territory. Bed entries are labeled
    by functional annotations associated with a gene.

    Ambiguities in territories are resolved by outputting 
    annotations for all genes within a territory.

    The output file contains annotations for both GO and GOSlim. These
    are prefixed by ``go:`` and ``goslim:``.
    '''

    territories_file = infiles[0]
    
    dbh = connect()

    PipelineGeneset.buildGenomicFunctionalAnnotation( territories_file, dbh, outfiles )

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
    | awk -v OFS="\\t" '{$4 = "ribosomal_coding"; print}'
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
                   | gzip 
                   > %(outfile)s;'''
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
                   | python %(scriptsdir)s/bed2bed.py 
                                  --method=filter-genome 
                                  --genome-file=%(genome_dir)s/%(genome)s 
                                  --log %(outfile)s.log 
                   | gzip
                   > %(outfile)s'''
    P.run()

############################################################
@follows( buildContigSizes )
@transform(buildGeneBed, regex(PARAMS['interface_genic_bed']), PARAMS['interface_upstream_flank_bed'] )
def buildUpstreamFlankBed( infile, outfile ):
    ''' build interval upstream of gene start for each entry in bed file'''

    statement = '''flankBed -i %(infile)s -g %(interface_contigs)s -l %(geneset_flank)i -r 0 -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome 
                             --genome-file=%(genome_dir)s/%(genome)s 
                             --log %(outfile)s.log 
                   | gzip
                   > %(outfile)s'''
    P.run()

############################################################
@follows( buildContigSizes )
@transform(buildGeneBed, regex(PARAMS['interface_genic_bed']), PARAMS['interface_downstream_flank_bed'] )
def buildDownstreamFlankBed( infile, outfile ):
    ''' build interval downstream of gene start for each entry in bed file'''

    statement = '''flankBed -i %(infile)s -g %(interface_contigs)s -l 0 -r %(geneset_flank)i -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome 
                                 --genome-file=%(genome_dir)s/%(genome)s 
                                 --log %(outfile)s.log 
                   | gzip
                   > %(outfile)s'''
    P.run()

############################################################
@follows( buildContigSizes )
@merge((buildGeneBed, buildUpstreamFlankBed, buildDownstreamFlankBed), PARAMS['interface_intergenic_bed'] )
def buildIntergenicBed( infiles, outfile ):
    ''' Genomic regions not associated with any other features'''
    inlist = " ".join(infiles)

    statement = '''zcat %(inlist)s | complementBed -i stdin -g %(interface_contigs)s 
                   | gzip 
                   > %(outfile)s'''
    P.run()

############################################################
@transform((buildUpstreamFlankBed,buildDownstreamFlankBed,buildIntergenicBed), suffix(".bed.gz"), ".gtf.gz" )
def convertBedToGtf( infile, outfile ):
    ''' convert bed files to GTF'''

    statement = '''zcat %(infile)s 
                  | python %(scriptsdir)s/bed2gff.py --as-gtf 
                  | grep -v "#" 
                  | gzip 
                  > %(outfile)s'''
    P.run()

##################################################################
##################################################################
##################################################################
## build G+C segmentation
##################################################################
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) + ".fasta",
        PARAMS["interface_gc_segmentation_bed"] )
def buildGenomeGCSegmentation( infile, outfile ):
    '''segment the genome into windows according to G+C content.'''

    to_cluster = True
    
    statement = '''
    python %(scriptsdir)s/fasta2bed.py 
        --method=fixed-width-windows-gc
        --window-size=%(segmentation_window_size)i 
        --log=%(outfile)s.log 
    < %(infile)s 
    | python %(scriptsdir)s/bed2bed.py 
        --method=bins 
        --num-bins=%(segmentation_num_bins)s 
        --binning-method=%(segmentation_method)s 
        --log=%(outfile)s.log 
    | bgzip
    > %(outfile)s'''

    P.run()

@files( os.path.join( PARAMS["genome_dir"], 
                      PARAMS["genome"]) + ".fasta",
        "gcprofile.bed.gz" )
def runGenomeGCProfile( infile, outfile ):
    '''segment the genome into windows according to G+C content.'''

    # on some cgat109 I got libstc++ error:
    # error while loading shared libraries: libstdc++.so.5
    # cannot open shared object file: No such file or directory
    to_cluster = False
    
    statement = '''
    cat %(infile)s
    | python %(scriptsdir)s/fasta2bed.py 
        --verbose=2
        --method=GCProfile
        --min-length=%(segmentation_min_length)i
        --halting-parameter=%(segmentation_halting_parameter)i
        --log=%(outfile)s.log
    | bgzip
    > %(outfile)s
    '''
    P.run()

@merge( runGenomeGCProfile, PARAMS["interface_gc_profile_bed"] )
def buildGenomeGCProfile( infile, outfile ):
    '''aggregate windows with similar G+C content into bins.
    '''
    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/bed2bed.py 
        --method=bins 
        --num-bins=%(segmentation_num_bins)s 
        --binning-method=%(segmentation_method)s 
        --log=%(outfile)s.log 
    | bgzip
    > %(outfile)s'''
    P.run()

##################################################################
##################################################################
##################################################################
## build a bed file with locations of CpGs
##################################################################
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta"), 
        PARAMS['interface_cpg_bed'] )
def buildCpGBed( infile, outfile ):
    '''bulid bed file with CpG locations.'''

    statement = ''' 
    python %(scriptsdir)s/fasta2bed.py 
        --method=cpg
        --log=%(outfile)s.log 
    < %(infile)s 
    | bgzip
    > %(outfile)s
    '''
    
    P.run()

    statement = '''
    tabix -p bed %(outfile)s
    '''
    P.run()

##################################################################
##################################################################
##################################################################
## download GWAS data
##################################################################
if PARAMS["genome"].startswith("hg"):
    @merge( None, "gwascatalog.txt" )
    def downloadGWASCatalog( infile, outfile ):
        '''download GWAS catalog.'''

        if os.path.exists( outfile ):
            os.path.remove(outfile)
        statement = '''wget http://www.genome.gov/admin/gwascatalog.txt'''
        P.run()
        
    @merge( downloadGWASCatalog, PARAMS["interface_gwas_catalog_bed"] )
    def buildGWASTracks( infile, outfile ):
        
        reader = csv.DictReader( IOTools.openFile( infile ), dialect = "excel-tab" )
        
        tracks = collections.defaultdict( lambda : collections.defaultdict( list ) )

        fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta" ) )
        contigsizes = fasta.getContigSizes()
        c = E.Counter()

        for row in reader:
            c.input += 1
            contig, pos, snp, disease = row['Chr_id'], row['Chr_pos'], row['SNPs'], row['Disease/Trait'] 
            if snp == "NR":             
                c.skipped += 1
                continue
            
            if pos == "":
                c.no_pos += 1
                continue

            # translate chr23 to X
            if contig == "23": contig = "X"

            contig = "chr%s" % contig
            
            try:
                tracks[disease][contig].append( int(pos) )
            except ValueError:
                print row
            c.output += 1
        
        E.info( c )
            
        extension = PARAMS["gwas_extension"]

        c = E.Counter()
        outf = IOTools.openFile(outfile, "w" )
        for disease, pp in tracks.iteritems():
            
            for contig, positions in pp.iteritems():
                contigsize = contigsizes[contig]
                regions = [ (max(0,x-extension), min( contigsize, x+extension)) for x in positions ]

                regions = Intervals.combine( regions )
                c[disease] += len(regions)

                for start,end in regions:
                    outf.write( "%s\t%i\t%i\t%s\n" % (contig, start, end, disease ) )
                
        outf.close()

        outf = IOTools.openFile(outfile + ".log", "w" )
        outf.write( "category\tcounts\n%s\n" % c.asTable() )
        outf.close()

    @merge( None, "gwas_distild.log" )
    def downloadDistiLD( infile, outfile ):
        '''download GWAS data from distild database.'''

        track = P.snip( outfile, ".log" )
        of = track + "_snps.tsv.gz"
        if os.path.exists( of ): os.path.remove(of)
        statement = '''wget http://distild.jensenlab.org/snps.tsv.gz -O %(of)s'''
        P.run()

        of = track + "_lds.tsv.gz"
        if os.path.exists( of ): os.path.remove(of)
        statement = '''wget http://distild.jensenlab.org/lds.tsv.gz -O %(of)s'''
        P.run()

        P.touch( outfile )
        
    
    @merge( downloadGWASCatalog, PARAMS["interface_gwas_distild_bed"] )
    def buildDistiLDTracks( infile, outfile ):
        '''build bed tracks from DistiLD database.'''
        
        track = P.snip( infile, ".log" )
        intervals = []
        c = E.Counter()
        for line in IOTools.openFile( track + "_snps.tsv.gz" ):
            pubmed_id, rs, pvalue, block, ensgenes, short, icd10 = line[:-1].split("\t")
            c.input += 1
            try:
                contig, start, end = re.match( "(\S+):(\d+)-(\d+)", block ).groups()
            except AttributeError:
                E.warn( "parsing error for %s" % block )
                c.errors += 1
                continue
            intervals.append( (contig, int(start), int(end), short) )
            c.parsed += 1

        intervals.sort()
        outf = IOTools.openFile( outfile, "w" )
        cc = E.Counter()
        for k, x in itertools.groupby( intervals, key = lambda x: x ):
            outf.write( "%s\t%i\t%i\t%s\n" % k )
            c.output += 1
            cc[k[3]] += 1
        outf.close()
        E.info( c )

        outf = IOTools.openFile(outfile + ".log", "w" )
        outf.write( "category\tcounts\n%s\n" % cc.asTable() )
        outf.close()

##################################################################
##################################################################
##################################################################
## build gff summary
##################################################################
@transform( (annotateGenome, 
             buildGeneTerritories, 
             importRNAAnnotationFromUCSC,
             importRepeatsFromUCSC ),
            suffix(".gff.gz"), ".summary.tsv.gz" )
def buildGFFSummary( infile, outfile ):
    '''summarize genomic coverage of gff file.'''
    tocluster = True
    statement = '''zcat %(infile)s 
                | python %(scriptsdir)s/gff2coverage.py 
                      --genome-file=%(genome_dir)s/%(genome)s
                | gzip > %(outfile)s
                '''
    P.run()

##################################################################
##################################################################
##################################################################
## build bed summary
##################################################################
@transform( (buildContigBed,
             buildPromotorRegions,
             buildTranscriptTSS,
             buildGeneTSS,
             buildMapableRegions,
             buildGenomicContext,
             buildGenomeGCSegmentation ),
            suffix(".bed.gz"), ".summary.tsv.gz" )
def buildBedSummary( infile, outfile ):
    '''summarize genomic coverage of bed file.'''
    tocluster = True
    statement = '''zcat %(infile)s 
                | python %(scriptsdir)s/bed2gff.py 
                | python %(scriptsdir)s/gff2coverage.py 
                      --genome-file=%(genome_dir)s/%(genome)s
                | gzip > %(outfile)s
                '''
    P.run()

@transform( (buildGFFSummary, buildBedSummary), suffix(".tsv.gz"), ".load" )
def loadIntervalSummary( infile, outfile ):
    P.load(infile, outfile)

##################################################################
##################################################################
##################################################################
## Primary targets
##################################################################
@follows( buildContigSizes,
          buildContigBed,
          loadGenomeInformation )
def genome():
    '''import information on geneset.'''
    pass

@follows( buildGeneSet,
          buildGeneTerritories,
          loadCDSTranscripts,
          loadTranscriptInformation,
          loadGeneStats,
          loadCDSStats,
          loadExonStats,
          loadGeneInformation,
          loadTranscriptSynonyms,
          buildExonTranscripts,
          buildCodingExonTranscripts,
          buildNonCodingExonTranscripts,
          buildPseudogenes,
          buildNUMTs,
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
          extendRegion, convertToGTF,
          buildTranscriptTTS, buildGeneTTS )
def tss():
    '''build promotors.'''
    pass
    
@follows( buildPromotorRegions )
def promotors():
    '''build promotors.'''
    pass

@follows( loadGOAssignments, 
          loadKEGGAssignments,
          buildGenomicFunctionalAnnotation)
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

@follows( loadIntervalSummary )
def summary():
    '''summary targets.'''
    pass


@follows( calculateMappability, countMappableBases,
          loadMappableBases, splitMappabiliyFileByContig,
          countMappableBasesPerContig, loadMappableBasesPerContig )
def gemMappability():
    '''Count mappable bases in genome'''
    pass

# taken out gemMappability as not fully configured
@follows( genome, geneset, repeats, fasta, tss, promotors, ontologies, GenicBedFiles )
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
    
