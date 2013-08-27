################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_rnaseq_geneset.py 2900 2012-03-38 14:38:00Z david $
#
#   Copyright (C) 2012 David Sims
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
========================
RNAseq Geneset Pipeline
========================

:Author: David Sims 
:Release: $Id: pipeline_rnaseq_genesst.py 2900 2012-03-28 14:38:00Z david $
:Date: |today|
:Tags: Python

The RNAseq geneset pipeline parses multiple GTF files derived from RNAseq experiments in different tissues and produces a consensus geneset

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>` documentation).
The configuration file is organized by section and the variables are documented within 
the file. In order to get a local configuration file in the current directory, type::

    python <codedir>/pipeline_rnaseq_genest.py config

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.


Input
-----

Input are GTF-formatted files. 

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
|BEDTools            |                   |interval comparison                             |
+--------------------+-------------------+------------------------------------------------+


Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.


Code
====

"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import gzip
import sqlite3
import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import cStringIO
import numpy
import CGAT.Masker as Masker
import fileinput
import CGAT.Experiment as E
import logging as L
from ruffus import *

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters(  ["pipeline.ini", ] )
PARAMS = P.PARAMS
#PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["geneset_dir"],"pipeline_annotations.py" )

###################################################
###################################################        
###################################################
## Parse transcripts from GTF file
@transform( "*.txt", regex(r"(\S+).txt"), r"\1.gtf" )
def convertTxtToGtf( infile, outfile ):
    '''Convert  txt to Gtf'''
    track = P.snip( os.path.basename(infile), ".txt" )
    statement = '''cat %(infile)s 
                   | awk 'OFS="\\t" {print "chr"$3,"rnaseq","exon",$4,$5,0,substr($6,0,1),".","transcript_id \\""$1"\\"; gene_id \\""$1"\\"; exon_number \\""substr($2,20,2)"\\"; "}' 
                   | sed s/0\\\\t1\\\\t\./0\\\\t+\\\\t\./g | sed s/chrUltra/Ultra/ | sed s/chrContig/Contig/ > %(outfile)s'''
    P.run()
    
###################################################
@transform( convertTxtToGtf, regex(r"(\S+).gtf"), r"\1.transcripts.gtf" )
def getGtfStrandedTranscripts( infile, outfile ):
    '''join exons to get transcripts from GTF file'''
    track = P.snip( os.path.basename(infile), ".gtf" )
    statement = '''cat %(infile)s | python %(scriptsdir)s/gtf2gtf.py --join-exons | sort -k1,1 -k4,4n > %(outfile)s'''
    P.run()
    
###################################################
@transform( getGtfStrandedTranscripts, suffix(".gtf"), ".coding.gtf")
def getEnsemblCodingGeneset( infile, outfile ):
    '''identify transcrpts that overlap an ensembl coding gene '''
    ensembl_genes = PARAMS["ensembl_genes"]
    ensembl_noncoding = PARAMS["ensembl_noncoding"]
    # need to remove transcripts that overlap 100% with noncoding transcripts
    statement = '''cat %(infile)s | intersectBed -a stdin -b %(ensembl_genes)s -u -s 
                   | intersectBed -a stdin -b %(ensembl_noncoding)s -v -s -f 1 -r
                   | sort -k1,1 -k4,4n > %(outfile)s;
                   echo "transcripts with ensembl coding overlap: " > %(outfile)s.count; 
                   cat %(outfile)s | wc -l >> %(outfile)s.count;'''
    P.run()

###################################################
@transform( getEnsemblCodingGeneset, suffix(".gtf"), ".rename.gtf")
def renameTranscripts( infile, outfile ):
    '''systematically rename transcripts to remove duplicate cuffdiff names '''
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \\"rnaseq_coding_transcript_"NR"\\"; "}' > %(outfile)s;'''
    P.run()
    
###################################################
@transform(renameTranscripts, suffix(".gtf"), ".ensg.gtf" )
def annotateTranscripts( infile, outfile ):
    ''' Add ensembl gene id to GTF file'''
    ensembl_genes = PARAMS["ensembl_genes"]
    statement = '''cat %(infile)s 
                   | intersectBed -a stdin -b %(ensembl_genes)s -wa -wb -s 
                   | awk -F "\\t" 'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \\""$13"\\"; "$9}'
                   | sort -k1,1 -k4,4n > %(outfile)s;'''
    P.run()
    
###################################################
@transform(annotateTranscripts, regex(r"(\S+).ensg.gtf"), r"all_transcripts.gtf" )
def addMissingEnsemblTranscripts( infile, outfile ):
    ''' Add ensembl gene id to GTF file'''
    ensembl_transcripts = PARAMS["ensembl_transcripts_gtf"]
    statement = '''cat %(infile)s 
                   | intersectBed -a %(ensembl_transcripts)s -b stdin -v -s -f 1 -r > missing_ensembl_transcripts.gtf;
                   cat %(infile)s missing_ensembl_transcripts.gtf | sort -k1,1 -k4,4n
                   > %(outfile)s;'''
    P.run()

###################################################
@follows( mkdir("geneset") )
@transform( addMissingEnsemblTranscripts, regex(r"all_transcripts.gtf"), "geneset/transcripts.gtf.gz" )
def renameTranscriptToExon( infile, outfile ):
    '''reformat transcript file for use in proj007 pipeline '''
    statement = """cat %(infile)s | sed s/\\\\ttranscript\\\\t/\\\\texon\\\\t/g
                 | gzip > %(outfile)s; """
    P.run()
    
###################################################
###################################################    
###################################################
## Non-coding
@transform( getGtfStrandedTranscripts, suffix(".gtf"), ".noncoding.gtf")
def getNoncodingGeneset( infile, outfile ):
    '''Assume that all transcripts the do not overlap with ensembl coding geneset are noncoding '''
    ensembl_transcripts = PARAMS["ensembl_transcripts"]
    statement = '''cat %(infile)s | intersectBed -a stdin -b %(ensembl_transcripts)s -v -s > %(outfile)s;
                   echo "transcripts without ensembl coding overlap: " > %(outfile)s.count; 
                   cat %(outfile)s | wc -l >> %(outfile)s.count;'''
    P.run()
    
###################################################
@transform( getNoncodingGeneset, suffix(".gtf"), ".rename.gtf")
def renameNoncodingTranscripts( infile, outfile ):
    '''systematically rename transcripts to remove duplicate cuffdiff names '''
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,"transcript_id \\"rnaseq_noncoding_transcript_"NR"\\"; "}' > %(outfile)s;'''
    P.run()
    
###################################################
@transform(renameNoncodingTranscripts, suffix(".gtf"), ".ensg.gtf" )
def annotateNoncodingTranscripts( infile, outfile ):
    ''' Add ensembl gene id to GTF file'''
    ensembl_noncoding = PARAMS["ensembl_noncoding_genes"]
    statement = '''intersectBed -a %(infile)s  -b %(ensembl_noncoding)s -wa -wb -s 
                   | awk -F "\\t" 'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \\""$13"\\"; "$9}' > known_noncoding_exons.gtf;  
                   intersectBed -a %(infile)s -b %(ensembl_noncoding)s -v -s
                   | awk -F "\\t" 'OFS="\\t" {print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \\"novel_gene_"NR"\\"; "$9}' > novel_noncoding_transcripts.gtf; 
                   cat known_noncoding_exons.gtf novel_noncoding_transcripts.gtf | sort -k1,1 -k4,4n
                   > %(outfile)s;'''
    P.run()
    
###################################################
@transform(annotateNoncodingTranscripts, regex(r"(\S+).ensg.gtf"), r"all_noncoding_transcripts.gtf" )
def addMissingNoncodingTranscripts( infile, outfile ):
    ''' Add ensembl gene id to GTF file'''
    ensembl_noncoding = PARAMS["ensembl_noncoding_gtf"]
    statement = '''intersectBed -a %(ensembl_noncoding)s -b %(infile)s  -v -s -f 1 -r > missing_ensembl_noncoding_transcripts.gtf;
                   cat %(infile)s missing_ensembl_noncoding_transcripts.gtf | sort -k1,1 -k4,4n
                   > %(outfile)s;'''
    P.run()

###################################################
###################################################
###################################################
## Load transcript_info table to database
@transform( addMissingEnsemblTranscripts, suffix(".gtf"), ".tab" )
def transcriptGtfToTab( infile, outfile ):
    '''Copy replicated Bed files generated by capseq pipline to geneset-specific output directory'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/gtf2tsv.py -f --log=%(outfile)s.log | awk '{if (NR==1) {print $0"\\tgene_biotype\\tgene_name"} else {print $0"\\tprotein_coding\\tunknown"}}' > %(outfile)s'''
    P.run()

###################################################
@transform( addMissingNoncodingTranscripts, suffix(".gtf"), ".tab" )
def noncodingGtfToTab( infile, outfile ):
    '''Copy replicated Bed files generated by capseq pipline to geneset-specific output directory'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/gtf2tsv.py -f --log=%(outfile)s.log | awk '{if (NR==1) {print $0"\\tgene_biotype\\tgene_name"} else {print $0"\\tnoncoding\\tunknown"}}' > %(outfile)s'''
    P.run()

###################################################
@follows(transcriptGtfToTab, noncodingGtfToTab)
@merge( "*.tab", "transcript_info.tab" )
def mergeTranscriptInfo( infiles, outfile ):
    '''Copy replicated Bed files generated by capseq pipline to geneset-specific output directory'''
    inlist = " ".join(infiles)
    statement = '''cat %(inlist)s | awk '{if (NR==1 || $1 != "contig") {print $0}}' > %(outfile)s'''
    P.run()
            
###################################################
@transform( mergeTranscriptInfo, suffix(".tab"), ".tab.load" )
def loadTranscripts( infile, outfile ):
    '''load GTF file into database '''
    headers = "contig,source,feature,start,end,score,strand,frame,gene_id,transcript_id,gene_biotype"
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --database=%(database)s
                         --table=transcript_info
                         --index=gene_id
                         --index=contig,start
                         --index=transcript_id
                 > %(outfile)s; """
    P.run()
        
###################################################
###################################################
###################################################
## create input files for CAPseq interval annotation pipeline
@follows( mkdir("geneset") )
@transform(getEnsemblCodingGeneset, regex(r"(\S+).coding.gtf"), r"geneset/\1.genes.bed" )
def buildGeneIntervals( infile, outfile ):
    ''' Merge all transcripts per gene (including utr) to get start and stop 
        coordinates for every protein-coding gene and store in a GTF file'''
    ensembl_genes = PARAMS["ensembl_genes_gtf"]
    statement = '''cat %(ensembl_genes)s %(infile)s | mergeBed -i stdin -s -n 
                   | awk 'OFS="\\t" {print $1,$2,$3,"gene"NR,$4,$5 }'
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
                   | sort -k1,1 -k2,2n > %(outfile)s;'''
    P.run()

###################################################
@transform(buildGeneIntervals, regex(r"geneset/(\S+).genes.bed"), r"geneset/genes.bed" )
def annotateGeneIntervals( infile, outfile ):
    ''' Add ensembl gene id to GTF file'''
    ensembl_genes = PARAMS["ensembl_genes"]
    statement = '''cat %(infile)s 
                   | intersectBed -a stdin -b %(ensembl_genes)s -wa -wb -s 
                   | awk 'OFS="\\t" {print $1,$2,$3,$10,$5,$6}' 
                   | sort -k1,1 -k2,2n > %(outfile)s;'''
    P.run()

###################################################
@transform(annotateGeneIntervals, regex(PARAMS['interface_genic_bed']), PARAMS['interface_upstream_flank_bed'] )
def buildUpstreamFlankBed( infile, outfile ):
    ''' build interval upstream of gene start for each entry in bed file'''
    window=PARAMS["geneset_flank"]
    faidx=PARAMS["faidx"]
    statement = '''flankBed -i %(infile)s -g %(faidx)s -l %(window)s -r 0 -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
@transform(annotateGeneIntervals, regex(PARAMS['interface_genic_bed']), PARAMS['interface_downstream_flank_bed'] )
def buildDownstreamFlankBed( infile, outfile ):
    ''' build interval downstream of gene start for each entry in bed file'''
    window=PARAMS["geneset_flank"]
    faidx=PARAMS["faidx"]
    statement = '''flankBed -i %(infile)s -g %(faidx)s -l 0 -r %(window)s -s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
@merge((annotateGeneIntervals, buildUpstreamFlankBed, buildDownstreamFlankBed), os.path.join("geneset",PARAMS['interface_intergenic_bed']) )
def buildIntergenicBed( infiles, outfile ):
    ''' Genomic regions not associated with any other features'''
    inlist = " ".join(infiles)
    statement = '''cat %(inlist)s | complementBed -i stdin -g %(faidx)s 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log > %(outfile)s'''
    P.run()

###################################################
@transform((annotateGeneIntervals,buildUpstreamFlankBed,buildDownstreamFlankBed,buildIntergenicBed), suffix(".bed"), ".gtf" )
def convertBedToGtf( infile, outfile ):
    ''' convert bed files to GTF'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/bed2gff.py --as-gtf | grep -v "#" > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
## TRANSCRIPTION START SITES
@transform(addMissingEnsemblTranscripts, regex(r"all_transcripts.gtf"), os.path.join("geneset",PARAMS["interface_tss_bed"]) )
def buildTranscriptTSS( infile, outfile ):
    '''annotate transcription start sites from reference gene set.
    Similar to promotors, except that the witdth is set to 1. '''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log  
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()

###################################################
@follows( convertBedToGtf )
@files( "geneset/genes.gtf", os.path.join("geneset", PARAMS["interface_tss_gene_bed"]) )
def buildGeneTSS( infile, outfile ):
    '''create a single TSS for each gene'''
    statement = """
        cat %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
###################################################
@transform( addMissingEnsemblTranscripts, regex(r"all_transcripts.gtf"), os.path.join("geneset",PARAMS["interface_tss_gene_interval_bed"]) )
def buildGeneTSSInterval( infile, outfile ):
    '''create a single interval that encompasses all annotated TSSs for a given gene'''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | sed s/\\\\ttranscript\\\\t/\\\\texon\\\\t/g 
        | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
###################################################
@transform( addMissingNoncodingTranscripts, regex(r"(\S+).gtf"), os.path.join("geneset",PARAMS["interface_tss_gene_noncoding_bed"]) )
def buildNoncodingGeneTSS( infile, outfile ):
    '''Assign a TSS for each non-coding gene'''
    statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome_dir)s/%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=gene_id --log=%(outfile)s.log 
        | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s --log %(outfile)s.log
        | gzip
        > %(outfile)s """
    P.run()
    
###################################################
@transform( (buildTranscriptTSS, buildGeneTSS, buildGeneTSSInterval, buildNoncodingGeneTSS), suffix(".bed.gz"), ".extended.bed.gz" )
def ExtendRegion( infile, outfile ):
    '''convert bed to gtf'''
    statement = """gunzip < %(infile)s 
                   | slopBed -i stdin -g %(faidx)s -b 1000  
                   | gzip
                   > %(outfile)s """
    P.run()
    
###################################################
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
## Copy relevant files from Ensembl build
@files( "pipeline.ini", "ensembl.copy.log" )
def copyEnsembl( infile, outfile ):
    '''convert bed to gtf'''
    ensembl_dir = PARAMS["ensembl_dir"]
    statement = """cp %(ensembl_dir)s/repeats.gtf geneset/ 2> %(outfile)s;
                   cp %(ensembl_dir)s/go.tsv.gz geneset/ 2> %(outfile)s;
                   cp %(ensembl_dir)s/goslim.tsv.gz geneset/ 2> %(outfile)s;
                   cp %(ensembl_dir)s/go_ontology.obo geneset/ 2> %(outfile)s;
                   cp %(ensembl_dir)s/goslim.obo geneset/ 2> %(outfile)s;"""
    P.run()

############################################################
@files( "pipeline.ini", "ensembldb.copy.log" )
def copyEnsemblDb( infile, outfile ):
    '''copy tables from ensembl database to rnaseq database'''
    table_list = P.asList(PARAMS["ensembl_tables"])
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    query = """ATTACH "%s" as ensembl;""" % PARAMS["ensembl_db"]
    cc.execute( query )
    for table in table_list:
        cc = dbhandle.cursor()
        query = """CREATE TABLE %s AS SELECT * FROM ensembl.%s;""" % (table, table)
        print query
        cc.execute( query )
    cc.close()
    statement = """touch %(outfile)s;"""
    P.run()
    
############################################################
@transform( loadTranscripts, regex(r"(\S+).load"), "gene_name.update.log" )
def updateGeneName( infile, outfile ):
    '''update gene name from ensembl database'''
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    query = """ATTACH "%s" as ensembl;""" % PARAMS["ensembl_db"]
    cc.execute( query )
    query = """UPDATE transcript_info SET gene_name=
               (SELECT gene_name from ensembl.transcript_info 
               WHERE transcript_info.gene_id=ensembl.transcript_info.gene_id);"""
    print query
    cc.execute( query )
    cc.close()
    statement = """touch %(outfile)s;"""
    P.run()
                
############################################################
############################################################
############################################################

@follows( convertTxtToGtf, getGtfStrandedTranscripts,
          getEnsemblCodingGeneset, renameTranscripts, annotateTranscripts,
          addMissingEnsemblTranscripts, renameTranscriptToExon )
def transcripts():
    '''build all targets.'''
    pass
    
@follows( getNoncodingGeneset, renameNoncodingTranscripts,
          annotateNoncodingTranscripts, addMissingNoncodingTranscripts )
def noncoding():
    '''build all targets.'''
    pass

@follows( transcriptGtfToTab, noncodingGtfToTab, mergeTranscriptInfo, loadTranscripts )
def loadTranscriptInfo():
    '''build all targets.'''
    pass
    
@follows( buildGeneIntervals, annotateGeneIntervals, buildUpstreamFlankBed,
          buildDownstreamFlankBed, buildIntergenicBed, convertBedToGtf )
def genes():
    '''build all targets.'''
    pass  

@follows( buildTranscriptTSS, buildGeneTSS, buildNoncodingGeneTSS,
          ExtendRegion, convertToGTF )
def tss():
    '''build all targets.'''
    pass 

@follows( copyEnsembl, copyEnsemblDb, updateGeneName )
def ensembl():
    '''build all targets.'''
    pass 

@follows( transcripts, noncoding, loadTranscriptInfo, genes, tss, ensembl )
def full():
    '''build all targets.'''
    pass 


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
    
