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
=========================
Benchmark: RNASeq-Mapping
=========================

:Author: Andreas Heger, Nicholas Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline takes one or more fastq formatted file from an RNASeq experiment
and applies various mappers on it. The focus is here on sensitivity
versus specificity. The pipeline will not

   * rigorously time the tools 
   * do an exhaustive parameter search
   * perform validation from simulated data sets (TODO)
 
.. note::

   This pipeline is currently configured for mapping in SOLiD colourspace, 
   single end reads and mapping against a mammalian (human) genome.

Overview
=========

Mapping
-------

The pipeline aligns reads to a genome using a variety of mappers

   * tophat
   * shrimp
   * bfast
   * bwa
   * bowtie

Mappers are configured to only report "best" matches', usually defined as those
with a minimum number of mismatches. Up to ten matches can be reported per read.
If there are more than ten matches for a read, the read is designated to be unmapped.

Tuning mappers depends on the input set, notably the length and quality of reads. 

For most mappers, the variables to tune are:

   * read trimming
   * the seed length
   * the number of mismatches in the seed
   * the number of mismatches in the full alignment

Validation
----------

Unless reads have come from a simulated data set, validation will be difficult as there
is no gold standard. Instead, the methods have to be compared against each other and with
respect to a reference gene set.

Validation without a gold standard
++++++++++++++++++++++++++++++++++

In the absence of a gold standard, a number of indirect criteria can be used to decide
which mapper is best. Some simple criteria are:

   * Number of reads mapped: More mapped reads is better as it increases coverage.
   
   * Number of uniquely mapped reads: A higher proportion of uniquely mapped reads
      is better.

These criteria are confounded. For example, if a mapper tolerates more mismatches,
the number of mapped reads will increase (specificity) while the number of uniquely 
mapped reads  might drop.

Comparison between mappers
++++++++++++++++++++++++++


Using strandedness
++++++++++++++++++

If the library is stranded, sense matches in protein coding transcripts should
vastly outnumber antisense matches. A measure of specificity is thus the number
of antisense reads in protein coding exons.

Using exon boundaries
+++++++++++++++++++++

For transcript assembly, exon boundaries should be delineated as accurately as possible.
Also, spliced reads are essential to chain reads together.



Configuration
=============

Input
=====

The pipeline expects one or more ``fastq.gz`` files in the
working directory.

Optional inputs
---------------

Notes
-----

Code
----

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time

import numpy, sqlite3
import GTF, IOTools
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro

import Stats

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
import PipelineMapping
import PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" )

USECLUSTER=True

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
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
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "coding_exons.gtf.gz" )
def buildCodingExons( infile, outfile ):
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

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "coding_regions.gtf.gz" )
def buildCodingRegions( infile, outfile ):
    to_cluster = True
    statement = '''
    zcat %(infile)s 
    | awk '$2 == "protein_coding" && $3 == "CDS"'
    | perl -p -e "s/CDS/exon/" 
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log 
    | gzip 
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
## Pipeline start
###################################################################
@transform( "data*.fastq.gz", 
            suffix(".fastq.gz" ),
            ".bfast.sam.gz" )
def mapReadsWithBFAST( infile, outfile ):
    '''map reads with bfast'''

    to_cluster = True
    job_options= "-pe dedicated %i -R y -l mem_free=16G" % PARAMS["bfast_threads"] 

    update = False
    if not os.path.exists( "%s.bmf" % outfile ):
        statement = '''
        bfast match -f %(bfast_genome_dir)s/%(genome)s.fa 
                    -A 1 
                    -n %(bfast_threads)i
                    -r <(gunzip < %(infile)s) 
                    -t 
                    %(bfast_match_options)s
                    > %(outfile)s.bmf 
                    2> %(outfile)s.bmf.log
        '''
        update = True
        P.run()

    if not os.path.exists( "%s.baf" % outfile ) or update:
        statement = '''
        bfast localalign -f %(bfast_genome_dir)s/%(genome)s.fa 
                    -A 1 
                    -n %(bfast_threads)i
                    -m %(outfile)s.bmf
                    %(bfast_align_options)s
                    > %(outfile)s.baf 
                    2> %(outfile)s.baf.log
        '''
        P.run()
    
    statement = '''
    bfast postprocess -f %(bfast_genome_dir)s/%(genome)s.fa 
                %(bfast_postprocess_options)s
               -A 1 
               -n %(bfast_threads)i
               -i %(outfile)s.baf 
    2> %(outfile)s.post.log 
    | python %(scriptsdir)s/bam2bam.py --sam --unset-unmapped-mapq --set-nh --log=%(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
@transform( "data*.fastq.gz", 
            suffix(".fastq.gz" ),
            ".shrimp.sam.gz" )
def mapReadsWithShrimp( infile, outfile ):
    '''map reads with shrimp'''

    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=64G" % PARAMS["shrimp_threads"] 

    statement = '''
    gmapper-cs --full-threshold 80%% --threads %(shrimp_threads)i --fastq --report 5 --sam
              --sam-unaligned
              %(shrimp_options)s
              %(infile)s 
              %(genome_dir)s/%(genome)s.fa 
    2> %(outfile)s.log
    | gzip 
    > %(outfile)s 
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( "data*.fastq.gz", 
            suffix(".fastq.gz" ),
            ".novoalign.sam.gz" )
def mapReadsWithNovoalign( infile, outfile ):
    '''map reads with shrimp'''

    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=64G" % PARAMS["novoalign_threads"] 

    statement = '''
    novoalignCS 
              -c %(novoalign_threads)s
              -d %(novoalign_genome_dir)s/%(genome)s_cs.ncx
              -f %(infile)s 
              -F BFASTQ 
              -o SAM
              %(novoalign_options)s
    2> %(outfile)s.log
    | gzip 
    > %(outfile)s 
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( "data*.bwafastq.gz", 
            suffix(".bwafastq.gz" ),
            ".bwa.sam.gz" )
def mapReadsWithBWA( infile, outfile ):
    '''map reads with shrimp'''

    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=64G" % PARAMS["bwa_threads"] 

    statement = '''
    bwa aln -t %(bwa_threads)s -c %(bwa_align_options)s %(bwa_genome_dir)s/%(genome)s_cs %(infile)s > %(outfile)s.sai
    '''
    P.run()

    statement = '''
    bwa samse %(bwa_samse_options)s %(bwa_genome_dir)s/%(genome)s_cs %(outfile)s.sai %(infile)s 
    | python %(scriptsdir)s/bam2bam.py --sam --set-nh --unset-unmapped-mapq --log=%(outfile)s.log
    | gzip 
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( "*.fastq.gz", 
            suffix(".fastq.gz" ),
            ".tophat.bam" )
def mapReadsWithTophat( infile, outfile ):
    '''map reads with tophat

    '''

    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=16G" % PARAMS["tophat_threads"] 

    tmpfile = P.getTempFilename( "." )

    #qualfile = P.snip(infile, "csfasta.gz" ) + "qual.gz"
    '''
    gunzip < %(infile)s > %(tmpfile)s.csfasta;
    checkpoint;
    gunzip < %(qualfile)s > %(tmpfile)s.qual;
    checkpoint;
    '''
    
    statement = '''
    zcat %(infile)s 
    | python %(scriptsdir)s/fastq2solid.py 
           --change-format=integer
           --pattern="%(tmpfile)s.%%s" >& %(outfile)s.log;
    checkpoint;
    tophat --output-dir %(outfile)s.dir                    
           --num-threads %(tophat_threads)s  
           --library-type %(tophat_library_type)s
           --color
           --quals
           --integer-quals
           %(tophat_options)s
           %(tophat_genome_dir)s/%(genome)s_cs
           %(tmpfile)s.csfasta %(tmpfile)s.qual
           >> %(outfile)s.log;
    checkpoint;
    mv %(outfile)s.dir/accepted_hits.bam %(outfile)s;
    checkpoint;
    samtools index %(outfile)s;
    checkpoint;
    rm -f %(tmpfile)s.csfasta %(tmpfile)s.qual
    '''

    P.run()

    os.unlink( tmpfile )

###################################################################
###################################################################
###################################################################
@transform( ("data*.fastq.gz", "strictdata.fastq.gz", "permissive.fastq.gz"),
            suffix(".fastq.gz" ),
            ".bowtie.sam.gz" )
def mapReadsWithBowtie( infile, outfile ):
    '''map reads with bowtie'''

    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=16G" % PARAMS["bowtie_threads"] 

    tmpfile = P.getTempFilename()
    
    statement = '''
    gunzip < %(infile)s > %(tmpfile)s;
    checkpoint;
    bowtie -q
           --sam 
           -C
           --threads %(bowtie_threads)s
           %(bowtie_options)s
           %(bowtie_genome_dir)s/%(genome)s_cs
           %(tmpfile)s
    | python %(scriptsdir)s/bam2bam.py --sam --set-nh --log=%(outfile)s.log
    | gzip
    > %(outfile)s;
    checkpoint;
    rm -f %(tmpfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( (mapReadsWithBowtie, 
             mapReadsWithShrimp,
             mapReadsWithNovoalign,
             mapReadsWithBFAST,
             mapReadsWithBWA ),
            suffix(".sam.gz"),
            add_inputs( os.path.join(PARAMS["annotations_dir"],
                                     PARAMS_ANNOTATIONS["interface_contigs"] ) ),
            ".bam" )
def buildBAMs( infiles, outfile ):
    '''build BAM files.'''

    infile, contigs = infiles
    to_cluster = USECLUSTER
    
    track = P.snip( outfile, ".bam" )

    # samtools sort does not set sorted flag, do it yourself
    statement = '''
    gunzip < %(infile)s
    | perl -p -e "if (/^\\@HD/) { s/\\bSO:\S+/\\bSO:coordinate/}"  
    | samtools import %(contigs)s - -
    | samtools sort - %(track)s;
    checkpoint;
    samtools index %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
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
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".fa")
def buildReferenceTranscriptome( infile, outfile ):
    '''build reference transcriptome. 

    The reference transcriptome contains all known 
    protein coding transcripts.
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

###################################################################
###################################################################
###################################################################
@transform( ("data*.fastq.gz" ),
            suffix(".fastq.gz" ),
            add_inputs( buildReferenceTranscriptome, 
                        os.path.join(PARAMS["annotations_dir"],
                                     PARAMS_ANNOTATIONS["interface_contigs"] ) ),
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
    to_cluster = USECLUSTER
    job_options= "-pe dedicated %i -R y -l mem_free=16G" % PARAMS["bowtie_threads"] 

    tmpfile = P.getTempFilename()

    infile, reffile, contigs = infiles
    track = P.snip( outfile, ".bam" )
    prefix = P.snip( reffile, ".fa" )
    
    statement = '''
    gunzip < %(infile)s > %(tmpfile)s;
    checkpoint;
    bowtie -q
           --sam 
           -C
           --threads %(bowtie_threads)s
           -v 2 --best --strata -a
           %(prefix)s_cs
           %(tmpfile)s
    | python %(scriptsdir)s/bam2bam.py --sam --set-nh --log=%(outfile)s.log
    | perl -p -e "if (/^\\@HD/) { s/\\bSO:\S+/\\bSO:coordinate/}"  
    | samtools import %(contigs)s - -
    | samtools sort - %(track)s;
    checkpoint;
    samtools index %(outfile)s
    checkpoint;
    rm -f %(tmpfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( (mapReadsWithTophat, buildBAMs),
            suffix(".bam"),
            ".bam" )
def mapping( infile, outfile): pass

###################################################################
###################################################################
###################################################################
@transform( mapping, 
            suffix(".bam"),
            add_inputs( buildCodingGeneSet, mapReadsWithBowtieAgainstTranscriptome),
            ".accepted.bam" )
def checkMappedReadsAgainstTranscriptome( infiles, outfile): 
    '''reconcile genomic and transcriptome matches.
    '''
    genome, reffile, transcriptome = infiles
    outfile_mismapped = P.snip(outfile, ".accepted.bam") + ".mismapped.bam"

    to_cluster = USECLUSTER

    options = []
    if "tophat_unique" in PARAMS and PARAMS["tophat_unique"]:
        options.append( "--unique" )

    if "tophat_remove_contigs" in PARAMS and PARAMS["tophat_remove_contigs"]:
        options.append( "--remove-contigs=%s" % PARAMS["tophat_remove_contigs"] )
        
    # if "tophat_protocol" in PARAMS and "stranded" in PARAMS["tophat_protocol"]:
    #     options.append( "--set-strand" )

    options = " ".join(options)

    # using mismatches to filter is inappropriate as mismatches are not handled in the same
    # way across all mappers:
    # * not all mappers have NM, for example (bwa does not) and 
    # * not all have CM (tophat does not). 

    statement = '''
    python %(scriptsdir)s/rnaseq_bams2bam.py 
       --force
       --filename-gtf=%(reffile)s
       --filename-mismapped=%(outfile_mismapped)s
       --ignore-mismatches
       %(options)s
       %(transcriptome)s %(genome)s %(outfile)s
    > %(outfile)s.log;
    checkpoint;
    samtools index %(outfile)s 2>&1 >> %(outfile)s.log;
    samtools index %(outfile_mismapped)s 2>&1 >> %(outfile)s.log;
    '''
    P.run()

@merge( checkMappedReadsAgainstTranscriptome, "transcriptome_validation.load" )
def loadTranscriptomeValidation( infiles, outfile ):
    '''load transcriptome validation data into database.'''

    to_cluster = USECLUSTER

    headers = ",".join( [P.quote( P.snip( x, ".accepted.bam")) for x in infiles ] )
    infiles = " ".join( ["%s.log" % x for x in infiles ] )
    
    tablename = P.toTable( outfile )

    statement = '''
    python %(scriptsdir)s/combine_tables.py 
         --headers=%(headers)s
         %(infiles)s
    | python %(scriptsdir)s/table2table.py --transpose
    | perl -p -e "s/bin/track/"
    | python %(scriptsdir)s/csv2db.py
         --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "fastqc" ) ) )
@transform( mapping, suffix( ".bam"), ".fastqc" )
def buildFastQCReport( infile, outfile ):
    '''run fastqc on aligned reads.'''
    
    to_cluster = USECLUSTER
    statement = '''fastqc --outdir=%(exportdir)s/fastqc %(infile)s >& %(outfile)s''' 
    P.run()

###################################################################
###################################################################
###################################################################
@transform("*.fastq.gz", suffix(".fastq.gz"), ".nreads")
def countReads( infile, outfile ):
    '''count number of reads in fastq file.'''
    statement = '''
    echo "`zcat %(infile)s | grep -v "#" | wc -l` / 4" | bc > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
@follows( countReads )
@transform( mapping,
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''
    to_cluster = USECLUSTER

    def __getReads( fn ):
        return int(open(fn).readlines()[0])

    track = infile[:infile.index(".")]
    nreads = __getReads( track + ".nreads" )

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --force-output
         --output-filename-pattern=%(outfile)s.%%s
         --input=%(nreads)i
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

    header = ",".join( [ P.quote( P.snip( x, ".readstats")) for x in infiles] )
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
                      --table=%(tname)s 
                >> %(outfile)s
                """
    
        P.run()

###################################################################
###################################################################
###################################################################
@transform( mapping,
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

###################################################################
###################################################################
###################################################################
@transform("*.fastq.gz",
           suffix(".fastq.gz"),
           ".readstats.tsv.gz" )
def buildReadStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/fastq2table.py
         --log=%(outfile)s.log
    | %(scriptsdir)s/hsort 1
    | gzip
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@merge( mapping,
        "read_correspondence.tsv.gz" )
def buildReadCorrespondence( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    
    headers = ",".join([ P.snip( x, ".bam") for x in infiles  ])
    sorters = " ".join([ "<( samtools view -h %s | %s/hsort 0 )" % (x, PARAMS["scriptsdir"]) for x in infiles ] )

    statement = '''
    python %(scriptsdir)s/rnaseq_bams_vs_bams.py
         --headers=%(headers)s
         --log=%(outfile)s.log
       %(sorters)s
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( buildReadStats )
@merge( (buildReadCorrespondence, "data.readstats.tsv.gz" ), "read_correspondence.load" )
def loadReadCorrespondence( infiles, outfile ):
    '''load read correspondence data into database.'''

    to_cluster = USECLUSTER

    infiles = " ".join(infiles)

    tablename = P.toTable( outfile )

    statement = '''
    python %(scriptsdir)s/combine_tables.py 
         %(infiles)s
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tablename)s 
                > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def mergeAndLoad( infiles, outfile, suffix ):
    '''load categorical tables (two columns) into a database.

    The tables are merged and entered row-wise.

    '''
    header = ",".join( [ P.quote( P.snip( x, suffix)) for x in infiles] )
    if suffix.endswith(".gz"):
        filenames = " ".join( [ "<( zcat %s | cut -f 1,2 )" % x for x in infiles ] )
    else:
        filenames = " ".join( [ "<( cat %s | cut -f 1,2 )" % x for x in infiles ] )

    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/" 
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
@merge( buildExonValidation, "exon_validation.load" )
def loadExonValidation( infiles, outfile ):
    '''merge alignment stats into single tables.'''
    suffix = suffix = ".exon.validation.tsv.gz" 
    mergeAndLoad( infiles, outfile, suffix = suffix )
    for infile in infiles:
        track = P.snip( infile, suffix )
        o = "%s_overrun.load" % track 
        P.load( infile + ".overrun.gz", o )

###################################################################
###################################################################
###################################################################
@transform( mapping,
            suffix(".bam"),
            add_inputs( buildCodingExons ),
            ".exon.coverage.tsv.gz" )
def buildExonCoverage( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    infile, exons = infiles
    statement = '''zcat < %(exons)s
    | python %(scriptsdir)s/gtf2table.py
         --reporter=genes
         --bam-file=%(infile)s
         --counter=read-coverage
         --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( mapping,
            suffix(".bam"),
            add_inputs( buildCodingRegions ),
            ".region.coverage.tsv.gz" )
def buildRegionCoverage( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = USECLUSTER
    infile, exons = infiles
    statement = '''zcat < %(exons)s
    | python %(scriptsdir)s/gtf2table.py
         --reporter=genes
         --bam-file=%(infile)s
         --counter=read-coverage
         --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
@transform( (buildExonCoverage, 
             buildRegionCoverage),
            suffix(".coverage.tsv.gz"), "_coverage.load" )
def loadCoverage( infile, outfile ):
    P.load( infile, outfile, options = "--index=gene_id" )

############################################################
############################################################
############################################################
@transform( mapping, 
            suffix(".bam" ), ".bam.stats")
def buildAlignmentStats( infile, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    to_cluster = USECLUSTER
    
    # replace the SO field from tophat/samtools with coordinate to indicate
    # that the file is sorted by coordinate.
    # naturally - the bam files should be sorted.

    # remove all columns without quality scores - picard chokes
    # on these (for some reason, many reads contain no quality scores
    # when mapping with tophat)
    statement = '''
    java -Xmx2g net.sf.picard.analysis.CollectMultipleMetrics
            I=<(samtools view -h %(infile)s | awk '$11 != "*"')
            O=%(outfile)s 
            R=%(genome_dir)s/%(genome)s.fasta
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

@follows( loadBAMStats, 
          loadCoverage, 
          buildFastQCReport, 
          loadAlignmentStats,
          loadTranscriptomeValidation,
          loadExonValidation,
          loadReadCorrespondence )
def validate(): pass

@follows( validate )
def full(): pass

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
    
    mapping_stats
    bam_stats
    '''

    tablename = P.toTable( outfile )
    # can not create views across multiple database, so use table
    view_type = "TABLE"
    
    dbhandle = connect()
    Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

    statement = '''
    CREATE %(view_type)s %(tablename)s AS
    SELECT *
    FROM bam_stats AS b
    '''

    Database.executewait( dbhandle, statement % locals() )

###################################################################
###################################################################
###################################################################
@follows( createViewMapping )
def views():
    pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows( views, mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( views, mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

@follows() 
def publish():
    '''publish files.'''
    P.publish_report()

###################################################################
###################################################################
###################################################################
## Epilog
###################################################################
if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
