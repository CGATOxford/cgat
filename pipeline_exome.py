################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
====================
readqc pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

The exome pipeline imports unmapped reads from one or more 
fastq or sra files and aligns them to the genome using BWA.
Post alignment quality control is performed using Picard. 
The pipeline then performs local realignment around indels 
and base quality score recalibration using GATK.
Next variants (SNVs and indels) are called using both GATK 
Unified Genotyper and Samtools. Variants are annotated and 
then filtered by both depth/rate and regions of interest.


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNVs & indels) using SAMtools and GATK
   5. Variant annotation using SNPeff
   6. Filter variants (SAMtools / BEDTools)
   7. Calculate variant statistics (vcf-tools)

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

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



Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|Stampy              |>=0.9.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|BWA                 |                   |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|SAMtools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |filtering                                       |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|vcf-tools           |                   | VCF filtering                                  |
+--------------------+-------------------+------------------------------------------------+
|GATK                | >2.0              | local realignment, BQSR, variant calling       |
+--------------------+-------------------+------------------------------------------------+
|SNPeff              |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is a single VCF file and an HTML quality control report.

Example
=======

ToDo: make exome sequencing example


Code
====

"""

# load modules
from ruffus import *
from rpy2.robjects import r as R

import Experiment as E
import logging as L
import Database
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
import numpy, sqlite3
import GTF, IOTools, IndexedFasta
import Tophat
import rpy2.robjects as ro
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import Pipeline as P

USECLUSTER = True

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s/pipeline.ini" % os.path.splitext(__file__)[0], "../exome.ini", "exome.ini" ] )
PARAMS = P.PARAMS

def getPicardOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

def getGATKOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"    
    
#########################################################################
#########################################################################
#########################################################################
## Load target and sample data
@files( PARAMS["roi_bed"], "roi.load" )
def loadROI( infiles, outfile ):
    '''Import regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    header = "chr,start,stop,feature"
    tablename = P.toTable( outfile )
    E.info( "loading regions of interest" )
    statement = '''cat %(infiles)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --header=%(header)s
              --index=feature
              --table=%(tablename)s 
            > %(outfile)s  '''      
    P.run()

#########################################################################
@files( PARAMS["roi_to_gene"], "roi2gene.load" )
def loadROI2Gene( infiles, outfile ):
    '''Import genes mapping to regions of interest bed file into SQLite.'''

    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable( outfile )
    E.info( "loading roi to gene mapping" )
    statement = '''cat %(infiles)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=feature
              --index=gene_symbol
              --table=%(tablename)s 
            > %(outfile)s  '''      
    P.run()

#########################################################################
@files( PARAMS["samples"], "samples.load" )
def loadSamples( infiles, outfile ):
    '''Import sample information into SQLite.'''

    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable( outfile )
    E.info( "loading samples" )
    statement = '''cat %(infiles)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track
              --index=category
              --table=%(tablename)s 
            > %(outfile)s  '''      
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Alignment to a reference genome
@follows(mkdir("bam"))
@transform( ("*.fastq.1.gz", "*.fastq.gz", "*.sra"),
            regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
            r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA (output=SAM), convert to BAM, sort and index BAM file '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(outfile), ".bam" )
    m = PipelineMapping.BWA( remove_unique = PARAMS["bwa_remove_non_unique"] )
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
## BAM file processing                
#########################################################################
@transform( mapReads, regex( r"bam/(\S+).bam"), r"bam/\1.dedup.bam")
def dedup(infile, outfile):
    '''Remove PE duplicate alignments from BAM files.'''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    track = P.snip( outfile, ".bam" )
    dedup_method = PARAMS["dedup_method"]
    if dedup_method == 'samtools':
        statement = '''samtools rmdup %(infile)s %(outfile)s; ''' % locals()    
    elif dedup_method == 'picard':
        statement = '''MarkDuplicates INPUT=%(infile)s ASSUME_SORTED=true OUTPUT=%(outfile)s METRICS_FILE=%(track)s.dupstats VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
@merge( dedup, "picard_duplicate_stats.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    tablename = P.toTable( outfile )
    outf = open('dupstats.txt','w')
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.bam" )
        statfile = P.snip(f, ".bam" )  + ".dupstats"
        if not os.path.exists( statfile ): 
            E.warn( "File %s missing" % statfile )
            continue
        lines = [ x for x in open( statfile, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        outf.write( "%s\t%s" % (track,lines[1] ))
    outf.close()
    tmpfilename = outf.name
    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()

#########################################################################
@transform( dedup, regex( r"bam/(\S+).dedup.bam"), r"bam/\1.reorder.bam")
def reorderBam(infile, outfile):
    '''Reorder BAM file using ordering of contigs in regference genome'''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    statement = '''ReorderSam INPUT=%(infile)s OUTPUT=%(outfile)s REFERENCE=%%(bwa_index_dir)s/%%(genome)s.fa  VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
@transform( reorderBam, regex( r"bam/(\S+).reorder.bam"), r"bam/\1.readgroups.bam")
def addReadGroups(infile, outfile):
    '''Add read groups to read names'''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    track = P.snip( os.path.basename(infile), ".reorder.bam" )
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    sample = PARAMS["readgroup_sample"]
    statement = '''AddOrReplaceReadGroups INPUT=%(infile)s OUTPUT=%(outfile)s RGLB=%(library)s RGPL=%(platform)s RGPU=%(platform_unit)s RGSM=%(sample)s VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Post-alignment QC
@transform( mapReads, regex( r"bam/(\S+).bam"),r"bam/\1.alignstats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

#########################################################################
@merge( buildPicardAlignStats, "picard_align_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    tablename = P.toTable( outfile )
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".alignstats" )
        if not os.path.exists( f ): 
            E.warn( "File %s missing" % f )
            continue
        lines = [ x for x in open( f, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        for i in range(1, len(lines)):
            outf.write( "%s\t%s" % (track,lines[i] ))
    outf.close()
    tmpfilename = outf.name
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s '''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( dedup, regex( r"bam/(\S+).bam"), r"bam/\1.isizestats" )
def buildPicardInsertSizeStats( infile, outfile ):
    '''Gather BAM file insert size statistics using Picard '''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectInsertSizeMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa ASSUME_SORTED=true OUTPUT=%(outfile)s HISTOGRAM_FILE=%(outfile)s.pdf VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

#########################################################################
@merge( buildPicardInsertSizeStats, "picard_isize_stats.load" )
def loadPicardInsertSizeStats( infiles, outfile ):
    '''Merge Picard insert size stats into single table and load into SQLite.'''
    tablename = P.toTable( outfile )
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.isizestats" )
        if not os.path.exists( f ): 
            E.warn( "File %s missing" % f )
            continue
        lines = [ x for x in open( f, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        outf.write( "%s\t%s" % (track,lines[1] ))
    outf.close()
    tmpfilename = outf.name
    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( dedup, regex( r"bam/(\S+).dedup.bam"), r"bam/\1.cov" )
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed file using BAMStats'''
    to_cluster = USECLUSTER
    job_options = getPicardOptions()
    filename = P.snip( os.path.basename(infile), ".dedup.bam")
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    statement = '''CalculateHsMetrics BAIT_INTERVALS=%(baits)s TARGET_INTERVALS=%(regions)s INPUT=%(infile)s OUTPUT=%(outfile)s VALIDATION_STRINGENCY=LENIENT''' % locals()
    print statement
    P.run()

#########################################################################
@merge( buildCoverageStats, "coverage_stats.load" )
def loadCoverageStats( infiles, outfile ):
    '''Import coverage statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    header = "track,feature,feature_length,cov_mean,cov_median,cov_sd,cov_q1,cov_q3,cov_2_5,cov_97_5,cov_min,cov_max"
    filenames = " ".join(infiles)
    tablename = P.toTable( outfile )
    E.info( "loading coverage stats..." )
    statement = '''cat %(filenames)s | sed -e /Track/D |  sed 's/[ \t]*$//' | sed 's/,//' | sed -e 's/[ \\t]\+/\\t/g' > covstats.txt;
                   cat covstats.txt  | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --allow-empty
                       --header=%(header)s
                       --index=track
                       --index=feature
                       --table=%(tablename)s 
                   > %(outfile)s; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
## GATK

#########################################################################
@transform( addReadGroups, regex( r"bam/(\S+).readgroups.bam"), r"bam/\1.indelrealignment.intervals")
def buildRealignmentTargets(infile, outfile):
    '''Identify regions of the genome that need to be realigned'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    track = P.snip( infile, ".readgroups.bam" )
    threads = PARAMS["gatk_threads"]
    statement = '''GenomeAnalysisTKLite -T RealignerTargetCreator -o %(outfile)s --num_threads %(threads)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s''' % locals()
    P.run()

#########################################################################
@follows( buildRealignmentTargets )
@transform( addReadGroups, regex( r"bam/(\S+).readgroups.bam"), add_inputs(r"bam/\1.indelrealignment.intervals"), r"bam/\1.indelrealigned.bam")
def localRealignmentAroundIndels(infiles, outfile):
    '''Perform local realignment around indels'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    infile, realignment_intervals = infiles
    threads = PARAMS["gatk_threads"]
    statement = '''GenomeAnalysisTKLite -T IndelRealigner -o %(outfile)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s -targetIntervals %(realignment_intervals)s''' % locals()
    P.run()

#########################################################################
@transform( localRealignmentAroundIndels, regex( r"bam/(\S+).indelrealigned.bam"), r"bam/\1.recal.grp")
def countCovariates(infile, outfile):
    '''Identify covariates for base quality score realignment'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTKLite -T BaseRecalibrator --out %(outfile)s --disable_indel_quals -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s --knownSites %(dbsnp)s''' % locals()
    P.run()

#########################################################################
@follows(countCovariates)
@transform( localRealignmentAroundIndels, regex( r"bam/(\S+).indelrealigned.bam"), add_inputs(r"bam/\1.recal.grp"), r"bam/\1.bqsr.bam")
def bqsr(infiles, outfile):
    '''base quality score realignment'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    threads = PARAMS["gatk_threads"]
    infile, recal = infiles
    statement = '''GenomeAnalysisTKLite -T PrintReads -o %(outfile)s -BQSR %(recal)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s ''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Variant Calling

#########################################################################
@follows( mkdir( "variants" ) )
@transform( bqsr, regex( r"bam/(\S+).bqsr.bam"), r"variants/\1.unifiedGenotyper.vcf")
def unifiedGenotyper(infile, outfile):
    '''Call SNVs and indels using GATK UnifiedGenotyper individually for each sample'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    track = P.snip( os.path.basename(infile), ".bqsr.bam" )
    metrics_file = track + ".unifiedGenotyper.metrics"
    threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTKLite -T UnifiedGenotyper --out %(outfile)s --num_threads %(threads)s --metrics_file %(metrics_file)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s --genotype_likelihoods_model BOTH --standard_min_confidence_threshold_for_calling 30.0 --standard_min_confidence_threshold_for_emitting 30.0 --dbsnp:dbsnp,vcf %(dbsnp)s''' % locals()
    P.run()

#########################################################################
@transform(unifiedGenotyper, regex(r"variants/(\S+).unifiedGenotyper.vcf"), add_inputs(r"bam/\1.bqsr.bam"), r"variants/\1.unifiedGenotyper.annotated.vcf")
def variantAnnotator( infiles, outfile ):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    infile, bamfile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTKLite -T VariantAnnotator -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamfile)s -o %(outfile)s
                   --variant %(infile)s -L %(infile)s --dbsnp %(dbsnp)s -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest''' % locals()
    P.run()
    
#########################################################################
@follows( mkdir( "variants" ) )
@transform( bqsr, regex( r"bam/(\S+).bqsr.bam"), r"variants/\1.samtools.vcf")
def callVariantsSAMtools(infile, outfile):
    '''Perform SNV and indel calling separately for each bam using SAMtools. '''
    to_cluster = USECLUSTER
    statement = '''samtools mpileup -ugf %%(genome_dir)s/%%(genome)s.fa %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()

#########################################################################
@transform(variantAnnotator, regex(r"variants/(\S+).unifiedGenotyper.annotated.vcf"), r"variants/\1.unifiedGenotyper.annotated.vqsr.recal")
def variantRecalibrator( infile, outfile ):
    '''Create variant recalibration file'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    track = P.snip( os.path.basename(outfile), ".recal" )
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTKLite -T VariantRecalibrator -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s 
                   -resource:omni,known=false,training=true,truth=false,prior=12.0 %(omni)s 
                   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %(dbsnp)s 
                   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ 
                   -mode SNP 
                   -recalFile %(outfile)s 
                   -tranchesFile %(track)s.tranches 
                   -rscriptFile %(track)s.plots.R ''' % locals()
    P.run()

#########################################################################
@follows( variantRecalibrator )
@transform(unifiedGenotyper, regex(r"variants/(\S+).unifiedGenotyper.vcf"), add_inputs(r"\1.unifiedGenotyper.vqsr.recal", r"\1.unifiedGenotyper.vqsr.tranches"), r"variants/\1.unifiedGenotyper.vqsr.vcf")
def applyVariantRecalibration( infiles, outfile ):
    '''Perform variant quality score recalibration using GATK '''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    infile, recal, tranches = infiles
    statement = '''GenomeAnalysisTKLite -T ApplyRecalibration -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   --ts_filter_level 99.0 
                   -tranchesFile %(tranches)s
                   -recalFile %(recal)s
                   -mode SNP 
                   -o %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Variants annotation - to be added
@transform( applyVariantRecalibration, regex( r"variants/(\S+).unifiedGenotyper.vqsr.vcf"), r"variants/\1.unifiedGenotyper.snpeff.vcf")
def annotateVariantsSNPeff( infile, outfile ):
    '''Annotate variants using SNPeff'''
    to_cluster = USECLUSTER
    genome = PARAMS["annotation_snpeff_genome"]
    statement = '''snpEff.sh eff -v -o vcf %(genome)s %(infile)s > %(outfile)s'''
    P.run()

#########################################################################
@transform( annotateVariantsSNPeff, regex( r"variants/(\S+).unifiedGenotyper.snpeff.vcf"), r"variants/\1.unifiedGenotyper.snpsift.vcf")
def annotateVariantsSNPsift( infile, outfile ):
    '''Annotate variants from SnpEff using SnpSift'''
    to_cluster = USECLUSTER
    db = PARAMS["annotation_snpsift_db"]
    statement = '''SnpSift.sh dbnsfp -v -o vcf %(db)s %(infile)s > %(outfile)s'''
    P.run()
    
#########################################################################
@transform( annotateVariantsSNPeff, regex( r"variants/(\S+).unifiedGenotyper.snpeff.vcf"), r"variants/\1.unifiedGenotyper.snpeff.vcf.load")
def loadVariantAnnotationSnpEff( infile, outfile ):
    '''Load VCF annotations into database'''
    statement = ''' '''
    P.run()
            
#########################################################################
#########################################################################
#########################################################################
## filter variants
@transform( unifiedGenotyper, regex( r"variants/(\S+).unifiedGenotyper.vcf"), r"variants/\1.unifiedGenotyper.roi.vcf")
def filterVariantsROI(infile, outfile):
    '''Filter variant calls in vcf format to regions of interest from a bed file'''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".vcf")
    statement  = '''intersectBed -u -a %(infile)s -b %%(roi_bed)s > variants/%(track)s.roi.tmp 2>>variants/roi.log; 
                    (zcat %(infiles)s | grep ^#; cat variants/%(track)s.roi.tmp;) > %(outfile)s 2>>variants/roi.log; 
                    rm variants/%(track)s.roi.tmp;''' % locals()
    P.run()

#########################################################################   
@transform( filterVariantsROI, regex( r"variants/(\S+).unifiedGenotyper.roi.vcf"), r"variants/\1.unifiedGenotyper.roi.qual.vcf")
def filterVariantsQuality(infile, outfile):    
    to_cluster = USECLUSTER
    statement = '''cat %(infile)s | vcfutils.pl varFilter %%(variant_filter)s > %(outfile)s 2>>%(outfile)s.log;  ''' % locals()
    P.run()

#########################################################################   
@merge( filterVariantsQuality, "variants/merged_variants.vcf")
def mergeVCFs(infiles, outfile):
    '''Merge multiple VCF files using VCF-tools. '''
    filenames = " ".join( infiles )
    statement = '''vcf-merge %(filenames)s > %(outfile)s 2>> %(outfile)s.log; ''' % locals()
    P.run()
    
#########################################################################   
@transform(filterVariantsQuality, regex( r"variants/(\S+).qual.vcf"), r"variants/\1.qual.vcf.gz")
def indexVariants(infile, outfile):
    '''Bzip vcf file and tabix index for random access'''
    to_cluster = USECLUSTER
    statement = '''cat %(infile)s | bgzip -c > %(outfile)s; tabix -p vcf %(outfile)s;  ''' % locals() 
    P.run()
    
#########################################################################
#########################################################################
#########################################################################
## Variant statistics
@transform(filterVariantsQuality, regex( r"variants/(\S+).unifiedGenotyper.qual.vcf"), r"variants/\1.unifiedGenotyper.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    to_cluster = USECLUSTER
    statement = '''vcf-stats %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()

#########################################################################
@merge( buildVCFstats, "vcf_stats.load" )
def loadVCFStats( infiles, outfile ):
    '''Import variant statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    filenames = " ".join(infiles)
    tablename = P.toTable( outfile )
    E.info( "Loading vcf stats..." )
    statement = '''python %(scriptsdir)s/vcfstats2db.py %(filenames)s >> %(outfile)s; '''
    statement += '''cat vcfstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=vcf_stats >> %(outfile)s; '''
    statement += '''cat sharedstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=vcf_shared_stats >> %(outfile)s; '''
    statement += '''cat indelstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=indel_stats >> %(outfile)s; '''
    statement += '''cat snpstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=snp_stats >> %(outfile)s; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( loadROI,
          loadROI2Gene,
          loadSamples)
def loadMetaData(): pass          

@follows( mapReads )
def mapping(): pass

@follows( dedup,
          reorderBam,
          addReadGroups,
          loadPicardDuplicateStats)

def processBAMs(): pass

@follows( buildPicardAlignStats,
          loadPicardAlignStats,
          buildPicardInsertSizeStats,
          loadPicardInsertSizeStats,
          buildCoverageStats)
def postMappingQC(): pass
          
@follows( buildRealignmentTargets,
          localRealignmentAroundIndels,
          countCovariates,
          bqsr )
def gatk(): pass

@follows( unifiedGenotyper,
          callVariantsSAMtools )
def callVariants(): pass

@follows( variantAnnotator,
          variantRecalibrator,
          applyVariantRecalibration )
def vqsr(): pass

@follows( annotateVariantsSNPeff )
def annotateVariants(): pass

@follows( filterVariantsROI,
          filterVariantsQuality )
def filterVariants(): pass

@follows( buildVCFstats,
          loadVCFStats )
def vcfstats(): pass
          
@follows( mapping,
          processBAMs,
          postMappingQC,
          gatk,
          callVariants,
          vqsr)  
def full(): pass

#########################################################################
#########################################################################
#########################################################################

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

