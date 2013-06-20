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
Exome pipeline
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
HaplotypeCaller and Samtools. Variants are then annotated 
and phased.


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNVs & indels) using SAMtools and GATK
   5. Variant annotation using SNPeff and GATK VariantAnnotator
   6. Phasing using GATK PhaseByTransmission and ReadBackedPhasing
   7. Further annotation using SnpSift

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <family>-<sample>-<condition>-<replicate>.<suffix>

``family`` = "single", "trio", or "multiplex" followed by numerical identifier.  ``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
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

If submitting trios or multiplex samples then a ped file must also be supplied.  This is a tab-delimited file with
no header and one individual per line according to the following pattern:

family_id sample_id father_id mother_id sex phenotype

family_id and sample_id should correspond to <family> and <sample> in the sra/fastq filenames, father_id and mother_id
should be '0' if unknown, sex should be '1' if male and '2' if female, and phenotype should be '0' if unaffected and
'1' if affected.

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
|GATK                | 2.5-2             | local realignment, BQSR, variant calling       |
+--------------------+-------------------+------------------------------------------------+
|SNPeff              | 3.3               |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is a csvdb containing quality control information by sample and variant information by family.

Example
=======

ToDo: make exome sequencing example


Code
====

"""

# load modules
from ruffus import *
from rpy2.robjects import r as R

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
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
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Tophat as Tophat
import rpy2.robjects as ro
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Stats as Stats
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.Pipeline as P

USECLUSTER = True

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s/pipeline.ini" % os.path.splitext(__file__)[0], "../exome.ini", "exome.ini" ] )
PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
## Load target and sample data
@files( PARAMS["roi_roi"], "roi.load" )
def loadROI( infile, outfile ):
    '''Import regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    header = "chr,start,stop,feature"
    tablename = P.toTable( outfile )
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --ignore-empty
              --retry
              --header=%(header)s
              --table=%(tablename)s 
            > %(outfile)s  '''      
    P.run()

#########################################################################
@files( PARAMS["roi_to_gene"], "roi2gene.load" )
def loadROI2Gene( infile, outfile ):
    '''Import genes mapping to regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable( outfile )
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s 
            > %(outfile)s  '''      
    P.run()

#########################################################################
@files( PARAMS["samples"], "samples.load" )
def loadSamples( infile, outfile ):
    '''Import sample information into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable( outfile )
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --ignore-empty
              --retry
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
    track = P.snip( outfile, ".bam" )
    dedup_method = PARAMS["dedup_method"]
    if dedup_method == 'samtools':
        statement = '''samtools rmdup %(infile)s %(outfile)s; ''' % locals()    
    elif dedup_method == 'picard':
        PipelineMappingQC.buildPicardDuplicateStats( infile, outfile )
    statement = '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
@merge( dedup, "picard_duplicate_stats.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats( infiles, outfile )

#########################################################################
@transform( dedup, regex( r"bam/(\S+).dedup.bam"), r"bam/\1.reorder.bam")
def reorderBam(infile, outfile):
    '''Reorder BAM file using ordering of contigs in regference genome'''
    to_cluster = USECLUSTER
    statement = '''ReorderSam INPUT=%(infile)s OUTPUT=%(outfile)s REFERENCE=%%(bwa_index_dir)s/%%(genome)s.fa  VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
@transform( reorderBam, regex( r"bam/(\S+).reorder.bam"), r"bam/\1.readgroups.bam")
def addReadGroups(infile, outfile):
    '''Add read groups to read names'''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".reorder.bam" )
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    statement = '''AddOrReplaceReadGroups INPUT=%(infile)s OUTPUT=%(outfile)s RGLB=%(library)s RGPL=%(platform)s RGPU=%(platform_unit)s RGSM=%(track)s VALIDATION_STRINGENCY=SILENT; ''' % locals()
    statement += '''samtools index %(outfile)s; ''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Post-alignment QC
@transform( mapReads, regex( r"bam/(\S+).bam"),r"bam/\1.picard_stats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    PipelineMappingQC.buildPicardAlignmentStats( infile, outfile, os.path.join( PARAMS["bwa_index_dir"], PARAMS["genome"] + ".fa") )

#########################################################################
@merge( buildPicardAlignStats, "picard_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )
    
#########################################################################
@transform( dedup, regex( r"bam/(\S+).dedup.bam"), r"bam/\1.cov" )
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed file using BAMStats'''
    to_cluster = USECLUSTER
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    statement = '''CalculateHsMetrics BAIT_INTERVALS=%(baits)s TARGET_INTERVALS=%(regions)s INPUT=%(infile)s OUTPUT=%(outfile)s VALIDATION_STRINGENCY=LENIENT''' % locals()
    P.run()

#########################################################################
@merge( buildCoverageStats, "coverage_stats.load" )
def loadCoverageStats( infiles, outfile ):
    '''Import coverage statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable( outfile )
    outf = open('coverage.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(  os.path.basename(f), ".cov" )
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
                      --ignore-empty
                      --retry 
                   > %(outfile)s '''
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
    infile, realignment_intervals = infiles
    threads = PARAMS["gatk_threads"]
    statement = '''GenomeAnalysisTKLite -T IndelRealigner -o %(outfile)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s -targetIntervals %(realignment_intervals)s''' % locals()
    P.run()

#########################################################################
@transform( localRealignmentAroundIndels, regex( r"bam/(\S+).indelrealigned.bam"), r"bam/\1.recal.grp")
def countCovariates(infile, outfile):
    '''Identify covariates for base quality score realignment'''
    to_cluster = USECLUSTER
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
    threads = PARAMS["gatk_threads"]
    infile, recal = infiles
    statement = '''GenomeAnalysisTKLite -T PrintReads -o %(outfile)s -BQSR %(recal)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s ''' % locals()
    P.run()

#########################################################################                                                                                                                                                                    
@collate( bqsr, regex( r"bam/(\S+?)-(\S+).bqsr.bam" ), r"bam/\1.list" )
def listOfBAMs( infiles, outfile ):
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile+'\n')

#########################################################################
#########################################################################
#########################################################################
## Variant Calling

#########################################################################
@follows( mkdir( "variants" ) )
@transform( listOfBAMs, regex( r"bam/(\S+).list"), r"variants/\1.haplotypeCaller.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a family together'''
    to_cluster = USECLUSTER
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    statement = '''GenomeAnalysisTK -T HaplotypeCaller -o %(outfile)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s --dbsnp %(dbsnp)s -L %(intervals)s -ip %(padding)s''' % locals()
    P.run()

##########################################################################
@follows( haplotypeCaller )
@transform( bqsr, regex( r"bam/(\S+).bqsr.bam"), r"variants/\1.samtools.vcf")
def callVariantsSAMtools(infile, outfile):
    '''Perform SNV and indel calling separately for each bam using SAMtools. '''
    to_cluster = USECLUSTER
    statement = '''samtools mpileup -ugf %%(bwa_index_dir)s/%%(genome)s.fa %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()

##########################################################################
##########################################################################
##########################################################################
## Variant Annotation
@transform( haplotypeCaller, regex( r"variants/(\S+).haplotypeCaller.vcf"), r"variants/\1.haplotypeCaller.snpeff.vcf")
def annotateVariantsSNPeff( infile, outfile ):
    '''Annotate variants using SNPeff'''
    to_cluster = USECLUSTER
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''snpEff.sh eff -c %(config)s -v %(snpeff_genome)s -o gatk %(infile)s > %(outfile)s''' % locals()
    P.run()

#########################################################################
@follows( annotateVariantsSNPeff )
@transform(haplotypeCaller, regex(r"variants/(\S+).haplotypeCaller.vcf"), add_inputs(r"bam/\1.list", r"variants/\1.haplotypeCaller.snpeff.vcf"), r"variants/\1.haplotypeCaller.annotated.vcf")
def variantAnnotator( infiles, outfile ):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK -T VariantAnnotator -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamlist)s -A SnpEff --snpEffFile %(effFile)s -o %(outfile)s
                   --variant %(infile)s -L %(infile)s --dbsnp %(dbsnp)s -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest''' % locals()
    P.run()

#########################################################################
@transform(variantAnnotator, regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"), r"variants/\1.haplotypeCaller.vqsr.recal")
def variantRecalibrator( infile, outfile ):
    '''Create variant recalibration file'''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(outfile), ".recal" )
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK -T VariantRecalibrator -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s 
                   -resource:omni,known=false,training=true,truth=false,prior=12.0 %(omni)s 
                   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %(dbsnp)s 
                   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --maxGaussians 4 -percentBad 0.05 
                   -mode SNP 
                   -recalFile %(outfile)s 
                   -tranchesFile %(track)s.tranches 
                   -rscriptFile %(track)s.plots.R ''' % locals()
    P.run()

#########################################################################
@follows( variantRecalibrator )
@transform(variantAnnotator, regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"), add_inputs(r"variants/\1.haplotypeCaller.vqsr.recal", r"\1.haplotypeCaller.vqsr.tranches"), r"variants/\1.haplotypeCaller.vqsr.vcf")
def applyVariantRecalibration( infiles, outfile ):
    '''Perform variant quality score recalibration using GATK '''
    to_cluster = USECLUSTER
    infile, recal, tranches = infiles
    statement = '''GenomeAnalysisTK -T ApplyRecalibration -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   --ts_filter_level 99.0 
                   -tranchesFile %(tranches)s
                   -recalFile %(recal)s
                   -mode SNP 
                   -o %(outfile)s ''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Phasing
@transform(applyVariantRecalibration, regex(r"variants/(trio\S+|multiplex\S+).haplotypeCaller.vqsr.vcf"), add_inputs(r"\1.ped"), r"variants/\1.haplotypeCaller.pbt.vcf")
def phaseByTransmission(infiles, outfile):
    '''Infer phase based on pedigree information'''
    to_cluster = USECLUSTER
    infile, pedfile = infiles
    statement = '''GenomeAnalysisTK -T PhaseByTransmission -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s -ped %(pedfile)s -o %(outfile)s''' % locals()
    P.run()

#########################################################################
@transform(phaseByTransmission, 
           regex(r"variants/(\S+).haplotypeCaller.pbt.vcf"), 
           add_inputs(r"bam/\1.list"), 
           r"variants/\1.haplotypeCaller.rbp.vcf")
def readBackedPhasing(infiles, outfile):
    '''Infer phase based on physical information'''
    to_cluster = USECLUSTER
    infile, bamlist = infiles
    statement = '''GenomeAnalysisTK -T ReadBackedPhasing -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamlist)s -V %(infile)s -o %(outfile)s --respectPhaseInInput''' % locals()
    P.run()

######################################################################### 
@transform( (readBackedPhasing, applyVariantRecalibration), regex( r"(variants/(\S+).haplotypeCaller.rbp|variants/(single\S+).haplotypeCaller.vqsr).vcf"), r"\1.snpsift.vcf")
def annotateVariantsSNPsift( infile, outfile ):
    '''Add annotations using SNPsift'''
    to_cluster = USECLUSTER
    dbNSFP = PARAMS["annotation_snpsift_dbnsfp"]
    statement = '''SnpSift.sh dbnsfp -v %(dbNSFP)s %(infile)s > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Tabulation
@transform( annotateVariantsSNPsift, regex( r"variants/(\S+).haplotypeCaller.(\S+).snpsift.vcf"), r"variants/\1.haplotypeCaller.snpsift.table")
def vcfToTable( infile, outfile ):
    '''Convert vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -o %(outfile)s''' % locals()
    P.run()

#########################################################################
@transform( vcfToTable, regex( r"variants/(\S+).haplotypeCaller.snpsift.table"), r"variants/\1.haplotypeCaller.snpsift.table.load")
def loadVariantAnnotation( infile, outfile ):
    '''Load VCF annotations into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty > %(outfile)s''' % locals()
    P.run()
            
#########################################################################
#########################################################################
#########################################################################
## vcf statistics
@transform(annotateVariantsSNPsift, regex( r"variants/(\S+).vcf"), r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    to_cluster = USECLUSTER
    statement = '''vcf-stats %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()

#########################################################################
@merge( buildVCFstats, "vcf_stats.load" )
def loadVCFstats( infiles, outfile ):
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
def loadMetadata(): pass          

@follows( mapReads )
def mapping(): pass

@follows( dedup,
          reorderBam,
          addReadGroups,
          loadPicardDuplicateStats)

def processBAMs(): pass

@follows( buildPicardAlignStats,
          loadPicardAlignStats,
          buildCoverageStats,
          loadCoverageStats)

def postMappingQC(): pass
          
@follows( buildRealignmentTargets,
          localRealignmentAroundIndels,
          countCovariates,
          bqsr,
          listOfBAMs)
def gatk(): pass

@follows( haplotypeCaller,
          callVariantsSAMtools )
def callVariants(): pass

@follows( annotateVariantsSNPeff,
          variantAnnotator)
def annotation(): pass

@follows( variantRecalibrator,
          applyVariantRecalibration )
def vqsr(): pass

@follows( phaseByTransmission,
          readBackedPhasing,
          annotateVariantsSNPsift )
def phasing(): pass

@follows( vcfToTable,
          loadVariantAnnotation )
def tabulation(): pass

@follows( buildVCFstats,
          loadVCFstats )
def vcfstats(): pass
          
@follows( loadMetadata,
          mapping,
          processBAMs,
          postMappingQC,
          gatk,
          callVariants,
          annotation,
          vqsr,
          phasing,
          tabulation,
          vcfstats)  
def full(): pass

#########################################################################
#########################################################################
#########################################################################

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

