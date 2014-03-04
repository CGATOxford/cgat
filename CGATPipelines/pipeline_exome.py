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
   4. Variant calling (SNVs & indels) using GATK
   5. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   6. For trios only, phasing genotypes using GATK and calling de novos using GATK and denovogear
   7. Flags variants within genes of interest (such as known disease genes) (GATK) (optional)

To do:
   1. Allow users to add other training sets for variant quality score recalibration
   2. Allow users to add annotations using SnpSift
   3. Add job options to all functions
   4. Add report functions

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

If you are submitting families then a .ped file for each family must be supplied within your working firectory. 
This is a tab-delimited file named <family>.ped (where <family> is the family ID in the title of the corresponding 
fastq files) with no header and one individual per line according to the following pattern:

family_id sample_id father_id mother_id sex phenotype

family_id and sample_id should correspond to <family> and <family>-<sample> in the sra/fastq filenames, father_id and 
mother_id should be '0' if unknown, sex should be '1' if male, '2' if female and '0' if unknown, and phenotype 
should be '1' if unaffected, '2' if affected and '0' if unknown.

If you are running the functions to look for compound heterozygotes in Multiplex families then there is a further
requirement for the .ped files.  The phasing tools expect a trio and therefore any other family members (other than
parents and one child) must be labelled as unrelated.  That is, the first additional family member could be labelled
"family0" in the family_id column, and subsequent additional family members could be "family1", "family2" and so on.

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf, make add_genes_of_interest=1 (default=0) and
provide a list of comma separated genes (without spaces) to the ini file.

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
import csv
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
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0], "../exome.ini", "exome.ini"])
PARAMS = P.PARAMS


def getPicardOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"


def getGATKOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

#########################################################################
#########################################################################
#########################################################################
# Load target and sample data


@files(PARAMS["roi_bed"], "roi.load")
def loadROI(infile, outfile):
    '''Import regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    header = "chr,start,stop,feature"
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --ignore-empty
              --retry
              --header=%(header)s
              --table=%(tablename)s 
            > %(outfile)s  '''
    P.run()

#########################################################################


@files(PARAMS["roi_to_gene"], "roi2gene.load")
def loadROI2Gene(infile, outfile):
    '''Import genes mapping to regions of interest bed file into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --ignore-empty
              --retry
              --table=%(tablename)s 
            > %(outfile)s  '''
    P.run()

#########################################################################


@files(PARAMS["samples"], "samples.load")
def loadSamples(infile, outfile):
    '''Import sample information into SQLite.'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
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
# Alignment to a reference genome


@collate(("*.fastq.1.gz", "*.fastq.2.gz", "*.fastq.gz", "*.sra"),
         regex(r"(\S+-\S+)-(\S+)-(\S+).(fastq.1.gz|fastq.2.gz|fastq.gz|sra)"),
         r"\1.\4")
def mergeFastqs(infiles, outfile):
    '''merge fastqs by library and read pair prior to mapping'''
    to_cluster = USECLUSTER
    inputfiles = " ".join(infiles)
    statement = '''zcat %(inputfiles)s | gzip > %(outfile)s '''
    P.run()

#########################################################################


@follows(mkdir("bam"))
@transform(mergeFastqs,
           regex(r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
           r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA (output=SAM), convert to BAM, sort and index BAM file '''
    to_cluster = USECLUSTER
    job_options = "-pe dedicated 2 -l mem_free=8G"
    track = P.snip(os.path.basename(outfile), ".bam")
    m = PipelineMapping.BWA(
        remove_unique=PARAMS["bwa_remove_non_unique"], align_stats=True, dedup=True)
    statement = m.build((infiles,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################
# BAM file processing
#########################################################################


@merge(mapReads, "picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)

#########################################################################
#########################################################################
#########################################################################
# Post-alignment QC
#########################################################################


@follows(mapReads)
@merge("bam/*.picard_stats", "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)

#########################################################################


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed file using BAMStats'''
    to_cluster = USECLUSTER
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    statement = '''CalculateHsMetrics BAIT_INTERVALS=%(baits)s TARGET_INTERVALS=%(regions)s INPUT=%(infile)s OUTPUT=%(outfile)s VALIDATION_STRINGENCY=LENIENT''' % locals(
    )
    P.run()

#########################################################################


@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    '''Import coverage statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    outf = open('coverage.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".cov")
        lines = [
            x for x in open(f, "r").readlines() if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
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
# GATK


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.bqsr.bam")
def GATKpreprocessing(infile, outfile):
    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    statement = '''ReorderSam INPUT=%(infile)s OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam REFERENCE=%%(bwa_index_dir)s/%%(genome)s.fa ALLOW_INCOMPLETE_DICT_CONCORDANCE=true VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()
    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ; checkpoint ;''' % locals()
    statement += '''AddOrReplaceReadGroups INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam OUTPUT=%(tmpdir_gatk)s/%(track)s.readgroups.bam RGLB=%(library)s RGPL=%(platform)s RGPU=%(platform_unit)s RGSM=%(track)s VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals(
    )
    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.readgroups.bam ; checkpoint ;''' % locals()
    statement += '''GenomeAnalysisTK -T RealignerTargetCreator -o %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals --num_threads %(threads)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(tmpdir_gatk)s/%(track)s.readgroups.bam ; checkpoint ;''' % locals(
    )
    statement += '''GenomeAnalysisTK -T IndelRealigner -o %(tmpdir_gatk)s/%(track)s.indelrealigned.bam -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(tmpdir_gatk)s/%(track)s.readgroups.bam -targetIntervals %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals ; checkpoint ;''' % locals(
    )
    statement += '''GenomeAnalysisTK -T BaseRecalibrator --out %(tmpdir_gatk)s/%(track)s.recal.grp -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(tmpdir_gatk)s/%(track)s.indelrealigned.bam --knownSites %(dbsnp)s %(solid_options)s ; checkpoint ;''' % locals(
    )
    statement += '''GenomeAnalysisTK -T PrintReads -o %(outfile)s -BQSR %(tmpdir_gatk)s/%(track)s.recal.grp -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(tmpdir_gatk)s/%(track)s.indelrealigned.bam ; checkpoint ;''' % locals(
    )
    statement += '''rm -rf %(tmpdir_gatk)s ;'''
    P.run()

#########################################################################


@collate(GATKpreprocessing, regex(r"bam/(\S+?)-(\S+).bqsr.bam"), r"bam/\1.list")
def listOfBAMs(infiles, outfile):
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')

#########################################################################
#########################################################################
#########################################################################
# Variant Calling

#########################################################################


@follows(mkdir("variants"))
@transform(listOfBAMs, regex(r"bam/(\S+).list"), r"variants/\1.haplotypeCaller.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a family together'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    hc_options = PARAMS["gatk_hc_options"]
    statement = '''GenomeAnalysisTK -T HaplotypeCaller -o %(outfile)s -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(infile)s --dbsnp %(dbsnp)s -L %(intervals)s -ip %(padding)s''' % locals(
    )
    P.run()

##########################################################################
##########################################################################
##########################################################################
# Variant Annotation and Recalibration


@transform(haplotypeCaller, regex(r"variants/(\S+).haplotypeCaller.vcf"), r"variants/\1.haplotypeCaller.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    to_cluster = USECLUSTER
    job_options = "-pe dedicated 4 -R y -l mem_free=6G"
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''snpEff.sh eff -c %(config)s -v %(snpeff_genome)s -o gatk %(infile)s > %(outfile)s''' % locals()
    P.run()

#########################################################################


@follows(annotateVariantsSNPeff)
@transform(haplotypeCaller, regex(r"variants/(\S+).haplotypeCaller.vcf"), add_inputs(r"bam/\1.list", r"variants/\1.haplotypeCaller.snpeff.vcf"), r"variants/\1.haplotypeCaller.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK -T VariantAnnotator -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamlist)s -A SnpEff --snpEffFile %(effFile)s -o %(outfile)s
                   --variant %(infile)s -L %(infile)s --dbsnp %(dbsnp)s -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest -A AlleleBalanceBySample''' % locals()
    P.run()

#########################################################################


@transform(variantAnnotator, regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"), r"variants/\1.haplotypeCaller.vqsr.recal")
def variantRecalibrator(infile, outfile):
    '''Create variant recalibration file'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    track = P.snip(os.path.basename(outfile), ".recal")
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK -T VariantRecalibrator -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s 
                   -resource:omni,known=false,training=true,truth=false,prior=12.0 %(omni)s 
                   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %(dbsnp)s 
                   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --maxGaussians 4 --numBadVariants 3000
                   -mode SNP 
                   -recalFile %(outfile)s 
                   -tranchesFile %(track)s.tranches 
                   -rscriptFile %(track)s.plots.R ''' % locals()
    P.run()

#########################################################################


@follows(variantRecalibrator)
@transform(variantAnnotator, regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"), add_inputs(r"variants/\1.haplotypeCaller.vqsr.recal", r"\1.haplotypeCaller.vqsr.tranches"), r"variants/\1.haplotypeCaller.vqsr.vcf")
def applyVariantRecalibration(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    infile, recal, tranches = infiles
    statement = '''GenomeAnalysisTK -T ApplyRecalibration -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
                   --ts_filter_level 99.0 
                   -tranchesFile %(tranches)s
                   -recalFile %(recal)s
                   -mode SNP 
                   -o %(outfile)s ''' % locals()
    P.run()

#########################################################################


@transform(applyVariantRecalibration, regex(r"variants/(\S+).haplotypeCaller.vqsr.vcf"), r"variants/\1.haplotypeCaller.snpsift.vcf")
def annotateVariantsSNPsift(infile, outfile):
    '''Add annotations using SNPsift'''
    to_cluster = USECLUSTER
    job_options = "-pe dedicated 4 -R y -l mem_free=6G"
    track = P.snip(os.path.basename(infile), ".vqsr.vcf")
    dbNSFP = PARAMS["annotation_snpsift_dbnsfp"]
#    statement = '''SnpSift.sh geneSets -v /ifs/projects/proj016/data/1000Genomes/msigdb.v4.0.symbols.gmt %(infile)s > variants/%(track)s_temp1.vcf; checkpoint;''' % locals()
    statement = '''SnpSift.sh dbnsfp -v %(dbNSFP)s %(infile)s > variants/%(track)s_temp1.vcf; checkpoint;''' % locals()
    statement += '''SnpSift.sh annotate /ifs/projects/proj016/data/1000Genomes/00-All.vcf variants/%(track)s_temp1.vcf > %(outfile)s ;''' % locals(
    )
#    statement += '''SnpSift.sh annotate /ifs/projects/proj016/data/1000Genomes/1000GENOMES-phase_1_EUR_new.vcf variants/%(track)s_temp.vcf > %(track)s_temp_new.vcf; checkpoint;''' % locals()
# statement += '''sed 's/##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">/##INFO=<ID=AF_EUR,Number=A,Type=Float,Description="Allele Frequency">/g' %(track)s_temp_new.vcf > %(outfile)s; checkpoint;''' % locals()
#    statement += '''rm -f variants/*temp*vcf;'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Genes of interest


@active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
@transform((annotateVariantsSNPsift), regex(r"variants/(\S+).haplotypeCaller.snpsift.vcf"), r"variants/\1.genes.vcf")
def findGenes(infile, outfile):
    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf if variant is within a gene of interest as defined in the ini file'''
    to_cluster = USECLUSTER
    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
    statement = '''GenomeAnalysisTK -T VariantFiltration -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'" --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Tabulation
CALIBRATION = {0: annotateVariantsSNPsift,
               1: findGenes}


@transform(CALIBRATION[PARAMS["annotation_add_genes_of_interest"]], regex(r"variants/((\S+).haplotypeCaller.snpsift|(\S+).genes).vcf"), r"variants/\1.table")
def vcfToTable(infile, outfile):
    '''Convert vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    # add AB and DP annotations to table
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(vcfToTable, regex(r"variants/(\S+).table"), r"variants/\1.table.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# De novos


@transform(annotateVariantsSNPsift, regex(r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"), add_inputs(r"\1.ped"), r"variants/\1.filtered.vcf")
def deNovoVariants(infiles, outfile):
    '''Filter variants based on provided jexl expression'''
    to_cluster = USECLUSTER
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    # need to filter on AB annotation but this is not present for all variants
    # - will try to fix this later
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -select 'vc.getGenotype("%(father)s").getDP()>=10&&vc.getGenotype("%(mother)s").getDP()>=10&&vc.getGenotype("%(child)s").getPL().0>20&&vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(child)s").getPL().2>0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(father)s").getPL().1>20&&vc.getGenotype("%(father)s").getPL().2>20&&vc.getGenotype("%(mother)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().1>20&&vc.getGenotype("%(mother)s").getPL().2>20&&vc.getGenotype("%(child)s").getAD().1>=3&&((vc.getGenotype("%(child)s").getAD().1)/(vc.getGenotype("%(child)s").getDP().floatValue()))>=0.25&&(vc.getGenotype("%(father)s").getAD().1==0||(vc.getGenotype("%(father)s").getAD().1>0&&((vc.getGenotype("%(father)s").getAD().1)/(vc.getGenotype("%(father)s").getDP().floatValue()))<0.05))&&(vc.getGenotype("%(mother)s").getAD().1==0||(vc.getGenotype("%(mother)s").getAD().1>0&&((vc.getGenotype("%(mother)s").getAD().1)/(vc.getGenotype("%(mother)s").getDP().floatValue()))<0.05))&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' > %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(deNovoVariants, regex(r"variants/(\S+).filtered.vcf"), r"variants/\1.filtered.table")
def tabulateDeNovos(infile, outfile):
    '''Tabulate filtered variants'''
    to_cluster = USECLUSTER
    # add AB and DP annotations to table
    statement = '''awk 'NR > 16' %(infile)s | sed '/^INFO/d' > %(infile)s.tmp ;'''
    statement += '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s.tmp --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s ;''' % locals(
    )
    statement += '''rm -f %(infile)s.tmp*'''
    P.run()

#########################################################################


@transform(tabulateDeNovos, regex(r"variants/(\S+).filtered.table"), r"variants/\1.filtered.table.load")
def loadDeNovos(infile, outfile):
    '''Load de novos into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################


@transform(annotateVariantsSNPsift, regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"), add_inputs(r"\1.ped"), r"variants/\1.denovos.vcf")
def lowerStringencyDeNovos(infiles, outfile):
    '''Filter variants based on provided jexl expression'''
    to_cluster = USECLUSTER
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    # need to filter on AB annotation but this is not present for all variants
    # - will try to fix this later
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -select 'vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().0==0&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' > %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(lowerStringencyDeNovos, regex(r"variants/(\S+).denovos.vcf"), r"variants/\1.denovos.table")
def tabulateLowerStringencyDeNovos(infile, outfile):
    '''Tabulate filtered variants'''
    to_cluster = USECLUSTER
    # add AB and DP annotations to table
    statement = '''awk 'NR > 16' %(infile)s | sed '/^INFO/d' > %(infile)s.tmp ;'''
    statement += '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s.tmp --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s ;''' % locals(
    )
    statement += '''rm -f %(infile)s.tmp*'''
    P.run()

#########################################################################


@transform(tabulateLowerStringencyDeNovos, regex(r"variants/(\S+).denovos.table"), r"variants/\1.denovos.table.load")
def loadLowerStringencyDeNovos(infile, outfile):
    '''Load de novos into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################
<<<<<<< HEAD
#########################################################################
#########################################################################
# Dominant


@transform(annotateVariantsSNPsift, regex(r"variants/(\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"), add_inputs(r"\1.ped"), r"variants/\1.dominant.vcf")
def dominantVariants(infiles, outfile):
    '''Filter variants according to autosomal dominant disease model'''
    to_cluster = USECLUSTER
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
        if row['status'] == '1':
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
            '").isHomRef()'
    # for some weird reason the 1000G filter doesn't work on these particular
    # files - will add later when I've figured out what's wrong
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -o %(outfile)s -select 'vc.getGenotype("%(affecteds_exp)s").getPL().1==0%(unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' ;''' % locals(
    )
    P.run()

#########################################################################


@transform(dominantVariants, regex(r"variants/(\S+).dominant.vcf"), r"variants/\1.dominant.table")
def tabulateDoms(infile, outfile):
    '''Tabulate filtered variants'''
    to_cluster = USECLUSTER
    # add AB and DP annotations to table
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant:VCF %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(tabulateDoms, regex(r"variants/(\S+).dominant.table"), r"variants/\1.dominant.table.load")
def loadDoms(infile, outfile):
    '''Load dominant disease candidates into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Recessive


@transform(annotateVariantsSNPsift, regex(r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"), add_inputs(r"\1.ped"), r"variants/\1.recessive.vcf")
def recessiveVariants(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    # need to add calling compound hets - will try using phased vcf and Gemini
    to_cluster = USECLUSTER
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    parents = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
            parents += [row['father'], row['mother']]
    affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
    if len(parents) == 0:
        parents_exp = ''
    else:
        parents_exp = '&&vc.getGenotype("' + \
            ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
            '").getPL().1==0'
    # need a way of specifying that other unaffecteds eg. sibs can't be
    # homozygous for alt allele
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -o %(outfile)s -select 'vc.getGenotype("%(affecteds_exp)s").getPL().2==0%(parents_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' ;''' % locals(
    )
    P.run()

#########################################################################


@transform(recessiveVariants, regex(r"variants/(\S+).recessive.vcf"), r"variants/\1.recessive.table")
def tabulateRecs(infile, outfile):
    '''Tabulate filtered variants'''
    to_cluster = USECLUSTER
    # add AB and DP annotations to table
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant:VCF %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(tabulateRecs, regex(r"variants/(\S+).recessive.table"), r"variants/\1.recessive.table.load")
def loadRecs(infile, outfile):
    '''Load dominant disease candidates into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################


@transform(annotateVariantsSNPeff, regex(r"variants/(\S*Multiplex\S+|\S*Trio\S+).haplotypeCaller.snpeff.vcf"), add_inputs(r"\1.ped", r"bam/\1.list"), r"variants/\1.compound_hets.table")
def compoundHets(infiles, outfile):
    '''Identify potentially pathogenic compound heterozygous variants'''
    to_cluster = USECLUSTER
    infile, pedfile, bamlist = infiles
    statement = '''GenomeAnalysisTK -T PhaseByTransmission -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s -ped %(pedfile)s -mvf %(infile)s.mvf -o %(infile)s.phased.vcf ;''' % locals(
    )
    statement += '''GenomeAnalysisTK -T ReadBackedPhasing -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamlist)s -V %(infile)s.phased.vcf -o %(infile)s.phased_rbp.vcf --respectPhaseInInput ;''' % locals(
    )
    statement += '''gemini load -v %(infile)s.phased_rbp.vcf -p %(pedfile)s -t snpEff %(infile)s.db ;'''
    statement += '''gemini comp_hets --only-affected --filter "(impact_severity = 'HIGH' OR impact_severity = 'MED') AND (in_esp = 0 OR aaf_esp_all < 0.01) AND (in_1kg = 0 OR aaf_1kg_all < 0.01)" %(infile)s.db > %(outfile)s'''
    P.run()

#########################################################################


@transform(compoundHets, regex(r"variants/(\S+).compound_hets.table"), r"variants/\1.compound_hets.table.load")
def loadCompoundHets(infile, outfile):
    '''Load compound heterozygous variants into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# vcf statistics


@transform((annotateVariantsSNPsift), regex(r"variants/(\S+).vcf"), r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    to_cluster = USECLUSTER
    statement = '''vcf-stats %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()

#########################################################################


@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    E.info("Loading vcf stats...")
    statement = '''python %(scriptsdir)s/vcfstats2db.py %(filenames)s >> %(outfile)s; '''
    statement += '''cat vcfstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=vcf_stats >> %(outfile)s; '''
    statement += '''cat sharedstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=vcf_shared_stats >> %(outfile)s; '''
    statement += '''cat indelstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=indel_stats >> %(outfile)s; '''
    statement += '''cat snpstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty --index=track --table=snp_stats >> %(outfile)s; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@follows(loadROI,
         loadROI2Gene,
         loadSamples)
def loadMetadata():
    pass


@follows(mergeFastqs,
         mapReads)
def mapping():
    pass


@follows(loadPicardDuplicateStats,
         loadPicardAlignStats,
         buildCoverageStats,
         loadCoverageStats)
def postMappingQC():
    pass


@follows(GATKpreprocessing,
         listOfBAMs)
def gatk():
    pass


@follows(haplotypeCaller)
def callVariants():
    pass


@follows(annotateVariantsSNPeff,
         variantAnnotator,
         variantRecalibrator,
         applyVariantRecalibration,
         annotateVariantsSNPsift)
def annotation():
    pass


@follows(findGenes)
def genesOfInterest():
    pass


@follows(vcfToTable,
         loadVariantAnnotation)
def tabulation():
    pass


@follows(deNovoVariants,
         tabulateDeNovos,
         loadDeNovos,
         lowerStringencyDeNovos,
         tabulateLowerStringencyDeNovos,
         loadLowerStringencyDeNovos,
         dominantVariants,
         tabulateDoms,
         loadDoms,
         recessiveVariants,
         tabulateRecs,
         loadRecs,
         compoundHets,
         loadCompoundHets)
def filtered():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(mapping,
         postMappingQC,
         gatk,
         callVariants,
         annotation,
         genesOfInterest,
         tabulation,
         filtered,
         vcfstats)
def full():
    pass

#########################################################################
#########################################################################
#########################################################################


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''
    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''
    E.info("updating documentation")
    P.run_report(clean=False)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
