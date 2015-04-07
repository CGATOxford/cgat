"""
====================
Exome pipeline
====================

:Author: David Sims & Kath Fawcett
:Release: $Id$
:Date: |today|
:Tags: Python

The exome pipeline imports unmapped reads from one or more fastq or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called, annotated,
and filtered according to various inheritance models (de novo,
dominant, and recessive).


   1. Align to genome using gapped alignment (BWA-MEM)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK and deduplication in Picard
   4. Calculate the ratio of reads on the X and Y chromosomes to assert sex
   5. Variant calling in families (SNVs & indels) using GATK HaplotypeCaller
   6. HapMap genotyping in individuals using GATK HaplotypeCaller
   7. Comparison of Hapmap genotypes to assess relatedness (VCFtools)
   8. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   9. Variant quality score recalibration (GATK)
   10. Flags variants within genes of interest (such as known disease genes)
      (GATK) (optional)
   11. Filters potential de novo variants
   12. Filters potential de novo variants using lower stringency criteria
   13. Filters potential dominant mutations
   14. Filters potential homozygous recessive mutations
   15. Filters potential compound heterozygous mutations
   16. Generates summary statistics for unfiltered vcf file
   17. Generates report

.. note::

   1. Great care should be taken in interpreting lower stringency de
      novo variants as it is expected almost all will be false
      positives. Users should examine them manually in the csvdb
      database and html report.

   2. Great care should be taken when interpreting compound
      heterozygous changes.  Gemini is very permissive and users
      should examine the genotypes in the csvdb table to make sure
      they are consistent with a recessive inheritance pattern.

To do:
   1. Allow users to add other training sets for variant quality score
      recalibration
   2. Allow users to add annotations using SnpSift

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <family>_<sample>-<condition>-<replicate>.<suffix>

``family`` = "Single", "Trio", or "Multiplex" followed by numerical
identifier.  ``sample`` and ``condition`` make up an
:term:`experiment`, while ``replicate`` denotes the :term:`replicate`
within an :term:`experiment`.  The ``suffix`` determines the file
type. The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted
   by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

If you are submitting families then a .ped file for each family must
be supplied within your working firectory.  This is a tab-delimited
file named <family>.ped (where <family> is the family ID in the title
of the corresponding fastq files) with no header and one individual
per line according to the following pattern:

family_id sample_id father_id mother_id sex phenotype

family_id and sample_id should correspond to <family> and
<family>-<sample> in the sra/fastq filenames, father_id and mother_id
should be '0' if unknown, sex should be '1' if male, '2' if female and
'0' if unknown, and phenotype should be '1' if unaffected, '2' if
affected and '0' if unknown.

If you are running the functions to look for compound heterozygotes in
Multiplex families then there is a further requirement for the .ped
files.  The phasing tools expect a trio and therefore any other family
members (other than parents and one child) must be labelled as
unrelated.  That is, the first additional family member could be
labelled "family0" in the family_id column, and subsequent additional
family members could be "family1", "family2" and so on.  For example,
a multiplex family called Multiplex1 may have two parents and two
affected children.  The .ped file would look like this:

Multiplex1 ID1 0 0 1 1
Multiplex1 ID2 0 0 2 1
Multiplex1 ID3 ID1 ID2 1 2
Family0 ID4 ID1 ID2 2 2

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

Pipeline output
===============

The major output is a csvdb containing quality control information by
sample and variant information by family and an html report with
similar information.

Example
=======

ToDo: make exome sequencing example

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Requirements:

* BWA >= 0.7.8
* picardtools >= 1.106
* samtools >= 1.1
* GATK >= 2.7
* snpEff >= 4.0
* Gemini >= ?
* VCFtools >= 0.1.8a

Code
====

"""

# load modules
from ruffus import *
from ruffus.combinatorics import *
from rpy2.robjects import r as R
import sys
import os
import csv
import glob
import re
import sqlite3
import shutil
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGATPipelines.PipelineExome as PipelineExome

###############################################################################
###############################################################################
###############################################################################
# load options from the config file


P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../exome.ini", "exome.ini"])

PARAMS = P.PARAMS
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")


def getGATKOptions():
    return "-l mem_free=1.4G -l picard=1"

###############################################################################
###############################################################################
###############################################################################
# Load target and sample data into database
# The following functions are designed to upload meta-data to the csvdb
# These haven't been fully implemented yet


@jobs_limit(1, "db")
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
              --header-names=%(header)s
              --table=%(tablename)s
            > %(outfile)s  '''
    P.run()

###############################################################################


@jobs_limit(1, "db")
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

###############################################################################


@jobs_limit(1, "db")
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

###############################################################################
###############################################################################
###############################################################################
# Alignment to a reference genome


@follows(mkdir("bam"))
@transform(INPUT_FORMATS, REGEX_FORMATS, r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA-MEM (output=SAM), convert to BAM,
    sort and index BAM file'''
    job_options = "-l mem_free=8G"
    job_threads = 2
    track = P.snip(os.path.basename(outfile), ".bam")
    m = PipelineMapping.BWAMEM(remove_unique=PARAMS["bwa_remove_non_unique"])
    statement = m.build((infiles,), outfile)
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Post-alignment QC


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.picard_stats")
def PicardAlignStats(infile, outfile):
    '''Run Picard CollectMultipleMetrics on each BAM file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    PipelineMappingQC.buildPicardAlignmentStats(infile, outfile, genome)

###############################################################################


@jobs_limit(1, "db")
@merge(PicardAlignStats, "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# GATK


@follows(mkdir("gatk"))
@transform(mapReads, regex(r"bam/(\S+).bam"), r"gatk/\1.readgroups.bam")
def GATKReadGroups(infile, outfile):
    '''Reorders BAM according to reference fasta and adds read groups using
    GATK'''
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = re.sub(r'-\w+-\w+\.bam', '', os.path.basename(infile))
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = 3

    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    PipelineExome.GATKReadGroups(infile, outfile, genome,
                                 library, platform,
                                 platform_unit, track)

###############################################################################
###############################################################################
###############################################################################
# Remove duplicates, realign and recalibrate lane-by-lane


@transform(GATKReadGroups, 
           regex(r"gatk/(\S+).readgroups.bam"), 
           r"gatk/\1.dedup.bam")
def RemoveDuplicatesLane(infile, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.buildPicardDuplicateStats(infile, outfile)
    IOTools.zapFile(infile)

###############################################################################


@jobs_limit(1, "db")
@merge(RemoveDuplicatesLane, "picard_duplicate_stats_lane.load")
def loadPicardDuplicateStatsLane(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)

###############################################################################


@transform(RemoveDuplicatesLane,
           regex(r"gatk/(\S+).dedup.bam"),
           r"gatk/\1.realigned.bam")
def GATKIndelRealignLane(infile, outfile):
    '''realigns around indels using GATK'''
    threads = PARAMS["gatk_threads"]
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    PipelineExome.GATKIndelRealign(infile, outfile, genome, threads)
    IOTools.zapFile(infile)

###############################################################################


@transform(GATKIndelRealignLane,
           regex(r"gatk/(\S+).realigned.bam"),
           r"gatk/\1.bqsr.bam")
def GATKBaseRecal(infile, outfile):
    '''recalibrates base quality scores using GATK'''
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    PipelineExome.GATKBaseRecal(infile, outfile, genome,
                                dbsnp, solid_options)
    IOTools.zapFile(infile)

###############################################################################
###############################################################################
###############################################################################
# Merge BAMs across different lanes for the same sample


@collate(GATKBaseRecal,
         regex(r"gatk/(\S+-\S+)-(\S+)-(\S+).bqsr.bam"),
         r"gatk/\1.merged.bam")
def mergeBAMs(infiles, outfile):
    '''merges BAMs for a single sample over multiple lanes'''
    inputfiles = " INPUT=".join(infiles)
    outf = open(outfile + ".count", "w")
    outf.write(str(len(infiles)))
    outf.close()
    statement = '''MergeSamFiles 
                   INPUT=%(inputfiles)s 
                   OUTPUT=%(outfile)s 
                   ASSUME_SORTED=true; '''
    statement += '''samtools index %(outfile)s ;''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Remove duplicates sample-by-sample


@transform(mergeBAMs, 
           regex(r"gatk/(\S+).merged.bam"),
           add_inputs(r"gatk/\1.merged.bam.count"),
           r"gatk/\1.dedup.bam")
def RemoveDuplicatesSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    infile, countfile = infiles
    countf = open(countfile, "r")
    if countf.read() != '1':
        PipelineMappingQC.buildPicardDuplicateStats(infile, outfile)
    else:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(infile + ".bai", outfile + ".bai")
    IOTools.zapFile(infile)

###############################################################################


@jobs_limit(1, "db")
@merge(RemoveDuplicatesSample, "picard_duplicate_stats_sample.load")
def loadPicardDuplicateStatsSample(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# Coverage of targetted area


@transform(RemoveDuplicatesSample, 
           regex(r"gatk/(\S+).dedup.bam"), 
           r"gatk/\1.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed
    file using Picard'''
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    PipelineMappingQC.buildPicardCoverageStats(infile, outfile,
                                               baits, regions)

###############################################################################


@jobs_limit(1, "db")
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    '''Import coverage statistics into SQLite'''
    PipelineMappingQC.loadPicardCoverageStats(infiles, outfile)

###############################################################################
###############################################################################
###############################################################################
# Realign sample-by-sample

@follows(buildCoverageStats)
@transform(RemoveDuplicatesSample,
           regex(r"gatk/(\S+).dedup.bam"),
           add_inputs(r"gatk/\1.merged.bam.count"),
           r"gatk/\1.realigned.bam")
def GATKIndelRealignSample(infiles, outfile):
    '''realigns around indels using GATK'''
    infile, countfile = infiles
    threads = PARAMS["gatk_threads"]
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    countf = open(countfile, "r")
    if countf.read() > '1':
        PipelineExome.GATKIndelRealign(infiles[0], outfile, genome, threads)
    else:
        shutil.copyfile(infile, outfile)
        shutil.copyfile(infile + ".bai", outfile + ".bai")
    IOTools.zapFile(infile)

###############################################################################
###############################################################################
###############################################################################
# Guess sex


@follows(mkdir("xy_ratio"))
@transform(GATKIndelRealignSample, 
           regex(r"gatk/(\S+).realigned.bam"),
           r"xy_ratio/\1.sex")
def calcXYratio(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads
    per megabase of sequence on X and Y'''
    PipelineExome.guessSex(infile, outfile)

###############################################################################


@merge(calcXYratio, "xy_ratio/xy_ratio.tsv")
def mergeXYRatio(infiles, outfile):
    '''merge XY ratios from all samples and load into database'''
    inlist = " ".join(infiles)
    statement = '''python %(scriptsdir)s/combine_tables.py
                   --add-file-prefix --regex-filename="xy_ratio/(\S+).sex"
                   --no-titles --missing-value=0 --ignore-empty
                   -L %(outfile)s.log -v 6
                   --cat=Track %(inlist)s
                   > %(outfile)s'''
    P.run()

###############################################################################


@transform(mergeXYRatio, regex(r"xy_ratio/xy_ratio.tsv"),
           r"xy_ratio/xy_ratio.load")
def loadXYRatio(infile, outfile):
    '''load into database'''
    P.load(infile, outfile, "--header-names=Track,X,Y,XY_ratio")

###############################################################################
###############################################################################
###############################################################################
# Variant Calling


@collate(GATKIndelRealignSample, 
         regex(r"gatk/(\S+?)-(\S+).realigned.bam"),
         r"gatk/\1.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for each family,
    for use in variant calling'''
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')

###############################################################################


@follows(mkdir("variants"))
@transform(listOfBAMs, regex(r"gatk/(\S+).list"),
           r"variants/\1.haplotypeCaller.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    job_options = getGATKOptions()
    job_threads = 3
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    options = PARAMS["gatk_hc_options"]
    PipelineExome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options)

###############################################################################
###############################################################################
###############################################################################
# HapMap Genotypes


@follows(mkdir("hapmap"))
@files(PARAMS["hapmap_vcf"], "hapmap/hapmap_exome.bed")
def SelectExonicHapmapVariants(infile, outfile):
    '''Select variants from HapMap project that are located within
    exome target regions. Assumes bgzipped & tabix indexed Hapmap VCF file.'''
    bed = PARAMS["roi_regions"]
    statement = '''tabix -B %(infile)s %(bed)s |
                   awk '{OFS="\\t"; if (!/^#/){print $1,$2-1,$2}}'
                   > %(outfile)s''' % locals()
    P.run()

###############################################################################


@follows(SelectExonicHapmapVariants)
@transform(GATKIndelRealignSample, 
           regex(r"gatk/(\S+).realigned.bam"),
           add_inputs("hapmap/hapmap_exome.bed"),
           r"hapmap/\1.hapmap.vcf")
def HapMapGenotype(infiles, outfile):
    '''Genotype HapMap SNPs using HaplotypeCaller in each individual'''
    infile, intervals = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    padding = PARAMS["hapmap_padding"]
    options = PARAMS["hapmap_hc_options"]
    PipelineExome.haplotypeCaller(infile, outfile, genome, dbsnp,
                                  intervals, padding, options)

###############################################################################


@transform(HapMapGenotype, regex(r"hapmap/(\S+).hapmap.vcf"),
           r"hapmap/\1.hapmap.vcf.gz")
def indexVCFs(infile, outfile):
    '''Genotype HapMap SNPs using HaplotypeCaller in each individual'''
    statement = '''bgzip -c %(infile)s > %(outfile)s;
                   tabix -p vcf %(outfile)s; '''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Compare Hapmap genotypes to assess relatedness


@follows(indexVCFs, mkdir("hapmap/vcfcompare"))
@permutations("hapmap/*.hapmap.vcf.gz",
              formatter("hapmap/(\S+).hapmap.vcf.gz$"),
              2,
              r"hapmap/vcfcompare/{basename[0][0]}_vs_{basename[1][0]}.vcfcompare")
def vcfCompare(infiles, outfile):
    '''Compare HapMap genotypes from each pair of calls '''
    sample1, sample2 = infiles
    name1 = P.snip(os.path.basename(sample1), ".hapmap.vcf.gz")
    name2 = P.snip(os.path.basename(sample2), ".hapmap.vcf.gz")
    statement = '''vcf-compare -g -m %(name1)s:%(name2)s
                   %(sample1)s %(sample2)s > %(outfile)s'''
    P.run()

###############################################################################


@transform(vcfCompare,
           regex(r"hapmap/vcfcompare/(\S+).vcfcompare"),
           r"hapmap/vcfcompare/\1.ndr")
def parseVcfCompare(infile, outfile):
    '''Extract non-reference discordance rate for each comparison'''
    statement = '''cat %(infile)s
                   | grep "Non-reference Discordance Rate (NDR):"
                   | cut -f 3
                   > %(outfile)s'''
    P.run()

###############################################################################


@merge(parseVcfCompare, "hapmap/vcfcompare/ndr.tsv")
def mergeNDR(infiles, outfile):
    '''Merge NDR values for all files'''
    outf = open(outfile, "w")
    outf.write("sample1\tsample2\tNDR\n")
    for infile in infiles:
        inbase = P.snip(os.path.basename(infile), ".ndr")
        sample1, sample2 = inbase.split("_vs_")
        inf = open(infile, "r")
        ndr = inf.readline().strip()
        result = "\t".join([sample1, sample2, ndr])
        outf.write(result + "\n")
        inf.close()
    outf.close()

###############################################################################


@transform(mergeNDR, regex(r"hapmap/vcfcompare/ndr.tsv"),
           r"hapmap/vcfcompare/ndr.load")
def loadNDR(infile, outfile):
    '''Load NDR into database'''
    P.load(infile, outfile)

###############################################################################
###############################################################################
###############################################################################
# Variant Annotation


@transform(haplotypeCaller,
           regex(r"variants/(\S+).haplotypeCaller.vcf"),
           r"variants/\1.haplotypeCaller.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    job_options = "-l mem_free=6G"
    job_threads = 4
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''/ifs/apps/bio/snpEff-3.3-dev/snpEff.sh eff
                    -c %(config)s
                    -v %(snpeff_genome)s
                    -o gatk %(infile)s > %(outfile)s''' % locals()
    P.run()

###############################################################################


@transform(annotateVariantsSNPeff,
           regex(r"variants/(\S+).haplotypeCaller.snpeff.vcf"),
           r"variants/\1.haplotypeCaller.snpeff.table")
def vcfToTableSnpEff(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["annotation_snpeff_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(1, "db")
@transform(vcfToTableSnpEff, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadTableSnpEff(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")


###############################################################################
# GATK Variant Annotator


@follows(annotateVariantsSNPeff)
@transform(haplotypeCaller,
           regex(r"variants/(\S+).haplotypeCaller.vcf"),
           add_inputs(r"gatk/\1.list",
                      r"variants/\1.haplotypeCaller.snpeff.vcf"),
           r"variants/\1.haplotypeCaller.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    vcffile, bamlist, snpeff_file = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    annotations = PARAMS["gatk_variant_annotations"]
    PipelineExome.variantAnnotator(vcffile, bamlist, outfile, genome,
                                   dbsnp, annotations, snpeff_file)


###############################################################################
# SNP Recalibration


@transform(variantAnnotator,
           regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"),
           r"variants/\1.haplotypeCaller.snp_vqsr.recal")
def variantRecalibratorSnps(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    dbsnp = PARAMS["gatk_dbsnp"]
    job_options = getGATKOptions()
    job_threads = 3
    track = P.snip(outfile, ".recal")
    kgenomes = PARAMS["gatk_kgenomes"]
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    mode = 'SNP'
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode, dbsnp,
                                      kgenomes, hapmap, omni)

###############################################################################


@follows(variantRecalibratorSnps)
@transform(variantAnnotator,
           regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"),
           add_inputs(r"variants/\1.haplotypeCaller.snp_vqsr.recal",
                      r"variants/\1.haplotypeCaller.snp_vqsr.tranches"),
           r"variants/\1.haplotypeCaller.snp_vqsr.vcf")
def applyVariantRecalibrationSnps(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    mode = 'SNP'
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches,
                                            outfile, genome, mode)

###############################################################################
# Indel recalibration

@transform(applyVariantRecalibrationSnps,
           regex(r"variants/(\S+).haplotypeCaller.snp_vqsr.vcf"),
           r"variants/\1.haplotypeCaller.vqsr.recal")
def variantRecalibratorIndels(infile, outfile):
    '''Create variant recalibration file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    job_options = getGATKOptions()
    job_threads = 3
    track = P.snip(outfile, ".recal")
    mills = PARAMS["gatk_mills"]
    dbsnp = PARAMS["gatk_dbsnp"]
    mode = 'INDEL'
    PipelineExome.variantRecalibrator(infile, outfile, genome, mode,
                                      mills=mills)

###############################################################################


@follows(variantRecalibratorIndels)
@transform(applyVariantRecalibrationSnps,
           regex(r"variants/(\S+).haplotypeCaller.snp_vqsr.vcf"),
           add_inputs(r"variants/\1.haplotypeCaller.vqsr.recal",
                      r"variants/\1.haplotypeCaller.vqsr.tranches"),
           r"variants/\1.haplotypeCaller.vqsr.vcf")
def applyVariantRecalibrationIndels(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    vcf, recal, tranches = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    mode = 'INDEL'
    PipelineExome.applyVariantRecalibration(vcf, recal, tranches,
                                            outfile, genome, mode)


###############################################################################
# SnpSift


@transform(applyVariantRecalibrationIndels,
           regex(r"variants/(\S+).haplotypeCaller.vqsr.vcf"),
           r"variants/\1.haplotypeCaller.snpsift.vcf")
def annotateVariantsSNPsift(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-l mem_free=6G"
    job_threads = 4
    track = P.snip(os.path.basename(infile), ".vqsr.vcf")
    dbNSFP = PARAMS["annotation_snpsift_dbnsfp"]
    thousand_genomes = PARAMS["annotation_thousand_genomes"]
    # The following statement is not fully implemented yet
    # statement = '''SnpSift.sh geneSets -v /ifs/projects/proj016/data/1000Genomes/msigdb.v4.0.symbols.gmt %(infile)s > variants/%(track)s_temp1.vcf; checkpoint;''' % locals()

    statement = '''SnpSift.sh dbnsfp -v -db %(dbNSFP)s %(infile)s
                    > variants/%(track)s_temp1.vcf; checkpoint;
                    SnpSift.sh annotate %(thousand_genomes)s
                    variants/%(track)s_temp1.vcf > %(outfile)s;
                    rm -f variants/%(track)s_temp1.vcf;''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Genes of interest


@active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
@transform((annotateVariantsSNPsift),
           regex(r"variants/(\S+).haplotypeCaller.snpsift.vcf"),
           r"variants/\1.genes.vcf")
def findGenes(infile, outfile):
    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
    if variant is within a gene of interest as defined in the ini
    file'''
    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
    statement = '''GenomeAnalysisTK -T VariantFiltration
                   -R %%(bwa_index_dir)s/%%(genome)s.fa
                   --variant %(infile)s
                   --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
                   --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Tabulation


TABULATION_INPUT = {0: annotateVariantsSNPsift, 1: findGenes}


@transform(TABULATION_INPUT[PARAMS["annotation_add_genes_of_interest"]],
           regex(r"variants/((\S+).haplotypeCaller.snpsift|(\S+).genes).vcf"),
           r"variants/\1.table")
def vcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(1, "db")
@transform(vcfToTable, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''
    P.load(infile, outfile, options="--retry --ignore-empty")

###############################################################################
###############################################################################
###############################################################################
# Confirm parentage (do novo trios only)

@follows(loadVariantAnnotation)
@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"),
           r"variants/\1.parentage")
def confirmParentage(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
        'family', 'sample', 'father', 'mother', 'sex', 'status'])
    trio = P.snip(os.path.basename(pedfile), ".ped")
    trio = trio.replace(".", "_").replace("-", "_")
    database = PARAMS["database"]
    proband = None
    mother = None
    father = None
    for row in pedigree:
        if row['status'] == '2':
            proband = row['sample'].replace(".", "_").replace("-", "_")
            mother = row['mother'].replace(".", "_").replace("-", "_")
            father = row['father'].replace(".", "_").replace("-", "_")
    E.info("proband:" + proband)
    E.info("mother:" + mother)
    E.info("father:" + father)

    query = '''SELECT '%(trio)s' as family,
               hom_nonref_trio/(hom_nonref_parents+0.0) as hom_nonref_conc,
               child_mat_het_nonref/(maternal_hom_nonref+0.0) as maternal_conc,
               child_pat_het_nonref/(paternal_hom_nonref+0.0) as paternal_conc,
                *
               FROM
               (SELECT count(*) as hom_nonref_trio
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as hom_nonref_parents
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as maternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_mat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND %(mother)s_PL LIKE '%%%%,0'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND  %(father)s_PL LIKE '0,%%%%'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as paternal_hom_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND cast(%(proband)s_DP as INTEGER) > 10),
               (SELECT count(*) as child_pat_het_nonref
               FROM %(trio)s_genes_table
               where  CHROM NOT IN ('X', 'Y')
               AND  %(father)s_PL LIKE '%%%%,0'
               AND cast(%(father)s_DP as INTEGER) > 10
               AND  %(mother)s_PL LIKE '0,%%%%'
               AND cast(%(mother)s_DP as INTEGER) > 10
               AND %(proband)s_PL LIKE '%%%%,0,%%%%'
               AND cast(%(proband)s_DP as INTEGER) > 10)''' % locals()
    statement = '''sqlite3 %(database)s "%(query)s" > %(outfile)s
                   2> %(outfile)s.log''' % locals()
    P.run()

###############################################################################
###############################################################################
###############################################################################
# De novos


@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"), r"variants/\1.filtered.vcf")
def deNovoVariants(infiles, outfile):
    '''Filter de novo variants based on provided jexl expression'''
    job_options = getGATKOptions()
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    infile, pedfile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(father)s").getDP()>=10&&vc.getGenotype("%(mother)s").getDP()>=10&&vc.getGenotype("%(child)s").getPL().0>20&&vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(child)s").getPL().2>0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(father)s").getPL().1>20&&vc.getGenotype("%(father)s").getPL().2>20&&vc.getGenotype("%(mother)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().1>20&&vc.getGenotype("%(mother)s").getPL().2>20&&vc.getGenotype("%(child)s").getAD().1>=3&&((vc.getGenotype("%(child)s").getAD().1)/(vc.getGenotype("%(child)s").getDP().floatValue()))>=0.25&&(vc.getGenotype("%(father)s").getAD().1==0||(vc.getGenotype("%(father)s").getAD().1>0&&((vc.getGenotype("%(father)s").getAD().1)/(vc.getGenotype("%(father)s").getDP().floatValue()))<0.05))&&(vc.getGenotype("%(mother)s").getAD().1==0||(vc.getGenotype("%(mother)s").getAD().1>0&&((vc.getGenotype("%(mother)s").getAD().1)/(vc.getGenotype("%(mother)s").getDP().floatValue()))<0.05))&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    PipelineExome.selectVariants(infile, outfile, genome, select)


@transform(deNovoVariants,
           regex(r"variants/(\S+).filtered.vcf"),
           r"variants/\1.filtered.table")
def tabulateDeNovos(infile, outfile):
    '''Tabulate de novo variants'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(1, "db")
@transform(tabulateDeNovos,
           regex(r"variants/(\S+).filtered.table"),
           r"variants/\1.filtered.table.load")
def loadDeNovos(infile, outfile):
    '''Load de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################


@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"), r"variants/\1.denovos.vcf")
def lowerStringencyDeNovos(infiles, outfile):
    '''Filter lower stringency de novo variants based on provided jexl
    expression'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    infile, pedfile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    select = '''vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().0==0&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    PipelineExome.selectVariants(infile, outfile, genome, select)


@transform(lowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.vcf"),
           r"variants/\1.denovos.table")
def tabulateLowerStringencyDeNovos(infile, outfile):
    '''Tabulate lower stringency de novo variants'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)


@jobs_limit(1, "db")
@transform(tabulateLowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.table"),
           r"variants/\1.denovos.table.load")
def loadLowerStringencyDeNovos(infile, outfile):
    '''Load lower stringency de novos into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Dominant


@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"),
           r"variants/\1.dominant.vcf")
def dominantVariants(infiles, outfile):
    '''Filter variants according to autosomal dominant disease model'''
    infile, pedfile = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
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
    # currently 1000G filter is performed at the report stage (not in csvdb)
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().1==0%(unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals()
    PipelineExome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(dominantVariants,
           regex(r"variants/(\S+).dominant.vcf"),
           r"variants/\1.dominant.table")
def tabulateDoms(infile, outfile):
    '''Tabulate dominant disease candidate variants'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(1, "db")
@transform(tabulateDoms, regex(r"variants/(\S+).dominant.table"),
           r"variants/\1.dominant.table.load")
def loadDoms(infile, outfile):
    '''Load dominant disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Recessive


@transform(annotateVariantsSNPsift,
           regex(
               r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"),
           r"variants/\1.recessive.vcf")
def recessiveVariants(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    parents = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
            if row['father'] != '0':
                parents += [row['father']]
            if row['mother'] != '0':
                parents += [row['mother']]
        elif row['status'] == '1' and row['sample'] not in parents:
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
    if len(unaffecteds) == 0:
        unaffecteds_exp = ''
    else:
        unaffecteds_exp = '&&vc.getGenotype("' + \
            ('").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)) + \
            '").getPL().2!=0'
    if len(parents) == 0:
        parents_exp = ''
    else:
        parents_exp = '&&vc.getGenotype("' + \
            ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
            '").getPL().1==0'
    # need a way of specifying that other unaffecteds eg. sibs can't be
    # homozygous for alt allele
    select = '''vc.getGenotype("%(affecteds_exp)s").getPL().2==0%(unaffecteds_exp)s%(parents_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")''' % locals(
    )
    PipelineExome.selectVariants(infile, outfile, genome, select)

###############################################################################


@transform(recessiveVariants,
           regex(r"variants/(\S+).recessive.vcf"),
           r"variants/\1.recessive.table")
def tabulateRecs(infile, outfile):
    '''Tabulate potential homozygous recessive disease variants'''
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    columns = PARAMS["gatk_vcf_to_table"]
    PipelineExome.vcfToTable(infile, outfile, genome, columns)

###############################################################################


@jobs_limit(1, "db")
@transform(tabulateRecs, regex(r"variants/(\S+).recessive.table"),
           r"variants/\1.recessive.table.load")
def loadRecs(infile, outfile):
    '''Load homozygous recessive disease candidates into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# Compound hets


# why does this not work from snpsift VCF file? DS
# @transform(annotateVariantsSNPsift,
#    regex(r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
@transform(annotateVariantsSNPeff,
           regex(
               r"variants/(\S*Multiplex\S+|\S*Trio\S+).haplotypeCaller.snpeff.vcf"),
           add_inputs(r"\1.ped", r"gatk/\1.list"),
           r"variants/\1.phased.vcf")
def phasing(infiles, outfile):
    '''phase variants with GATK'''
    infile, pedfile, bamlist = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa" 
    statement = '''GenomeAnalysisTK -T PhaseByTransmission
                   -R %(genome)s
                   -V %(infile)s
                   -ped %(pedfile)s
                   -mvf %(infile)s.mvf
                   -o %(outfile)s ;''' % locals()
    P.run()

###############################################################################


@transform(phasing, regex(r"variants/(\S+).phased.vcf"),
                          add_inputs(r"gatk/\1.list"), r"variants/\1.rbp.vcf")
def readbackedphasing(infiles, outfile):
    '''phase variants with ReadBackedPhasing'''
    job_options = getGATKOptions()
    infile, bamlist = infiles
    genome = PARAMS["bwa_index_dir"] + "/" + PARAMS["genome"] + ".fa"
    statement = '''GenomeAnalysisTK -T ReadBackedPhasing
                   -R %(genome)s
                   -I %(bamlist)s
                   -V %(infile)s
                   -o %(outfile)s; ''' % locals()
    P.run()

###############################################################################


@transform(readbackedphasing,
           regex(r"variants/(\S*Multiplex\S+|\S*Trio\S+).rbp.vcf"),
           add_inputs(r"\1.ped"),
           r"variants/\1.compound_hets.table")
def compoundHets(infiles, outfile):
    '''Identify potentially pathogenic compound heterozygous variants
    (phasing with GATK followed by compound het calling using Gemini'''
    infile, pedfile = infiles
    statement = '''gemini load -v %(infile)s
                   -p %(pedfile)s -t snpEff %(infile)s.db ;
                   gemini comp_hets
                   --only-affected
                   --filter
                   "(impact_severity = 'HIGH' OR impact_severity = 'MED')
                   AND (in_esp = 0 OR aaf_esp_all < 0.01)
                   AND (in_1kg = 0 OR aaf_1kg_all < 0.01)"
                   %(infile)s.db > %(outfile)s'''
    P.run()

###############################################################################


@jobs_limit(1, "db")
@transform(compoundHets, regex(r"variants/(\S+).compound_hets.table"),
           r"variants/\1.compound_hets.table.load")
def loadCompoundHets(infile, outfile):
    '''Load compound heterozygous variants into database'''
    P.load(infile, outfile,
           options="--retry --ignore-empty --allow-empty-file")

###############################################################################
###############################################################################
###############################################################################
# vcf statistics


@transform((annotateVariantsSNPsift), regex(
    r"variants/(\S+).vcf"), r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    statement = '''vcf-stats %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals(
    )
    P.run()

###############################################################################


@jobs_limit(1, "db")
@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    scriptsdir = PARAMS["general_scriptsdir"]
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    E.info("Loading vcf stats...")
    statement = '''python %(scriptsdir)s/vcfstats2db.py %(filenames)s >> %(outfile)s; '''
    statement += '''cat vcfstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=vcf_stats >> %(outfile)s; '''
    statement += '''cat sharedstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=vcf_shared_stats >> %(outfile)s; '''
    statement += '''cat indelstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=indel_stats >> %(outfile)s; '''
    statement += '''cat snpstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=snp_stats >> %(outfile)s; '''
    P.run()

###############################################################################
###############################################################################
###############################################################################
# Targets


@follows(loadROI,
         loadROI2Gene,
         loadSamples)
def loadMetadata():
    pass


@follows(mapReads,
         loadPicardAlignStats)
def mapping():
    pass


@follows(GATKBaseRecal,
         loadPicardDuplicateStatsLane,
         loadCoverageStats)
def gatk():
    pass


@follows(loadXYRatio,
         indexVCFs,
         loadNDR)
def sampleFeatures():
    pass


@follows(haplotypeCaller)
def callVariants():
    pass


@follows(loadTableSnpEff,
         loadVariantAnnotation)
def annotation():
    pass


@follows(loadDeNovos)
def denovo():
    pass


@follows(loadLowerStringencyDeNovos)
def denovo2():
    pass


@follows(loadDoms)
def dominant():
    pass


@follows(loadRecs)
def recessive():
    pass


@follows(loadCompoundHets)
def compoundHet():
    pass


@follows(denovo,
         denovo2,
         dominant,
         recessive)
def filtering():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
    pass


@follows(mapping,
         gatk,
         # sampleFeatures,
         callVariants,
         annotation,
         confirmParentage,
         filtering)
def full():
    pass

###############################################################################
###############################################################################
###############################################################################
# Reports


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
