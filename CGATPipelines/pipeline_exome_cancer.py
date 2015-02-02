
"""====================
Exome Cancer pipeline
====================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

###########
To do:

Document fully
make phone home/key option work - GATK public key?
Summarise Indel calling (size of indels called)
Example

######


The exome cancer pipeline imports unmapped reads from matched sample fastqs or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called and filtered


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNPs) on control samples using muTect to generate
      a "panel of normal" variants
   5a. Variant calling (SNPs) with tumour samples using muTect including
      filtering
   5b. Variant calling (indels) using Strelka
   6. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   7. Generates report

.. note::

   An optional downsampling analysis can also be performed to assess how
   coverage a control sample affects the called variants

   1. Currently the pipeline is not able to deal with replicates, i.e
      replicates will be treated seperately.



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

   <patientID>-<tissue>-<replicate>.<suffix>

``patientID`` and ``tissue`` make up an :term:`experiment`, while ``replicate``
denotes the :term:`replicate` within an :term:`experiment`.
The ``suffix`` determines the file type.
The following suffixes/file types are possible:

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

Documentation
-------------

If you would like the genes of interest to be flagged in your vcf,
make add_genes_of_interest=1 (default=0) and provide a list of comma
separated genes (without spaces) in the ini file.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------------+------------+-------------------------------------------+
|*Program*           |*Version*   |*Purpose*                                  |
+--------------------+------------+-------------------------------------------+
|Stampy              |>=0.9.0     |read mapping                               |
+--------------------+------------+-------------------------------------------+
|BWA                 |            |read mapping                               |
+--------------------+------------+-------------------------------------------+
|SAMtools            |            |filtering, SNV / indel calling             |
+--------------------+------------+-------------------------------------------+
|BEDTools            |            |filtering                                  |
+--------------------+------------+-------------------------------------------+
|sra-tools           |            |extracting reads from .sra files           |
+--------------------+------------+-------------------------------------------+
|picard              |>=1.38      |bam/sam files. The .jar files need to be in|
|                    |            |your CLASSPATH environment variable.       |
+--------------------+------------+-------------------------------------------+
|vcf-tools           |            |VCF filtering                              |
+--------------------+------------+-------------------------------------------+
|GATK                | 2.5-2      |local realignment, BQSR, variant calling   |
+--------------------+------------+-------------------------------------------+
|SNPeff              | 3.3        |                                           |
+--------------------+------------+-------------------------------------------+

Pipeline output
===============

The major output is a csvdb containing quality control information
and variant information by patientID and an html report with
similar information.

Example
=======


Code
====

"""

# load modules
from ruffus import *
# from rpy2.robjects import r as R

import numpy
import CGAT.Experiment as E
import sys
import os
import csv
import sqlite3
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Pipeline as P
import glob
import pandas as pd
import itertools
import re
import CGATPipelines.PipelineExome as Exome
USECLUSTER = True

#########################################################################
#########################################################################


def connect():
    '''connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''
    dbh = sqlite3.connect(PARAMS["database"])

    return dbh


#########################################################################
# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS


def getGATKOptions():
    # removed picard=1, surely not neccessary?
    return "-l mem_free=4G"


def getMuTectOptions():
    return "-l mem_free=6G"


#########################################################################
#########################################################################
#########################################################################
# Load target and sample data
# The following functions are designed to upload meta-data to the csvdb
# These haven't been fully implemented yet


# @files(PARAMS["roi_bed"], "roi.load")
# def loadROI(infile, outfile):
#    '''Import regions of interest bed file into SQLite.'''
#    header = "chr,start,stop,feature"
#    tablename = P.toTable(outfile)
#    statement = '''cat %(infile)s
#            | python %%(scriptsdir)s/csv2db.py %(csv2db_options)s
#              --ignore-empty
#              --retry
#              --header-names=%(header)s
#              --table=%(tablename)s
#            > %(outfile)s  '''
#    P.run()

#########################################################################


# @files(PARAMS["roi_to_gene"], "roi2gene.load")
# def loadROI2Gene(infile, outfile):
#    '''Import genes mapping to regions of interest bed file into SQLite.'''
#    tablename = P.toTable(outfile)
#    statement = '''cat %(infile)s
#            | python %%(scriptsdir)s/csv2db.py %(csv2db_options)s
#              --ignore-empty
#              --retry
#              --table=%(tablename)s
#            > %(outfile)s  '''
#    P.run()

#########################################################################

# not currently implemented
# @files(PARAMS["samples"], "samples.load")
# def loadSamples(infile, outfile):
#    '''Import sample information into SQLite.'''
#    tablename = P.toTable(outfile)
#    statement = '''cat %(infile)s
#            | python %%(scriptsdir)s/csv2db.py %(csv2db_options)s
#              --ignore-empty
#              --retry
#              --table=%(tablename)s
#            > %(outfile)s  '''
#    P.run()

#########################################################################
#########################################################################
#########################################################################
# Alignment to a reference genome


@follows(mkdir("bam"))
@transform(("*.fastq.1.gz", "*.fastq.gz", "*.sra"),
           regex(r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
           r"bam/\1.bam")
def mapReads(infile, outfile):
    '''Map reads to the genome using BWA, sort and index BAM file,
    generate alignment statistics and deduplicate using Picard'''

    job_threads = PARAMS["bwa_threads"]
    job_options = "-l mem_free=8G"
    job_threads = 2

    if PARAMS["bwa_algorithm"] == "aln":
        m = PipelineMapping.BWA(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False, align_stats=True, dedup=True)

    elif PARAMS["bwa_algorithm"] == "mem":
        m = PipelineMapping.BWAMEM(
            remove_non_unique=PARAMS["bwa_remove_non_unique"],
            strip_sequence=False, align_stats=True, dedup=True)
    else:
        raise ValueError("bwa algorithm '%s' not known" % algorithm)

    statement = m.build((infile,), outfile)
    P.run()


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
    '''Generate coverage statistics for regions of interest from a
       bed file using Picard'''
    # baits file requires modification to make picard accept it
    # this is performed before CalculateHsMetrics
    to_cluster = USECLUSTER
    baits = PARAMS["roi_baits"]
    modified_baits = infile + "_temp_baits_final.bed"
    regions = PARAMS["roi_regions"]
    statement = '''samtools view -H %(infile)s > %(infile)s_temp_header.txt;
                awk 'NR>2' %(baits)s |
                awk -F '\\t' 'BEGIN { OFS="\\t" } {print $1,$2,$3,"+",$4;}'
                > %(infile)s_temp_baits.bed;
                cat  %(infile)s_temp_header.txt %(infile)s_temp_baits.bed
                > %(modified_baits)s; checkpoint ;
                rm -rf %(infile)s_temp_baits.bed %(infile)s_temp_header.txt
                ''' % locals()
    P.run()

    PipelineMappingQC.loadPicardAlignmentStats(
        infile, outfile, modified_baits, modified_baits)

    IOTools.zapFile(modified_baits)


@follows(buildCoverageStats)
@merge(buildCoverageStats, "coverage_stats.load")
def loadCoverageStats(infiles, outfile):
    '''Import coverage statistics into SQLite'''
    tablename = P.toTable(outfile)
    outf = open('coverage.txt', 'w')
    first = True
    for f in infiles:
        track = P.snip(os.path.basename(f), ".cov")
        lines = [x for x in open(f, "r").readlines()
                 if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
    outf.close()
    tmpfilename = outf.name
    statement = '''cat %(tmpfilename)s
                   | python %%(scriptsdir)s/csv2db.py
                      --add-index=track
                      --table=%(tablename)s
                      --ignore-empty
                      --retry
                   > %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
# GATK realign bams
#########################################################################


@transform(mapReads,
           regex(r"bam/(\S+).bam"),
           r"bam/\1.bqsr.bam")
def GATKpreprocessing(infile, outfile):
    '''Reorders BAM according to reference fasta and add read groups using
       SAMtools, realigns around indels and recalibrates base quality scores
       using GATK'''

    to_cluster = USECLUSTER
    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('/ifs/scratch')
    job_options = getGATKOptions()
    job_threads = 6
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    gatk_threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]
    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    outfile1 = re.sub(".bqsr", ".readgroups.bqsr", outfile)
    outfile2 = re.sub(".bqsr", ".realign.bqsr", outfile)

    Exome.GATKreadGroups(infile, outfile1, genome, library, platform,
                         platform_unit)

    Exome.GATKrealign(outfile1, outfile2, genome, gatk_threads)

    Exome.GATKrescore(outfile2, outfile, genome, dbsnp, solid_options)

    IOTools.zapFile(outfile1 + ".bam")
    IOTools.zapFile(outfile1 + ".bai")
    IOTools.zapFile(outfile2 + ".bam")
    IOTools.zapFile(outfile2 + ".bai")


@transform(GATKpreprocessing,
           regex("bam/(\S+)-Control-(\d+).bqsr.bam"),
           r"bam/\1-Control-\2.merged.bam")
def mergeSampleBams(infile, outfile):
    '''merge control and tumor bams'''
    # Note: need to change readgroup headers for merge and subsequent
    # splitting of bam files
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    job_threads = 6
    # tmpdir_gatk = P.getTempDir('tmpbam')
    tmpdir_gatk = P.getTempDir('/ifs/scratch')
    # threads = PARAMS["gatk_threads"]

    outfile_tumor = outfile.replace("Control", PARAMS["mutect_tumour"])
    infile_tumor = infile.replace("Control", PARAMS["mutect_tumour"])

    infile_base = os.path.basename(infile)
    infile_tumor_base = infile_base.replace("Control", PARAMS["mutect_tumour"])

    track = P.snip(os.path.basename(infile), ".bam")
    track_tumor = track.replace("Control", PARAMS["mutect_tumour"])

    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]

    control_id = "Control.bam"
    tumor_id = control_id.replace("Control", PARAMS["mutect_tumour"])
    tmpdir_gatk = P.getTempDir('.')

    statement = '''AddOrReplaceReadGroups
                    INPUT=%(infile)s
                    OUTPUT=%(tmpdir_gatk)s/%(infile_base)s
                    RGLB=%(library)s RGPL=%(platform)s
                    RGPU=%(platform_unit)s RGSM=%(track)s
                    ID=%(track)s
                    VALIDATION_STRINGENCY=SILENT ;
                    checkpoint ;''' % locals()
    statement += '''AddOrReplaceReadGroups
                    INPUT=%(infile_tumor)s
                    OUTPUT=%(tmpdir_gatk)s/%(infile_tumor_base)s
                    RGLB=%(library)s RGPL=%(platform)s
                    RGPU=%(platform_unit)s RGSM=%(track_tumor)s
                    ID=%(track_tumor)s
                    VALIDATION_STRINGENCY=SILENT ;
                    checkpoint ;''' % locals()
    statement += '''samtools merge -rf
                    %(outfile)s
                    %(tmpdir_gatk)s/%(infile_base)s
                    %(tmpdir_gatk)s/%(infile_tumor_base)s
                    ; checkpoint ;''' % locals()
    statement += '''samtools index %(outfile)s ;
                    checkpoint ;'''
    statement += '''rm -rf %(tmpdir_gatk)s ;
                    checkpoint ; ''' % locals()
    P.run()
    IOTools.zapFile(infile)
    IOTools.zapFile(infile_tumor)


@transform(mergeSampleBams,
           regex("bam/(\S+)-Control-(\d+).merged.bam"),
           r"bam/\1-Control-\2.realigned.bqsr.bam")
def realignMatchedSample(infile, outfile):
    ''' repeat realignments with merged bam of control and tumor
        this should help avoid problems with sample-specific realignments'''

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    Exome.GATKrealign(infile, outfile, genome)

    IOTools.zapFile(infile)


@transform(realignMatchedSample,
           regex("bam/(\S+)-Control-(\d+).realigned.bqsr.bam"),
           r"bam/\1-Control-\2.realigned.split.bqsr.bam")
def splitMergedRealigned(infile, outfile):
    ''' split realignment file and truncate intermediate bams'''

    track = P.snip(os.path.basename(infile), "realigned.bqsr.bam") + ".bqsr"
    track_tumor = track.replace("Control", PARAMS["mutect_tumour"])
    outfile_tumor = outfile.replace("Control", PARAMS["mutect_tumour"])

    statement = '''samtools view -hb %(infile)s
                   -r %(track)s > %(outfile)s;
                   samtools view -hb %(infile)s
                   -r %(track_tumor)s > %(outfile_tumor)s; checkpoint ;
                   samtools index %(outfile)s;
                   samtools index %(outfile_tumor)s; checkpoint;''' % locals()
    P.run()
    IOTools.zapFile(infile)


@transform(splitMergedRealigned,
           regex("bam/(\S+)-Control-(\S+).realigned.split.bqsr.bam"),
           r"bam/\1-Control-\2.realigned.picard_stats")
def runPicardOnRealigned(infile, outfile):
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    job_threads = 6
    tmpdir_gatk = P.getTempDir('/ifs/scratch')
    # threads = PARAMS["gatk_threads"]

    outfile_tumor = outfile.replace("Control", PARAMS["mutect_tumour"])
    infile_tumor = infile.replace("Control", PARAMS["mutect_tumour"])

    track = P.snip(os.path.basename(infile), ".bam")
    track_tumor = track.replace("Control", PARAMS["mutect_tumour"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    PipelineMappingQC.buildPicardAlignmentStats(infile, outfile, genome)
    PipelineMappingQC.buildPicardAlignmentStats(infile_tumour,
                                                outfile_tumour, genome)

    # check above functions then remove statement
    statement = '''
    cat %(infile)s
    | python %%(scriptsdir)s/bam2bam.py -v 0 --method=set-sequence
    | CollectMultipleMetrics
    INPUT=/dev/stdin
    REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa
    ASSUME_SORTED=true
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile)s;
    cat %(infile_tumor)s
    | python %%(scriptsdir)s/bam2bam.py -v 0
    --method=set-sequence --output-sam
    | CollectMultipleMetrics
    INPUT=/dev/stdin
    REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa
    ASSUME_SORTED=true
    OUTPUT=%(outfile_tumor)s
    VALIDATION_STRINGENCY=SILENT
    >& %(outfile_tumor)s;''' % locals()


@follows(runPicardOnRealigned)
@merge("bam/*.realigned.picard_stats", "realigned_picard_stats.load")
def loadPicardRealigenedAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)

#########################################################################
#########################################################################
#########################################################################
# Variant Calling
#########################################################################


@follows(mkdir("normal_panel_variants"))
@transform(splitMergedRealigned,
           regex(r"bam/(\S+)-Control-(\S).realigned.split.bqsr.bam"),
           r"normal_panel_variants/\1_normal_mutect.vcf")
def callControlVariants(infile, outfile):
    '''run mutect to call snps in tumor sample'''
    cluster_options = getMuTectOptions()
    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"
    gatk_key = PARAMS["mutect_key"]
    cosmic, dbsnp, = (PARAMS["mutect_cosmic"],
                      PARAMS["gatk_dbsnp"])
    genome = "%(bwa_index_dir)s/%(genome)s.fa" % locals()

    Exome.mutectSNPCaller(infile, outfile, mutect_log, genome, cosmic,
                          dbsnp, call_stats_out, cluster_options)
    P.run()


@follows(mkdir("normal_panel_variants"))
@transform(realignMatchedSample,
           regex(r"bam/(\S+)-Control-(\S).realigned.bqsr.bam"),
           r"normal_panel_variants/\1_normal_mutect.vcf")
def callControlVariants2(infile, outfile):
    '''run mutect to call snps in tumor sample'''
    job_options = getMuTectOptions()
    job_threads = 2
    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"
    # mutect repeatedly hangs-up with multithreading
    # furthermore, multithreading doesn't speed up even nearly linearly
    # threads = PARAMS["gatk_threads"]

    if PARAMS["mutect_key"]:
        key = "-et NO_ET -K %s" % PARAMS["mutect_key_path"]
    else:
        key = ""

    cosmic, dbsnp, = (
        PARAMS["mutect_cosmic"],
        PARAMS["gatk_dbsnp"])

    statement = '''java -Xmx4g -jar
    /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
    --analysis_type MuTect
    --reference_sequence %%(bwa_index_dir)s/%%(genome)s.fa
    --cosmic %(cosmic)s  --dbsnp %(dbsnp)s
    --input_file:tumor %(infile)s
    --out %(call_stats_out)s
    --vcf %(outfile)s --artifact_detection_mode
    %(key)s > %(mutect_log)s
    ''' % locals()

    P.run()


@transform(callControlVariants,
           suffix(".vcf"),
           "_slim.vcf.gz")
def indexControlVariants(infile, outfile):
    '''index control vcf for intersection by vcftools'''

    outfile = P.snip(outfile, ".gz")

    statement = '''cut -f1-8 %(infile)s > %(outfile)s;
                   bgzip -f %(outfile)s;
                   tabix -f %(outfile)s.gz''' % locals()
    P.run()


# paramaterise vcf intersection (number of req. observations - currently 1)
@merge(indexControlVariants,
       "normal_panel_variants/combined.vcf")
def mergeControlVariants(infiles, outfile):
    ''' intersect control vcfs to generate a panel of normals for mutect'''
    infiles = " ".join(infiles)

    statement = '''vcf-isec -o -n +1 %(infiles)s
                   > %(outfile)s
                   ''' % locals()
    P.run()


@follows(mkdir("variants"), callControlVariants)
@transform(splitMergedRealigned,
           regex(r"bam/(\S+)-Control-(\S).realigned.split.bqsr.bam"),
           add_inputs(mergeControlVariants),
           r"variants/\1.mutect.snp.vcf")
def runMutect(infiles, outfile):
    '''calls somatic SNPs using MuTect'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        "Control", PARAMS["mutect_tumour"])
    cluster_options = getMuTectOptions()
    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, gatk_key) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"], PARAMS["mutect_key"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    Exome.mutectSNPCaller(infile_tumour, outfile, mutect_log, genome,
                          cosmic, dbsnp, call_stats_out,
                          cluster_options, quality, max_alt_qual,
                          max_alt, max_fraction, tumor_LOD,
                          normal_panel, infile_matched=infile)
    P.run()


# delete once above function checked
@follows(mkdir("variants"), callControlVariants)
@transform(realignMatchedSample,
           regex(r"bam/(\S+)-Control-(\S).realigned.bqsr.bam"),
           add_inputs(mergeControlVariants),
           r"variants/\1.mutect.snp.vcf")
#            r"variants/\1_call_stats.out"])
def runMutect2(infiles, outfile):
    '''calls somatic SNPs using MuTect'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        "Control", PARAMS["mutect_tumour"])
    # mutect repeatedly hangs-up with multithreading
    # furthermore, multithreading doesn't speed up even nearly linearly
    # threads = PARAMS["gatk_threads"]
    job_options = getMuTectOptions()
    job_threads = 2
    # outfile, extended_out = outfiles
    basename = P.snip(outfile, ".mutect.snp.vcf")
    call_stats_out = basename + "_call_stats.out"
    coverage_wig_out = basename + "_coverage.wig"
    mutect_log = basename + ".log"

    cosmic, dbsnp, = (
        PARAMS["mutect_cosmic"],
        PARAMS["gatk_dbsnp"])
    tumor_LOD = PARAMS["mutect_lod"]

    # problems with public key for GATK so this is not yet implemented
    if PARAMS["mutect_key"]:
        key = "-et NO_ET -K %s" % PARAMS["mutect_key_path"]
    else:
        key = ""

    statement = '''java -Xmx4g -jar
    /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
    --analysis_type MuTect
    --reference_sequence %%(bwa_index_dir)s/%%(genome)s.fa
    --cosmic %(cosmic)s
    --dbsnp %(dbsnp)s
    --input_file:normal %(infile)s
    --input_file:tumor %(infile_tumour)s
    --out %(call_stats_out)s
    --coverage_file %(coverage_wig_out)s
    --vcf %(outfile)s
    --min_qscore 20
    --max_alt_alleles_in_normal_qscore_sum 150
    --max_alt_alleles_in_normal_count 5
    --max_alt_allele_in_normal_fraction 0.05
    --gap_events_threshold 2
    --tumor_lod %(tumor_LOD)s
    --enable_extended_output
    --normal_panel %(normal_panel)s
    %(key)s > %(mutect_log)s''' % locals()

    P.run()


# @transform(realignMatchedSample,
#           regex(r"bam/(\S+)-Control-(\S+).bqsr.bam"),
@transform(splitMergedRealigned,
           regex(r"bam/(\S+)-Control-(\S).realigned.split.bqsr.bam"),
           r"variants/\1/results/all.somatic.indels.vcf")
def indelCaller(infile, outfile):
    '''Call somatic indels using Strelka'''
    infile_tumour = infile.replace("Control", PARAMS["mutect_tumour"])
    outdir = "/".join(outfile.split("/")[0:2])
    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])
    cluster_options = "-pe dedicated 12 -R y -l mem_free=1.9G"
    config = "config.ini"
    job_options = "-l mem_free=1.9G"
    job_threads = 12

    PipelineExome.strelkaINDELCaller(infile, infile_tumour, outfile,
                                     genome, config, outdir, cluster_options)

##########################################################################
##########################################################################
##########################################################################
# repeat mutect in reverse and on subsampled control bam as quality control
##########################################################################
# this analysis should be part of an optional check of mutect parameters
# mutect paramters should be identical to the runMutect function above


@follows(mergeControlVariants)
@transform(splitMergedRealigned,
           regex(r"bam/(\S+)-Control-(\S).realigned.split.bqsr.bam"),
           add_inputs(mergeControlVariants),
           r"variants/\1.mutect.reverse.snp.vcf")
def runMutectReverse(infiles, outfile):
    '''Use control as tumor and vis versa to estimate false positive rate'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        "Control", PARAMS["mutect_tumour"])

    cluster_options = getMuTectOptions()
    basename = P.snip(outfile, "_normal_mutect.vcf")
    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    job_options = getMuTectOptions()
    job_threads = 2
    basename = P.snip(outfile, ".mutect.reverse.snp.vcf")
    call_stats_out = basename + "_call_stats.reverse.out"
    coverage_wig_out = basename + "_coverage.reverse.wig"
    mutect_log = basename + ".reverse.log"

    if PARAMS["mutect_key"]:
        key = "-et NO_ET -K %s" % PARAMS["mutect_key_path"]
    else:
        key = ""

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, gatk_key) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"], PARAMS["mutect_key"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    Exome.mutectSNPCaller(infile, outfile, mutect_log, genome,
                          cosmic, dbsnp, call_stats_out,
                          cluster_options, quality, max_alt_qual,
                          max_alt, max_fraction, tumor_LOD,
                          normal_panel, infile_matched=infile_tumour)


# generalise the functions below
# 1. identify sample with highest coverage in control
# - should this check coverage in tumour also?
# 2. subset control bam
# 3. run mutect calling function with subset against unsubsetted tumour
# 4. summary table

adeno_bam = "bam/NU16C-Control-1.realigned.bqsr.bam"


@subdivide(adeno_bam,
           regex("(\S+).bqsr.bam"),
           [r"\1.0.1.bqsr.bam",
            r"\1.0.2.bqsr.bam",
            r"\1.0.3.bqsr.bam",
            r"\1.0.4.bqsr.bam",
            r"\1.0.5.bqsr.bam",
            r"\1.0.6.bqsr.bam",
            r"\1.0.7.bqsr.bam",
            r"\1.0.8.bqsr.bam",
            r"\1.0.9.bqsr.bam",
            r"\1.1.0.bqsr.bam"])
def subsetControlBam(infile, outfiles):
    statements = []
    n = 0
    for fraction in numpy.arange(0.1, 1.1, 0.1):
        outfile = outfiles[n]
        n += 1
        statement = ('''samtools view -s %(fraction)s -b %(infile)s
                     > %(outfile)s''' % locals())
        P.run()


@transform(subsetControlBam,
           suffix(".bam"),
           ".bam.bai")
def indexSubsets(infile, outfile):
    statement = '''samtools index %(infile)s''' % locals()
    P.run()


@follows(indexSubsets)
@transform(subsetControlBam,
           regex(r"bam/(\S+)-Control-1.realigned.(\S+).bqsr.bam"),
           add_inputs(mergeControlVariants),
           r"variants/\1-downsampled-\2.mutect.snp.vcf")
def runMutectOnDownsampled(infiles, outfile):
    '''call somatic SNPs using MuTect on downsampled bams'''
    infile, normal_panel = infiles
    infile_tumour = infile.replace(
        "Control", PARAMS["mutect_tumour"])
    cluster_options = getMuTectOptions()
    basename = P.snip(outfile, "_normal_mutect.vcf")

    call_stats_out = basename + "_call_stats.out"
    mutect_log = basename + ".log"

    (cosmic, dbsnp, quality, max_alt_qual, max_alt,
     max_fraction, tumor_LOD, gatk_key) = (
         PARAMS["mutect_cosmic"], PARAMS["gatk_dbsnp"],
         PARAMS["mutect_quality"], PARAMS["mutect_max_alt_qual"],
         PARAMS["mutect_max_alt"], PARAMS["mutect_max_fraction"],
         PARAMS["mutect_LOD"], PARAMS["mutect_key"])

    genome = "%s/%s.fa" % (PARAMS["bwa_index_dir"],
                           PARAMS["genome"])

    Exome.mutectSNPCaller(infile_tumour, outfile, mutect_log, genome,
                          cosmic, dbsnp, call_stats_out,
                          cluster_options, quality, max_alt_qual,
                          max_alt, max_fraction, tumor_LOD,
                          normal_panel, infile_matched=infile)

##############################################################################
##############################################################################
##############################################################################
# Variant Annotation and Recalibration
##############################################################################


@collate(realignMatchedSample,
         regex(r"bam/(\S+)-(\S+)-(\S+).realigned.bqsr.bam"),
         r"bam/\1.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for each patient,
       for use in variant calling'''
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')


@transform(runMutect,
           regex(r"variants/(\S+).mutect.snp.vcf"),
           r"variants/\1.mutect.snp.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate SNP variants using SNPeff'''
    to_cluster = USECLUSTER
    job_options = "-pe dedicated 2 -R y -l mem_free=4G"

    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''java -Xmx4G -jar /ifs/apps/bio/snpEff-3.3-dev/snpEff.jar
                   -c %(config)s -v %(snpeff_genome)s -o gatk
                   %(infile)s > %(outfile)s''' % locals()
    P.run()


@transform(indelCaller,
           regex("variants/(\S+)/results/all.somatic.indels.vcf"),
           r"variants/\1.indels.snpeff.vcf")
def annotateVariantsINDELsSNPeff(infile, outfile):
    '''Annotate INDEL variants using SNPeff'''
    to_cluster = USECLUSTER
    job_options = "-pe dedicated 2 -R y -l mem_free=4G"

    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''java -Xmx4G -jar /ifs/apps/bio/snpEff-3.3-dev/snpEff.jar
                   -c %(config)s -v %(snpeff_genome)s -o gatk
                   %(infile)s > %(outfile)s''' % locals()
    P.run()


#########################################################################
# Annotate SNP and INDEL variants
#########################################################################

# note: these annotations are currently not used.
# Need to check whether variant annotatot is using both bams
# from a single patient?
# should just be the tumour bam or else scores will be wrong!

@follows(annotateVariantsSNPeff, listOfBAMs)
@transform(runMutect,
           regex(r"variants/(\S+).mutect.snp.vcf"),
           add_inputs(r"bam/\1.list",
                      r"variants/\1.mutect.snp.snpeff.vcf"),
           r"variants/\1.mutect.snp.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''module unload apps/java/jre1.6.0_26;
                   java -Xmx2g -jar
                    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
                   -T VariantAnnotator
                   -R %%(bwa_index_dir)s/%%(genome)s.fa
                   -I %(bamlist)s
                   -A SnpEff --snpEffFile %(effFile)s
                   -o %(outfile)s
                   --variant %(infile)s
                   -L %(infile)s
                   --dbsnp %(dbsnp)s
                   -A HaplotypeScore
                   -A MappingQualityRankSumTest
                   -A ReadPosRankSumTest
                   -A AlleleBalanceBySample''' % locals()
    P.run()


@follows(annotateVariantsINDELsSNPeff, listOfBAMs)
@transform(indelCaller,
           regex("variants/(\S+)/results/all.somatic.indels.vcf"),
           add_inputs(r"bam/\1.list", r"variants/\1.indels.snpeff.vcf"),
           r"variants/\1.indels.annotated.vcf")
def variantAnnotatorIndels(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    to_cluster = USECLUSTER
    infile, bamlist, effFile = infiles
    statement = '''module unload apps/java/jre1.6.0_26;
                   java -Xmx2g -jar
                    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
                   -T VariantAnnotator
                   -R %%(bwa_index_dir)s/%%(genome)s.fa
                   -I %(bamlist)s
                   -A SnpEff --snpEffFile %(effFile)s
                   -o %(outfile)s
                   --variant %(infile)s
                   -L %(infile)s
                   -A Coverage
                   -A FisherStrand
                   -A HaplotypeScore
                   -A MappingQualityRankSumTest
                   -A ReadPosRankSumTest
                   -A AlleleBalanceBySample
                   -A RMSMappingQuality''' % locals()
    P.run()


######################################################################

# this does not work - insufficient number of indels in mills+
@transform(variantAnnotatorIndels,
           suffix(".annotated.vcf"),
           ".annotated.recalibrated.vcf")
def variantRecalibrator(infile, outfile):
    '''Create variant recalibration file for indels'''
    to_cluster = USECLUSTER
    job_options = getGATKOptions()
    job_threads = 6
    track = P.snip(os.path.basename(outfile), ".annotated.recalibrated.vcf")
    mills = PARAMS["gatk_mills"]

    statement = '''module unload apps/java/jre1.6.0_26;
                   java -Xmx4g -jar
                   /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
                   -T VariantRecalibrator
                   -R %%(bwa_index_dir)s/%%(genome)s.fa
                   -input %(infile)s
                   -resource:mills,known=true,training=true,truth=true,prior=12.0
                   %(mills)s
                   -an DP -an MQRankSum -an ReadPosRankSum
                   -mode INDEL
                   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
                   --maxGaussians 4
                   -recalFile %(outfile)s
                   -tranchesFile variants/%(track)s.tranches
                   -rscriptFile variants/%(track)s.plots.R''' % locals()
    P.run()

#########################################################################


@transform(variantAnnotatorIndels,
           suffix(".annotated.vcf"),
           ".passed.annotated.vcf")
def filterIndels(infile, outfile):
    ''' use SnpSift to filter INDELS using VCF fields'''
    statement = '''cat %(infile)s |
                   java -Xmx2g -jar /ifs/apps/bio/snpEff-3.1/SnpSift.jar filter
                   "(QSI_NT>20 & IHP<12 & RC<12 & IC<12) "
                   > %(outfile)s ''' % locals()
    P.run()


#########################################################################
#########################################################################
# convert vcf to tsv files and load into database


@transform(variantAnnotator,
           regex("variants/(\S+).annotated.vcf"),
           r"variants/\1.annotated.tsv")
def snpvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    statement = '''module unload apps/java/jre1.6.0_26;
                   java -Xmx2g -jar
                    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
                   -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -GF GT -GF AD -GF SS -GF FA -GF AB -GF DP
                   -o %(outfile)s''' % locals()
    P.run()


@transform(variantAnnotatorIndels,
           regex("variants/(\S+).annotated.vcf"),
           r"variants/\1.annotated.tsv")
def indelvcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    to_cluster = USECLUSTER
    statement = '''module unload apps/java/jre1.6.0_26;
                   java -Xmx2g -jar
                    /ifs/apps/bio/GATK-2.7-2/GenomeAnalysisTK.jar
                   -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa
                   -V %(infile)s --showFiltered --allowMissingData
                   -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER
                   -F INFO -F BaseQRankSum
                   -F HaplotypeScore -F MQRankSum -F -F ReadPosRankSum
                   -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
                   -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE
                   -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE
                   -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
                   -F TQSI -F TSQI_NT -F DP -F IC -F IHP -F NT
                   -F QSI -F QSI_NT -F RC -F RU -F SGT
                   -GF DP -GF DP2 -GF DP50 -GF SUBDP50 -GF TAR -GF TIR -GF TOR
                   -o %(outfile)s''' % locals()
    P.run()


@transform([snpvcfToTable,
            indelvcfToTable],
           regex(r"variants/(\S+).(?P<suffix>annotated.tsv|call_stats.out)"),
           r"variants/\1.\g<suffix>.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''

    if infile.endswith("indels.annotated.tsv"):
        index = "contig"
    elif infile.endswith("mutect.snp.annotated.tsv"):
        index = "CHROM"
    print "index: " + index

    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   python %%(scriptsdir)s/csv2db.py
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s''' % locals()
    P.run()


@follows(runMutect)
@transform("variants/*call_stats.out",
           regex(r"variants/(\S+)_call_stats.out"),
           r"variants/\1_call_stats.out.load")
def loadMutectExtendedOutput(infile, outfile):
    '''Load mutect extended output into database'''

    index = "CHROM, POS"

    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   python %%(scriptsdir)s/csv2db.py
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s''' % locals()
    P.run()


#########################################################################
# Genes of interest
# check this will run in the correct position if option selected

# @active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
# @transform((annotateVariantsSNPsift),
#           regex(r"variants/(\S+).haplotypeCaller.snpsift.vcf"),
#           r"variants/\1.genes.vcf")
# def findGenes(infile, outfile):
#    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
#    if variant is within a gene of interest as defined in the ini
#    file'''
#
#    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
#    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
#    statement = '''GenomeAnalysisTK -T VariantFiltration
#    -R %%(bwa_index_dir)s/%%(genome)s.fa
#    --variant %(infile)s
#    --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
#    --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
#    P.run()

#########################################################################
#########################################################################
#########################################################################
# vcf statistics -   this only summarises the nucleotide changes
# this currently does not provide useful output!


@transform((variantAnnotator,
            variantAnnotatorIndels),
           regex(r"variants/(\S+).vcf"),
           r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    to_cluster = USECLUSTER
    statement = '''vcf-stats %(infile)s
                   > %(outfile)s 2>>%(outfile)s.log;''' % locals()
    P.run()


@merge(buildVCFstats, "vcf_stats.load")
def loadVCFstats(infiles, outfile):
    '''Import variant statistics into SQLite'''
    filenames = " ".join(infiles)
    tablename = P.toTable(outfile)
    csv2db_options = PARAMS["csv2db_options"]
    E.info("Loading vcf stats...")
    statement = '''python %%(scriptsdir)s/vcfstats2db.py
                   %(filenames)s >> %(outfile)s; ''' % locals()
    statement += '''cat vcfstats.txt |
                    python %%(scriptsdir)s/csv2db.py %(csv2db_options)s
                    --allow-empty-file --add-index=track --table=vcf_stats
                    >> %(outfile)s; ''' % locals()
    P.run()

#########################################################################


@transform(runMutect,
           suffix(".mutect.snp.vcf"),
           "_mutect_filtering_summary.tsv")
def summariseFiltering(infile, outfile):
    PipelineExome.parseMutectCallStats(infile, outfile)


@transform(summariseFiltering,
           regex(r"variants/(\S+)_mutect_filtering_summary.tsv"),
           r"variants/\1_mutect_filtering_summary.load")
def loadMutectFilteringSummary(infile, outfile):
    '''Load mutect extended output into database'''

    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   python %%(scriptsdir)s/csv2db.py
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# load Network of Cancer Genes table


# parameterise file location:
@originate("cancergenes.load")
def loadNCG(outfile):
    '''Load NCG into database'''

    # infile = PARAMS["cancergenes_table"]
    infile = "../backup/NCG/cancergenes.tsv"
    index = "symbol"
    dbh = connect()
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s |
                   python %%(scriptsdir)s/csv2db.py
                   --table %(tablename)s --retry --ignore-empty
                   > %(outfile)s''' % locals()
    P.run()


#########################################################################
#########################################################################
#########################################################################
# analyse mutational siganture of filtered variants

@merge(runMutect,
       ["variants/mutational_signature.tsv",
        "variants/mutational_signature_table.tsv"])
def mutationalSignature(infiles, outfiles):

    min_t_alt = PARAMS["filter_minimum_tumor_allele"]
    min_t_alt_freq = PARAMS["filter_minimum_tumor_allele_frequency"]
    min_n_depth = PARAMS["filter_minimum_normal_depth"]
    max_n_alt_freq = PARAMS["filter_maximim_normal_allele_frequency"]
    tumour = PARAMS["mutect_tumour"]

    ExomeCancer.compileMutationalSignature(infiles, outfiles,
                                           min_t_alt, min_n_depth,
                                           max_n_alt_freq, min_t_alt_freq,
                                           tumour, submit=True)


@transform(mutationalSignature,
           suffix(".tsv"),
           ".load")
def loadMutationalSignature(infiles, outfile):
    outfile2 = re.sub(".load", "_table.load", outfile)
    P.load(infiles[0], outfile)
    P.load(infiles[1], outfile2)


#########################################################################
#########################################################################
#########################################################################

@follows(loadMutectFilteringSummary,
         loadMutectExtendedOutput,
         loadVariantAnnotation,
         loadCoverageStats,
         loadPicardRealigenedAlignStats,
         loadPicardAlignStats,
         loadNCG,
         loadMutationalSignature)
def full():
    pass


@follows(loadMutectFilteringSummary)
def test():
    pass


@follows(runMutectOnDownsampled,
         runMutectReverse)
def TestMutect():
    '''This target runs function which can be used to assess the chosen
    mutect parameters'''
    pass


# @follows(loadROI,
#         loadROI2Gene)
# def loadMetadata():
#    pass


@follows(mapReads)
def mapping():
    pass


@follows(loadPicardDuplicateStats,
         loadPicardAlignStats,
         buildCoverageStats,
         loadCoverageStats)
def postMappingQC():
    pass


@follows(GATKpreprocessing,
         runPicardOnRealigned)
def gatk():
    pass


@follows(runMutect,
         indelCaller)
def callVariants():
    pass


@follows(loadVariantAnnotation)
def tabulation():
    pass


@follows(buildVCFstats,
         loadVCFstats)
def vcfstats():
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
