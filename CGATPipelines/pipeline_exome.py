"""
====================
Exome pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

.. note::

   This pipeline needs refactoring:

   1. It is very repetetive, such as all the tabulateXYZ functions
   2. Use P.load() for loading data
   3. annotateVariantsSNPsift contains a hardcoded data path

The exome pipeline imports unmapped reads from one or more fastq or
sra files and aligns them to the genome using BWA.  Post alignment
quality control is performed using Picard.  The pipeline then performs
local realignment around indels and base quality score recalibration
using GATK.  Next variants (SNVs and indels) are called, annotated,
and filtered according to various inheritance models (de novo,
dominant, and recessive).


   1. Align to genome using gapped alignment (BWA)
   2. Check alignment quality and target region coverage (Picard)
   3. Local realignment and BQSR in GATK
   4. Variant calling (SNVs & indels) using GATK HaplotypeCaller
   5. Variant annotation using SNPeff, GATK VariantAnnotator, and SnpSift
   6. Variant quality score recalibration (GATK)
   7. Flags variants within genes of interest (such as known disease genes)
      (GATK) (optional)
   8. Filters potential de novo variants
   9. Filters potential de novo variants using lower stringency criteria
   10. Filters potential dominant mutations
   11. Filters potential homozygous recessive mutations
   12. Filters potential compound heterozygous mutations
   13. Generates summary statistics for unfiltered vcf file
   14. Generates report

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

Requirements:

* picardtools >= 1.106
* samtools >= 1.1
* GATK >= 2.7
* snpEff >= 4.0
* Gemini >= ?

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Code
====

"""

# load modules
from ruffus import *
from rpy2.robjects import r as R
import sys
import os
import csv
import sqlite3
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC


#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../exome.ini", "exome.ini"])

PARAMS = P.PARAMS
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")
USECLUSTER = True


def getPicardOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"


def getGATKOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

#########################################################################
#########################################################################
#########################################################################
# Load target and sample data
# The following functions are designed to upload meta-data to the csvdb
# These haven't been fully implemented yet


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


@follows(mkdir("bam"))
@transform(INPUT_FORMATS, REGEX_FORMATS, r"bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA (output=SAM), convert to BAM,
    sort and index BAM file, generate alignment statistics and
    deduplicate using Picard
    '''
    job_options = "-pe dedicated 2 -l mem_free=8G"
    track = P.snip(os.path.basename(outfile), ".bam")
    m = PipelineMapping.BWA(
        remove_unique=PARAMS["bwa_remove_non_unique"], align_stats=True,
        dedup=True)
    statement = m.build((infiles,), outfile)
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Post-alignment QC
#########################################################################


@merge(mapReads, "picard_duplicate_stats.load")
def loadPicardDuplicateStats(infiles, outfile):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardDuplicateStats(infiles, outfile)


@follows(mapReads)
@merge("bam/*.picard_stats", "picard_stats.load")
def loadPicardAlignStats(infiles, outfile):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    PipelineMappingQC.loadPicardAlignmentStats(infiles, outfile)

#########################################################################


@transform(mapReads, regex(r"bam/(\S+).bam"), r"bam/\1.cov")
def buildCoverageStats(infile, outfile):
    '''Generate coverage statistics for regions of interest from a bed
    file using Picard'''
    baits = PARAMS["roi_baits"]
    regions = PARAMS["roi_regions"]
    statement = '''CalculateHsMetrics BAIT_INTERVALS=%(baits)s
    TARGET_INTERVALS=%(regions)s
    INPUT=%(infile)s
    OUTPUT=%(outfile)s
    VALIDATION_STRINGENCY=LENIENT''' % locals()
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
        lines = [x for x in open(f, "r").readlines()
                 if not x.startswith("#") and x.strip()]
        if first:
            outf.write("%s\t%s" % ("track", lines[0]))
        first = False
        outf.write("%s\t%s" % (track, lines[1]))
    outf.close()
    tmpfilename = outf.name
    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --add-index=track
                      --table=%(tablename)s
                      --ignore-empty
                      --retry
                   > %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
# GATK


@follows(mkdir("gatk"))
@transform(mapReads, regex(r"bam/(\S+).bam"), r"gatk/\1.bqsr.bam")
def GATKpreprocessing(infile, outfile):
    '''Reorders BAM according to reference fasta and add read groups using
    SAMtools, realigns around indels and recalibrates base quality
    scores using GATK

    '''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    library = PARAMS["readgroup_library"]
    platform = PARAMS["readgroup_platform"]
    platform_unit = PARAMS["readgroup_platform_unit"]
    threads = PARAMS["gatk_threads"]
    dbsnp = PARAMS["gatk_dbsnp"]
    solid_options = PARAMS["gatk_solid_options"]

    statement = '''ReorderSam
    INPUT=%(infile)s
    OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
    REFERENCE=%%(bwa_index_dir)s/%%(genome)s.fa
    ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ;
    checkpoint ;''' % locals()

    statement += '''AddOrReplaceReadGroups
    INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
    OUTPUT=%(tmpdir_gatk)s/%(track)s.readgroups.bam
    RGLB=%(library)s
    RGPL=%(platform)s
    RGPU=%(platform_unit)s
    RGSM=%(track)s
    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.readgroups.bam ;
    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
    -T RealignerTargetCreator
    -o %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals
    --num_threads %(threads)s
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(tmpdir_gatk)s/%(track)s.readgroups.bam ; checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
    -T IndelRealigner
    -o %(tmpdir_gatk)s/%(track)s.indelrealigned.bam
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(tmpdir_gatk)s/%(track)s.readgroups.bam
    -targetIntervals %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals ;
    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
    -T BaseRecalibrator
    --out %(tmpdir_gatk)s/%(track)s.recal.grp
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(tmpdir_gatk)s/%(track)s.indelrealigned.bam
    --knownSites %(dbsnp)s %(solid_options)s ;
    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
    -T PrintReads -o %(outfile)s
    -BQSR %(tmpdir_gatk)s/%(track)s.recal.grp
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(tmpdir_gatk)s/%(track)s.indelrealigned.bam ;
    checkpoint ;''' % locals()

    statement += '''rm -rf %(tmpdir_gatk)s ;'''
    P.run()


@collate(GATKpreprocessing, regex(r"gatk/(\S+?)-(\S+).bqsr.bam"),
         r"gatk/\1.list")
def listOfBAMs(infiles, outfile):
    '''generates a file containing a list of BAMs for each family, for use
    in variant calling

    '''
    with IOTools.openFile(outfile, "w") as outf:
        for infile in infiles:
            outf.write(infile + '\n')

#########################################################################
#########################################################################
#########################################################################
# Variant Calling

#########################################################################


@follows(mkdir("variants"))
@transform(listOfBAMs, regex(r"gatk/(\S+).list"),
           r"variants/\1.haplotypeCaller.vcf")
def haplotypeCaller(infile, outfile):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    job_options = getGATKOptions()
    dbsnp = PARAMS["gatk_dbsnp"]
    intervals = PARAMS["roi_intervals"]
    padding = PARAMS["roi_padding"]
    hc_options = PARAMS["gatk_hc_options"]
    statement = '''GenomeAnalysisTK
    -T HaplotypeCaller
    -o %(outfile)s
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(infile)s
    --dbsnp %(dbsnp)s
    -L %(intervals)s
    -ip %(padding)s''' % locals()
    P.run()

##########################################################################
##########################################################################
##########################################################################
# Variant Annotation and Recalibration


@transform(haplotypeCaller,
           regex(r"variants/(\S+).haplotypeCaller.vcf"),
           r"variants/\1.haplotypeCaller.snpeff.vcf")
def annotateVariantsSNPeff(infile, outfile):
    '''Annotate variants using SNPeff'''
    job_options = "-pe dedicated 4 -R y -l mem_free=6G"
    snpeff_genome = PARAMS["annotation_snpeff_genome"]
    config = PARAMS["annotation_snpeff_config"]
    statement = '''snpEff.sh eff
    -c %(config)s
    -v %(snpeff_genome)s
    -o gatk %(infile)s > %(outfile)s''' % locals()
    P.run()


@follows(annotateVariantsSNPeff)
@transform(haplotypeCaller,
           regex(r"variants/(\S+).haplotypeCaller.vcf"),
           add_inputs(r"bam/\1.list",
                      r"variants/\1.haplotypeCaller.snpeff.vcf"),
           r"variants/\1.haplotypeCaller.annotated.vcf")
def variantAnnotator(infiles, outfile):
    '''Annotate variant file using GATK VariantAnnotator'''
    infile, bamlist, effFile = infiles
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK
    -T VariantAnnotator
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    -I %(bamlist)s
    -A SnpEff
    --snpEffFile %(effFile)s
    -o %(outfile)s
    --variant %(infile)s
    -L %(infile)s
    --dbsnp %(dbsnp)s
    -A HaplotypeScore
    -A MappingQualityRankSumTest
    -A ReadPosRankSumTest
    -A AlleleBalanceBySample''' % locals()
    P.run()


@transform(variantAnnotator,
           regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"),
           r"variants/\1.haplotypeCaller.vqsr.recal")
def variantRecalibrator(infile, outfile):
    '''Create variant recalibration file'''
    job_options = getGATKOptions()
    track = P.snip(os.path.basename(outfile), ".recal")
    hapmap = PARAMS["gatk_hapmap"]
    omni = PARAMS["gatk_omni"]
    dbsnp = PARAMS["gatk_dbsnp"]
    statement = '''GenomeAnalysisTK
    -T VariantRecalibrator
    -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s
    -resource:omni,known=false,training=true,truth=false,prior=12.0 %(omni)s
    -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %(dbsnp)s
    -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ
    --maxGaussians 4 --numBadVariants 3000
    -mode SNP
    -recalFile %(outfile)s
    -tranchesFile %(track)s.tranches
    -rscriptFile %(track)s.plots.R ''' % locals()
    P.run()

#########################################################################


@follows(variantRecalibrator)
@transform(variantAnnotator,
           regex(r"variants/(\S+).haplotypeCaller.annotated.vcf"),
           add_inputs(r"variants/\1.haplotypeCaller.vqsr.recal",
                      r"\1.haplotypeCaller.vqsr.tranches"),
           r"variants/\1.haplotypeCaller.vqsr.vcf")
def applyVariantRecalibration(infiles, outfile):
    '''Perform variant quality score recalibration using GATK '''
    job_options = getGATKOptions()
    infile, recal, tranches = infiles
    statement = '''GenomeAnalysisTK -T ApplyRecalibration
    -R %%(bwa_index_dir)s/%%(genome)s.fa -input %(infile)s
    --ts_filter_level 99.0
    -tranchesFile %(tranches)s
    -recalFile %(recal)s
    -mode SNP
    -o %(outfile)s ''' % locals()
    P.run()

#########################################################################


@transform(applyVariantRecalibration,
           regex(r"variants/(\S+).haplotypeCaller.vqsr.vcf"),
           r"variants/\1.haplotypeCaller.snpsift.vcf")
def annotateVariantsSNPsift(infile, outfile):
    '''Add annotations using SNPsift'''
    job_options = "-pe dedicated 4 -R y -l mem_free=6G"
    track = P.snip(os.path.basename(infile), ".vqsr.vcf")
    dbNSFP = PARAMS["annotation_snpsift_dbnsfp"]
# The following statement is not fully implemented yet
#    statement = '''SnpSift.sh geneSets -v /ifs/projects/proj016/data/1000Genomes/msigdb.v4.0.symbols.gmt %(infile)s > variants/%(track)s_temp1.vcf; checkpoint;''' % locals()

    statement = '''SnpSift.sh dbnsfp -v %(dbNSFP)s %(infile)s
    > variants/%(track)s_temp1.vcf; checkpoint;''' % locals()

    raise NotImplementedError(
        "this method contains hardcoded path to data!")
    statement += '''SnpSift.sh annotate
    /ifs/projects/proj016/data/1000Genomes/00-All.vcf
    variants/%(track)s_temp1.vcf > %(outfile)s ;''' % locals()
#    statement += '''rm -f variants/*temp*vcf;'''

    P.run()

#########################################################################
#########################################################################
#########################################################################
# Genes of interest


@active_if(PARAMS["annotation_add_genes_of_interest"] == 1)
@transform((annotateVariantsSNPsift),
           regex(r"variants/(\S+).haplotypeCaller.snpsift.vcf"),
           r"variants/\1.genes.vcf")
def findGenes(infile, outfile):
    '''Adds expression "GENE_OF_INTEREST" to the FILTER column of the vcf
    if variant is within a gene of interest as defined in the ini
    file

    '''
    geneList = P.asList(PARAMS["annotation_genes_of_interest"])
    expression = '\'||SNPEFF_GENE_NAME==\''.join(geneList)
    statement = '''GenomeAnalysisTK -T VariantFiltration
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    --variant %(infile)s
    --filterExpression "SNPEFF_GENE_NAME=='%(expression)s'"
    --filterName "GENE_OF_INTEREST" -o %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Tabulation
CALIBRATION = {0: annotateVariantsSNPsift,
               1: findGenes}


@transform(CALIBRATION[PARAMS["annotation_add_genes_of_interest"]],
           regex(r"variants/((\S+).haplotypeCaller.snpsift|(\S+).genes).vcf"),
           r"variants/\1.table")
def vcfToTable(infile, outfile):
    '''Converts vcf to tab-delimited file'''
    # add AB and DP annotations to table
    statement = '''GenomeAnalysisTK
    -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa
    -V %(infile)s --showFiltered --allowMissingData
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN
    -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC
    -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB
    -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
    -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME
    -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
    -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid
    -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score
    -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF
    -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF
    -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF
    -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD
    -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD
    -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals()
    P.run()


@transform(vcfToTable, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadVariantAnnotation(infile, outfile):
    '''Load VCF annotations into database'''
    # AH: why not use P.load()
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
    --table %(tablename)s --retry --ignore-empty > %(outfile)s''' % locals()
    P.run()


@transform(annotateVariantsSNPeff,
           regex(r"variants/(\S+).haplotypeCaller.snpeff.vcf"),
           r"variants/\1.haplotypeCaller.snpeff.table")
def snpeffToTable(infile, outfile):
    '''Converts snpeff file to table as this supplies all the annotations
    for a given variant - used for reports'''
    statement = '''GenomeAnalysisTK -T VariantsToTable
    -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s
    --showFiltered --allowMissingData
    -F CHROM -F POS -F ID -F REF -F ALT -F EFF
    -o %(outfile)s''' % locals()
    P.run()


@transform(snpeffToTable, regex(r"variants/(\S+).table"),
           r"variants/\1.table.load")
def loadSnpeffAnnotation(infile, outfile):
    '''Load snpeff annotations into database'''
    # AH: why not P.load()
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
    | python %(scriptsdir)s/csv2db.py --table %(tablename)s
    --retry --ignore-empty > %(outfile)s''' % locals()
    P.run()


# the following function is not working - not sure why yet
@follows(loadVariantAnnotation, loadSnpeffAnnotation)
@transform(loadSnpeffAnnotation,
           regex(r"variants/(\S+).haplotypeCaller.snpeff.table.load"),
           r"variants/\1.snpeff_snpsift.table.load")
def createAnnotationsTable(infile, outfile):
    '''Create new annotations table from snpeff and snpsift tables'''
    dbh = sqlite3.connect(PARAMS["database"])
    table = P.toTable(outfile)
    track = P.snip(os.path.basename(outfile), ".snpeff_snpsift.table.load")
    setrack = track.replace('.', '_')
    if PARAMS["annotation_add_genes_of_interest"] == 1:
        sstrack = setrack + '_genes'
    else:
        sstrack = setrack + 'haplotypeCaller_snpsift'
    cc = dbh.cursor()
    cc.execute("""DROP TABLE IF EXISTS %(table)s """ % locals())
    cc.execute("""CREATE TABLE %(table)s AS SELECT *
    FROM %(sstrack)s_table
    INNER JOIN %(setrack)s_haplotypeCaller_snpeff_table
    ON %(sstrack)s_table.CHROM = %(setrack)s_haplotypeCaller_snpeff_table.CHROM
    AND %(sstrack)s_table.POS = %(setrack)s_haplotypeCaller_snpeff_table.POS
    AND %(sstrack)s_table.REF = %(setrack)s_haplotypeCaller_snpeff_table.REF
    AND %(sstrack)s_table.ALT = %(setrack)s_haplotypeCaller_snpeff_table.ALT
    """ % locals())
    cc.close()
    P.touch(outfile)


#########################################################################
#########################################################################
#########################################################################
# De novos


@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"), r"variants/\1.filtered.vcf")
def deNovoVariants(infiles, outfile):
    '''Filter de novo variants based on provided jexl expression'''
    infile, pedfile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    statement = '''GenomeAnalysisTK
    -T SelectVariants
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    --variant %(infile)s
    -select 'vc.getGenotype("%(father)s").getDP()>=10&&vc.getGenotype("%(mother)s").getDP()>=10&&vc.getGenotype("%(child)s").getPL().0>20&&vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(child)s").getPL().2>0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(father)s").getPL().1>20&&vc.getGenotype("%(father)s").getPL().2>20&&vc.getGenotype("%(mother)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().1>20&&vc.getGenotype("%(mother)s").getPL().2>20&&vc.getGenotype("%(child)s").getAD().1>=3&&((vc.getGenotype("%(child)s").getAD().1)/(vc.getGenotype("%(child)s").getDP().floatValue()))>=0.25&&(vc.getGenotype("%(father)s").getAD().1==0||(vc.getGenotype("%(father)s").getAD().1>0&&((vc.getGenotype("%(father)s").getAD().1)/(vc.getGenotype("%(father)s").getDP().floatValue()))<0.05))&&(vc.getGenotype("%(mother)s").getAD().1==0||(vc.getGenotype("%(mother)s").getAD().1>0&&((vc.getGenotype("%(mother)s").getAD().1)/(vc.getGenotype("%(mother)s").getDP().floatValue()))<0.05))&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")'
    > %(outfile)s''' % locals()

    P.run()

#########################################################################


@transform(deNovoVariants,
           regex(r"variants/(\S+).filtered.vcf"),
           r"variants/\1.filtered.table")
def tabulateDeNovos(infile, outfile):
    '''Tabulate de novo variants'''
    statement = '''awk 'NR > 16' %(infile)s | sed '/^INFO/d' > %(infile)s.tmp ;'''
    statement += '''GenomeAnalysisTK
    -T VariantsToTable
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    --variant %(infile)s.tmp
    --showFiltered
    --allowMissingData
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN
    -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC
    -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB
    -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
    -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME
    -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
    -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid
    -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score
    -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF
    -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF
    -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF
    -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD
    -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD
    -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s ;''' % locals()
    statement += '''rm -f %(infile)s.tmp*'''
    P.run()


@transform(tabulateDeNovos,
           regex(r"variants/(\S+).filtered.table"),
           r"variants/\1.filtered.table.load")
def loadDeNovos(infile, outfile):
    '''Load de novos into database'''
    # AH: why not P.load()
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
    | python %(scriptsdir)s/csv2db.py --table %(tablename)s
    --retry --ignore-empty --allow-empty-file > %(outfile)s''' % locals()
    P.run()


@transform(annotateVariantsSNPsift,
           regex(r"variants/(\S*Trio\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"), r"variants/\1.denovos.vcf")
def lowerStringencyDeNovos(infiles, outfile):
    '''Filter lower stringency de novo variants based on provided jexl
    expression'''
    infile, pedfile = infiles
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    for row in pedigree:
        if row['status'] == '2':
            father = row['father']
            mother = row['mother']
            child = row['sample']
    statement = '''GenomeAnalysisTK
    -T SelectVariants
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    --variant %(infile)s
    -select 'vc.getGenotype("%(child)s").getPL().1==0&&vc.getGenotype("%(father)s").getPL().0==0&&vc.getGenotype("%(mother)s").getPL().0==0&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' > %(outfile)s''' % locals()
    P.run()

#########################################################################


@transform(lowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.vcf"),
           r"variants/\1.denovos.table")
def tabulateLowerStringencyDeNovos(infile, outfile):
    '''Tabulate lower stringency de novo variants'''
    statement = '''awk 'NR > 16' %(infile)s | sed '/^INFO/d' > %(infile)s.tmp ;'''
    statement += '''GenomeAnalysisTK
    -T VariantsToTable
    -R %%(bwa_index_dir)s/%%(genome)s.fa
    --variant %(infile)s.tmp
    --showFiltered
    --allowMissingData
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN
    -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC
    -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB
    -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS
    -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME
    -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID
    -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid
    -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score
    -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF
    -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF
    -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF
    -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD
    -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD
    -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s ;''' % locals()
    statement += '''rm -f %(infile)s.tmp*'''
    P.run()


@transform(tabulateLowerStringencyDeNovos,
           regex(r"variants/(\S+).denovos.table"),
           r"variants/\1.denovos.table.load")
def loadLowerStringencyDeNovos(infile, outfile):
    '''Load lower stringency de novos into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s
    | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry
    --ignore-empty --allow-empty-file > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Dominant


@transform(annotateVariantsSNPsift, regex(r"variants/(\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
           add_inputs(r"\1.ped"), r"variants/\1.dominant.vcf")
def dominantVariants(infiles, outfile):
    '''Filter variants according to autosomal dominant disease model'''
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
    # currently 1000G filter is performed at the report stage (not in csvdb)
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -o %(outfile)s -select 'vc.getGenotype("%(affecteds_exp)s").getPL().1==0%(unaffecteds_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' ;''' % locals(
    )
    P.run()

#########################################################################


@transform(dominantVariants, regex(
    r"variants/(\S+).dominant.vcf"), r"variants/\1.dominant.table")
def tabulateDoms(infile, outfile):
    '''Tabulate dominant disease candidate variants'''
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant:VCF %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(tabulateDoms, regex(r"variants/(\S+).dominant.table"),
           r"variants/\1.dominant.table.load")
def loadDoms(infile, outfile):
    '''Load dominant disease candidates into database'''
    # AH: Why not P.load()
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty-file > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# Recessive


@transform(
    annotateVariantsSNPsift,
    regex(
        r"variants/(\S*Trio\S+|\S*Multiplex\S+).haplotypeCaller.snpsift.vcf"),
    add_inputs(r"\1.ped"), r"variants/\1.recessive.vcf")
def recessiveVariants(infiles, outfile):
    '''Filter variants according to autosomal recessive disease model'''
    infile, pedfile = infiles
    pedigree = csv.DictReader(open(pedfile), delimiter='\t', fieldnames=[
                              'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    parents = []
    unaffecteds = []
    for row in pedigree:
        if row['status'] == '2':
            affecteds += [row['sample']]
            parents += [row['father'], row['mother']]
        elif row['status'] == '1' and row['sample'] not in parents:
            unaffecteds += [row['sample']]
    affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
    unaffecteds_exp = '").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)
    if len(parents) == 0:
        parents_exp = ''
    else:
        parents_exp = '&&vc.getGenotype("' + \
            ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
            '").getPL().1==0'
    # need a way of specifying that other unaffecteds eg. sibs can't be
    # homozygous for alt allele
    statement = '''GenomeAnalysisTK -T SelectVariants -R %%(bwa_index_dir)s/%%(genome)s.fa --variant %(infile)s -o %(outfile)s -select 'vc.getGenotype("%(affecteds_exp)s").getPL().2==0&&vc.getGenotype("%(unaffecteds_exp)s").getPL().2!=0%(parents_exp)s&&(SNPEFF_IMPACT=="HIGH"||SNPEFF_IMPACT=="MODERATE")' ;''' % locals(
    )
    P.run()

#########################################################################


@transform(recessiveVariants, regex(
    r"variants/(\S+).recessive.vcf"), r"variants/\1.recessive.table")
def tabulateRecs(infile, outfile):
    '''Tabulate potential homozygous recessive disease variants'''
    statement = '''GenomeAnalysisTK -T VariantsToTable -R %%(bwa_index_dir)s/%%(genome)s.fa --variant:VCF %(infile)s --showFiltered --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AF -F AN -F BaseQRankSum -F DB -F DP -F Dels -F FS -F HaplotypeScore -F MLEAC -F MLEAF -F MQ -F MQ0 -F MQRankSum -F QD -F ReadPosRankSum -F SB -F SNPEFF_EFFECT -F SNPEFF_IMPACT -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_CODON_CHANGE -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_GENE_NAME -F SNPEFF_GENE_BIOTYPE -F SNPEFF_TRANSCRIPT_ID -F SNPEFF_EXON_ID -F dbNSFP_GERP++_RS -F dbNSFP_GERP++_NR -F dbNSFP_Ensembl_transcriptid -F dbNSFP_Uniprot_acc -F dbNSFP_Interpro_domain -F dbNSFP_SIFT_score -F dbNSFP_Polyphen2_HVAR_pred -F dbNSFP_29way_logOdds -F dbNSFP_1000Gp1_AF -F dbNSFP_1000Gp1_AFR_AF -F dbNSFP_1000Gp1_EUR_AF -F dbNSFP_1000Gp1_AMR_AF -F dbNSFP_1000Gp1_ASN_AF -F dbNSFP_ESP6500_AA_AF -F dbNSFP_ESP6500_EA_AF -F RSPOS -F SSR -F SAO -F VP -F VC -F PM -F TPA -F PMC -F MUT -F VLD -F OTHERKG -F PH3 -F CDA -F MTP -F OM -F CAF -F COMMON -GF GT -GF AD -GF GQ -GF PL -GF PQ -GF TP -GF AB -GF DP -o %(outfile)s''' % locals(
    )
    P.run()

#########################################################################


@transform(tabulateRecs, regex(r"variants/(\S+).recessive.table"),
           r"variants/\1.recessive.table.load")
def loadRecs(infile, outfile):
    '''Load homozygous recessive disease candidates into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty-file > %(outfile)s''' % locals()
    P.run()

#########################################################################


@transform(annotateVariantsSNPeff, regex(r"variants/(\S*Multiplex\S+|\S*Trio\S+).haplotypeCaller.snpeff.vcf"),
           add_inputs(r"\1.ped", r"bam/\1.list"), r"variants/\1.compound_hets.table")
def compoundHets(infiles, outfile):
    '''Identify potentially pathogenic compound heterozygous variants (phasing with GATK followed by compound het
    calling using Gemini'''
    infile, pedfile, bamlist = infiles
    statement = '''GenomeAnalysisTK -T PhaseByTransmission -R %%(bwa_index_dir)s/%%(genome)s.fa -V %(infile)s -ped %(pedfile)s -mvf %(infile)s.mvf -o %(infile)s.phased.vcf ;''' % locals(
    )
    statement += '''GenomeAnalysisTK -T ReadBackedPhasing -R %%(bwa_index_dir)s/%%(genome)s.fa -I %(bamlist)s -V %(infile)s.phased.vcf -o %(infile)s.phased_rbp.vcf --respectPhaseInInput ;''' % locals(
    )
    statement += '''gemini load -v %(infile)s.phased_rbp.vcf -p %(pedfile)s -t snpEff %(infile)s.db ;'''
    statement += '''gemini comp_hets --only-affected --method=filter --filter-method "(impact_severity = 'HIGH' OR impact_severity = 'MED') AND (in_esp = 0 OR aaf_esp_all < 0.01) AND (in_1kg = 0 OR aaf_1kg_all < 0.01)" %(infile)s.db > %(outfile)s'''
    P.run()

#########################################################################


@transform(compoundHets, regex(r"variants/(\S+).compound_hets.table"),
           r"variants/\1.compound_hets.table.load")
def loadCompoundHets(infile, outfile):
    '''Load compound heterozygous variants into database'''
    scriptsdir = PARAMS["general_scriptsdir"]
    tablename = P.toTable(outfile)
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py --table %(tablename)s --retry --ignore-empty --allow-empty-file > %(outfile)s''' % locals()
    P.run()

#########################################################################
#########################################################################
#########################################################################
# vcf statistics


@transform((annotateVariantsSNPsift), regex(
    r"variants/(\S+).vcf"), r"variants/\1.vcfstats")
def buildVCFstats(infile, outfile):
    '''Calculate statistics on VCF file'''
    statement = '''vcf-stats %(infile)s > %(outfile)s 2>>%(outfile)s.log;''' % locals(
    )
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
    statement += '''cat vcfstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=vcf_stats >> %(outfile)s; '''
    statement += '''cat sharedstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=vcf_shared_stats >> %(outfile)s; '''
    statement += '''cat indelstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=indel_stats >> %(outfile)s; '''
    statement += '''cat snpstats.txt | python %(scriptsdir)s/csv2db.py %(csv2db_options)s --allow-empty-file --add-index=track --table=snp_stats >> %(outfile)s; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@follows(loadROI,
         loadROI2Gene,
         loadSamples)
def loadMetadata():
    pass


@follows(mapReads,
         loadPicardDuplicateStats,
         loadPicardAlignStats,
         buildCoverageStats,
         loadCoverageStats)
def mapping():
    pass


def postMappingQC(mapping):
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
         loadVariantAnnotation,
         snpeffToTable,
         loadSnpeffAnnotation,
         createAnnotationsTable)
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
