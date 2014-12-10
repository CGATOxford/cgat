"""
======================================================
PipelineExome.py - common tasks for Variant Calling
======================================================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----


"""
# Import modules
import os
import re
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
import CGAT.Experiment as E

# Set PARAMS in calling module
PARAMS = {}


def getGATKOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

#########################################################################


def GATKreadGroups(infile, outfile, genome,
                   library="unknown", platform="Illumina",
                   platform_unit="1", threads=4):
    '''Reorders BAM according to reference fasta and adds read groups'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()

    statement = '''ReorderSam
                    INPUT=%(infile)s
                    OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    REFERENCE=%(genome)s
                    ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ;
                    checkpoint ;''' % locals()

    statement += '''AddOrReplaceReadGroups
                    INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    OUTPUT=%(outfile)s
                    RGLB=%(library)s
                    RGPL=%(platform)s
                    RGPU=%(platform_unit)s
                    RGSM=%(track)s
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(outfile)s ;
                    checkpoint ;''' % locals()
    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()

    P.run()


def GATKrealign(infile, outfile, genome):
    '''Realigns BAMs around indels using GATK'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()

    statement = '''GenomeAnalysisTK
                    -T RealignerTargetCreator
                    -o %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals
                    --num_threads %(threads)s
                    -R %(genome)s
                    -I %(infile)s ;
                    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
                    -T IndelRealigner
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(tmpdir_gatk)s/%(track)s.readgroups.bam
                    -targetIntervals
                    %(tmpdir_gatk)s/%(track)s.indelrealignment.intervals ;
                    checkpoint ;''' % locals()
    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()
    P.run()


def GATKrescore(infile, outfile, genome, dbsnp, solid_options=""):
    '''Recalibrates base quality scores using GATK'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()

    statement = '''GenomeAnalysisTK
                    -T BaseRecalibrator
                    --out %(tmpdir_gatk)s/%(track)s.recal.grp
                    -R %(genome)s
                    -I %(infile)s
                    --knownSites %(dbsnp)s %(solid_options)s ;
                    checkpoint ;''' % locals()

    statement += '''GenomeAnalysisTK
                    -T PrintReads -o %(outfile)s
                    -BQSR %(tmpdir_gatk)s/%(track)s.recal.grp
                    -R %(genome)s
                    -I %(tmpdir_gatk)s/%(track)s.indelrealigned.bam ;
                    checkpoint ;''' % locals()

    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()
    P.run()

#########################################################################


def haplotypeCaller(infile, outfile, genome,
                    dbsnp, intervals, padding, options):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    job_options = getGATKOptions()
    statement =  '''GenomeAnalysisTK
                    -T HaplotypeCaller
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    --dbsnp %(dbsnp)s
                    -L %(intervals)s
                    -ip %(padding)s''' % locals()
    P.run()

#########################################################################


def mutectSNPCaller(infile, outfile, mutect_log, genome, cosmic,
                    dbsnp, call_stats_out, cluster_options,
                    quality=20, max_alt_qual=150, max_alt=5,
                    max_fraction=0.05, tumor_LOD=6.3,
                    normal_panel=None, gatk_key=None,
                    infile_matched=None):
    '''Call SNVs using Broad's muTect'''
    # get mutect to work from module file without full path.

    job_options = cluster_options
    statement = '''java -Xmx4g -jar
                   /ifs/apps/bio/muTect-1.1.4/muTect-1.1.4.jar
                   --analysis_type MuTect
                   --reference_sequence %(genome)s
                   --cosmic %(cosmic)s
                   --dbsnp %(dbsnp)s
                   --input_file:tumor %(infile)s
                   --out %(call_stats_out)s
                   --enable_extended_output
                   --vcf %(outfile)s --artifact_detection_mode
                ''' % locals()
    if infile_matched:
        statement += '''--min_qscore %(quality)s
                        --gap_events_threshold 2
                        --max_alt_alleles_in_normal_qscore_sum %(max_alt_qual)s
                        --max_alt_alleles_in_normal_count %(max_alt)s
                        --max_alt_allele_in_normal_fraction %(max_fraction)s
                        --tumor_lod %(tumor_LOD)s
                        --input_file:normal %(infile_matched)s ''' % locals()
    if normal_panel:
        statement += ''' --normal_panel %(normal_panel)s ''' % locals()

    if gatk_key:
        statement += " -et NO_ET -K %(gatk_key)s " % locals()

    statement += " > %(mutect_log)s " % locals()

    P.run()


#########################################################################


def variantAnnotator(vcffile, bamlist, outfile, genome, 
                     dbsnp, annotations, snpeff_file=""):
    '''Annotate variant file using GATK VariantAnnotator'''
    job_options = getGATKOptions()
    anno = annotations.split(",")
    anno = " -A " + " -A ".join(anno)
    statement = '''GenomeAnalysisTK -T VariantAnnotator
                    -R %(genome)s
                    -I %(bamlist)s
                    -A SnpEff
                    --snpEffFile %(snpeff_file)s
                    -o %(outfile)s
                    --variant %(vcffile)s
                    -L %(vcffile)s
                    --dbsnp %(dbsnp)s
                    %(anno)s''' % locals()
    P.run()

#########################################################################


def variantRecalibrator(infile, outfile, genome,
                        dbsnp, hapmap, omni):
    '''Create variant recalibration file'''
    job_options = getGATKOptions()
    track = P.snip(outfile, ".recal")
    statement = '''GenomeAnalysisTK -T VariantRecalibrator
                    -R %(genome)s
                    -input %(infile)s
                    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s
                    -resource:omni,known=false,training=true,truth=false,prior=12.0 %(omni)s
                    -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 %(dbsnp)s
                    -an QD -an HaplotypeScore -an MQRankSum 
                    -an ReadPosRankSum -an FS -an MQ
                    --maxGaussians 4 
                    --numBadVariants 3000
                    -mode SNP
                    -recalFile %(outfile)s
                    -tranchesFile %(track)s.tranches
                    -rscriptFile %(track)s.plots.R ''' % locals()
    P.run()

#########################################################################


def applyVariantRecalibration(vcf, recal, tranches, outfile, genome):
    '''Perform variant quality score recalibration using GATK '''
    job_options = getGATKOptions()
    statement = '''GenomeAnalysisTK -T ApplyRecalibration
                    -R %(genome)s 
                    -input %(vcf)s
                    -recalFile %(recal)s
                    -tranchesFile %(tranches)s
                    --ts_filter_level 99.0
                    -mode SNP
                    -o %(outfile)s ''' % locals()
    P.run()

#########################################################################


def vcfToTable(infile, outfile, genome, columns):
    '''Converts vcf to tab-delimited file'''
    job_options = getGATKOptions()
    statement = '''GenomeAnalysisTK -T VariantsToTable 
                   -R %(genome)s
                   -V %(infile)s 
                   --showFiltered 
                   --allowMissingData
                   %(columns)s
                   -o %(outfile)s''' % locals()
    P.run()


#########################################################################


def selectVariants(infile, outfile, genome, select):
    '''Filter de novo variants based on provided jexl expression'''
    statement = '''GenomeAnalysisTK -T SelectVariants
                    -R %(genome)s
                    --variant %(infile)s
                    -select '%(select)s'
                    -log %(outfile)s.log
                    -o %(outfile)s''' % locals()
    P.run()

#########################################################################

def buildSelectStatementfromPed(filter_type, pedfile, template):
    '''Build a select statement from a template and a pedigree file'''
    pedigree = csv.DictReader(
        IOTools.openFile(pedfile), delimiter='\t', fieldnames=[
            'family', 'sample', 'father', 'mother', 'sex', 'status'])
    affecteds = []
    unaffecteds = []
    parents = []
    select = None
    # loop over pedigree file and establish relationships
    for row in pedigree:
        if row['status'] == '2':
            if filter_type == "denovo":
                father = row['father']
                mother = row['mother']
                proband = row['sample']
            elif filter_type == "dominant" or filter_type == "recessive":
                affecteds += [row['sample']]
            if filter_type == "recessive":
                parents += [row['father'], row['mother']]
        if row['status'] == '1':
            if filter_type == "dominant":
                unaffecteds += [row['sample']]
            elif filter_type == "recessive":
                if row['sample'] not in parents:
                    unaffecteds += [row['sample']]
    
    # Build select statement from template
    if filter_type == "denovo":
        select = template.replace("father", father)
        select = select.replace("mother",mother)
        select = select.replace("proband",proband)
    elif filter_type == "dominant":
        affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
        if len(unaffecteds) == 0:
            unaffecteds_exp = ''
        else:
            unaffecteds_exp = '&&vc.getGenotype("' + \
                ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
                '").isHomRef()'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp",unaffecteds_exp)
    elif filter_type == "recessive":
        affecteds_exp = '").getPL().2==0&&vc.getGenotype("'.join(affecteds)
        unaffecteds_exp = '").getPL().2!=0&&vc.getGenotype("'.join(unaffecteds)
        if len(parents) == 0:
            parents_exp = ''
        else:
            parents_exp = '&&vc.getGenotype("' + \
                ('").getPL().1==0&&vc.getGenotype("'.join(parents)) + \
                '").getPL().1==0'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp",unaffecteds_exp)
        select = select.replace("parents_exp",parents_exp)
    
    return select

#########################################################################


def guessSex(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads 
    per megabase of sequence on X and Y'''
    statement = '''calc `samtools idxstats %(infile)s | grep 'X' 
                    | awk '{print $3/($2/1000000)}'`
                    /`samtools idxstats %(infile)s | grep 'Y' 
                    | awk '{print $3/($2/1000000)}'` 
                    | tr -d " " | tr "=" "\\t" | tr "/" "\\t"
                    > %(outfile)s'''
    P.run()
