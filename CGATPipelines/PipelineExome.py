"""
======================================================
PipelineExome.py - common tasks for Variant Calling
======================================================

:Author: David Sims & Tom Smith
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
import CGAT.IOTools as IOTools
import CGAT.Pipeline as P
from CGAT.Pipeline import cluster_runnable
import numpy as np
import pandas as pd
import CGAT.CSV as csv
# Set PARAMS in calling module
PARAMS = {}


def getGATKOptions():
    return "-l mem_free=1.4G -l picard=1"

##############################################################################


def GATKReadGroups(infile, outfile, genome,
                   library="unknown", platform="Illumina",
                   platform_unit="1", track = "unknown"):
    '''Reorders BAM according to reference fasta and adds read groups'''

    if track == 'unknown':
        track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = 3

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

##############################################################################


def GATKIndelRealign(infile, outfile, genome, threads=4):
    '''Realigns BAMs around indels using GATK'''

    intervalfile = outfile.replace(".bam", ".intervals")
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T RealignerTargetCreator
                    -o %(intervalfile)s
                    --num_threads %(threads)s
                    -R %(genome)s
                    -I %(infile)s; ''' % locals()

    statement += '''GenomeAnalysisTK
                    -T IndelRealigner
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    -targetIntervals %(intervalfile)s;''' % locals()
    P.run()

##############################################################################


def GATKBaseRecal(infile, outfile, genome, dbsnp, solid_options=""):
    '''Recalibrates base quality scores using GATK'''

    track = P.snip(os.path.basename(infile), ".bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = getGATKOptions()
    job_threads = 3

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
                    -I %(infile)s ;
                    checkpoint ;''' % locals()

    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()
    P.run()

##############################################################################


def haplotypeCaller(infile, outfile, genome,
                    dbsnp, intervals, padding, options):
    '''Call SNVs and indels using GATK HaplotypeCaller in all members of a
    family together'''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK
                    -T HaplotypeCaller
                    -o %(outfile)s
                    -R %(genome)s
                    -I %(infile)s
                    --dbsnp %(dbsnp)s
                    -L %(intervals)s
                    -ip %(padding)s
                    %(options)s''' % locals()
    P.run()

##############################################################################


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

##############################################################################


def strelkaINDELCaller(infile_control, infile_tumour, outfile, genome, config,
                       outdir, cluster_options):
    '''Call INDELs using Strelka'''

    job_options = cluster_options
    statement = '''
    rm -rf %(outdir)s;
    /ifs/apps/bio/strelka-1.0.14/bin/configureStrelkaWorkflow.pl
    --normal=%(infile)s  --tumor=%(infile_tumor)s
    --ref=%(genome)s  --config=%(config)s  --output-dir=%(outdir)s;
    checkpoint ; make -j 12 -C %(outdir)s''' % locals()

    P.run()

##############################################################################


def variantAnnotator(vcffile, bamlist, outfile, genome,
                     dbsnp, annotations="", snpeff_file=""):
    '''Annotate variant file using GATK VariantAnnotator'''
    job_options = getGATKOptions()
    job_threads = 3

    if annotations!="":
        anno = annotations.split(",")
        anno = " -A " + " -A ".join(anno)
    else:
        anno = ""
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

##############################################################################


def variantRecalibrator(infile, outfile, genome, mode, dbsnp=None, 
                        kgenomes=None, hapmap=None, omni=None, mills=None):
    '''Create variant recalibration file'''
    job_options = getGATKOptions()
    job_threads = 3

    track = P.snip(outfile, ".recal")
    if mode == 'SNP':
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %(hapmap)s
        -resource:omni,known=false,training=true,truth=true,prior=12.0 %(omni)s
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %(dbsnp)s
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 %(kgenomes)s
        -an QD -an SOR -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
        P.run()
    elif mode == 'INDEL':
        statement = '''GenomeAnalysisTK -T VariantRecalibrator
        -R %(genome)s
        -input %(infile)s
        -resource:mills,known=true,training=true,truth=true,prior=12.0 %(mills)s
        -an QD -an MQRankSum
        -an ReadPosRankSum -an FS -an MQ
        --maxGaussians 4
        --minNumBadVariants 5000
        -mode %(mode)s
        -recalFile %(outfile)s
        -tranchesFile %(track)s.tranches
        -rscriptFile %(track)s.plots.R ''' % locals()
        P.run()

##############################################################################


def applyVariantRecalibration(vcf, recal, tranches, outfile, genome, mode):
    '''Perform variant quality score recalibration using GATK '''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK -T ApplyRecalibration
    -R %(genome)s
    -input:VCF %(vcf)s
    -recalFile %(recal)s
    -tranchesFile %(tranches)s
    --ts_filter_level 99.0
    -mode %(mode)s
    -o %(outfile)s ''' % locals()
    P.run()

##############################################################################


def vcfToTable(infile, outfile, genome, columns):
    '''Converts vcf to tab-delimited file'''
    job_options = getGATKOptions()
    job_threads = 3

    statement = '''GenomeAnalysisTK -T VariantsToTable
                   -R %(genome)s
                   -V %(infile)s
                   --showFiltered
                   --allowMissingData
                   %(columns)s
                   -o %(outfile)s''' % locals()
    P.run()

##############################################################################


def selectVariants(infile, outfile, genome, select):
    '''Filter de novo variants based on provided jexl expression'''
    statement = '''GenomeAnalysisTK -T SelectVariants
                    -R %(genome)s
                    --variant %(infile)s
                    -select '%(select)s'
                    -log %(outfile)s.log
                    -o %(outfile)s''' % locals()
    P.run()

##############################################################################


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
        select = select.replace("mother", mother)
        select = select.replace("proband", proband)
    elif filter_type == "dominant":
        affecteds_exp = '").getPL().1==0&&vc.getGenotype("'.join(affecteds)
        if len(unaffecteds) == 0:
            unaffecteds_exp = ''
        else:
            unaffecteds_exp = '&&vc.getGenotype("' + \
                ('").isHomRef()&&vc.getGenotype("'.join(unaffecteds)) + \
                '").isHomRef()'
        select = template.replace("affecteds_exp", affecteds_exp)
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
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
        select = select.replace("unaffecteds_exp", unaffecteds_exp)
        select = select.replace("parents_exp", parents_exp)

    return select

##############################################################################


def guessSex(infile, outfile):
    '''Guess the sex of a sample based on ratio of reads
    per megabase of sequence on X and Y'''
    statement = '''calc `samtools idxstats %(infile)s
                    | grep 'X'
                    | awk '{print $3/($2/1000000)}'`
                    /`samtools idxstats %(infile)s | grep 'Y'
                    | awk '{print $3/($2/1000000)}'`
                    | tr -d " " | tr "=" "\\t" | tr "/" "\\t"
                    > %(outfile)s'''
    P.run()

##############################################################################


# the following two functions should be generalised
# currently they operate only on mutect output
@cluster_runnable
def compileMutationalSignature(infiles, outfiles, min_t_alt, min_n_depth,
                               max_n_alt_freq, min_t_alt_freq, tumour):
    '''takes a list of mutect output files and compiles per sample mutation
    signatures'''

    delim = ":"

    def lookup(b1, b2):
        '''return lookup key for a pair of bases'''
        return(b1 + delim + b2)

    def breakKey(key):
        '''take a lookup key and return the elements'''
        return key.split(delim)

    def comp(base):
        '''return complementary base'''
        comp_dict = {"C": "G", "G": "C", "A": "T", "T": "A"}
        return comp_dict[base]

    outfile1 = open(outfiles[0], "w")
    mutations = ["C:T", "C:A", "C:G", "A:C", "A:T", "A:G"]

    outfile1.write("%s\t%s\t%s\t%s\t%s\n" % ("patient_id", "base_change",
                                             "ref", "alt", "frequency"))
    patient_freq = {}

    for infile in infiles:
        patient_id = P.snip(os.path.basename(infile), ".mutect.snp.vcf")
        mut_dict = {}
        for comb in mutations:
            mut_dict[comb] = 0

        with open(infile, "r") as f:
            for line in f.readlines():
                # need to find location of control and tumor columns
                if line.startswith('#CHROM'):
                    columns = line.split("\t")
                    for x in range(0, len(columns)):
                        if "Control" in columns[x]:
                            control_col = x
                        elif tumour in columns[x]:
                            tumor_col = x
                if not line.startswith('#'):
                    values = line.split("\t")
                    if values[6] == "PASS":
                        t_values = values[tumor_col].split(":")
                        t_ref, t_alt = map(float, (t_values[1].split(",")))
                        t_depth = t_alt + t_ref
                        n_values = values[control_col].split(":")
                        n_ref, n_alt = map(float, (n_values[1].split(",")))
                        n_depth = n_alt + n_ref
                        np.seterr(divide='ignore')
                        if (t_alt > min_t_alt and n_depth >= min_n_depth and
                            np.divide(n_alt, n_depth) <= max_n_alt_freq and
                            (((np.divide(t_alt, t_depth) /
                               np.divide(n_alt, n_depth)) >= min_t_alt_freq) or
                             ((np.divide(t_alt, t_depth) /
                               np.divide(n_alt, n_depth)) == 0))):
                                key = lookup(values[3], values[4])
                                if key in mut_dict:
                                    mut_dict[key] += 1
                                else:
                                    comp_key = lookup(
                                        comp(values[3]), comp(values[4]))
                                    mut_dict[comp_key] += 1
                        np.seterr(divide='warn')
        patient_freq[patient_id] = mut_dict

    for mutation in mutations:
        base1, base2 = breakKey(mutation)
        for infile in infiles:
            patient_id = P.snip(os.path.basename(infile), ".mutect.snp.vcf")
            outfile1.write("%s\t%s\t%s\t%s\t%s\n" % (patient_id, mutation,
                                                     base1, base2,
                                                     patient_freq[patient_id]
                                                     [mutation]))
    outfile1.close()

    outfile2 = open(outfiles[1], "w")
    outfile2.write("%s\t%s\n" % ("patient_id",
                                 "\t".join(mutations)))
    for infile in infiles:
        patient_id = P.snip(os.path.basename(infile), ".mutect.snp.vcf")
        frequencies = "\t".join(map(str, [patient_freq[patient_id][x]
                                          for x in mutations]))
        outfile2.write("%s\t%s\n" % (patient_id, frequencies))
    outfile2.close()

##############################################################################


@cluster_runnable
def parseMutectCallStats(infile, outfile):
    '''take the call stats outfile from mutect and summarise the
    reasons for variant rejection'''
    single_dict = {}
    combinations_dict = {}

    def updateDict(dictionary, key):
        if key in dictionary:
            dictionary[key] += 1
        else:
            dictionary[key] = 1

    with open(infile, "rb") as infile:
        lines = infile.readlines()
        for i, line in enumerate(lines):
            if i < 2:
                continue
                values = line.strip().split("\t")
                judgement, justification = (values[-1], values[-2])
                if judgement == "REJECT":
                    reasons = justification.split(",")
                    if len(reasons) == 1:
                        updateDict(single_dict, reasons[0])
                    else:
                        for reason in reasons:
                            updateDict(combinations_dict, reasons[0])

    df = pd.DataFrame([single_dict, combinations_dict])
    df = df.transpose()
    df.columns = ["single", "combination"]
    df = df.sort("single", ascending=False)
    df.to_csv(df, outfile, header=True, index=False)
