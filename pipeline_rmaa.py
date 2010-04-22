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

:Author: Tildon Grant Belgard
:Release: $Id: pipeline_kamilah.py 2869 2010-03-03 10:20:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Rmaa mapping pipeline.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Notes
-----

The original pipeline was implemented by Tildon Grant Belgard. The code below
has been adapted for the MRC by Andreas Heger.

Additional features added:
   * fastq files are compressed

.. todo::

   * There are separate execution paths for ``nonmulti`` and ``multi``. These
   could probably be merged by appropriate parameterizing.
   * Use scratch disk for fastq files

Code
----

"""


# load modules
from ruffus import *

import sys, os, re, shutil, itertools, math, glob, logging
from ConfigParser import *

# for plotting
import matplotlib
import pylab
import numpy

# from basic_headers import run_cmd, run_cmd_output, run_cmd_qsub
# from commands import getoutput

import Pipeline as P
import Experiment as E

PARAMS = P.getParameters()

# load options from the config file
configdict=ConfigParser()
configdict.read('./rmaa.conf')

#TOPHAT_PATH=configdict.defaults()['tophat_path']
REMOVE_BASES_FROM_RIGHT=int(configdict.get('GENERAL','remove_bases_from_right')) # how many bases to cut off from the export file (for example base 51 used for phasing in a 50 bp PE read)
LIBRARY_CODES=dict(configdict.items('LIBRARIES'))
SEQUENCE_START_DATE=str(configdict.get('GENERAL','sequence_start_date'))
NEW_QUALS=configdict.getboolean('GENERAL','new_quals')
MAX_INTRON=configdict.get('GENERAL','max_intron')
JUNCTIONS_FILE=configdict.get('GENERAL','junctions_file')
GENOME=configdict.get('GENERAL','genome')
GENES_FILE=configdict.get('GENERAL','genes_file')
RMAAPATH=os.getenv("RMAAPATH")
DATABASES=dict(configdict.items('DATABASES'))
FOLD_DIFFS=str(configdict.get('GENERAL','fold_diffs')).split(',')
MIN_READS_TO_CLASSIFY=configdict.get('GENERAL','min_reads_to_classify')
K_NEIGHB=str(configdict.get('GENERAL','k_neighb')).split(',')
REF_GTF=configdict.get('GENERAL','ref_gtf')
SEQUENCE_DIR=configdict.get('GENERAL','sequence_dir')
CHROM_SIZES=configdict.get('GENERAL','chrom_sizes')
MIN_RPKM_TREE=str(configdict.get('GENERAL','min_rpkm_tree')).split(',')
MIN_GENE_LEN=configdict.get('GENERAL','min_gene_len')
MIN_RPKM=str(configdict.get('GENERAL','min_rpkm')).split(',')
MIN_DIFF=str(configdict.get('GENERAL','min_diff')).split(',')
MIN_RPKM_VAR=str(configdict.get('GENERAL','min_rpkm_var'))
MIN_RPKM_TERM_VAR=str(configdict.get('GENERAL','min_rpkm_term_var'))
MIN_GENES_TERM_VAR=str(configdict.get('GENERAL','min_genes_term_var'))
CUFF_MIN_ISOFORM=str(configdict.get('GENERAL','cuff_min_isoform'))
CUFF_PRE_MRNA=str(configdict.get('GENERAL','cuff_pre_mrna'))
CLUSTER_LVL=str(configdict.get('GENERAL','cluster_lvl')).split(',')
FDR=configdict.get('GENERAL','fdr')
MIN_GENES_TO_CLUSTER=configdict.get('GENERAL','min_genes_to_cluster')   # min # genes required to cluster the results

# there can be several samples per tissue
parameters = ( ["reads/tissue1/sample1.export.gz", 
                ("reads/tissue1/sample1.1.fq.gz", 
                 "reads/tissue1/sample1.2.fq.gz") ],
               ["reads/tissue2/sample1.export.gz", 
                ("reads/tissue2/sample1.1.fq.gz", 
                 "reads/tissue2/sample1.2.fq.gz") ],
               )
parameters2 = ( ( ["reads/tissue1/sample1.export.gz",], "reads/tissue1/insert_sizes"),
                ( ["reads/tissue2/sample1.export.gz",], "reads/tissue2/insert_sizes"), )

USECLUSTER=False

@files(parameters)
def exportToFastQ( infile, outfiles):
    """
    Creates fastq files of paired-ended reads from export files.
 
    """

    to_cluster = USECLUSTER

    outfile1, outfile2 = outfiles
    statement = '''
    python %(rmaadir)s/exports_to_fq.py 
          %(infile)s 
          %(outfile1)s 
          %(outfile2)s 
          %(remove_bases_from_right)i 
          %(new_quals)s
    ''' 
   
    P.run()

@files(parameters2)
def getInsertSize( infiles, outfile):
    """
    Plots the internal insert size distribution and calculates the average and standard deviation based on the FWHM
    """
    
    infiles = " ".join(infiles)

    to_cluster = USECLUSTER

    statement = '''
    zcat %(infiles)s | python %(rmaadir)s/return_insert_sizes.py > %(outfile)s
    '''
    P.run()

    ins_sizes_array=numpy.array( [map(int, x[:-1].split("\t")) for x in open(outfile, "r")] )
    max_freq=ins_sizes_array[:,1].max()
    half_max=float(max_freq)/2.0
    E.info( "maximum frequency=%i, halfwidth=%i" % (max_freq, half_max))

    # get half width coordinates
    for bin, value in ins_sizes_array:
        if value < half_max: continue
        FWHMmin=bin
        break

    for bin, value in ins_sizes_array[::-1]:
        if value < half_max: continue
        FWHMmax=bin
        break

    FWHM=FWHMmax-FWHMmin
    std_dev=int(float(FWHM)/2.3548)
    ins_size=int(FWHMmin+float(FWHM)/2.0)-REMOVE_BASES_FROM_RIGHT

    my_str= "".join(["For ", ",".join(infiles), " FWHM is ", str(FWHM), " ranging from ", str(FWHMmin), " to ", str(FWHMmax), ". std dev ", str(std_dev), " and ins size ", str(ins_size)])
    loglibraries.write(my_str)

    E.info( "".join(["For ", infiles, " FWHM is ", str(FWHM), " ranging from ", str(FWHMmin), " to ", str(FWHMmax), ". std dev ", 
                     str(std_dev), " and ins size ", str(ins_size)] ) )

    x, y= [], []
    
    for bin,value in ins_sizes_array:
        if FWHMmin - 2 * std_dev < bin < FWHMmax + 2 * std_dev:
            x.append(bin)
            y.append(value)

    pylab.title("Insert size")
    pylab.xlabel('inner distance between sequenced ends')
    pylab.ylabel('frequency based on unique eland mappings')
    pylab.scatter(x,y)
    pylab.savefig(outfile + ".png")

    fwhm_file=open(outfile + ".txt", 'w')
    my_str='%s\t%s\n' % (ins_size, std_dev)
    fwhm_file.write(my_str)
    fwhm_file.close()

def getInsertSizes( dirname ):

    with open( '%s/insert_sizes.txt' % dirname, 'r') as tmp_f:
       ins_size, std_dev = map( int, tmp_f.readline().rstrip('\n').split('\t'))
       
    return ins_size, std_dev

@follows( exportToFastQ, getInsertSize, mkdir("logs/juncs"))
@collate( exportToFastQ, regex(r"reads\/(.+)\/(.+)\..\.fq.gz$"), r'reads/\1/\2.junctions' )
def findJunctions(infiles, outfile):
    '''map reads using all known junctions in order to identify new possible junctions - 
    cat the junctions together and delete the tophat output directories
    '''

    ins_size, std_dev = getInsertSizes( os.path.dirname( outfile[:-len(".junctions") ] ) )

    nslots = 4
    fastq1, fastq2 = infiles[0]

    tmpfilename = P.getTempFilename()
    
    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

    job_options= "-pe dedicated 4-8 -l mem_free=3G -R y"
    
    to_cluster = USECLUSTER

    # tophat does a seek operation on the fq files, hence they
    # need to unpacked into real files
    statement = '''
    gunzip < %(fastq1)s > %(tmpfilename)s.1.fq;
    gunzip < %(fastq2)s > %(tmpfilename)s.2.fq;
    tophat --output-dir %(tmpfilename)s
           --butterfly-search 
           --min-anchor-length 5 
           --min-isoform-fraction 0.0 
           --mate-inner-dist %(ins_size)i 
           --mate-std-dev %(std_dev)i 
           --max-intron-length %(max_intron)i 
           --raw-juncs %(junctions_file)s 
           --closure-search 
           --microexon-search 
           -p %(nslots)i 
           %(bowtiedir)s/%(genome)s
           %(tmpfilename)s.1.fq
           %(tmpfilename)s.2.fq
    >& %(outfile)s.log;
    mv %(tmpfilename)s/junctions.bed %(outfile)s >& %(outfile)s.log2;
    mv %(tmpfilename)s/logs %(outfile)s.logs >& %(outfile)s.log3;
    rm -rf %(tmpfilename)s %(tmpfilename)s.1.fq %(tmpfilename)s.2.fq >& %(outfile)s.log4
    ''' 
    P.run()

@merge(findJunctions, "reads/all.junctions")
def combineJunctions(infiles, outfile):
    '''collate all junctions found with tophat together.'''
    
    infiles = " ".join(infiles)
    statement = '''
    cat %(infiles)s
    | grep -v description 
    | python %(rmaadir)s/combine_junctions.py 
    | sort 
    | uniq 
    > %(outfile)s'''
    P.run()

@follows( combineJunctions)
@collate( exportToFastQ, 
          regex(r"reads\/(.+)\/(.+)\..\.fq.gz$"), 
          r'reads/\1/\2.bam' )
def mapReads(infiles, outfile):
    '''map reads using all known junctions and all junctions found before.

    This method requires the explicit genome in bowtiedir together with the
    samtools index. Using a flattened genome file will not work due to 
    the limit of a line length of 65536 in samtools.
    '''

    if not os.path.exists( "%(bowtiedir)s/%(genome)s.fa" % PARAMS ):
       raise ValueError( "genome %(bowtiedir)s/%(genome)s.fa does not exist - create with bowtie-inspect first" % PARAMS)

    ins_size, std_dev = getInsertSizes( os.path.dirname( outfile[:-len(".bam") ] ) )
    
    nslots = 4
    fastq1, fastq2 = infiles[0]

    tmpfilename = P.getTempFilename()
    
    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

    job_options= "-pe dedicated 4-8 -l mem_free=3G -R y"

    to_cluster = USECLUSTER

    junctions_file = "reads/all.junctions"

    statement = '''
    gunzip < %(fastq1)s > %(tmpfilename)s.1.fq;
    gunzip < %(fastq2)s > %(tmpfilename)s.2.fq;
    tophat --output-dir %(tmpfilename)s
           --min-isoform-fraction 0.0 
           --mate-inner-dist %(ins_size)i 
           --mate-std-dev %(std_dev)i 
           --raw-juncs %(junctions_file)s 
           -p %(nslots)i 
           %(bowtiedir)s/%(genome)s
           %(tmpfilename)s.1.fq
           %(tmpfilename)s.2.fq
    >& %(outfile)s.log;
    samtools view -bT %(bowtiedir)s/%(genome)s.fa %(tmpfilename)s/accepted_hits.sam 2>%(outfile)s.log1 > %(outfile)s; 
    rm -f %(tmpfilename)s/accepted_hits.sam >& %(outfile)s.log2;
    mv %(tmpfilename)s %(outfile)s.final >& %(outfile)s.log3
    ''' 

    P.run()

@collate( mapReads,
          regex(r"reads/(.+)/(.+).bam$"),
          r'mappings/\1.multi.bam' )
def combineBams(infiles, outfile):
    '''collate all resultant bams together and index.

    This method assumes that BAM files have been sorted consistently by bowtie.
    '''

    to_cluster = USECLUSTER

    infile = infiles[0]
    if len(infiles) > 1:
       infiles = " ".join(infiles)
       statement = '''samtools merge -h %(infile)s %(outfile)s %(infiles)s >& %(outfile)s.log'''
       P.run()
    else:
       shutil.copyfile( infile, outfile )

    # assume that files are sorted
    # statement = '''samtools sort %(outfile)s''' 
    # P.run()

    statement = '''samtools index %(outfile)s >& %(outfile)s.log''' 
    P.run()

@transform(combineBams, suffix(".multi.bam"), ".unique.bam")
def uniquifyBams( infile, outfile ):
    '''extract unique hits'''

    to_cluster = USECLUSTER

    statement = '''python %(rmaadir)s/uniquify_bam.py %(infile)s %(outfile)s'''
    P.run()

    statement = '''samtools index %(outfile)s''' 
    P.run()

@transform(combineBams, suffix(".multi.bam"), ".gtf")
def buildGeneModels(infile, outfile):
    '''build transcript models - run cufflinks on each region seperately'''

    to_cluster = USECLUSTER    

    track = os.path.basename( outfile[:-len(".gtf")] )
    ins_size, std_dev = getInsertSizes( "reads/%s" % track )

    tmpfilename = P.getTempFilename()
    nslots = 4

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    infile = os.path.abspath( infile )
    outfile = os.path.abspath( outfile )

    statement = '''mkdir %(tmpfilename)s; 
    samtools view %(infile)s | sort -k3,3 -k4,4n 2> %(outfile)s.log1 > %(tmpfilename)s/temp.sam;
    cd %(tmpfilename)s; 
    cufflinks --inner-dist-mean %(ins_size)i
              --inner-dist-stddev %(std_dev)i
              --label %(track)s           
              --num-threads %(nslots)i 
              --min-isoform-fraction %(cuff_min_isoform)f
              --pre-mrna-fraction %(cuff_pre_mrna)f 
               %(tmpfilename)s/temp.sam >& %(outfile)s.log2;
    mv transcripts.gtf %(outfile)s >& %(outfile)s.log3;
    rm -rf %(tmpfilename)s >& %(outfile)s.log4 
    '''

    P.run()

@follows( mkdir("transcripts"))
@collate( buildGeneModels,
          regex(r"mappings/(.+).gtf$"), 
          'transcripts/summary.txt' )
def compareGeneModels(infiles, outfile):
    '''compare transcript models, using a reference GTF'''
    to_cluster = USECLUSTER

    infiles = " ".join(infiles)
    statement = '''
              cuffcompare 
                    -o %(outfile)s
                    -r %(ref_gtf)s 
                    -s %(sequence_dir)s 
                    %(infiles)s >& %(outfile)s.log 
    '''
    P.run()

@follows( mkdir("transcripts"))
@merge( buildGeneModels, 
        ["transcripts/all.shortreads.gtf", "transcripts/all.combined.gtf"] )
def combineGeneModels(infiles, outfiles):
    '''combine Cufflinks gene models together, and also combine with the given reference GTF.'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=10G"
    infiles = " ".join(infiles)
    outfile1, outfile2 = outfiles
    statement = '''cat %(infiles)s 
                 | awk '$3 == "exon"'
                 | python %(scriptsdir)s/gtf2gtf.py --merge-genes --log=%(outfile1)s.log
                 | python %(scriptsdir)s/gtf2gtf.py --renumber-genes="SR%%010i" --log=%(outfile1)s.log
                 > %(outfile1)s'''
    P.run()

    statement = '''cat %(infiles)s %(ref_gtf)s 
                 | awk '$3 == "exon"'
                 | python %(scriptsdir)s/gtf2gtf.py --merge-genes --log=%(outfile2)s.log
                 | python %(scriptsdir)s/gtf2gtf.py --renumber-genes="ALL%%010i" --log=%(outfile2)s.log
                 > %(outfile2)s'''
    P.run()

@transform(uniquifyBams, suffix(".bam"), ".counts")
def countReads(infile, outfile):
    '''count reads in Ensembl protein-coding models.'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=10G"
    
    statement = '''
    python %(rmaadir)s/count_reads_in_genes.py %(genes_file)s %(infile)s > %(outfile)s 
    '''
    P.run()

@follows(combineBams)
@transform(combineBams, suffix(".bam"), ".counts")
def countReadsMulti(infile, outfile):
    '''count MULTI reads in Ensembl protein-coding models'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=10G"
    
    statement = '''
    python %(rmaadir)s/count_reads_in_genes.py %(genes_file)s %(infile)s > %(outfile)s 
    '''
    P.run()

@merge( countReads, 
        ["mappings/unique.counts.all", "mappings/unique.ratios.all"] )
def adjustCounts(infiles, outfiles):
    '''normalize raw read counts to adjusted counts'''
    to_cluster = USECLUSTER
    job_options = "-l mem_free=20G"

    infiles = " ".join(infiles)
    outfile0, outfile1 = outfiles
    statement = '''
          python %(rmaadir)s/build_adjusted_counts.py 
              %(outfile1)s %(infiles)s 
          > %(outfile0)s'''
    P.run()

@merge( countReadsMulti, 
        ["mappings/multi.counts.all", "mappings/multi.ratios.all"] )
def adjustCountsMulti(infiles, outfiles):
    '''normalize raw read counts to adjusted counts for MULTI reads'''
    to_cluster = USECLUSTER
    job_options = "-l mem_free=20G"

    infiles = " ".join(infiles)
    outfile0, outfile1 = outfiles
    statement = '''
          python %(rmaadir)s/build_adjusted_counts.py 
              %(outfile1)s %(infiles)s 
          > %(outfile0)s'''
    P.run()

@follows( mkdir("graphs"))
@transform(adjustCounts, suffix(".counts.all"), ".rpkm.all")
def calculateRPKMs(infiles, outfile):
    '''calculate RPKMs from adjusted read counts'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=5G"

    infile = infiles[0]
    statement = '''
    python %(rmaadir)s/calculate_rpkms.py %(infile)s 
    > %(outfile)s
    '''
    P.run()

    statement = '''
    python %(rmaadir)s/sort_genes_by_rel_std_dev.py 
                     %(outfile)s %(min_rpkm_var)i 
    > %(outfile)s.rel_std_dev.%(min_rpkm_var)i'''
    
    P.run()

    statement = '''
    python %(rmaadir)s/sort_genes_by_expn.py %(outfile)s 
    > %(outfile)s.sorted_by_expression
    '''
    P.run()

def generate_calculate_term_params():
    for dbname, location in DATABASES.iteritems():
        yield [ "mappings/unique.rpkm.all", "mappings/unique.term_patterns.%s" % dbname, location ]

@follows(calculateRPKMs)
@files( generate_calculate_term_params )
def calculate_term_patterns(input, output, params):
    my_cmd = 'qrsh -now n -cwd -V -l mem_free=5G $RMAAPATH/sort_patterns_by_rel_std_dev.py %s %s %s %s > %s' % ( input, params, MIN_RPKM_TERM_VAR, MIN_GENES_TERM_VAR, output ); run_cmd(my_cmd)

# 
@follows( mkdir("graphs"))
@transform(adjustCountsMulti, suffix(".counts.all"), ".rpkm.all")
def calculateRPKMsMulti(infiles, outfile):
    '''calculate RPKMs from MULTI adjusted read counts'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=5G"
    
    infile = infiles[0]
    statement = '''
    python %(rmaadir)s/calculate_rpkms.py %(infile)s 
    > %(outfile)s
    '''
    P.run()

def _makeBigwig( infile, outfile, normfile ):

    with open(normfile, "r" ) as f:
        headers=f.readline().rstrip('\n').split('\t')
        ratios_list=f.readline().rstrip('\n').split('\t')

    ratios = dict( zip(headers, ratios_list) )
    ratio = float( ratios[ infile.rsplit('/',1)[1].split('.')[0] ] )

    to_cluster = USECLUSTER
    job_options = "-l mem_free=25G"
    
    outfile2 = outfile.rsplit('.',1)[0] + '.wig'
    statement = '''
    samtools pileup %(infile)s
    | awk '{ print $1 "\\t" $2 "\\t" $4 * %(ratio)f }' 
    | python %(rmaadir)s/pileup_to_wig.py > %(outfile2)s
    '''

    P.run()

    statement = '''
    wigToBigWig %(outfile2)s %(chrom_sizes)s %(outfile)s
    '''
    
    P.run()

@follows(adjustCounts, uniquifyBams)
@transform(uniquifyBams, suffix(".bam"), ".bw")
def makeBigwigs(infile, outfile):
    '''make normalized bigwigs.'''
    return _makeBigwig( infile, outfile, "mappings/unique.ratios.all" )

@follows(adjustCountsMulti, combineBams)
@transform(combineBams, suffix(".bam"), ".bw")
def makeBigwigsMulti(infile, outfile):
    '''make NORMALIZED bigwigs for MULTI counts'''
    return _makeBigwig( infile, outfile, "mappings/multi.ratios.all" )

@follows(calculateRPKMs, mkdir("trees") )
@files( [ ('mappings/unique.rpkm.all', 'trees/genes.%s.tree' % min_rpkm, min_rpkm) \
              for min_rpkm in PARAMS["min_rpkm_tree"] ])
def makeTreesAllGenes(infile, outfile, rpkmfile):
    '''build region relatedness trees for all genes'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=5G"

    statement = '''
    python %(rmaadir)s/make_trees_weight_genes.py %(infile)s %(rpkmfile)s %(outfile)s
    '''
    
    P.run()

############################################################################
############################################################################
############################################################################
## Clustering
############################################################################
@follows( calculateRPKMs, mkdir("fuzzy_k") )
@files( [ ( "mappings/unique.rpkm.all", 
            ("fuzzy_k/all-%s-%s.cdt" % (x,y),
             "fuzzy_k/background-%s-%s" % (x,y) ),
            x, y ) \
             for x,y in itertools.product( PARAMS["min_rpkm"], PARAMS["min_diff"] ) ] )

def buildCDTFiles( infile, outfiles, min_rpkm, min_diff ):
    '''build cdt files for fuzzy k clustering.'''
   
    cdt_filename, background_filename = outfiles
    # create cdt file
    labels = open(infile).readline().rstrip('\n').split('\t')[2::]

    min_diff_threshold = math.log( min_diff, 2)
    min_gene_len = PARAMS["min_gene_len"]
    background_file = open( background_filename, "w" )

    with open( cdt_filename, "w" ) as cdt_file:
       cdt_file.write( "UID\tNAME\tGWEIGHT\t%s\n" % ("\t".join( labels ) ) )
       cdt_file.write( "EWEIGT\t\t\t%s\n" % ( "\t".join( ["1"] * len(labels))))

       counts = E.Counter()

       for line in open(infile,"r"):

           if line.startswith("Name"): continue

           data = line[:-1].split("\t")
           counts.input += 1

           # exclude genes that are too short
           if int(data[1]) < min_gene_len:    
               counts.skipped_length += 1
               continue

           name = data[0]
           la = map(float, data[2:])

           # exclude lowly expressed genes
           if max(la) < min_rpkm: 
               counts.skipped_rpkm += 1
               continue

           background_file.write(name + "\n")

           # to handle any zero values, add 0.01 to every RPKM
           la = map(lambda x: x + 0.01, la)    
           avg_rpkm = float(sum(la) ) / len(la)
           ratios = [ math.log( x/avg_rpkm, 2) for x in la]

           if max(ratios) < min_diff_threshold:
               counts.skipped_diff += 1
               continue

           cdt_file.write( "%s\t%s\t%i\t%s\n" % (name, name, 1, 
                                                 "\t".join(map(str,ratios)) ) )

           counts.output += 1

    background_file.close()
   
    E.info( "%s\n" % str(counts) )

    # if we have too few genes to cluster anyway, mark it in a bad file so things downstream don't happen    
    l = len( open(cdt_filename).readlines())
    if l - 2 < PARAMS["min_genes_to_cluster"]:    
        bad_clusters = open('fuzzy_k/bad_clusters','a')
        bad_clusters.write( cdt_filename + '\n')
        bad_clusters.close()

@transform( buildCDTFiles, 
            regex(r"all-(.*).cdt"), 
            (r"instructions-\1", r"centroids-\1", r"membership-\1" ) )
def buildClusters( infiles, outfiles ):
    '''run c-means clustering on expression level data.'''

    to_cluster = USECLUSTER
    job_options = "-l mem_free=10G"

    # ignore the background file (why is it included in infiles?)
    infile, _ = infiles
    instructions_filename, centroid_filename, membership_filename = outfiles

    instructions_filename = os.path.abspath( instructions_filename )
    cdt_filename = os.path.abspath( infile )
    kmeans_clusters = PARAMS["kmeans_clusters"]

    # run aerie in a temporary directory
    tmpdir = P.getTempDir(".")

    with open( instructions_filename, "w" ) as outf:
       outf.write( '''load %(cdt_filename)s
fuzzy %(kmeans_clusters)i
%(tmpdir)s/all
exit
''' % locals())
       
    statement = '''
    aerie < %(instructions_filename)s >& %(instructions_filename)s.log
    '''
    P.run()

    shutil.move( os.path.join( tmpdir, "all.fct"), centroid_filename )
    shutil.move( os.path.join( tmpdir, "all.mb"), membership_filename )

    shutil.rmtree( tmpdir )

# extract & calculate enrichments in fuzzyK clusters over the background sets; 
# extract for several different cutoffs of "membership" 
# note some clusters may be degenerate (NEED ANOTHER SCRIPT TO PREPROCESS!)...
if not os.path.isfile('fuzzy_k/bad_clusters'):
    P.touch( 'fuzzy_k/bad_clusters')

def generate_fuzzy_clusters_params():
    bad_clusters=[]
    bad_clustersf = open('fuzzy_k/bad_clusters','r')
    for line in bad_clustersf:
        bad_clusters.append(line.rstrip('\n'))
    bad_clustersf.close()
    for min_rpkm in MIN_RPKM:
        for min_diff in MIN_DIFF:
            for cluster_lvl in CLUSTER_LVL:
                if not "fuzzy_k/all-%s-%s.pcl" % (min_rpkm, min_diff) in bad_clusters:
                    yield [ "fuzzy_k/membership-%s-%s" % (min_rpkm, min_diff), "fuzzy_k/cluster-%s-%s-%s.0" % (min_rpkm, min_diff, cluster_lvl), cluster_lvl ]

@follows( buildClusters )
@files( generate_fuzzy_clusters_params )
def extractClusters(infile, outfile, param0 ):
    '''build fuzzy clusters.'''

    to_cluster = USECLUSTER
    outp = outfile.rsplit('.',1)[0]

    statement = '''
    python %(rmaadir)s/extract_clusters.py %(infile)s %(outp)s %(param0)s
    '''
    
    P.run()

def generate_fuzzy_enrich_params():
    '''find enrichments.'''

    bad_clusters=[ x[:-1] for x in open('fuzzy_k/bad_clusters','r').readlines()]

    for min_rpkm, min_diff, cluster_lvl in \
           itertools.product( PARAMS["min_rpkm"], PARAMS["min_diff"], PARAMS["cluster_lvl" ] ):
       for dbname, location in DATABASES.iteritems():
          if "fuzzy_k/all-%s-%s.pcl" % (min_rpkm, min_diff) in bad_clusters: continue
          yield [ ["fuzzy_k/background-%s-%s" % (min_rpkm, min_diff), 
                   glob.glob("fuzzy_k/cluster-%s-%s-%s.*" % (min_rpkm, min_diff, cluster_lvl)), location], 
                  ["fuzzy_k/%s-summary-cluster-%s-%s-%s.0" % (dbname, min_rpkm, min_diff, cluster_lvl), 
                   "fuzzy_k/%s-expanded-cluster-%s-%s-%s.0" % (dbname, min_rpkm, min_diff, cluster_lvl)] ]

@follows( extractClusters )
@files( generate_fuzzy_enrich_params )
def computeEnrichments(infiles, outfiles):
    '''extract enrichments.'''
    to_cluster = USECLUSTER

    background_filename, foreground_filenames, association_filename = infiles
    
    for foreground_filename in foreground_filenames:
        cluster = foreground_filename.rsplit('.',1)[1]
        summary_filename = outfiles[0].rsplit('.',1)[0]+'.'+cluster
        expanded_filename = outfiles[1].rsplit('.',1)[0]+'.'+cluster
 
        statement = '''
        python %(rmaadir)s/compute_enrichments.py 
            %(foreground_filename)s 
            %(background_filename)s 
            %(association_filename)s 
            %(summary_filename)s 
            %(expanded_filename)s
        '''
        
        P.run()

def generate_enrichments_high_expn_params():
    for min_rpkm in PARAMS["min_rpkm_tree"]:
        for dbname, location in DATABASES.iteritems():
            yield [ ['mappings/unique.rpkm.all', location], 
                    "high_expn/high_expn.%s.%s" % (dbname, min_rpkm), min_rpkm ]

@follows( calculateRPKMs, mkdir("high_expn") )
@files(generate_enrichments_high_expn_params)
def computeEnrichmentsHighlyExpressedGenes( infiles, outfile, min_rpkm):

    to_cluster = USECLUSTER
    rpkm_filename, association_filename = infiles
    statement = '''cat %(rpkm_filename)s 
                   | sed 1d 
                   | awk '{ print $1 }' 
                   > %(outfile)s.background
    '''
    P.run()

    statement = '''cat %(rpkm_filename)s
                   | awk '{ sum=0; for(i=3; i<=NF; i++) { sum += $i}; sum/=NF; if (sum> %(min_rpkm)i) print $1 }' 
                   > %(outfile)s.foreground
    '''
    P.run()

    statement = '''
        python %(rmaadir)s/compute_enrichments.py 
            %(outfile)s.foreground
            %(outfile)s.background
            %(association_filename)s 
            %(outfile)s 
            %(outfile)s.expanded
    '''
    P.run()

# Calculate results of all enrichments.  For each database, find the best combination of parameters and report this file to a final analysis directory along with the relevant reference files.
def generate_params_report():
    for dbname, location in DATABASES.iteritems():
        yield [glob.glob("fuzzy_k/%s-summary-cluster-*-*-*.*" % dbname), 
               "best_conditions/clustering_summary_%s.txt" % dbname]

@follows( computeEnrichments, mkdir("best_conditions") )
@files( generate_params_report )
def reportBestParams(infiles, outfile):
    '''report the best parameters.'''

    bestfile = ''
    bestnumb = 0
    outf = open(outfile, 'w')
    numbs = {}
    fdr = PARAMS["fdr"]
    for infile in infiles:
        pref = infile.rsplit('.',1)[0]
        numbs[pref] = numbs.get( pref, 0 )
        numb = 0
        with open(infile,"r") as inf:
            for line in inf:
                if line.startswith("Name"): continue
                l = line[:-1].split("\t")
                # exit if fdr is above threshold
                if float(l[3]) >= fdr: break
                numb += 1
              
        numbs[pref] += numb

        outf.write( infile + '\t' + str(numb) + '\n' )

        if numbs[pref] > bestnumb:
            bestnumb = numbs[pref]
            bestfile = infile

    prefix = bestfile.rsplit('.',1)[0]
    for x in glob.glob( "%s.*"% prefix ):
        shutil.copyfile( x, "best_conditions/%s" % os.path.basename(x) )

    a = bestfile.rsplit('.',1)[0].rsplit('-',3)[1]
    b = bestfile.rsplit('.',1)[0].rsplit('-',3)[2] 

    for x in glob.glob( "fuzzy_k/centroids-%s-%s" % (a,b) ):
        shutil.copyfile( x, "best_conditions/%s" % os.path.basename(x) )
       
def generate_params_totalRNAfunctions():
    for dbname, location in DATABASES.iteritems():
        yield [['mappings/unique.rpkm.all', location], 
               ["overall_annotations/totalRNA.%s" % dbname, 
                "overall_annotations/totalRNA_diffs.%s" % dbname]]
@follows( calculateRPKMs, mkdir("overall_annotations") )
@files( generate_params_totalRNAfunctions )
def reportTotalRNAFunctions(infiles, outfiles):
    '''report total RNA functions.'''

    to_cluster = USECLUSTER
   
    rpkm_filename, annotations_filename = infiles
    expression_filename, diff_filename = outfiles
    statement = '''
    python %(rmaadir)s/report_totalRNA_annotations.py 
           %(rpkm_filename)s 
           %(annotations_filename)s 
           %(expression_filename)s 
           %(diff_filename)s
    '''
    
    P.run()

@follows( reportBestParams, 
          makeTreesAllGenes, 
          combineJunctions, 
          combineBams, 
          uniquifyBams, 
          combineGeneModels, 
          compareGeneModels, 
          calculateRPKMs, 
          calculateRPKMsMulti, 
          makeBigwigs, 
          makeBigwigsMulti, 
          reportTotalRNAFunctions, 
          computeEnrichmentsHighlyExpressedGenes,
          mkdir("web") )
def full(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
