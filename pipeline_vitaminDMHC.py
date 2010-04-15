################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_vitaminDMHC.py 2820 2009-11-24 16:06:58Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
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

:Author: Andreas Heger
:Release: $Id: pipeline_vitaminDMHC.py 2820 2009-11-24 16:06:58Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

vitamin D project pipeline

* Subproblem: mapping reads within the MHC region
  on chr6:28000000-34000000

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""

CONTIG,START,END = "chr6", 28000000,34000000

import sys, tempfile, optparse, shutil, itertools, random

from ruffus import *
import Experiment as E
import Pipeline as P
import IndexedFasta, IndexedGenome
import pysam

import pipeline_vitaminD
PARAMS = pipeline_vitaminD.PARAMS

############################################################
############################################################
############################################################
@files( (("genome.fasta", 'subsequence.fasta'),) )
def extractSequence( infile, outfile ):
    '''extract genomic sequence to be aligned against.'''
    
    fasta = IndexedFasta.IndexedFasta( infile[:-len(".fasta")] )
    outs = open( outfile,"w")
    outs.write( ">%s\n%s\n" % (CONTIG, fasta.getSequence( CONTIG, "+", START, END) ))
    outs.close()

############################################################
############################################################
############################################################
@files_re( extractSequence, '(.*).fasta', r'\1.ebwt')
def buildBowtieIndex( infile, outfile ):
    statement = '''bowtie-build %(infile)s %(outfile)s > %(outfile)s'''
    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@follows( extractSequence, buildBowtieIndex)
@files_re( '../*.bam', '../(.*).bam', (r"../\1.bam", "subsequence.ebwt"), r"\1.bam" )
def remapWithBowtie( infiles, outfile ):
    '''re-map unaligned reads.

    Select those reads that have not been mapped from a bam file (flag-value = 4)
    and map again with Bowtie.
    '''

    to_cluster = True

    tmpfilename = P.getTempFilename()

    prefix = outfile[:-len(".bam")]

    infile, subsequence = infiles
    start = START 
    statement = '''
    samtools view %(infile)s |\
    awk '$2 == 4 {printf("@%%s\\n%%s\\n+\\n%%s\\n", $1,$10,$11);}' |\
    bowtie --sam -n 3 %(subsequence)s - 2>%(outfile)s.log |\
    awk -v OFS="\\t" '/^@/ {print;next;} {if ($4 > 0) { $4 += %(start)s } print; }' |\
    samtools import %(genome)s - %(tmpfilename)s >& %(outfile)s.log;
    samtools sort %(tmpfilename)s %(prefix)s;
    samtools index %(outfile)s;
    rm -f %(tmpfilename)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

@files_re( remapWithBowtie,
           "(.*).bam",
           r"\1.bigwig" )
def buildBigwig( infile, outfile ):
    pipeline_vitaminD.buildBigwig( infile, outfile )

@files_re( remapWithBowtie,
           "(.*).bam",
           r"\1.readstats" )
def buildBAMStats( infile, outfile ):
    pipeline_vitaminD.buildBAMStats( infile, outfile )

def getMappedReads( infile ):
    '''return number of reads mapped.
    '''
    for lines in open(infile,"r"):
        data = lines[:-1].split("\t")
        if data[1].startswith( "mapped"):
            return int(data[0])
    return


def getMinimumMappedReads( infiles ):
    '''find the minimum number of mapped reads in infiles.'''
    v = []
    for infile in infiles:
        x = getMappedReads( infile )
        if x: v.append( x )
    return min(v)
    
@follows( buildBAMStats )
@files_re( remapWithBowtie,
           "(.*).bam",
           (r"\1.bam", r"\1.readstats" ),
           r"\1.normbam" )
def buildNormalizedBAM( infiles, outfile ):
    '''run MACS.'''
    
    min_reads = getMinimumMappedReads( glob.glob("*.readstats") )
    infile, statsfile = infiles
    num_reads = getMappedReads( statsfile )
    
    pysam_in = pysam.Samfile( infile, "rb" )
    pysam_out = pysam.Samfile( outfile, "wb", template = pysam_in )

    ninput, noutput = 0, 0

    take = [1] * min_reads + [0] * (num_reads-min_reads)
    random.shuffle( take )

    # iterate over mapped reads
    for read in pysam_in.fetch():
        if take[ninput]:
            pysam_out.write( read )
            noutput += 1
        ninput += 1

    pysam_in.close()
    pysam_out.close()

    P.info( "buildNormalizedBam: %i input, %i output (%5.2f%%), should be %i" % (ninput, noutput, 100.0*noutput/ninput, min_reads ))

@files_re( buildNormalizedBAM,
           "(.*).normbam",
           r"\1.macs" )
def runMACS( infile, outfile ):

    to_cluster = False

    track = infile[:-len("normbam")]
    try:
        control = pipeline_vitaminD.getControl( track ) + ".bam"
    except AssertionError:
        return

    statement = '''
    macs -t %(infile)s -c %(control)s \
          --name=%(outfile)s \
          --format=bam --tsize=35 --bw=110 --mfold=8 --gsize=6000000 >& %(outfile)s''' 

    P.run( **dict( locals().items() + PARAMS.items() ) )

@follows( remapWithBowtie, buildBigwig, runMACS )
def full():
    pass

if __name__== "__main__":
    P.checkFiles( ("genome.fasta", "genome.idx" ) )
    P.checkExecutables( ("liftOver",) )
    sys.exit( P.main(sys.argv) )

