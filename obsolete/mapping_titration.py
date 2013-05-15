################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_cpg.py 2900 2011-05-24 14:38:00Z david $
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
==========================================
Read Mapping parameter titration pipeline
==========================================

:Author: David Sims 
:Release: $Id: mapping_titration.py 2900 2011-05-24 14:38:00Z david $
:Date: |today|
:Tags: Python

   * align reads to the genome using a range of different parameters
   * calculate alignment statistics

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie_             |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+


Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Glossary
========

.. glossary::

   bowtie
      bowtie_ - a read mapper

.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml

Code
====

"""
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import CGAT.Experiment as E
import logging as L
from ruffus import *
import csv
import sqlite3
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import cStringIO
import pysam
import numpy
import gzip
import fileinput
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMapping as PipelineMapping
from bein.util import *

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters(  ["%s.ini" % __file__[:-len(".py")],  "../pipeline.ini", "pipeline.ini" ] )
PARAMS = P.PARAMS

bowtie_options = {'n0m1':"-n 0 -a --best --strata -m 1 -3 1",'n1m1':"-n 1 -a --best --strata -m 1 -3 1",'n2m1':"-n 2 -a --best --strata -m 1 -3 1",'n3m1':"-n 3 -a --best --strata -m 1 -3 1",
                  'n0m2':"-n 0 -a --best --strata -m 2 -3 1",'n1m2':"-n 1 -a --best --strata -m 2 -3 1",'n2m2':"-n 2 -a --best --strata -m 2 -3 1",'n3m2':"-n 3 -a --best --strata -m 2 -3 1",
                  'n0m3':"-n 0 -a --best --strata -m 3 -3 1",'n1m3':"-n 1 -a --best --strata -m 3 -3 1",'n2m3':"-n 2 -a --best --strata -m 3 -3 1",'n3m3':"-n 3 -a --best --strata -m 3 -3 1",
                  'n0m4':"-n 0 -a --best --strata -m 4 -3 1",'n1m4':"-n 1 -a --best --strata -m 4 -3 1",'n2m4':"-n 2 -a --best --strata -m 4 -3 1",'n3m4':"-n 3 -a --best --strata -m 4 -3 1",
                  'n0m5':"-n 0 -a --best --strata -m 5 -3 1",'n1m5':"-n 1 -a --best --strata -m 5 -3 1",'n2m5':"-n 2 -a --best --strata -m 5 -3 1",'n3m5':"-n 3 -a --best --strata -m 5 -3 1",
                  'v0m1':"-v 0 -a --best --strata -m 1 -3 1",'v1m1':"-v 1 -a --best --strata -m 1 -3 1",'v2m1':"-v 2 -a --best --strata -m 1 -3 1",'v3m1':"-v 3 -a --best --strata -m 1 -3 1",
                  'v0m2':"-v 0 -a --best --strata -m 2 -3 1",'v1m2':"-v 1 -a --best --strata -m 2 -3 1",'v2m2':"-v 2 -a --best --strata -m 2 -3 1",'v3m2':"-v 3 -a --best --strata -m 2 -3 1",
                  'v0m3':"-v 0 -a --best --strata -m 3 -3 1",'v1m3':"-v 1 -a --best --strata -m 3 -3 1",'v2m3':"-v 2 -a --best --strata -m 3 -3 1",'v3m3':"-v 3 -a --best --strata -m 3 -3 1",
                  'v0m4':"-v 0 -a --best --strata -m 4 -3 1",'v1m4':"-v 1 -a --best --strata -m 4 -3 1",'v2m4':"-v 2 -a --best --strata -m 4 -3 1",'v3m4':"-v 3 -a --best --strata -m 4 -3 1",
                  'v0m5':"-v 0 -a --best --strata -m 5 -3 1",'v1m5':"-v 1 -a --best --strata -m 5 -3 1",'v2m5':"-v 2 -a --best --strata -m 5 -3 1",'v3m5':"-v 3 -a --best --strata -m 5 -3 1"}

###################################################################
###################################################################
###################################################################
## MAP READS
@files( [ (PARAMS["test_file"], "%s.bam" % x, bowtie_options.get(x))  for x in bowtie_options.keys() ] )
def buildBAM( infile, outfile, options ):
    '''map reads with bowtie'''
    to_cluster = True
    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    m = PipelineMapping.Bowtie()
    reffile = PARAMS["samtools_genome"]
    bowtie_options = options
    statement = m.build( (infile,), outfile ) 
    #print(statement)
    P.run()

#########################################################################
@transform( buildBAM, 
            regex(r"(\S+).bam"),
            r"\1.nsrt.bam" )
def sortByName( infile, outfile ):
    '''Add number of hits tags to sam file'''
    to_cluster = USECLUSTER
    track = P.snip(outfile, ".bam")
    statement = '''samtools sort -n %(infile)s %(track)s;'''
    P.run()

#########################################################################
@transform( sortByName, 
            regex(r"(\S+).nsrt.bam"),
            r"\1.nh.bam" )
def addNHTag( infile, outfile ):
    '''Add number of hits tags to sam file'''
    to_cluster = USECLUSTER

    inf = pysam.Samfile(infile, "rb")
    outf = pysam.Samfile(outfile, "wb", template=inf)
    for readset in read_sets(inf,keep_unmapped=True):
        nh = len(readset)
        for read in readset:
            if (read.is_unmapped):
                nh = 0
            read.tags = read.tags+[("NH",nh)]
            outf.write(read)
    inf.close()
    outf.close()

#########################################################################
@transform( addNHTag, 
            regex(r"(\S+).bam"),
            r"\1.srt.bam" )
def sortByPosition( infile, outfile ):
    '''Add number of hits tags to sam file'''
    to_cluster = USECLUSTER
    track = P.snip(outfile, ".bam")
    statement = '''samtools sort %(infile)s %(track)s;'''
    P.run()

#########################################################################
@transform( sortByPosition,
            regex( r"(\S+).nh.srt.bam"),
            r"\1.dedup.bam")
def dedup(infiles, outfile):
        '''Remove duplicate alignments from BAM files.'''
        to_cluster = USECLUSTER
        track = P.snip( outfile, ".bam" )
        statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s METRICS_FILE=%(track)s.dupstats VALIDATION_STRINGENCY=SILENT; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
        #print statement
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
                > %(outfile)s
               '''
    P.run()

#########################################################################
@transform( dedup, 
            regex(r"(\S+).dedup.bam"),
            r"\1.readstats" )
def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. '''
    to_cluster = USECLUSTER
    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''python %(scriptsdir)s/bam2stats.py --force 
                   --output-filename-pattern=%(outfile)s.%%s < %(infile)s > %(outfile)s'''
    P.run()

#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''Import bam statistics into SQLite'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".readstats") for x in infiles] )
    filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s"""
    P.run()

    for suffix in ("nm", "nh"):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()



#########################################################################
@transform( dedup, 
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.alignstats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
@merge( buildPicardAlignStats, "picard_align_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )

    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.alignstats" )
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
                > %(outfile)s
               '''
    P.run()

    os.unlink( tmpfilename )


############################################################
############################################################
############################################################
## Pipeline organisation
@follows( buildBAM, sortByName, addNHTag, sortByPosition, dedup, 
          loadPicardDuplicateStats, buildBAMStats, loadBAMStats)
def mapReads():
    '''Align reads to target genome.'''
    pass

@follows( mapReads )
def full():
    '''run the full pipeline.'''
    pass

############################################################
############################################################
############################################################
## REPORTS
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

