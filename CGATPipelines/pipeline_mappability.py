################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_capseq.py 2900 2011-05-24 14:38:00Z david $
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
=====================
Mappability pipeline
=====================

Pipeline to count mappable bases in a given genome

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
import gzip
import sqlite3
import CGAT.Experiment as E
import logging as L
from ruffus import *
import CGAT.Pipeline as P

P.getParameters(  ["%s.ini" % os.path.splitext(__file__)[0],  "pipeline.ini" ] )
PARAMS = P.PARAMS
USECLUSTER = True

###################################################################
###################################################################
###################################################################
@files( os.path.join(PARAMS["gem_dir"],PARAMS["genome"]+".gem"),  PARAMS["genome"]+".mappability" )
def calculateMappability( infile, outfile ):
    '''Calculate mappability using GEM '''
    index = P.snip(infile, ".gem")
    to_cluster = True
    window_size = PARAMS("gem_window_size")
    threads = PARAMS("gem_threads")
    mismatches = PARAMS("gem_mismatches")
    max_indel_length = PARAMS("gem_max_indel_length")
    statement = '''gem-mappability -t %(threads)s -m %(mismatches)s --max-indel-length %(max_indel_length)s -l %(window_size)s -I %(index)s -o %(outfile)s ''' % locals()
    P.run()

###################################################################
@transform( calculateMappability, suffix(".mappability"), ".mappability.count" )
def countMappableBases( infile, outfile ):
    '''Count mappable bases in genome'''
    to_cluster = True
    statement = '''cat %(infile)s | tr -cd ! | wc -c > %(outfile)s''' % locals()
    P.run()
    
###################################################################
@transform( countMappableBases, suffix(".count"), ".count.load" )
def loadMappableBases( infile, outfile ):
    '''load count of mappable bases in genome'''
    to_cluster = True
    header = "total_mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=total_mappable_bases
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

###################################################################
@transform( calculateMappability, suffix(".mappability"), ".split.log" )
def splitMappabiliyFileByContig( infile, outfile ):
    '''Count mappable bases in genome'''
    to_cluster = True
    track = P.snip( os.path.basename(infile), ".mappability" )
    statement = '''mkdir contigs; 
                   csplit -k -f contigs/contig %(infile)s '/^~[a-zA-Z]/' {100000} > %(outfile)s;
                   rm contigs/contig00;''' % locals()
    P.run()

###################################################################
@follows( splitMappabiliyFileByContig )
@merge( "contigs/contig*", PARAMS["genome"]+"_mappability_per_contig.tsv" )
def countMappableBasesPerContig( infiles, outfile ):
    '''Count mappable bases for each contig'''
    for infile in infiles:
        statement = '''grep '~' %(infile)s | sed s/~//g >> %(outfile)s; cat %(infile)s | tr -cd ! | wc -c >> %(outfile)s'''
        P.run()
    
    statement = '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s;'''
    P.run()

###################################################################
@transform( countMappableBasesPerContig, suffix(".tsv"), ".tsv.load" )
def loadMappableBasesPerContig( infile, outfile ):
    '''load count of mappable bases per contig '''
    to_cluster = True
    header = "contig,mappable_bases"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=mappable_bases_per_contig
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

###################################################################
###################################################################
###################################################################
@follows( calculateMappability, countMappableBases,
          loadMappableBases, splitMappabiliyFileByContig,
          countMappableBasesPerContig, loadMappableBasesPerContig )
def full():
    '''Count mappable bases in genome'''
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
    
