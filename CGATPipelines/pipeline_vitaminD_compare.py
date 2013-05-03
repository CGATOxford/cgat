################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_vitaminD_compare.py 2861 2010-02-23 17:36:32Z andreas $
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
:Release: $Id: pipeline_vitaminD_compare.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

vitamin D project pipeline

Input

   experimental tracks should be called run<cellline><condition><replicate>
   control tracks should be called control<cellline>

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

TODO: currently the bed files and the intervals are inconsistent 
    (due to filtering, there are more intervals in the bed files than
     in the table. The ids do correspond).

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *
import csv
import sqlite3
import IOTools
import pysam
import numpy
import gzip

PARAMS=P.getParameters()

VERSIONS = ("version1", "version2", "version3", "version4", "version5")

@transform( "version1_*.bed", suffix(".bed"), ".compare")
def compare( infile, outfile ):
    '''compare several bed-files.'''
    pattern = re.match("version\d+_(.*).bed", infile).groups()[0]

    files = " ".join( sorted(glob.glob( "version*_%s.bed" % pattern )) )

    statement = '''
    python %(scriptsdir)s/diff_bed.py %(files)s > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )


@transform( "version1_*.bed", suffix(".bed"), ".rest")
def difference_to1( infile, outfile ):
    '''compare several bed-files. List the number of intervals that are not present in the other versions
    compared to version 1.'''
    track = re.match("version\d+_(.*).bed", infile).groups()[0]

    tmpfile = P.getTempFilename()
    
    for version in VERSIONS:
        t = tmpfile + "%s" % version
        if version == "version1":
            statement = '''cut -f 5 < %(version)s_%(track)s.bed |\
                    python %(toolsdir)s/data2histogram.py --headers=%(version)s --bin-size=1 --min-value=1 > %(t)s
                    '''
        else:
            statement = '''
            intersectBed -v -a version1_%(track)s.bed -b %(version)s_%(track)s.bed | cut -f 5 |\
            python %(toolsdir)s/data2histogram.py --headers=%(version)s --bin-size=1 --min-value=1 > %(t)s
            '''
        P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
    python %(toolsdir)s/combine_tables.py --sort-keys=numeric %(tmpfile)s* > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

@transform( "version2_*.bed", suffix(".bed"), ".rest2")
def difference_to2( infile, outfile ):
    '''compare several bed-files.'''
    track = re.match("version\d+_(.*).bed", infile).groups()[0]

    tmpfile = P.getTempFilename()
    
    for version in VERSIONS:
        t = tmpfile + "%s" % version
        if version == "version2":
            statement = '''cut -f 5 < %(version)s_%(track)s.bed |\
                    python %(toolsdir)s/data2histogram.py --headers=%(version)s --bin-size=1 --min-value=1 > %(t)s
                    '''
        else:
            statement = '''
            intersectBed -v -a version2_%(track)s.bed -b %(version)s_%(track)s.bed | cut -f 5 |\
            python %(toolsdir)s/data2histogram.py --headers=%(version)s --bin-size=1 --min-value=1 > %(t)s
            '''
        P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
    python %(toolsdir)s/combine_tables.py --sort-keys=numeric %(tmpfile)s* > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
