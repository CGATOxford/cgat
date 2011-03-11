################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
'''
rnaseq_countreads.py - count context that reads map to
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes as input a :term:`BAM` file from an RNASeq experiment 
and a :term:`bed` formatted file.

It counts the number of alignments overlapping between the :term:`bam`
file and the :term:`bed` file. Annotations in the :term:`bed` file can
be overlapping - they are counted independently per name.

This scripts requires bedtools to be installed.

Usage
-----

Example::

   python script_template.py in.bam in.bed.gz

Type::

   python script_template.py --help

for command line help.

Documentation
-------------


Code
----

'''

import os, sys, re, optparse, time, subprocess, tempfile, collections

import Experiment as E
import IOTools

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--min-overlap", dest="min_overlap", type="float",
                       help = "minimum overlap [%default]" )

    parser.set_defaults(
        min_overlap = 0.5,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError( "please supply a bam and a bed file." )
    
    bamfile, bedfile = args
    
    E.info( "intersecting the two files" )

    tmpfile = tempfile.NamedTemporaryFile( delete=False )
    tmpfile.close()
    tmpfilename = tmpfile.name
    min_overlap = options.min_overlap
    
    statement = """intersectBed -abam %(bamfile)s -b %(bedfile)s -bed -wo -f %(min_overlap)f | groupBy -i stdin -g 4 -c 10 -o collapse > %(tmpfilename)s""" % locals()

    E.info( "running %s" % statement )
    E.run( statement )

    infile = open( tmpfilename, "r")
    counts_per_alignment = collections.defaultdict(int)

    E.info( "counting" )

    for line in infile:
        read, annotations = line[:-1].split("\t")
        annotations = annotations.split(",")[:-1]
        for anno in annotations:
            counts_per_alignment[anno] += 1
    infile.close()

    options.stdout.write( "category\talignments\n" )
    for key, counts  in counts_per_alignment.iteritems():
        options.stdout.write( "%s\t%i\n" % (key, counts) )

    os.unlink( tmpfilename )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

