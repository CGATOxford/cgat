################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
rnaseq_bam_vs_bed.py - count context that reads map to
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

   python cgat_script_template.py in.bam in.bed.gz

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------


Code
----

'''

import os
import sys
import re
import optparse
import time
import subprocess
import tempfile
import collections
import itertools

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.Bed as Bed

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--min-overlap", dest="min_overlap", type="float",
                       help = "minimum overlap [%default]" )

    parser.add_option( "-k", "--keep-temp", dest="keep_temp", action="store_true",
                       help = "do not delete temporary files [%default]" )

    parser.set_defaults(
        min_overlap = 0.5,
        keep_temp = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError( "please supply a bam and a bed file or two bed-files." )
    
    bamfile, bedfile = args
    
    E.info( "intersecting the two files" )

    tmpfile = tempfile.NamedTemporaryFile( delete=False )
    tmpfile.close()
    tmpfilename = tmpfile.name

    min_overlap = options.min_overlap

    options.stdout.write( "category\talignments\n" )

    # get number of columns of reference bed file
    for bed in Bed.iterator(IOTools.openFile( bedfile )):
        ncolumns_bed = bed.columns
        break
    E.info( "assuming %s is bed%i format" % (bedfile, ncolumns_bed))

    if ncolumns_bed < 4:
        raise ValueError("please supply a name attribute in the bed file")

    # get information about 
    if bamfile.endswith(".bam"):
        format = "-abam"
        samfile = pysam.Samfile( bamfile, "rb" )
        total = samfile.mapped
        # latest bedtools uses bed12 format when bam is input
        ncolumns_bam = 12
        # count per read
        sort_key = lambda x: x.name
    else: 
        format = "-a"
        total = IOTools.getNumLines( bamfile )
        # get bed format
        ncolumns_bam = 0
        for bed in Bed.iterator(IOTools.openFile( bamfile )):
            ncolumns_bam = bed.columns
            break

        if ncolumns_bam > 0:
            E.info( "assuming %s is bed%i fomat" % (bamfile, ncolumns_bam ))
            if ncolumns_bam == 3:
                # count per interval
                sort_key = lambda x: (x.contig, x.start, x.end)
            else:
                # count per interval category
                sort_key = lambda x: x.name

    # use fields for bam/bed file (regions to count with)
    data_fields = [ "contig", "start", "end", "name",
                    "score", "strand", "thickstart", "thickend", "rgb",
                    "blockcount", "blockstarts", "blockends" ][:ncolumns_bam]

    # add fields for second bed (regions to count in)
    data_fields.extend( [ "contig2", "start2", "end2", "name2",
                          "score2", "strand2", "thickstart2", "thickend2", "rgb2",
                          "blockcount2", "blockstarts2", "blockends2" ][:ncolumns_bed] )

    # add bases overlap
    data_fields.append( "bases_overlap"  )

    data = collections.namedtuple( "data", data_fields )

    options.stdout.write( "total\t%i\n" % total )

    if total == 0:
        E.warn( "no data in %s" % bamfile )
        return

    statement = """intersectBed %(format)s %(bamfile)s -b %(bedfile)s -bed -wo -f %(min_overlap)f > %(tmpfilename)s""" % locals()

    E.info( "running %s" % statement )
    retcode = E.run( statement )

    if retcode != 0:
        raise ValueError( "error while executing statement %s" % statement )

    infile = open( tmpfilename, "r")
    counts_per_alignment = collections.defaultdict(int)

    E.info( "counting" )

    take_columns = len(data._fields)

    def iter( infile ):
        for line in infile:
            if not line.strip(): continue
            yield data._make( line[:-1].split()[:take_columns] )

    for read, overlaps in itertools.groupby( iter(infile), key = sort_key ):
        annotations = [x.name2 for x in overlaps ]
        for anno in annotations:
            counts_per_alignment[anno] += 1
    infile.close()

    for key, counts  in counts_per_alignment.iteritems():
        options.stdout.write( "%s\t%i\n" % (key, counts) )

    if not options.keep_temp:
        os.unlink( tmpfilename )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

