################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Convert BAM into BED.

The fragment length of paired-end reads is the read lengths of both pairs
plus the insert size. 

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse
import pysam

import Experiment as E

import pyximport
pyximport.install(build_in_temp=False)
import _bam2bed

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )


    parser.add_option("-r", "--region", dest="region", type="string",
                      help="samtools region string [default=%default]."  )

    parser.add_option("-m", "--merge-pairs", dest="merge_pairs", action="store_true",
                      help="merge paired-ended reads into a single bed interval [default=%default]. " )

    parser.add_option( "--max-insert-size", dest="max_insert_size", type = "int",
                      help = "only merge if insert size less that # bases. 0 turns of this filter [default=%default]."  )

    parser.add_option( "--min-insert-size", dest="min_insert_size", type = "int",
                       help = "only merge paired-end reads if they are at least # bases apart. "
                              " 0 turns of this filter. [default=%default]" )


    parser.set_defaults(
        region = None,
        merge_pairs = None,
        call_peaks = None,
        min_insert_size = 0,
        max_insert_size = 0,
        )

    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 1:
        raise ValueError( "no samfile specified - see --help for usage" )
    
    samfile = pysam.Samfile( args[0], "rb" )

    if options.merge_pairs != None:
        counter = _bam2bed.merge_pairs( samfile, 
                                        options.stdout,
                                        min_insert_size = options.min_insert_size,
                                        max_insert_size = options.max_insert_size )

        options.stdlog.write( "category\tcounts\n%s\n" % counter.asTable() )

    else:

        if options.region != None:
            it = samfile.fetch( region = options.region )
        else:
            it = samfile.fetch()

        # more comfortable cigar parsing will
        # come with the next pysam elease
        BAM_CMATCH = 0
        BAM_CDEL = 2
        BAM_CREF_SKIP = 3
        take = (BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP)
        outfile = sys.stdout

        for read in it:
            if read.is_unmapped: continue

            t = 0
            for op, l in read.cigar:
                if op in take: t += l

            start = read.pos
            if read.is_reverse: strand = "-"
            else: strand = "+"
            outfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %\
                          ( read.rname,
                            read.pos,
                            read.pos+t,
                            read.qname,
                            read.mapq,
                            strand) )            

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

