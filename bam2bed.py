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

Convert BAM into BED. Example for samtools-devel.

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

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )


    parser.add_option("-r", "--region", dest="region", type="string",
                      help="samtools region string [default=%default]."  )

    parser.set_defaults(
        region = None,
        )

    (options, args) = parser.parse_args( argv[1:] )

    if len(args) != 1:
        raise ValueError( "no samfile specified - see --help for usage" )
    
    samfile = pysam.Samfile( args[0], "rb" )

    if options.region != None:
        it = samfile.fetch( options.region )
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

