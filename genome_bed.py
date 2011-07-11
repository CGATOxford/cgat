################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2011 David Sims
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
====================
Create a bed file tiling a genome from a fai file
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

This program takes one input file: a Samtools genome idex file.
It also takes two input parameters: the window/tile size (bases) and the shift size.
By default the shift size is equal to the window size.
The default window size is 1000.


Usage
=====

genome_bed -g <genome.fai> -o <output.bed> -w window size -s shift size


Code
====

"""

# load modules

import logging as L
import Experiment as E
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random, optparse
import GTF, IOTools, IndexedFasta

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: snp2table.py 2861 2010-02-23 17:36:32Z davids $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with Samtools indexed genome [default=%default]."  )
    parser.add_option("-w", "--window-size", dest="window", type="string",
                      help="Window size for tiling [default=%default]."  )
    parser.add_option("-s", "--shift-size", dest="shift", type="string",
                      help="Window shift for tiling [default=%default]."  )
    parser.add_option("-o", "--output", dest="output_file", type="string",
                      help="output filename  [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        window=1000,
        shift=1000,
        output_file = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ninput, nunchanged, nchanged = 0, 0, 0

    # Open input file
    E.info("Opening input file: %s" % options.genome_file)
    faidx = open(options.genome_file, "r" )

    # Open output file
    bed = open(options.output_file, "w")

    # Loop over input files and convert to soft clipped
    nwindows = 0
    ncontigs = 0
    for line in faidx:
        contig, stop = line.split()[0:2]
        stop = int(stop)
        ncontigs = ncontigs + 1
        i=1
        while (i < stop):
            if stop >= i+options.window:
                j=i+options.window-1
            else:
                j=stop
            #bed.write( """%(contig)s\t%(i)i\t%(j)i\t%(contig)s:%(i)i..%(j)i\n""" % locals() )
            bed.write( """%(contig)s\t%(i)i\t%(j)i\n""" % locals() )
            nwindows = nwindows + 1
            i = i+options.shift
        
    # Close all files
    faidx.close()
    bed.close()
                        
    # Report statistics
    E.info( "ncontigs=%i, nwindows=%i" % (ncontigs,nwindows) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

