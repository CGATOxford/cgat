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
==================================================================
Combine hard masked + unmasked genomes to give soft masked genome
==================================================================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

This program takes two input files: a hard masked genome and a complementary unmasked genome ans outputs a soft masked genome.


Usage
=====

softmask -u <unmasked.fasta> -h <hard_masked.fasta> -o <outfile.fasta>


Code
====

"""

# load modules

import logging as L
import CGAT.Experiment as E
import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random
import optparse
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-u", "--unmasked", dest="unmasked_genome_file", type="string",
                      help="filename with unmasked genome [default=%default]."  )
    parser.add_option("-r", "--hard-masked", dest="hard_masked_genome_file", type="string",
                      help="filename with hard masked genome [default=%default]."  )
    parser.add_option("-o", "--output", dest="output_file", type="string",
                      help="output filename  [default=%default]."  )

    parser.set_defaults(
        unmasked_genome_file = None,
        hard_masked_genome_file = None,
        output_file = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ninput, nunchanged, nchanged = 0, 0, 0

    # Open input files
    E.info("Opening input file: %s" % options.unmasked_genome_file)
    unmask = open(options.unmasked_genome_file, "r" )
    E.info("Opening input file: %s" % options.hard_masked_genome_file)
    hardmask = open(options.hard_masked_genome_file, "r" )

    # Open output file
    softmask = open(options.output_file, "w")

    # Loop over input files and convert to soft clipped
    for line in hardmask:
        unmask_line = unmask.readline()
        ninput = ninput + 1

        if line.startswith( ">" ) or line==unmask_line:
            softmask.write( line )
            nunchanged = nunchanged + 1
        else:
            softmask_line = ""
            for i, c in enumerate(line):
                if c=="N":
                    softmask_line = softmask_line + unmask_line[i].lower()
                else:
                    softmask_line = softmask_line + c
            softmask.write( softmask_line )
            nchanged = nchanged + 1
        

    # Close all files
    unmask.close()
    hardmask.close()
    softmask.close()
                        
    # Report statistics
    E.info( "ninput=%i, nunchanged=%i, nchanged=%i" % (ninput,nunchanged,nchanged) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

