 ################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $
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
bed2graph.py - compute the overlap graph between two bed files
==============================================================

:Author: Andreas Heger
:Release: $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This script ouputs a list of the names of all overlapping intervals 
between two bed files.

Usage
-----

Type::

   python bed2graph.py A.bed.gz B.bed.gz > graph.out

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import numpy

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-o", "--output", dest="output", type="choice",
                      choices = ("full", "name"),
                      help="output either ``full`` overlapping entries, only the ``name``s. [default=%default]." )

    parser.set_defaults(
        output = "full",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError( "two arguments required" )

    if args[0] == "-":
        infile1 = options.stdin
    else:
        infile1 = IOTools.openFile( args[0], "r")

    infile2 = IOTools.openFile(args[1], "r")        
    
    idx = Bed.readAndIndex( infile2, with_values = True )

    output = options.output
    outfile = options.stdout
    
    if output == "name":
        outfile.write( "name1\tname2\n" )
        outf = lambda x: x.fields[0]
    else:
        outf = str
        
    for bed in Bed.iterator( infile1 ):
        try:
            overlaps = idx[bed.contig].find( bed.start, bed.end )
        except (KeyError, IndexError):
            # ignore missing contig and zero length intervals
            continue

        for o in overlaps:
            outfile.write( "\t".join( (outf(bed), outf(o[2]))) + "\n" )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )
