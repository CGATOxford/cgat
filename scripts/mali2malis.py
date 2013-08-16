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
'''
mali2malis.py - split a multiple alignment into smaller multiple alignments
===========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

split a multiple alignment into its constituent parts.

This is the inverse operation performed by :doc:`malis2mali`.

Usage
-----

Example::

   python mali2malis.py --help

Type::

   python mali2malis.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import math
import time
import random

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: mali2malis.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-a", "--pattern-mali", dest="pattern_mali", type="string",
                      help="filename pattern for multiple alignment files."  )

    parser.add_option( "--filename-coordinates", dest="filename_coordinates", type="string",
                      help="filename of coordinates that constitute the multiple alignment."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal" ),
                      help="input format of multiple alignment"  )
    
    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=( "fasta", "codeml", "phylip" ),
                      help="output format of multiple alignment"  )

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        filename_coordinates = None,
        pattern_mali= "%s.fasta",
        )

    (options, args) = E.Start( parser )

    ## read coordinates
    if options.filename_coordinates:
        coordinates = []
        for line in open(options.filename_coordinates, "r" ):
            if line[0] == "#": continue
            id, length, position = line[:-1].split("\t")
            if id == "component": continue
            coordinates.append( (id, int(length), int(position) ) )
            
    mali = Mali.Mali()
    mali.readFromFile( sys.stdin, format=options.input_format )

    ids = mali.getIdentifiers()

    ninput, noutput = 0, 0
    
    for id, length, position in coordinates:
        
        ninput += 1
        part_mali = Mali.Mali()

        for x in ids:
            part_mali.addSequence( x, 0, length, mali[x][position:position+length] )

        outfile_name = options.pattern_mali % id

        outfile = open(outfile_name, "w" )

        part_mali.writeToFile( outfile, format = options.output_format )

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# input=%i, output=%i\n" % (ninput, noutput) )

    E.Stop()
    
