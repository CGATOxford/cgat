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
mali2summary.py - compute column stats for a multiple alignment
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a multiple alignment from stdin and outputs
a :term:`tsv` formatted table on stdout. Each row contains
summary information on a multiple alignment column.

Usage
-----

Example::

   python mali2table.py --help

Type::

   python mali2table.py --help

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

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali
import scipy
import scipy.stats

from CGAT.SequenceProperties import *

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id$",
                                    usage = globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm" ),
                      help="input format of multiple alignment"  )
    
    parser.add_option("-a", "--alphabet", dest="alphabet", type="choice",
                      choices=("aa", "na"),
                      help="alphabet to use [default=%default].", )

    parser.add_option("-s", "--sections", dest="sections", type="choice", action="append",
                      choices = ("length", "composition", "entropy", "all"),
                      help="which sections to output" )

    parser.add_option( "-u", "--allow-duplicates", dest="allow_duplicates", action="store_true",
                      help="permit duplicate entries [default=%default]."  )

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        mask_chars = "nN",
        gap_chars = "-.",
        alphabet="na",
        sections = [],
        allow_duplicates = False,
        )

    (options, args) = E.Start( parser )

    if len(options.sections) == 0:
        raise ValueError("please supply at least one method." )

    if "all" in options.sections:
        options.sections = ["length", "composition", "entropy" ]

    counters = []

    def getCounter( section ):

        if options.alphabet == "na":
            if section == "length":
                s = SequencePropertiesLength()
            elif section == "composition":
                s = SequencePropertiesNA()
            elif section == "entropy":
                s = SequencePropertiesEntropy( "ACGT" )
            else:
                raise ValueError("unknown section %s" % section)
        elif options.alphabet == "aa":
            if section == "length":
                s = SequencePropertiesLength()
            elif section == "composition":
                s = SequencePropertiesAminoAcids()
            elif section == "entropy":
                s = SequencePropertiesEntropy( "ACDEFGHIKLMNPQRSTVWY" )
            else:
                raise ValueError("unknown section %s" % section)
        return s

    # read multiple alignment in various formats
    ## 1. read multiple alignment in various formats
    if options.allow_duplicates:
        mali = Mali.SequenceCollection()
    else:
        mali = Mali.Mali()

    mali.readFromFile( options.stdin, format = options.input_format )

    # do not use column, as it is a reserved word in sql
    options.stdout.write ("col" )
    for section in options.sections:
        options.stdout.write( "\t" + "\t".join(getCounter(section).getHeaders()))
    options.stdout.write( "\n" )

    columns = mali.getColumns()
    counter = E.Counter()

    for x, column in enumerate( columns ):
        counter.input += 1
        sequence = "".join(column)
        options.stdout.write( "%i" % x )

        for section in options.sections:
            s = getCounter( section )
            s.loadSequence( sequence )
            options.stdout.write( "\t" + "\t".join(s.getFields()) )

        options.stdout.write("\n")
        counter.output += 1

    E.info( "%s" % str(counter ) )
        
    E.Stop()
    
##------------------------------------------------------------
if __name__ == '__main__':
    sys.exit( main( sys.argv ) )

