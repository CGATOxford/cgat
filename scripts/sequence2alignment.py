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
sequence2alignment.py - convert aligned sequences to an alignment
=================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

map residue positions in sequences with gaps to positions in
the aligned string.

This script is useful to map aligned strings to their multiple
alignment positions.

Usage
-----

Example::

   python sequence2alignment.py --help

Type::

   python sequence2alignment.py --help

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
import types

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import alignlib
import CGAT.FastaIterator as FastaIterator

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: sequence2alignment.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.set_defaults(
        )

    (options, args) = E.Start( parser )

    iterator = FastaIterator.FastaIterator( sys.stdin )

    ninput, noutput, nskipped = 0, 0, 0
    
    options.stdout.write( "query\tsbjct\tquery_from\tquery_to\tsbjct_from\tsbjct_to\tquery_starts\tsbjct_starts\tblock_sizes\n" )

    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        ninput += 1
        
        sequence = re.sub( " ", "", cur_record.sequence)
        l = len(sequence)

        map_sequence2mali = alignlib.makeAlignmentVector()        

        alignlib.AlignmentFormatExplicit( 0, sequence,
                                          0, "X" * l ).copy( map_sequence2mali )

        options.stdout.write( "\t".join( (
                cur_record.title,
                "ref",
                str( alignlib.AlignmentFormatBlocks( map_sequence2mali ) ) ) ) + "\n" )
        
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i.\n" % (ninput, noutput, nskipped))
        
    E.Stop()
    
