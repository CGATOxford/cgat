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
concatenate_sequences.py - concatenate sequences from multiple fasta files
==========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads sequences from two or more :term:`fasta` formatted
files and outputs a new file with the sequences concatenated per
entry. 

All files must have the same number of sequences and the id of
the first file is output.

Usage
-----

Example::

   python concatenate_sequences.py a.fasta b.fasta > c.fasta

Type::

   python concatenate_sequences.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import re

USAGE="""python %s [OPTIONS] in1 in2 [...]



""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: concatenate_sequences.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    (options, args) = E.Start( parser )

    if len(args) < 2:
        raise ValueError("please supply at least two filenames to concatenate.")

    iterators = []
    for a in args:
        iterators.append(FastaIterator.FastaIterator( IOTools.openFile(a,"r") ) )

    ninput, noutput, nerrors = 0, 0, 0
    
    while 1:

        sequences = []
        ids = []

        for iterator in iterators:
            try:
                cur_record = iterator.next()
            except StopIteration:
                break

            sequences.append( re.sub( " ", "", cur_record.sequence) )
            ids.append( cur_record.title )

        if not sequences: break
        ninput += 1

        if len(sequences) != len(iterators):
            raise "unequal number of sequences in files."

        noutput += 1

        options.stdout.write( ">%s\n%s\n" % (ids[0],
                                             "".join( sequences )))

    E.info( "ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput, nerrors) )
    
    E.Stop()
    
    
