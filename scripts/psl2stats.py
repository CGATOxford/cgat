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
psl2stats.py - output genomic coverage from psl formatted alignments
====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python psl2stats.py --help

Type::

   python psl2stats.py --help

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
import tempfile
import subprocess
import optparse
import math


import CGAT.Experiment as E
import CGAT.Blat as Blat
import CGAT.IOTools as IOTools

# import psyco_full
import sys
import bx.bitset
import bx.bitset_builders

##---------------------------------------------------------------------------------------------
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: psl2stats.py 2781 2009-09-10 11:33:14Z andreas $",
                                    usage = globals()["__doc__"])

    parser.set_defaults(
        )
    
    (options, args) = E.Start( parser )

    query_bitsets, target_bitsets = {}, {}

    def addRange( bitset, id, size, iterator ):
        
        if id not in bitset: bitset[id] = bx.bitset.BinnedBitSet( size )
        b = bitset[id]

        for start, end in iterator:
            b.set_range( start, end-start )

    for psl in Blat.iterator( options.stdin ):

        addRange( query_bitsets, 
                  psl.mQueryId, 
                  psl.mQueryLength,
                  psl.iterator_query_exons() )

        addRange( target_bitsets, 
                  psl.mSbjctId, 
                  psl.mSbjctLength,
                  psl.iterator_sbjct_exons() )
        
    def printBitset( outfile, bitsets ):


        outfile.write( "contig\tcovered\tsize\tpcovered\n" )
        total, total_len = 0, 0
        for chrom in sorted(bitsets):
            
            l = bitsets[chrom].size 
            s = bitsets[chrom].count_range( 0, l )
            if l > 0:
                outfile.write( "%s\t%i\t%i\t%6.4f\n" % (chrom, s,l,100.0 * s / l) )
            total += s
            total_len += l

        if total_len > 0:
            outfile.write("total\t%i\t%i\t%6.4f\n" % (total,total_len, 100.0 * total / total_len))        
        
    options.stdout.write("# query\n" )
    printBitset( options.stdout, query_bitsets )
    options.stdout.write("# target\n" )
    printBitset( options.stdout, target_bitsets )

    E.Stop()
