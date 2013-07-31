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
blast2fasta.py - 
======================================================

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

   python blast2fasta.py --help

Type::

   python blast2fasta.py --help

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
import getopt
import tempfile
import time
import optparse
import math

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.BlastAlignments as BlastAlignments
import alignlib

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: blast2fasta.py 2782 2009-09-10 11:40:29Z andreas $

Convert a blast graph into a pairwise alignment graph

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --peptides=                 filename with peptide sequences
""" % sys.argv[0]

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: blast2fasta.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-s", "--sequences", dest="filename_sequences", type="string",
                      help="filename with sequences."  )
    parser.add_option("-f", "--format", dest="format", type="string",
                      help="output format."  )
    
    parser.set_defaults(
        filename_sequences = None,
        format = "fasta",
        )

    (options, args) = E.Start( parser )

    if not options.filename_sequences:
        raise "please supply filename with sequences."

    sequences = Genomics.ReadPeptideSequences( open(options.filename_sequences, "r") )

    if options.loglevel >= 1:
        print "# read %i sequences" % len(sequences)
        
    for k in sequences.keys():
        sequences[k] = alignlib.makeSequence( sequences[k] )

    if options.loglevel >= 2:
        print "# converted %i sequences" % len(sequences)
    
    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0
    link = BlastAlignments.Link()

    ali = alignlib.makeAlignataVector()
    
    for line in sys.stdin:
        
        if line[0] == "#": continue

        link.Read( line )
        ninput += 1

        if link.mQueryToken not in sequences or link.mSbjctToken not in sequences:
            nskipped += 1
            continue
        
        ali.Clear()
        alignlib.fillAlignataCompressed( ali, link.mQueryFrom, link.mQueryAli, link.mSbjctFrom, link.mSbjctAli )


        result = alignlib.writePairAlignment( sequences[link.mQueryToken], sequences[link.mSbjctToken], ali ).split("\n")

        if len(result) != 3:
            nfailed += 1

        if options.format == "fasta":
            print ">%s %i-%i\n%s\n>%s %i-%i\n%s\n" %\
                  (link.mQueryToken, link.mQueryFrom, link.mQueryTo, result[0].split("\t")[1],
                   link.mSbjctToken, link.mSbjctFrom, link.mSbjctTo, result[1].split("\t")[1] )
            
        noutput += 1
        
    print "# ninput=%i, noutput=%i, nskipped=%i, nfailed=%i" % (ninput, noutput, nskipped, nfailed)
    E.Stop()

            
    
