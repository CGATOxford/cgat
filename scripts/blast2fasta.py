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

Command line options
--------------------

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
        sequences[k] = alignlib.py_makeSequence( sequences[k] )

    if options.loglevel >= 2:
        print "# converted %i sequences" % len(sequences)
    
    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0
    link = BlastAlignments.Link()

    ali = alignlib.py_makeAlignataVector()
    
    for line in sys.stdin:
        
        if line[0] == "#": continue

        link.Read( line )
        ninput += 1

        if link.mQueryToken not in sequences or link.mSbjctToken not in sequences:
            nskipped += 1
            continue
        
        ali.Clear()
        alignlib.py_fillAlignataCompressed( ali, link.mQueryFrom, link.mQueryAli, link.mSbjctFrom, link.mSbjctAli )


        result = alignlib.py_writePairAlignment( sequences[link.mQueryToken], sequences[link.mSbjctToken], ali ).split("\n")

        if len(result) != 3:
            nfailed += 1

        if options.format == "fasta":
            print ">%s %i-%i\n%s\n>%s %i-%i\n%s\n" %\
                  (link.mQueryToken, link.mQueryFrom, link.mQueryTo, result[0].split("\t")[1],
                   link.mSbjctToken, link.mSbjctFrom, link.mSbjctTo, result[1].split("\t")[1] )
            
        noutput += 1
        
    print "# ninput=%i, noutput=%i, nskipped=%i, nfailed=%i" % (ninput, noutput, nskipped, nfailed)
    E.Stop()

            
    
