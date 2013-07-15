################################################################################
#   Gene prediction pipeline 
#
#   $Id: fasta2counts.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
fasta2counts.py - basic stats from collection of sequences
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes same basic counts from a fasta-formatted file. For each sequence in
the :term:`fasta` formatted and indexed input file it will output in a :term:`tsv` formatted table:

contig
   Sequence identifier
nresidues
   Number of residues
ngaps
   Number of gaps
nseqregions
   Number of ungapped regions
ngapregions
   Number of gapped regions
nA, nC, nG, nT, nN, nX, nO
   Number of A,C,T,G,N,X,O characters in sequence

Usage
-----

Example::

   python fasta2counts.py h19

Type::

   python fasta2counts.py --help

for command line help.

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
import glob

import CGAT.Experiment as Experiment
import CGAT.IndexedFasta as IndexedFasta

def writeHeader( outfile ):
    outfile.write( "\t".join( ("contig",
                               "nresidues",
                               "ngaps",
                               "nseqregions",                                      
                               "ngapregions", 
                               "nA", "nC", "nG", "nT",
                               "nN", "nX", "nO" ) ) + "\n" )
    

    

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: fasta2counts.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome."  )

    parser.set_defaults(
        genome_file = None,
        gap_chars = 'xXnN',
        filename_total = None,
        )

    (options, args) = Experiment.Start( parser )

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    contigs = fasta.getContigSizes()

    totals = { 'A' : 0, 'C': 0, 'G' : 0, 'T' : 0, 'X' : 0, 'N' : 0 }
    nothers_total = 0
    nresidues_total = 0
    ngap_regions_total = 0
    nseq_regions_total = 0

    writeHeader( options.stdout )

    for contig in contigs.keys():
        
        subtotals = {}
        for x in totals.keys(): subtotals[x] = 0
        
        nothers = 0
        
        s = fasta.getSequence( contig, "+", 0, 0 )

        ngaps = 0
        
        was_gap = not (s[0] in options.gap_chars)
        
        ngap_regions = 0
        nseq_regions = 0

        x = 0
        xx = len(s)
        while x < xx:
            c = s[x].upper() 
            if c in subtotals:
                subtotals[c] += 1
            else:
                nothers += 1

            is_gap = c in options.gap_chars
            if is_gap:
                if not was_gap: ngap_regions += 1
            else:
                if was_gap: nseq_regions += 1
            was_gap = is_gap

            x += 1
            
        ngaps = subtotals['N'] + subtotals['X']
        
        options.stdout.write( "\t".join( map(str, (contig, len(s), ngaps,
                                                   nseq_regions,
                                                   ngap_regions,
                                                   subtotals['A'], subtotals['C'], subtotals['G'], subtotals['T'],
                                                   subtotals['N'], subtotals['X'], nothers ) ) ) + "\n" )
        

        for x in subtotals.keys(): totals[x] += subtotals[x]
        nothers_total += nothers
        nresidues_total += len(s)
        ngap_regions_total += ngap_regions
        nseq_regions_total += nseq_regions
        
    ngaps = totals['N'] + totals['X']

    if options.filename_total:
        outfile = open(options.filename_total, "w" )
        writeHeader( outfile )
        outfile.write( "\t".join( map(str, (len(contigs), nresidues_total, ngaps,
                                               nseq_regions_total,
                                               ngap_regions_total,
                                               totals['A'], totals['C'], totals['G'], totals['T'],
                                               totals['N'], totals['X'], nothers_total ) ) ) + "\n" )
        
        outfile.close()
        
    Experiment.Stop()
