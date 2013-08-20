################################################################################
#   MRC FGU Computational Genomics Group
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
Antonio changes up to 15/08/13:
- edits to documentation, modified header, added tags after Python, added under Usage, added under Example
- added code (commented)
 


fasta2counts.py - basic stats from collection of sequences
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python, fasta, counts, contigs, residues, gaps

Purpose
-------

This script computes some basic counts from a fasta-formatted file. For each sequence in
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

-g, --genome-file is a positional argument which requires an indexed fasta file created with index_fasta.py


Example::

#Download the example data here:

   wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Pseudomonas_aeruginosa_PAO1_uid57945/NC_002516.fna

#Index your fasta file:

   python /ifs/devel/antoniob/cgat/scripts/index_fasta.py NC_002516_db NC_002516.fna 

#which outputs:

NC_002516_db.fasta
NC_002516_db.idx

#Run fasta2counts.py:

   python /ifs/devel/antoniob/cgat/scripts/fasta2counts.py --stdin=NC_002516.fna --log=fasta2counts_NC_002516.log --error=fasta2counts_NC_002516.error --stdout=fasta2counts_NC_002516.output --genome-file=NC_002516_db.fasta 

#Results should look like:

   cat fasta2counts_NC_002516.output 
contig	nresidues	ngaps	nseqregions	ngapregions	nA	nC	nG	nT	nN	nX	nO
gi|110645304|ref|NC_002516.2|	6264404	0	1	0	1056134	2102687	2066633	1038950	0	0	0

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

import CGAT.Experiment as E
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

#added ' usage = globals()["__doc__"]) '
    parser = E.OptionParser( version = "%prog version: $Id: fasta2counts.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

#added ' ..., requires indexed fasta file" ' 
    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome, requires indexed fasta file"  )

    parser.set_defaults(
        genome_file = None,
        gap_chars = 'xXnN',
        filename_total = None,
        )

    (options, args) = E.Start( parser )

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    contigs = fasta.getContigSizes( with_synonyms = False )

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
        
    E.Stop()

