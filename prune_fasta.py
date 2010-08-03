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
prune_fasta.py - prune fasta sequences
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Given input like::

    >BcDNA:AT01047 AY069026 [start_codon:116]
    CATACATAGTTCCCAGCAACTCAGGCCGGCAGCTACATAAGACCTTCGTA
    CCAAATCCAAAACGAAAACCCATTTTTCGCTCCGTTTCGATTCGAACCGT
    ATTCCAGTCCGCCCAGATGGATATGAACTACCAGTACAAGAAGGACCACT
    CGTTCGACAAGCGCCGCAACGAAGGCGACAAGATCCGGCGCAAGTATCCG
    GACCGTGTGCCCGTCATCGTGGAAAAGGCGCCGAAGACGCGTTACGCGGA
    GCTGGACAAGAAGAAGTACCTGGTGCCGGCGGACCTGACAGTGGGCCAGT
    TCTACTTTCTCATCCGCAAGCGTATCAATCTGCGTCCCGACGACGCCCTC
    TTCTTCTTCGTAAACAATGTGATCCCACCGACATCGGCCACCATGGGTGC
    ACTGTACCAGGAGCACTTCGACAAGGACTACTTCCTCTACATTTCCTATA
    CCGATGAGAACGTCTATGGACGGCAGTAGACGCGGGCTTGACTCGCGTAA
    TCGTTTTGGCGGTCGGTTGAGTTTCAATTGTCATTTTTTGCCATTGGCAC
    TACTCGCTCAATAAAGAATGACTGCTCCTAGCCAGCTAACCAAAAAAAAA
    AAAAAAAAA

this scripts outputs sequence from start codon to stop codon.

Usage
-----

Example::

   python prune_fasta.py --help

Type::

   python prune_fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, re, string, os

"""

Version: $Id: prune_fasta.py 18 2005-08-09 15:32:24Z andreas $

"""

param_stop_codons = ("TAG", "TAA", "TGA")
param_remove_errors = 1

def Write( description, sequence ):


    try:
        start = int(re.search( "\[start_codon:(\d+)\]", description).groups()[0])
    except AttributeError:
        if param_remove_errors: return
        description += " :Error, start=0"
        start = 0

    print ">" + description
    
    for stop in range(start, len(sequence), 3):
        if sequence[stop:stop+3] in param_stop_codons:
            break
        
    fragment = sequence[start:stop]
    print fragment
    
if __name__ == "__main__":

    sequence = ""
    description = None
    
    for line in sys.stdin:
        
        if line[0] == ">":
            if description:
                Write( description, sequence )
            description = line[1:-1]
            sequence = ""
            continue

        sequence += line[:-1]

    Write( description, sequence )
    

    
            
    
