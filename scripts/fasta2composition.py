################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
fasta2composition.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes as input a multi-fasta file from stdin and computes a k-nucleotide
composition for each contig in the set. The output is a tab-delimited file of kmer comunts:

     contig1  contig2  contig3  contig4  
n1
n2
n3


where n is the kmer and contig is the fasta entry

The user can specify the nucleotides that are to be searched for example
tetra, pentamer etc in which case all tetramers from the ATCG will be computed
and tested.


Usage
-----

Example::

   python fasta2composition.py --help

Type::

   python fasta2composition.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import CGAT.FastaIterator as FastaIterator
import itertools
import collections
import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )
    parser.add_option("-k", "--kmer", dest="kmer", type="int",
                      help="supply kmer length")
    parser.add_option("-p", dest = "proportion", action="store_true",
                      help="output proportions")
    

    # add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # do not allow greater than octonucleotide
    assert options.kmer <= 8, "cannot handle kmer of length %i" % options.kmer
    
    # how we deal with the nucleotides depends on the kmer length
    nucleotides = []
    for nucleotide in ["A", "C", "T", "G"]:
        nucleotides = nucleotides + [x for x in itertools.repeat(nucleotide, options.kmer)]
 
    E.info("retrieving %imer sequences" % options.kmer)
    # get all kmer sequences to query 
    kmers = set()
    for kmer in itertools.permutations(nucleotides, options.kmer):
        kmers.add(kmer)
        
    E.info("matching %imers in file" % options.kmer)
    # count the number of kmers in each sequence
    
    result = {}

    # NB assume that non fasta files are caught by FastaIterator
    for fasta in FastaIterator.iterate(options.stdin):
        result[fasta.title] = {}
        for kmer in kmers:
            counts = [m.start() for m in re.finditer("".join(kmer), fasta.sequence)]
            result[fasta.title][kmer] = len(counts) 
            
    E.info("writing results")
    # write out the results
    headers = result.keys()
    rows = set()
    for kmer_counts in result.values():
        for kmer, count in kmer_counts.iteritems():
            rows.add("".join(kmer))
    
    # write header row
    options.stdout.write("kmer\t" + "\t".join(headers) + "\n")

    # output proportions if required - normalises by 
    # sequence length
    E.info("computing total counts")
    totals = {}
    for header in headers:
        totals[header] = sum([result[header][tuple(row)] for row in rows])

    for row in rows:
        if options.proportion:
            options.stdout.write("\t".join([row] + [str(float(result[header][tuple(row)])/totals[header]) for header in headers]) + "\n")
        else:
            options.stdout.write("\t".join([row] + [str(result[header][tuple(row)]) for header in headers]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
