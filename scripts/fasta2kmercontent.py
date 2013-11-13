'''
fasta2kmercontent.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Genomics Sequences FASTA Summary

Purpose
-------

This script takes as input a :term:`fasta` file from stdin and computes a k-nucleotide
content for each contig in the file. The output is a tab-delimited file of kmer counts:

     contig1  contig2  contig3  contig4  
n1
n2
n3


where n is the kmer and contig is the fasta entry.

The user specifies the kmer that is to be searched. Note that the longer the kmer, the
longer the script will take to run.


Usage
-----

Example::

   zcat in.fasta.gz | python fasta2kmercontent.py --kmer 4 > tetranucleotide_counts.tsv

In this example, for each contig in in.fasta.gz we count the occurrence of each four base
combination.


Alternative example::

   zcat in.fasta.gz | python fasta2kmercontent.py --kmer 4 --proportion > tetranucleotide_proportions.tsv

In this example, for each contig in in.fasta.gz we return the proportion of each four base
combination out of the total tetranucleotide occurences. --proportion overides the count
output. 


Type::

   python fasta2composition.py --help

for command line help.


Command line options
--------------------

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
    parser = E.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                             usage = globals()["__doc__"] )

    parser.add_option("-k", "--kmer", dest="kmer", type="int",
                      help="supply kmer length")

    parser.add_option("-p", "--proportion", dest = "proportion", action="store_true",
                      help="output proportions - overides the default output")
    

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
    total_entries = 0
    for fasta in FastaIterator.iterate(options.stdin):
        total_entries += 1
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

    E.info("written kmer counts for %i contigs" % total_entries)
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
