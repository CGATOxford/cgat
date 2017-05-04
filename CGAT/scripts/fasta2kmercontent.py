'''
fasta2kmercontent.py
=============================================

:Tags: Genomics Sequences FASTA Summary

Purpose
-------

This script takes an input :term:`fasta` file from stdin and computes a
k-nucleotide content for each contig in the file. The output is a
tab-delimited file of kmer counts::

         contig1  contig2  contig3  contig4
    n1
    n2
    n3

where n is the kmer and contig is the fasta entry.

The user specifies the kmer length that is to be searched. Note that the longer
the kmer, the longer the script will take to run.

Note the order of output will not necessarily be the same order as the input.

Usage
-----
Example::

   zcat in.fasta.gz | head::

    >NODE_1_length_120_cov_4.233333
    TCACGAGCACCGCTATTATCAGCAACTTTTAAGCGACTTTCTTGTTGAATCATTTCAATT
    GTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCTACAAGACGGAAG
    CGTTTTGTAGCTGAAAGTGGGCGAGTTTCCATGATACGAACGATATCGCC

    >NODE_3_length_51_cov_33.000000
    CGAGTTTCCATGATACGAACGATATCGCCTTCTTTAGCAACGTTGTTTTCGTCATGTGCT
    TTATATTTTTTAGAATAGTTGATACGTTTACCATAGACTGG

   zcat in.fasta.gz | python fasta2kmercontent.py
                      --kmer-size 4
                      > tetranucleotide_counts.tsv

   head tetranucleotide_counts.tsv::

     kmer NODE_228_length_74_cov_506.432434 NODE_167_length_57_cov_138.438599
     GTAC 0                                 0
     TGCT 0                                 0
     GTAA 2                                 0
     CGAA 1                                 1
     AAAT 1                                 0
     CGAC 0                                 0

In this example, for each contig in in.fasta.gz the occurrence of each four
nucleotide combination is counted.

Alternative example::

   zcat in.fasta.gz | python fasta2kmercontent.py
                      --kmer-size 4
                      --output-proportion
                      > tetranucleotide_proportions.tsv

In this example, for each contig in in.fasta.gz we return the proportion of
each four base combination out of the total tetranucleotide occurences.
``--output-proportion`` overides the count output.

Options
-------
Two options control the behaviour of fasta2kmercontent.py; ``--kmer-size`` and
``--output-proportion``.

``--kmer-size``::
  The kmer length to count over in the input fasta file

``--output-proportion``::
  The output values are proportions rather than absolute counts


Type::

   python fasta2composition.py --help

for command line help.


Command line options
--------------------

'''

import sys
import re
import CGAT.FastaIterator as FastaIterator
import itertools
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-k", "--kmer-size", dest="kmer", type="int",
                      help="supply kmer length")

    parser.add_option(
        "-p", "--output-proportion", dest="proportion", action="store_true",
        help="output proportions - overides the default output")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # do not allow greater than octonucleotide
    assert options.kmer <= 8, "cannot handle kmer of length %i" % options.kmer

    # how we deal with the nucleotides depends on the kmer length
    nucleotides = []
    for nucleotide in ["A", "C", "T", "G"]:
        nucleotides = nucleotides + \
            [x for x in itertools.repeat(nucleotide, options.kmer)]

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
            counts = [m.start()
                      for m in re.finditer("".join(kmer), fasta.sequence)]
            result[fasta.title][kmer] = len(counts)

    E.info("writing results")
    # write out the results
    headers = sorted(result.keys())
    rows = set()
    for kmer_counts in list(result.values()):
        for kmer, count in kmer_counts.items():
            rows.add("".join(kmer))

    # write header row
    options.stdout.write("kmer\t" + "\t".join(headers) + "\n")

    # output proportions if required - normalises by
    # sequence length
    E.info("computing total counts")
    totals = {}
    for header in headers:
        totals[header] = sum([result[header][tuple(row)] for row in rows])

    for row in sorted(rows):
        if options.proportion:
            options.stdout.write("\t".join(
                [row] + [str(float(result[header][tuple(row)]) / totals[header]) for header in headers]) + "\n")
        else:
            options.stdout.write(
                "\t".join([row] + [str(result[header][tuple(row)]) for header in headers]) + "\n")

    E.info("written kmer counts for %i contigs" % total_entries)
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
