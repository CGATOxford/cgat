'''
contigs2stats.py
====================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Compute basic statistics on a set of contigs in a multifasta file.

Usage
-----

Example::

   zcat infile.fasta.gz | python contigs2fasta.py -n 50 > out.stats

Type::

   python contigs2fasta.py --help

for command line help.

Command line options
--------------------

'''

import sys
import numpy as np
import CGAT.FastaIterator as FastaIterator

import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-n", dest="N", type="int",
                      help="e.g N50 - the length at which 50% of contigs are "
                      "equal or above")
    parser.add_option("-f", "--filter-length", dest="filter_length",
                      type="int",
                      help="calculate stats on contigs longer than -f")
    parser.add_option("--use-length", dest="use_length", action="store_true",
                      help="use a predefined length on which to calculate "
                      "the NX")
    parser.add_option("-l", "--length", dest="length", type="int",
                      help="if use-length is set then provide the "
                      "length to calculate N50")

    parser.set_defaults(N=50,
                        filter_length=0,
                        use_length=False,
                        length=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    f = options.filter_length

    # iterate over the contigs/scaffolds and return stats
    number_of_contigs = 0

    N = options.N
    contig_lengths = []

    for record in FastaIterator.iterate(options.stdin):
        contig_length = len(list(record.sequence))
        if contig_length >= f:
            number_of_contigs += 1
            contig_lengths.append(contig_length)

    # mean, median and max contig/scaffold lengths
    mean_length = np.mean(contig_lengths)
    median_length = np.median(contig_lengths)
    max_length = max(contig_lengths)

    # iterate over contigs/scaffolds sorted by longest
    # and caculate the NX
    index = 0
    cum_length = 0
    if options.use_length:
        assert options.length, "must supply a length when --use-length is set"
        total_length = options.length
    else:
        total_length = sum(contig_lengths)

    for length in sorted(contig_lengths, reverse=True):
        while cum_length <= total_length * (float(N) / 100):
            index += 1
            cum_length += length

    # output the results
    options.stdout.write(
        "nscaffolds\tscaffold_length\tN%i\tmedian_length\tmean_length\tmax_length\n" % N)
    options.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
        number_of_contigs,
        total_length,
        sorted(
            contig_lengths, reverse=True)[index],
        str(median_length),
        str(mean_length),
        str(max_length)))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
