'''fastqs2fastq.py - merge reads in fastq files
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS FASTQ FASTQ Manipulation

Purpose
-------

This script takes two paired-ended fastq files and outputs a single
fastq file in which reads have merged.

The two files must be sorted by read identifier.

Note that this script is currently a proof-of-principle implementation
and has not been optimized for speed or functionality.

Usage
-----

Example::

   python fastqs2fastq.py myReads.1.fastq.gz myReads.2.fastq.gz
          --method=join
          > join.fastq

In this example we take a pair of fastq files, join the reads and save
the output in :file:`join.fastq`.

Type::

   python fastqs2fastq.py --help

for command line help.

Command line options
--------------------

'''

import sys
import collections
import copy

import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import CGAT.Genomics as Genomics


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=('join', ),
                      help="method to apply [default=%default].")

    parser.set_defaults(
        method="join",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError(
            "please supply at least two fastq files on the commandline")

    fn1, fn2 = args
    c = E.Counter()
    outfile = options.stdout

    if options.method == "join":
        # merge based on diagonals in dotplot
        iter1 = Fastq.iterate(IOTools.openFile(fn1))
        iter2 = Fastq.iterate(IOTools.openFile(fn2))
        tuple_size = 2
        for left, right in zip(iter1, iter2):
            c.input += 1

            # build dictionary of tuples
            s1, q1 = left.seq, left.quals
            d = collections.defaultdict(list)
            for x in range(len(s1) - tuple_size):
                d[s1[x:x + tuple_size]].append(x)

            s2, q2 = right.seq, right.quals
            # reverse complement
            s2 = Genomics.complement(s2)
            q2 = q2[::-1]

            # compute list of offsets/diagonals
            offsets = collections.defaultdict(int)
            for x in range(len(s2) - tuple_size):
                c = s2[x:x + tuple_size]
                for y in d[c]:
                    offsets[x - y] += 1

            # find maximum diagonal
            sorted = sorted([(y, x) for x, y in offsets.items()])
            max_count, max_offset = sorted[-1]

            E.debug('%s: maximum offset at %i' % (left.identifier,
                                                  max_offset))

            # simple merge sequence
            take = len(s2) - max_offset
            merged_seq = s1 + s2[take:]

            # simple merge quality scores
            merged_quals = q1 + q2[take:]

            new_entry = copy.copy(left)
            new_entry.seq = merged_seq
            new_entry.quals = merged_quals
            outfile.write(new_entry)
            c.output += 1

    # write footer and output benchmark information.
    E.info("%s" % str(c))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
