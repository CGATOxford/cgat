'''
psl2stats.py - output genomic coverage from psl formatted alignments
====================================================================

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

   python psl2stats.py --help

Type::

   python psl2stats.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E
import CGAT.Blat as Blat
import bx.bitset


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: psl2stats.py 2781 2009-09-10 11:33:14Z andreas $",
                            usage=globals()["__doc__"])

    parser.set_defaults(
    )

    (options, args) = E.Start(parser)

    query_bitsets, target_bitsets = {}, {}

    def addRange(bitset, id, size, iterator):

        if id not in bitset:
            bitset[id] = bx.bitset.BinnedBitSet(size)
        b = bitset[id]

        for start, end in iterator:
            b.set_range(start, end - start)

    for psl in Blat.iterator(options.stdin):

        addRange(query_bitsets,
                 psl.mQueryId,
                 psl.mQueryLength,
                 psl.iterator_query_exons())

        addRange(target_bitsets,
                 psl.mSbjctId,
                 psl.mSbjctLength,
                 psl.iterator_sbjct_exons())

    def printBitset(outfile, bitsets):

        outfile.write("contig\tcovered\tsize\tpcovered\n")
        total, total_len = 0, 0
        for chrom in sorted(bitsets):

            l = bitsets[chrom].size
            s = bitsets[chrom].count_range(0, l)
            if l > 0:
                outfile.write("%s\t%i\t%i\t%6.4f\n" %
                              (chrom, s, l, 100.0 * s / l))
            total += s
            total_len += l

        if total_len > 0:
            outfile.write("total\t%i\t%i\t%6.4f\n" %
                          (total, total_len, 100.0 * total / total_len))

    options.stdout.write("# query\n")
    printBitset(options.stdout, query_bitsets)
    options.stdout.write("# target\n")
    printBitset(options.stdout, target_bitsets)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
