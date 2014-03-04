'''
align_all_vs_all.py - all-vs-all pairwise alignment
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes all-vs-all alignments between
sequences in a :term:`fasta` formatted file.

Currently only Smith-Waterman protein alignment is
implemented.

Usage
-----

Example::

   python align_all_vs_all.py --help

Type::

   python align_all_vs_all.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import getopt
import time
import optparse
import math
import tempfile

import CGAT.Experiment as E

import alignlib_lite
import CGAT.FastaIterator as FastaIterator


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: align_all_vs_all.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--sequences", dest="filename_sequences", type="string",
                      help="input file with sequences")

    parser.set_defaults(
        filename_sequences=None,
        gop=-10.0,
        gep=-1.0,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if options.filename_sequences:
        infile = open(options.filename_sequences, "r")
    else:
        infile = sys.stdin

    parser = FastaIterator.FastaIterator(infile)

    sequences = []
    while 1:
        cur_record = iterator.next()

        if cur_record is None:
            break
        sequences.append((cur_record.title, alignlib_lite.py_makeSequence(
            re.sub(" ", "", cur_record.sequence))))

    if options.filename_sequences:
        infile.close()

    alignator = alignlib_lite.py_makeAlignatorFullDP(options.gop, options.gep)
    map_a2b = alignlib_lite.py_makeAlignataVector()
    nsequences = len(sequences)

    for x in range(0, nsequences - 1):
        for y in range(x + 1, nsequences):
            alignator.Align(sequences[x][1], sequences[y][1], map_a2b)

            row_ali, col_ali = alignlib_lite.py_writeAlignataCompressed(
                map_a2b)

            options.stdout.write("%s\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%s\t%i\t%i\t%i\t%i\n" % (
                sequences[x][0], sequences[y][0],
                map_a2b.getScore(),
                map_a2b.getRowFrom(),
                map_a2b.getRowTo(),
                row_ali,
                map_a2b.getColFrom(),
                map_a2b.getColTo(),
                col_ali,
                map_a2b.getScore(),
                100 *
                alignlib_lite.py_calculatePercentIdentity(
                    map_a2b, sequences[x][1], sequences[y][1]),
                sequences[x][1].getLength(),
                sequences[y][1].getLength()))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
