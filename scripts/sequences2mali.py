'''
sequences2mali.py - build pileup mali from a set of sequences
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a set of sequences into a multiple alignment

Usage
-----

Example::

   python sequences2mali.py --help

Type::

   python sequences2mali.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import optparse
import math
import time

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali
import alignlib_lite
import CGAT.FastaIterator as FastaIterator


def convertMali2Mali(mali):
    """convert a mali to a profile."""

    new_mali = alignlib_lite.py_makeMultipleAlignment()
    for id in mali.getIdentifiers():
        s = alignlib_lite.py_makeAlignatumFromString(mali[id])
        s.thisown = 0
        new_mali.addAlignatum(s)

    return new_mali

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: sequences2mali.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=(
                          "plain", "fasta", "clustal", "stockholm", "phylip"),
                      help="input format of multiple alignment")

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("plain", "fasta", "stockholm", "phylip"),
                      help="output format of multiple alignment")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("add",),
                      help="""method to use to build multiple alignment.""")

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameter stack for methods that require one.")

    parser.add_option("-a", "--alignment-method", dest="alignment_method", type="choice",
                      choices=("sw", "nw"),
                      help="alignment_method [%default].")

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        method=None,
        parameters="",
        gop=-10.0,
        gep=-1.0,
        alignment_method="sw",
    )

    (options, args) = E.Start(parser)

    options.parameters = options.parameters.split(",")

    iterator = FastaIterator.iterate(sys.stdin)

    if options.method == "add":

        mali = Mali.Mali()

        mali.readFromFile(
            open(options.parameters[0], "r"), format=options.input_format)
        del options.parameters[0]

        old_length = mali.getLength()

        new_mali = convertMali2Mali(mali)

        if options.alignment_method == "sw":
            alignator = alignlib_lite.py_makeAlignatorFullDP(
                options.gop, options.gep)
        else:
            alignator = alignlib_lite.py_makeAlignatorFullDPGlobal(
                options.gop, options.gep)

        while 1:
            cur_record = iterator.next()
            if cur_record is None:
                break

            map_mali2seq = alignlib_lite.py_makeAlignataVector()

            sequence = alignlib_lite.py_makeSequence(cur_record.sequence)
            profile = alignlib_lite.py_makeProfileFromMali(new_mali)

            if options.loglevel >= 4:
                options.stdlog.write(profile.Write())

            alignator.Align(profile, sequence, map_mali2seq)

            if options.loglevel >= 3:
                options.stdlog.write(map_mali2seq.Write())

            # add sequence to mali
            a = alignlib_lite.py_makeAlignatumFromString(cur_record.sequence)
            a.thisown = 0

            new_mali.addAlignatum(a, map_mali2seq, 1, 1, 1, 1, 1)

            id = cur_record.title
            mali.mIdentifiers.append(id)
            mali.mMali[id] = Mali.AlignedString(id, 0, len(
                cur_record.sequence), new_mali.getRow(new_mali.getWidth() - 1).getString())

        # substitute
        for x in range(old_length):
            mali.mMali[mali.mIdentifiers[x]].mString = new_mali.getRow(
                x).getString()

        mali.writeToFile(sys.stdout, format=options.output_format)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
