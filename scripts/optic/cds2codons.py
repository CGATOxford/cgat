##########################################################################
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
##########################################################################
'''
optic/cds2codons.py - remove frameshifts from cdna sequences
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Given a map between peptide to cds sequences, print out
sequences with only codons.

Usage
-----

Example::

   python optic/cds2codons.py --help

Type::

   python optic/cds2codons.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import time
import sets
import optparse
import math
import tempfile

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import alignlib_lite


class Map:

    def __init__(self):
        pass

    def Read(self, line):
        (self.mToken,
         self.mOldFrom, self.mOldTo, self.mOldAli,
         self.mNewFrom, self.mNewTo, self.mNewAli,
         self.mOldLength, self.mNewLength) = line[:-1].split("\t")

        (self.mOldFrom, self.mOldTo, self.mNewFrom, self.mNewTo) = \
            map(int, (self.mOldFrom, self.mOldTo, self.mNewFrom, self.mNewTo))
        self.mMapOld2New = None

    def Expand(self):
        self.mMapOld2New = alignlib_lite.makeAlignmentVector()
        alignlib_lite.AlignmentFormatEmissions(
            self.mOldFrom, self.mOldAli,
            self.mNewFrom, self.mNewAli).copy(self.mMapOld2New)

    def Clear(self):
        if self.mMapOld2New:
            self.mMapOld2New.clear()
        self.mMapOld2New = None

    def __str__(self):
        return string.join(map(str, (self.mToken,
                                     self.mOldFrom, self.mOldTo, self.mOldAli,
                                     self.mNewFrom, self.mNewTo, self.mNewAli,
                                     self.mOldLength, self.mNewLength)), "\t")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/cds2codons.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--map", dest="filename_map", type="string",
                      help="filename with mapping information.")
    parser.add_option("-f", "--format", dest="format", type="string",
                      help="output file format [fasta-codons].")
    parser.add_option("-c", "--codons", dest="codons", action="store_true",
                      help="print codons separated by spaces.")

    parser.set_defaults(
        filename_cds=None,
        codons=False,
        format="fasta",
        filename_map=None,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if not options.filename_map:
        raise "please supply filename with map between peptide to cds."

    if options.filename_map:
        map_old2new = {}
        for line in open(options.filename_map, "r"):
            if line[0] == "#":
                continue
            m = Map()
            m.Read(line)
            map_old2new[m.mToken] = m
    else:
        map_old2new = {}

    if options.filename_cds:
        sequences = Genomics.ReadPeptideSequences(
            open(options.filename_cds, "r"))
    else:
        sequences = Genomics.ReadPeptideSequences(sys.stdin)

    if options.loglevel >= 1:
        print "# read %i sequences" % len(sequences)
        sys.stdout.flush()

    ninput, nskipped, noutput, nerrors, nstops = 0, 0, 0, 0, 0

    for key, s in sequences.items():

        ninput += 1

        if key not in map_old2new:
            nskipped += 1
            continue

        out_seq = []

        m = map_old2new[key]
        m.Expand()
        mm = m.mMapOld2New

        if mm.getColTo() > len(s):
            options.stderr.write("# error for %s: sequence shorter than alignment: %i < %i\n" % (
                key, len(s), mm.getColTo()))
            nerrors += 1
            continue

        for x in range(mm.getRowFrom(), mm.getRowTo() + 1):

            y = mm.mapRowToCol(x)
            if y > 0:
                out_seq.append(s[y - 1])

        m.Clear()

        out_seq = "".join(out_seq)
        translation = Genomics.TranslateDNA2Protein(out_seq)

        if "X" in translation:
            nstops += 1

        if options.codons:
            out_seq = " ".join([out_seq[x:x + 3]
                               for x in range(0, len(out_seq), 3)])

        noutput += 1
        options.stdout.write(">%s\n%s\n" % (key, out_seq))

    options.stderr.write("# input=%i, output=%i, errors=%i, stops=%i\n" % (
        ninput, noutput, nerrors, nstops))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
