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
gpipe/filter_fasta.py - select sequences in a fasta file
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts implements some more complex filtering. The sequence
identifiers need to correspond to the ```gpipe id`` scheme.

Usage
-----

Example::

   python gpipe/filter_fasta.py --help

Type::

   python gpipe/filter_fasta.py --help

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
import optparse
import math
import time
import tempfile
import subprocess


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator

# Class for calling masking programs.


class Masker:

    mLogLevel = 0
    mExecutable = None
    mOptions = ""

    def __init__(self):
        pass

    def __call__(self, peptide_sequence):
        """mask peptide sequence
        """
        Masker.__init__(self)

        outfile, filename_peptide = tempfile.mkstemp()
        os.write(outfile, ">test\n%s\n" % (peptide_sequence))
        os.close(outfile)

        statement = string.join(map(str, (
            self.mExecutable,
            filename_peptide,
            self.mOptions
        )), " ")

        if self.mLogLevel >= 3:
            print "# statement: %s" % statement
            sys.stdout.flush()

        s = subprocess.Popen(statement,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=True)

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise "Error in running %s \n%s\nTemporary directory" % (
                self.mExecutable, err)

        os.remove(filename_peptide)

        masked_sequence = re.sub(
            "\s", "", string.join(out.split("\n")[1:], ""))

        return masked_sequence


class MaskerBias (Masker):

    mLogLevel = 0
    mExecutable = "biasdb.pl"
    mOptions = ""


class MaskerSeg (Masker):

    mLogLevel = 0
    mExecutable = "seg"
    mOptions = "12 2.2 2.5 -x"

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/filter_fasta.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("longest-transcript", "ids", "quality"),
                      help="""method to apply to sequences ["longest-transcript", "ids", "quality"]."""  )

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameter stack for methods that require one.")

    parser.add_option("-t", "--sequence-type", dest="type", type="choice",
                      choices=("aa", "na"),
                      help="sequence type (aa or na).")

    parser.set_defaults(
        methods="",
        parameters="",
        type="na",
        aa_mask_chars="xX",
        aa_mask_char="x",
        na_mask_chars="nN",
        na_mask_char="n",
        gap_chars="-.",
        gap_char="-",
        template_identifier="ID%06i",
        separator="|",
    )

    (options, args) = E.Start(parser)
    options.parameters = options.parameters.split(",")

    iterator = FastaIterator.FastaIterator(sys.stdin)

    if options.method == "quality":
        filter_quality = set(options.parameters)
    else:
        filter_quality = None

    sequences = []
    ninput, noutput, nskipped = 0, 0, 0

    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        ninput += 1

        if filter_quality:
            id = re.split(" ", cur_record.title)[0]
            species, transcript, gene, quality = id.split(options.separator)

            if quality not in filter_quality:
                nskipped += 1
                continue

        sequences.append(cur_record)

    take = None

    if options.method == "longest-transcript":

        take = []
        lengths = []
        for x in range(len(sequences)):
            l = len(re.sub(" ", "", sequences[x].sequence))
            id = re.split(" ", sequences[x].title)[0]
            species, transcript, gene = id.split(options.separator)[:3]
            lengths.append((species, gene, -l, x))

        lengths.sort()

        last_species = None
        last_gene = None

        for species, gene, l, x in lengths:
            if last_species == species and last_gene == gene:
                continue
            take.append(x)
            last_species, last_gene = species, gene

    elif options.method == "ids":

        take = []
        ids, nerrors = IOTools.ReadList(open(options.parameters[0], "r"))
        del options.parameters[0]

        ids = set(ids)

        for x in range(len(sequences)):
            id = re.split(" ", sequences[x].title)[0]
            if id in ids:
                take.append(x)

    if take is not None:
        sequences = map(lambda x: sequences[x], take)

    noutput = len(sequences)

    for sequence in sequences:
        options.stdout.write(
            ">%s\n%s\n" % (sequence.title,  sequence.sequence))

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
