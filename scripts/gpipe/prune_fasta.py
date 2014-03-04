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
gpipe/prune_fasta.py - prune fasta sequences
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Given input like::

    >BcDNA:AT01047 AY069026 [start_codon:116]
    CATACATAGTTCCCAGCAACTCAGGCCGGCAGCTACATAAGACCTTCGTA
    CCAAATCCAAAACGAAAACCCATTTTTCGCTCCGTTTCGATTCGAACCGT
    ATTCCAGTCCGCCCAGATGGATATGAACTACCAGTACAAGAAGGACCACT
    CGTTCGACAAGCGCCGCAACGAAGGCGACAAGATCCGGCGCAAGTATCCG
    GACCGTGTGCCCGTCATCGTGGAAAAGGCGCCGAAGACGCGTTACGCGGA
    GCTGGACAAGAAGAAGTACCTGGTGCCGGCGGACCTGACAGTGGGCCAGT
    TCTACTTTCTCATCCGCAAGCGTATCAATCTGCGTCCCGACGACGCCCTC
    TTCTTCTTCGTAAACAATGTGATCCCACCGACATCGGCCACCATGGGTGC
    ACTGTACCAGGAGCACTTCGACAAGGACTACTTCCTCTACATTTCCTATA
    CCGATGAGAACGTCTATGGACGGCAGTAGACGCGGGCTTGACTCGCGTAA
    TCGTTTTGGCGGTCGGTTGAGTTTCAATTGTCATTTTTTGCCATTGGCAC
    TACTCGCTCAATAAAGAATGACTGCTCCTAGCCAGCTAACCAAAAAAAAA
    AAAAAAAAA

this scripts outputs sequence from start codon to stop codon.

Usage
-----

Example::

   python gpipe/prune_fasta.py --help

Type::

   python gpipe/prune_fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import string
import os
import CGAT.Experiment as E
import optparse

param_stop_codons = ("TAG", "TAA", "TGA")
param_remove_errors = 1


def Write(description, sequence):

    try:
        start = int(
            re.search("\[start_codon:(\d+)\]", description).groups()[0])
    except AttributeError:
        if param_remove_errors:
            return
        description += " :Error, start=0"
        start = 0

    print ">" + description

    for stop in range(start, len(sequence), 3):
        if sequence[stop:stop + 3] in param_stop_codons:
            break

    fragment = sequence[start:stop]
    print fragment


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    sequence = ""
    description = None

    for line in options.stdin:

        if line[0] == ">":
            if description:
                Write(description, sequence)
            description = line[1:-1]
            sequence = ""
            continue

        sequence += line[:-1]

    Write(description, sequence)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
