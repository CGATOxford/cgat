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
optic/update_blast.py - take a blast output file and rerun with BLAST
===============================================================

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

   python optic/update_blast.py --help

Type::

   python optic/update_blast.py --help

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
import CGAT.BlastAlignments as BlastAlignments


def WriteUnique(file, clusters):

    for x, y in clusters.items():
        file.write(";".join(y) + "\n")


def GetMapSequences(sequences):

    map_cluster2sequence = {}
    map_sequence2cluster = {}
    # cluster sequences by identity
    # (clumsy sort, use hashes for bigger sets)
    for key, sequence in sequences.items():
        h = Genomics.GetHID(sequence)

        if h not in map_cluster2sequence:
            map_cluster2sequence[h] = []

        map_sequence2cluster[key] = h
        map_cluster2sequence[h].append(key)

    return map_cluster2sequence, map_sequence2cluster


def CompareSequenceSets(seqs1, seqs2):
    """compare two sequence sets.

    """

    map_cluster2sequence1, map_sequence2cluster1 = GetMapSequences(seqs1)
    map_cluster2sequence2, map_sequence2cluster2 = GetMapSequences(seqs2)

    unique1 = {}
    unique2 = {}
    common = {}

    for x, y in map_cluster2sequence1.items():
        if x in map_cluster2sequence2:
            common[x] = [y, ]
        else:
            unique1[x] = y
    for x, y in map_cluster2sequence2.items():
        if x in map_cluster2sequence1:
            common[x].append(y)
        else:
            unique2[x] = y

    return common, unique1, unique2


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/update_blast.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--new-query-sequences", dest="filename_new_query_sequences", type="string",
                      help="new sequences.")
    parser.add_option("--old-query-sequences", dest="filename_old_query_sequences", type="string",
                      help="old sequences.")
    parser.add_option("-o", "--output-filename-pattern", dest="filename_pattern", type="string",
                      help="pattern for output files.")
    parser.add_option("-f", "--write-fasta", dest="write_fasta", action="store_true",
                      help="write fasta files for runs.")
    parser.add_option("-a", "--fasta-pattern", dest="filename_pattern_fasta", type="string",
                      help="pattern for fasta files.")

    parser.set_defaults(
        filename_new_query_sequences=None,
        filename_old_query_sequences=None,
        filename_pattern="diff_%s",
        filename_pattern_fasta="%s.fasta",
        write_fasta=False,
    )

    (options, args) = E.Start(parser)

    if not options.filename_new_query_sequences:
        raise "please supply filename with new query sequences."
    if not options.filename_old_query_sequences:
        raise "please supply filename with old query sequences."

    new_query_sequences = Genomics.ReadPeptideSequences(
        open(options.filename_new_query_sequences, "r"))
    old_query_sequences = Genomics.ReadPeptideSequences(
        open(options.filename_old_query_sequences, "r"))

    if options.loglevel >= 1:
        print "# read %i/%i sequences for old/new query" % (len(old_query_sequences), len(new_query_sequences))

    common, unique_old, unique_new = CompareSequenceSets(
        old_query_sequences, new_query_sequences)

    filename = options.filename_pattern % "old2new.map"
    outfile = open(filename, "w")
    for x, y in common.items():
        outfile.write(";".join(y[0]) + "\t" + ";".join(y[1]) + "\n")
    outfile.close()

    WriteUnique(
        open(options.filename_pattern % "unique_query_old", "w"), unique_old)
    WriteUnique(
        open(options.filename_pattern % "unique_query_new", "w"), unique_new)

    if options.write_fasta:

        # write new query sequences
        outfile = open(options.filename_pattern_fasta % "query_new", "w")
        for x, yy in unique_new.items():
            for y in yy:
                outfile.write(">%s\n%s\n" % (y, new_query_sequences[y]))
        outfile.close()

    if options.loglevel >= 1:
        print "# query: common=%i, old=%i, new=%i" % (len(common), len(unique_old), len(unique_new))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
