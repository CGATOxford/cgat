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
gpipe/predictions2pseudogenes.py - 
======================================================

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

   python gpipe/predictions2pseudogenes.py --help

Type::

   python gpipe/predictions2pseudogenes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import getopt
import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser


USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/predictions2pseudogenes.py 18 2005-08-09 15:32:24Z andreas $

Compile collinear predictions into new predictions.

The entries have to be sorted by sbjct_token.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-f, --format=                   input format [exonerate|predictions]
-i, --max-intron                maximum intron length
-d, --max-difference            maximum difference between peptide and genomic gap
-c, --conserved-frame           require conserved frames
-a, --assembly                  do contig assembly
""" % sys.argv[0]

param_loglevel = 1

# maximum intron size
param_max_intron = 50000

param_format = "exonerate"

param_long_options = ["verbose=", "help",
                      "format=", "max-intron=",
                      "conserve-frame", "max-difference=",
                      "version"]


param_short_options = "v:hf:i:d:ca"

param_max_difference = 10

param_conserve_frame = 0

# ------------------------------------------------------------


def ProcessSegments(segments):
    """process a set of segments for a given query.

    1. Resolve exon permutations

    Exon permutations are streches, were the peptide fragment is
    not aligned in the right order to genomic DNA. This is not
    crucial here, as we are interested only in the genomic region.

    However, we do not want to extend genomic stretches due
    to spurious matches. Thus, delete exon permutations at
    the beginning and the end and only take the core.

    """

    if param_loglevel >= 3:
        print "## processing %i segments" % ntotal_segments

    # combine segments

    new_entries = []

    for x in range(len(segments) - 1):
        for y in range(x, len(segments)):

            # check for no overlap on genome
            if (min(segments[x].mSbjctGenomeTo, segments[y].mSbjctGenomeTo) -
                    max(segments[x].mSbjctGenomeFrom,
                        segments[y].mSbjctGenomeFrom)) > 0:
                continue

            # check for no overlap of sbjct
            if (min(segments[x].mQueryTo, segments[y].mQueryTo) -
                    max(segments[x].mQueryFrom, segments[y].mQueryFrom)) > 0:
                continue

            # check for collinearity
            d_aa = segments[y].mQueryFrom - segments[x].mQueryTo + 1
            d_na = segments[y].mSbjctGenomeFrom - segments[x].mSbjctGenomeTo

            if abs(d_aa * 3 - d_na) < param_max_difference:

                dframe = d_na % 3

                if param_loglevel >= 2:
                    print "# collinear sequences with d_aa=%i, d_na=%i, delta=%i, dframe=%i" % \
                          (d_aa, d_na, d_aa * 3 - d_na, dframe)

                if param_loglevel >= 3:
                    print "# part1:", str(segments[x])
                    print "# part2:", str(segments[y])

                if param_conserve_frame and dframe:
                    continue

                new_entry = segments[x].GetCopy()
                new_entry.Add(segments[y])
                new_entries.append(new_entry)

    return new_entries

# ------------------------------------------------------------


def ProcessChunk(entries):

    if param_loglevel >= 2:
        print "# received %i entries." % (len(entries))

    # array with predictions after segments have been merged
    new_entries = []

    if len(entries) > 0:

        # sort entries by query and genomic region
        entries.sort(lambda x, y: cmp((x.mQueryToken, x.mSbjctToken,
                                       x.mSbjctStrand, x.mSbjctGenomeFrom),
                                      (y.mQueryToken, y.mSbjctToken,
                                       y.mSbjctStrand, y.mSbjctGenomeFrom)))

        # array with distinct segmental regions
        segments = []

        last_entry = entries[0]
        segments.append(last_entry)

        for entry in entries[1:]:

            is_new_chunk = 0
            # check, if we are within the same "gene"
            # same gene is:
            # * same query, same chromosome, same strand
            # * gap not longer than param_max_intron
            if last_entry.mSbjctToken != entry.mSbjctToken or \
               last_entry.mSbjctStrand != entry.mSbjctStrand or \
               last_entry.mQueryToken != entry.mQueryToken or \
               (entry.mSbjctGenomeFrom - last_entry.mSbjctGenomeTo) > param_max_intron:

                new_entries += ProcessSegments(segments)
                segments = []

            segments.append(entry)
            last_entry = entry

        new_entries += ProcessSegments(segments)

    if param_loglevel >= 2:
        print "# number of predictions: %i" % len(new_entries)

    return new_entries

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-f", "--format"):
            param_format = a
        elif o in ("-i", "--max-intron"):
            param_max_intron = int(a)
        elif o in ("-d", "--max-difference"):
            param_max_difference = int(a)
        elif o in ("-f", "--conserve-frame"):
            param_conserve_frame = 1

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()

    parser = PredictionParser.PredictionParserExonerate()

    new_entries = []
    ninput = 0
    max_id = 0

    if param_format == "exonerate":
        for line in sys.stdin:

            if line[0] == "#":
                continue
            if line[:3] != "diy":
                continue

            data = string.split(line[:-1], "\t")
            query_token = data[1]

            ninput += 1

            # parser has to go inside, because GetBestMatch returns reference
            # copy
            result = parser.Parse([line, ])
            if not result:
                print "# ERROR: parsing line", line[:-1]
                continue
            entry = result.GetBestMatch()

            max_id = max(entry.mPredictionId, max_id)

            query_key = re.match("(\S+)", entry.mQueryToken).groups()[0]

            if param_correct_offset:
                data = string.split(entry.mSbjctToken, "_")
                if len(data) >= 3:
                    # truncate sbjct_token
                    entry.mSbjctToken = string.join(data[:-2], "_")
                    sbjct_offset_positive_from, sbjct_offset_negative_from = map(
                        int, data[-2:])

                    if entry.mSbjctStrand == "+":
                        sbjct_offset_from = sbjct_offset_positive_from
                    else:
                        sbjct_offset_from = sbjct_offset_negative_from
                else:
                    raise "parsing error for offset: key = %s" % sbjct_token

                entry.mSbjctGenomeFrom += sbjct_offset_from
                entry.mSbjctGenomeTo += sbjct_offset_from

            if param_loglevel >= 1:
                print "# received\t%s" % str(entry)

            entries.append(entry)

        new_entries += ProcessChunk(entries)

    elif param_format == "predictions":

        last_entry = None
        entries = []
        for line in sys.stdin:
            if line[0] == "#":
                continue

            entry = PredictionParser.PredictionParserEntry(expand=1)

            try:
                entry.Read(line)
            except ValueError:
                print "# warning: parsing error in line %s" % line[:-1]
                continue

            ninput += 1
            max_id = max(entry.mPredictionId, max_id)
            if last_entry:
                if last_entry.mSbjctToken != entry.mSbjctToken:
                    if entries:
                        new_entries += ProcessChunk(entries)
                    entries = []

            entries.append(entry)
            last_entry = entry

        new_entries += ProcessChunk(entries)

    max_id += 1
    first_pseudo_id = max_id

    for p in new_entries:
        p.mPredictionId = max_id
        print str(p)
        max_id += 1

    if param_loglevel >= 1:
        print "# nread=%i, new=%i, first_id=%i" % (ninput, len(new_entries), first_pseudo_id)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
