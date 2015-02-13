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
gpipe/exonerate2regions.py - 
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

   python gpipe/exonerate2regions.py --help

Type::

   python gpipe/exonerate2regions.py --help

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
import CGAT.Genomics as Genomics
import CGAT.Intervalls as Intervalls
import CGAT.PredictionParser as PredictionParser

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/exonerate2regions.py 1799 2008-03-28 11:44:19Z andreas $

Analyse exonerate matches and calculate summary statistics for an alignment.

Resolve conflicts of overlapping segments.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-g, --genome=                   file with the genomic DNA (FASTA).
-p, --peptides-fasta-file=                 file with peptide sequences (FASTA).
-o, --correct-offset                use offset information in sbjct_token to correct genomic location.
-f, --format=                   input format [exonerate|predictions]
-c, --force-contiguous          force contiguous blocks
-d, --dump                      just dump results
""" % sys.argv[0]

HEADER = """# QUERY:        1  query
# QFROM:        2  query first residue
# QTO:          3  query last residue
# QSTR:         4  query strand
# QLEN:         5  query length
# SBJCT:        6  sbjct
# SFROM:        7  sbjct first residue
# STO:          8  sbjct last residue
# SSTR:         9  sbjct strand
# SLEN:         10 sbjct length
# SCORE:        11 alignment score
# CQUERY:       12 coverage of query (in percent)
# NGAPS:        13 number of gaps in alignment 
# NFR:          14 number of frame-shifts
# NINTRON:      15 number of introns
# NSPLIT:       16 number of split codons
# PIDE:         17 percent identity
# PSIM:         18 percent similarity
# NSEGS:        19 number of segments
# NPER:         20 number of removed permuted exons
# EXONS:        21 list of exons"""

SHORT_HEADER_SUMMARY = """#s\tQUERY\tQFROM\tQTO\tQSTR\tQLEN\tSBJCT\tSFROM\tSTO\tSSTR\tSLEN\tSCORE\tCQUERY\tNGAPS\tNFR\tNINTRON\tNSPLIT\tPIDE\tPSIM\tNSEGS\tNPER\tEXONS"""
SHORT_HEADER_ENTRY = """#e\tQUERY\tQFROM\tQTO\tQSTR\tQLEN\tSBJCT\tSFROM\tSTO\tSSTR\tSLEN\tSCORE\tCQUERY\tNGAPS\tNFR\tNINTRON\tNSPLIT\tPIDE\tPSIM"""

param_stop_codons = ("TAG", "TAA", "TGA")

param_loglevel = 1

# maximum intron size
param_max_intron = 50000

# check for stop codons within exons
param_border_stop_codon = 3

# correct genomic offset
param_correct_offset = None

param_format = "exonerate"

param_long_options = ["verbose=", "help", "genome=", "peptides=",
                      "min-score=", "correct-offset", "format=",
                      "max-intron=", "force-contiguous", "dump",
                      "version"]

param_short_options = "v:hg:p:e:o:c:s:of:i:cd"

param_filename_peptides = None
param_filename_genome = None
param_force_contiguous = 0
param_just_dump = 0


def RemoveExonPermutationsFromFront(segments):
    """remove exon permutations from the front.

    Only permutations are removed, that are completely out of sync
    with query. Overlapping segments are retained, as they might
    correspond to conflicting starting points and both should be
    checked by genewise.
    """

    if len(segments) <= 1:
        return segments

    first_index = 0

    while first_index + 1 < len(segments):
        if segments[first_index].mQueryFrom < segments[first_index + 1].mQueryTo:
            break

        first_index += 1

    return segments[first_index:]


def RemoveExonPermutationsFromBack(segments):
    """remove exon permutations from the back

    Only permutations are removed, that are completely out of sync
    with query. Overlapping segments are retained, as they might
    correspond to conflicting starting points and both should be
    checked by genewise.
    """

    if len(segments) <= 1:
        return segments

    first_index = len(segments) - 1

    while first_index - 1 > 0:
        if segments[first_index].mQueryTo > segments[first_index].mQueryFrom:
            break

        first_index -= 1

    return segments[:first_index + 1]

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

    # resolve exon permutations from the front
    ntotal_segments = len(segments)

    if len(segments) > 1:
        segments = RemoveExonPermutationsFromFront(segments)
        segments = RemoveExonPermutationsFromBack(segments)

    # combine segments
    s = segments[0]

    parts = [(s.mQueryFrom, s.mQueryTo)]

    for x in segments[1:]:
        parts.append((x.mQueryFrom, x.mQueryTo))
        s.mQueryFrom = min(x.mQueryFrom, s.mQueryFrom)
        s.mQueryTo = max(x.mQueryTo, s.mQueryTo)
        s.mSbjctGenomeFrom = min(x.mSbjctGenomeFrom, s.mSbjctGenomeFrom)
        s.mSbjctGenomeTo = max(x.mSbjctGenomeTo, s.mSbjctGenomeTo)
        s.mNIntrons += x.mNIntrons + 1
        s.mPercentIdentity += x.mPercentIdentity
        s.mPercentSimilarity += x.mPercentSimilarity

    s.mSbjctFrom = 0
    s.mSbjctTo = 0
    s.mAlignmentString = ""
    s.mPercentIdentity /= len(segments)
    s.mPercentSimilarity /= len(segments)

    new_intervalls = Intervalls.CombineIntervallsLarge(parts)
    # calculate covered region
    covered = reduce(
        lambda x, y: x + y, map(lambda x: x[1] - x[0] + 1, new_intervalls))

    s.mQueryCoverage = 100 * covered / s.mQueryLength

    return [s]

# ------------------------------------------------------------


def ProcessChunk(entries):

    if param_loglevel >= 1:
        print "# received %i entries." % (len(entries))

    # array with predictions after segments have been merged
    predictions = []

    if len(entries) > 0:

        # sort entries by genomic region
        entries.sort(lambda x, y: cmp((x.mQueryToken, x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
                                      (y.mQueryToken, y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

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

                is_new_chunk = 1

            else:
                if last_entry.mQueryTo > entry.mQueryFrom:
                    if param_force_contiguous:
                        is_new_chunk = 1
                    else:
                        if param_loglevel >= 1:
                            print "# WARNING: exon permutation in alignment between %s and %s" %\
                                  (entry.mQueryToken, entry.mSbjctToken)

            if is_new_chunk:
                if param_loglevel >= 3:
                    print SHORT_HEADER_ENTRY

                predictions += ProcessSegments(segments)
                segments = []

            segments.append(entry)
            last_entry = entry

        predictions += ProcessSegments(segments)

    if param_loglevel >= 1:
        print "# number of predictions: %i" % len(predictions)

    for prediction in predictions:
        prediction.Write()

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
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-g", "--genome"):
            param_filename_genome = a
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-o", "--correct-offset"):
            param_correct_offset = 1
        elif o in ("-f", "--format"):
            param_format = a
        elif o in ("-i", "--max-intron"):
            param_max_intron = int(a)
        elif o in ("-c", "--force-contiguous"):
            param_force_contiguous = 1
        elif o in ("-d", "--dump"):
            param_just_dump = 1

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()

    # read complete genomic sequence
    if param_filename_genome:
        forward_sequences, reverse_sequences = Genomics.ReadGenomicSequences(
            open(param_filename_genome, "r"))
    else:
        forward_sequences = None
        reverse_sequences = None

    # read peptide sequences
    if param_filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(param_filename_peptides, "r"))
    else:
        peptide_sequences = {}

    # print HEADER

    if param_loglevel >= 2:
        print SHORT_HEADER_SUMMARY

    # aligned entries from exonerate
    entries = []

    parser = PredictionParser.PredictionParserExonerate()

    if param_format == "exonerate":
        for line in sys.stdin:

            if line[0] == "#":
                continue
            if line[:3] != "diy":
                continue

            data = string.split(line[:-1], "\t")

            query_token = data[1]

            # parser has to go inside, because GetBestMatch returns reference
            # copy
            result = parser.Parse([line, ])
            if not result:
                print "# ERROR: parsing line", line[:-1]
                continue
            entry = result.GetBestMatch()

            query_key = re.match("(\S+)", entry.mQueryToken).groups()[0]

            if forward_sequences and forward_sequences.has_key(entry.mSbjctToken):
                if entry.mSbjctStrand == "+":
                    genomic_sequence = forward_sequences[entry.mSbjctToken]
                else:
                    genomic_sequence = reverse_sequences[entry.mSbjctToken]
            else:
                genomic_sequence = None

            if genomic_sequence:
                entry.SetTranslation(genomic_sequence)

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
                    raise "parsing error for offset: key = %s" % entry.mSbjctToken

                entry.mSbjctGenomeFrom += sbjct_offset_from
                entry.mSbjctGenomeTo += sbjct_offset_from

            if param_just_dump:
                print str(entry)
            else:
                if param_loglevel >= 1:
                    print "# received\t%s" % str(entry)

                entries.append(entry)

        ProcessChunk(entries)

    elif param_format == "predictions":

        last_sbjct_token = None

        for line in sys.stdin:
            if line[0] == "#":
                continue

            entry = PredictionParser.PredictionParserEntry(expand=0)
            try:
                entry.Read(line)
            except ValueError:
                print "# warning: parsing error in line %s" % line[:-1]
                continue

            if last_sbjct_token != entry.mSbjctToken:
                if entries:
                    ProcessChunk(entries)
                entries = []
                last_sbjct_token = entry.mSbjctToken

            if param_just_dump:
                print str(entry)
            else:
                entries.append(entry)

        ProcessChunk(entries)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
