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
gpipe/exonerate_combine_regions.py - 
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

   python gpipe/exonerate_combine_regions.py --help

Type::

   python gpipe/exonerate_combine_regions.py --help

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

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/exonerate_combine_regions.py 18 2005-08-09 15:32:24Z andreas $

Analyse exonerate matches and calculate summary statistics
for an alignment.

Resolve conflicts of overlapping segments.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-g, --genome=                   file with the genomic DNA (FASTA).
-p, --peptides-fasta-file=                 file with peptide sequences (FASTA).
-i, --max-intron=               maximum intron size
-o, --max-percent-overlap=      maximal percentage overlap for conflict resolution.
-c, --min-coverage-query=       minimum coverage of query
-s, --min-score=                minimum total alignment score
""" % sys.argv[0]


HEADER = """# QUERY:        query
# QFROM:        query first residue
# QTO:          query last residue
# QSTR:         query strand
# QLEN:         query length
# SBJCT:        sbjct
# SFROM:        sbjct first residue
# STO:          sbjct last residue
# SSTR:         sbjct strand
# SLEN:         sbjct length
# SCORE:        alignment score
# CQUERY:       coverage of query (in percent)
# NGAPS:        number of gaps in alignment 
# NFR:          number of frame-shifts
# NINTRON:      number of introns
# NSPLIT:       number of split codons
# PIDE:         percent identity
# PSIM:         percent similarity
# NSEGS:        number of segments
# NPER:         number of removed permuted exons"""

SHORT_HEADER_SUMMARY = """#s\tQUERY\tQFROM\tQTO\tQSTR\tQLEN\tSBJCT\tSFROM\tSTO\tSSTR\tSLEN\tSCORE\tCQUERY\tNGAPS\tNFR\tNINTRON\tNSPLIT\tPIDE\tPSIM\tNSEGS\tNPER"""
SHORT_HEADER_ENTRY = """#e\tQUERY\tQFROM\tQTO\tQSTR\tQLEN\tSBJCT\tSFROM\tSTO\tSSTR\tSLEN\tSCORE\tCQUERY\tNGAPS\tNFR\tNINTRON\tNSPLIT\tPIDE\tPSIM"""

param_stop_codons = ("TAG", "TAA", "TGA")

param_loglevel = 1

# overlap allowed for matches on genomic region
param_max_percent_overlap = 0.2

# maximum intron size
param_max_intron = 100000

# minimum score
param_min_score = 80

# minimum coverage of query
param_min_coverage_query = 10

# check for stop codons within exons
param_border_stop_codon = 3

param_long_options = ["verbose=", "help", "genome=", "peptides=",
                      "max-percent-overlap=", "min-coverage-query=", "min-score=",
                      "version"]
param_short_options = "v:hg:p:e:fo:c:s:"

param_filename_peptides = None
param_filename_genome = None


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
        this_query_from, this_query_to = segments[first_index][1:3]
        next_query_from, next_query_to = segments[first_index + 1][1:3]

        if this_query_from < next_query_to:
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
        this_query_from, this_query_to = segments[first_index][1:3]
        next_query_from, next_query_to = segments[first_index - 1][1:3]

        if this_query_to > next_query_from:
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

    However, we do not want to align to genomic stretches due
    to spurious matches. Thus, delete exon permutations at
    the beginning and the end and only take the core.

    """

    # resolve exon permutations from the front
    ntotal_segments = len(segments)

    if len(segments) > 1:
        segments = RemoveExonPermutationsFromFront(segments)
        segments = RemoveExonPermutationsFromBack(segments)

    ncleaned_segments = len(segments)
    # calculate summary statistics
    summary = list(reduce(lambda x, y: (x[0],
                                        min(x[1], y[1]),
                                        max(x[2], y[2]),
                                        x[3],
                                        x[4],              # query_length
                                        x[5],
                                        min(x[6], y[6]),   # sbjct_from
                                        max(x[7], y[7]),   # sbjct_to
                                        x[8],              # sbjct_strand
                                        x[9],              # sbjct_length
                                        x[10] + y[10],       # score
                                        x[11] + y[11],       # query_coverage
                                        x[12] + y[12],
                                        x[13] + y[13],
                                        x[14] + y[14],
                                        x[15] + y[15],
                                        x[16] + y[16],
                                        x[17] + y[17],
                                        ), segments))

    summary[16] /= ncleaned_segments
    summary[17] /= ncleaned_segments

    exon_boundaries = map(lambda x: (x[6], x[7]), segments)

    (query_token, min_query_from, max_query_to, query_strand, query_length,
     sbjct_token, min_sbjct_from, max_sbjct_to, sbjct_strand, sbjct_length,
     total_score,
     total_coverage_query,
     total_ngaps, total_nframeshifts, total_nintrons, total_nsplits,
     avg_percent_identity, avg_percent_similarity) = summary

    summary = tuple(summary)

    if param_loglevel >= 2:
        print "#s\t" + string.join(map(str, summary + (ntotal_segments, ntotal_segments - ncleaned_segments)), "\t")

    return [summary + (ntotal_segments, ntotal_segments - ncleaned_segments, exon_boundaries)]

# ------------------------------------------------------------


def CheckPrediction(prediction):
    """check if prediction is sensible.

    Apply thresholds to check validity of prediction:

    1. minimum score
    2. minimum coverage
    """

    (query_token, min_query_from, max_query_to, query_strand, query_length,
     sbjct_token, min_sbjct_from, max_sbjct_to, sbjct_strand, sbjct_length,
     total_score, total_coverage_query,
     total_ngaps, total_nframeshifts, total_nintrons, total_nsplits,
     avg_percent_identity, avg_percent_similarity,
     nsegments, npermutations, exon_boundaries) = prediction

    if total_score < param_min_score or \
            total_coverage_query < param_min_coverage_query:
        return 0

    return 1

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
        elif o in ("-g", "--genome"):
            param_filename_genome = a
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-e", "--extend"):
            param_extend = int(a)
        elif o in ("-i", "--max-intron"):
            param_max_intron = int(a)
        elif o in ("-c", "--min-coverage-query"):
            param_min_coverage_query = float(a)
        elif o in ("-s", "--min-score"):
            param_min_total_score = float(a)

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

    print HEADER

    if param_loglevel >= 2:
        print SHORT_HEADER_SUMMARY

    last_query_token = None
    last_sbjct_token = None
    last_sbjct_strand = None
    last_query_to = 0
    last_sbjct_to = 0

    # array with final predictions
    predictions = []

    # aligned segments from exonerate
    segments = []

    for line in sys.stdin:

        if line[0] == "#":
            continue

        data = string.split(line[:-1], "\t")

        (query_token, query_length, query_from, query_to, query_lali,
         sbjct_token, sbjct_length, sbjct_from, sbjct_to, sbjct_lali, sbjct_strand,
         score, rank, percent_identity, percent_similarity,
         equivalent_total, equivalent_identity, equivalent_similarity, equivalent_mismatches) = data[0:19]

        # skip over next entries, just take query_strand
        (xquery_token, xquery_from, xquery_to, query_strand,
         xsbjct_token, xsbjct_from, xsbjct_to, xsbjct_strand,
         xscore) = data[20:29]

        (query_from, query_to, query_length,
         sbjct_from, sbjct_to, sbjct_length,
         score, rank) =\
            map(int, (query_from, query_to, query_length,
                      sbjct_from, sbjct_to, sbjct_length,
                      score, rank))

        percent_identity, percent_similarity = map(
            float, (percent_identity, percent_similarity))

        query_key = re.match("(\S+)", query_token).groups()[0]

        alignment = data[29:]

        # number of frameshifts
        nframeshifts = 0

        # number of split codons
        nsplitcodons = 0

        # number of stop codons not within exon boundary
        nstopcodons = 0

        # number of introns
        nintrons = 0

        # number of split codons
        nsplits = 0

        # number of query
        coverage_query = (query_to - query_from) * 100 / query_length

        # number of gaps
        ngaps = 0

        current_pos_protein = query_from
        current_pos_genome = sbjct_from

        # setup up sequence to use
        if sbjct_strand == "+":
            genomic_sequence = forward_sequences[sbjct_token]
        else:
            genomic_sequence = reverse_sequences[sbjct_token]

        for x in range(0, len(alignment), 3):

            state, l_protein, l_genome = alignment[x:x + 3]
            l_protein = int(l_protein)
            l_genome = int(l_genome)
            # matching state
            if state == "M":
                # check for stop codons
                if genomic_sequence:
                    for y in range(current_pos_genome + param_border_stop_codon, current_pos_genome + l_genome - param_border_stop_codon, 3):
                        if genomic_sequence[y:y + 3] in param_stop_codons:
                            nstopcodons += 1

            elif state == "I":
                nintrons += 1

            elif state == "F":
                nframeshifts += 1

            elif state == "G":
                ngaps += 1

            elif state == "S":
                nsplits += 1

        # check, if we are within the same "gene"
        # same gene is:
        # * same query, same chromosome, same strand
        # * gap not longer than param_max_intron
        if last_query_token != query_token or \
                last_sbjct_token != sbjct_token or \
                last_sbjct_strand != sbjct_strand or \
                (sbjct_from - last_sbjct_to) > param_max_intron:

            if last_query_token:
                predictions += ProcessSegments(segments)

            segments = []

            if param_loglevel >= 3:
                print SHORT_HEADER_ENTRY

        else:
            if last_query_to > query_from:
                if param_loglevel >= 1:
                    print "# WARNING: exon permutation in alignment between %s and %s" % (query_token, sbjct_token)

        last_query_token = query_token
        last_sbjct_token = sbjct_token
        last_sbjct_strand = sbjct_strand
        last_sbjct_to = sbjct_to
        last_query_to = query_to

        segment = (query_token, query_from, query_to, query_strand, query_length,
                   sbjct_token, sbjct_from, sbjct_to, sbjct_strand, sbjct_length,
                   score, coverage_query,
                   ngaps, nframeshifts, nintrons, nsplits,
                   percent_identity, percent_similarity)

        if param_loglevel >= 3:
            print "#e\t" + string.join(map(str, segment), "\t")

        segments.append(segment)

    if last_query_token:
        predictions += ProcessSegments(segments)

    if param_loglevel >= 1:
        print "# number of predictions: %i" % len(predictions)

    for prediction in predictions:
        print string.join(map(str, prediction), "\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
