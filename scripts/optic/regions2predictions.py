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
optic/regions2predictions.py - 
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

   python optic/regions2predictions.py --help

Type::

   python optic/regions2predictions.py --help

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
import tempfile
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Intervalls as Intervalls
import CGAT.PredictionParser as PredictionParser
import CGAT.PredictionFile as PredictionFile
import alignlib_lite

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: optic/regions2predictions.py 1799 2008-03-28 11:44:19Z andreas $

Resolve conflicts of overlapping exonerate segments.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --peptides-fasta-file=                 file with peptide sequences (FASTA).
-b, --benchmark=                benchmarking output. File in #
-y, --benchmark-synonyms=       list of synonymous benchmarking ids.
-o, --max-percent-overlap=      maximal percentage overlap for conflict resolution.
-c, --min-coverage-query=       minimum coverage of query
-s, --min-score=                minimum total alignment score
-i, --min-percent-identity=     minimum percent identity
-m, --max-matches=              maximum number of matches per query to show [0: all]
-j, --join-regions=             join regions of same query token and less than # residues apart
-x, --conserve-memory           conserve memory.
--disable-conflict              turn of resolution of conflicts
--disable-overlap               turn of resolution of overlaps
--disable-suboptimal            turn of elimination of suboptimal predictions
--disable-activation            turn of reactivation of eliminated queries
""" % sys.argv[0]

param_loglevel = 2

# overlap allowed for matches on genomic region
param_max_percent_overlap = 20
param_gop = -10.0
param_gep = -2.0
param_min_score_overlap = 200

# threshold for filtering bad predictions:

# minimum score
param_min_total_score = 80

# joining regions
param_join_regions = 0

# minimum coverage of query
param_min_coverage_query = 10

# conserve memory
param_conserve_memory = 0

# minimum percent identity
param_min_percent_identity = 0

# minimum length
param_min_length = 0

param_max_matches = 0

param_long_options = ["verbose=", "help", "max-percent-overlap=",
                      "min-coverage-query=", "min-score=", "min-percent-identity=",
                      "max-matches=", "peptides=", "min-length=",
                      "disable-conflict", "disable-overlap", "disable-suboptimal", "disable-activation",
                      "join-regions=", "conserve-memory",
                      "benchmark=", "benchmark-synonyms=", "test=", "version"]
param_short_options = "v:ho:c:s:i:m:p:j:xb:s:t:"

param_filename_peptides = None

param_filter_overlaps = 1
param_filter_conflicts = 1
param_filter_suboptimal = 1
param_reactivate_missed = 1

param_min_relative_coverage = 0.5
param_min_relative_score = 0.5
param_min_relative_percent_identity = 0.5

# minimum difference between non-correlated conflicts to keep them both.
param_conflicts_min_difference = 0.1

# benchmarking data
param_benchmarks = None
param_benchmark_synonyms = None
param_filename_benchmark = None
param_filename_benchmark_synonyms = None

param_test = None

# ------------------------------------------------------------


def CheckBenchmark(p, kept_p=None):

    key = "%s;%s" % (p.mSbjctToken, p.mSbjctStrand)

    if not param_benchmarks.has_key(key):
        return

    matches = param_benchmarks[key]
    x = 0
    while x < len(matches) and matches[x][1] < p.mSbjctGenomeFrom:
        x += 1

    found = 0
    while x < len(matches) and matches[x][0] < p.mSbjctGenomeTo:
        if matches[x][2] == p.mQueryToken:
            found = 1
            break
        x += 1

    status = "UNK"
    if kept_p and param_benchmark_synonyms:
        if param_benchmark_synonyms.has_key(p.mQueryToken) and \
                param_benchmark_synonyms.has_key(kept_p.mQueryToken):
            if param_benchmark_synonyms[p.mQueryToken] == param_benchmark_synonyms[kept_p.mQueryToken]:
                status = "SYN "
            else:
                status = "DIF "
    if found:
        print "# %s BENCHMARK REMOVED: %s (%i-%i):" % (status, matches[x][2], matches[x][0], matches[x][1]), str(p)

    return found
# ------------------------------------------------------------


def EvaluateBenchmark(predictions, loglevel=2):
    """evaluate.
    """

    if not param_benchmarks:
        return

    # build hash of all benchmarking regions
    # hash is indexed by
    # sbjct_token;sbjct_strand;sbjct_genome_from;sbjct_genome_to;query_token
    benchmark_regions = {}
    for key in param_benchmarks.keys():
        for sbjct_genome_from, sbjct_genome_to, query_token in param_benchmarks[key]:
            benchmark_regions[
                "%s;%i;%i;%s" % (key, sbjct_genome_from, sbjct_genome_to, query_token)] = []

    found = 0
    nchecked = 0
    # count matches to benchmarking regions
    for p in predictions:

        nchecked += 1

        key = "%s;%s" % (p.mSbjctToken, p.mSbjctStrand)

        if not param_benchmarks.has_key(key):
            continue

        # iterate over matches to chromosome and strand
        # and record matches
        matches = param_benchmarks[key]

        x = 0
        while x < len(matches) and matches[x][1] < p.mSbjctGenomeFrom:
            x += 1

        while x < len(matches) and matches[x][0] < p.mSbjctGenomeTo:
            big_key = "%s;%i;%i;%s" % ((key,) + matches[x])
            benchmark_regions[big_key].append(p.mQueryToken)
            if loglevel >= 3:
                print "# BENCHMARK FOUND: %s (%i-%i):" % (matches[x][2], matches[x][0], matches[x][1]), str(p)
            x += 1

    # evaluate matches
    total_nfound = 0
    total_nredundant = 0
    total_nmissed = 0
    total_nnonsyn = 0
    total_nloci = 0
    for big_key in benchmark_regions.keys():
        total_nloci += 1
        n = len(benchmark_regions[big_key])

        sbjct_token, sbjct_strand, bench_from, bench_to, bench_token = big_key.split(
            ";")

        bench_synonym = param_benchmark_synonyms[bench_token]
        nidentical = 0
        nsyn = 0
        nnonsyn = 0
        for found_token in benchmark_regions[big_key]:

            if found_token == bench_token:
                nidentical += 1
            else:
                if bench_synonym == param_benchmark_synonyms[found_token]:
                    nsyn += 1
                else:
                    nnonsyn += 1

        nfound = nsyn + nidentical

        if loglevel >= 2:
            if n > 0:
                print "# BENCHMARK %s: found=%i, identical=%i, synonymous=%i, mismatches=%i:" %\
                      (big_key, nfound, nidentical, nsyn, nnonsyn), string.join(
                          benchmark_regions[big_key], ",")

        if nfound:
            total_nfound += 1
        else:
            total_nmissed += 1

        if nfound > 1:
            total_nredundant += 1
        if nnonsyn > 0:
            total_nnonsyn += 1

    if loglevel >= 1:
        print "# BENCHMARK SUMMARY: nchecked=%i, nloci=%i, nfound=%i, nmissed=%i, nredundant=%i, mismatches=%i" % \
              (nchecked, total_nloci, total_nfound,
               total_nmissed, total_nredundant, total_nnonsyn)

# ------------------------------------------------------------


def CombineSegments(predictions):
    """combine predictions to a new prediction.
    """

    # combine segments
    s = segments[0]
    s.mMapPeptide2Translation.clear()

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

    return s

# ------------------------------------------------------------


def JoinRegions(predictions):
    # sort entries by genomic region
    predictions.sort(lambda x, y: cmp((x.mQueryToken, x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
                                      (y.mQueryToken, y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

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

        if is_new_chunk:
            predictions.append(CombinePredictions(segments))
            segments = []

        segments.append(entry)
        last_entry = entry

    predictions.append(CombinePredictions(segments))

    return predictions

# ------------------------------------------------------------


def FilterBadPredictions(predictions):

    new_predictions = []
    removed_predictions = []

    for p in predictions:

        # check if prediction looks ok:
        if p.score < param_min_total_score:
            if param_loglevel >= 2:
                print "# PRUNING: reason: score below minimum: removing: %s" % str(p)
            removed_predictions.append(p)
            continue
        elif p.mQueryCoverage < param_min_coverage_query:
            if param_loglevel >= 2:
                print "# PRUNING: reason: coverage below minimum: removing: %s" % str(p)
            removed_predictions.append(p)
            continue
        elif p.mPercentIdentity < param_min_percent_identity:
            if param_loglevel >= 2:
                print "# PRUNING: reason: percent identity below minimum: removing: %s" % str(p)
            removed_predictions.append(p)
            continue

        new_predictions.append(p)

    return new_predictions, removed_predictions

# ------------------------------------------------------------


class MyBestEntry:
    mQueryCoverage = 0
    score = 0
    mPercentIdentity = 0

    def __init__(self):
        pass

# ------------------------------------------------------------


def FilterSuboptimal(old_predictions,
                     new_predictions,
                     removed_predictions,
                     min_relative_coverage=0.0,
                     min_relative_score=0.0,
                     min_relative_pide=0.0):
    """remove suboptimal alignments.

    """

    best_predictions = {}

    for p in old_predictions:
        if not best_predictions.has_key(p.mQueryToken):
            best_predictions[p.mQueryToken] = MyBestEntry()

        x = best_predictions[p.mQueryToken]
        x.mQueryCoverage = max(x.mQueryCoverage, p.mQueryCoverage)
        x.score = max(x.score, p.score)
        x.mPercentIdentity = max(x.mPercentIdentity, p.mPercentIdentity)

    nnew = 0
    for p in old_predictions:
        x = best_predictions[p.mQueryToken]

        if p.mQueryCoverage / x.mQueryCoverage < min_relative_coverage:
            if param_loglevel >= 2:
                print "# PRUNING: reason: coverage below best: removing %s" % str(p)
            if param_benchmarks:
                CheckBenchmark(p)
            removed_predictions.append(p)
            continue

        if p.score / x.score < min_relative_score:
            if param_loglevel >= 2:
                print "# PRUNING: reason: score below best: removing %s" % str(p)
            if param_benchmarks:
                CheckBenchmark(p)
            removed_predictions.append(p)
            continue

        if p.mPercentIdentity / x.mPercentIdentity < min_relative_pide:
            if param_loglevel >= 2:
                print "# PRUNING: reason: percent identity below best: removing %s" % str(p)
            if param_benchmarks:
                CheckBenchmark(p)
            removed_predictions.append(p)
            continue

        new_predictions.append(p)
        nnew += 1

    return nnew

# ------------------------------------------------------------


def FilterMaxMatches(predictions, max_matches):

    predictions.sort(lambda x, y: cmp((x.mQueryToken, -x.mQueryCoverage, -x.score),
                                      (y.mQueryToken, -y.mQueryCoverage, -y.score)))

    new_predictions = []
    removed_predictions = []

    last_token = None
    n = 0

    for p in predictions:
        if last_token != p.mQueryToken:
            n = 0
            last_token = p.mQueryToken
        n += 1
        if n <= param_max_matches:
            new_predictions.append(p)
        else:
            removed_predictions.append(p)
            if param_benchmarks:
                CheckBenchmark(p)

    return new_predictions, removed_predictions

# ------------------------------------------------------------


def FilterOverlaps(old_predictions, new_predictions, removed_predictions):
    """remove overlapping entries.

    Overlapping entries align to the same query. This usually means repeats of
    some sort.
    """

    # sort predictions by query and genomic region
    if isinstance(old_predictions, PredictionFile.PredictionFile):
        old_predictions.sort(
            ('mQueryToken', 'mSbjctToken', 'mSbjctStrand', 'mSbjctGenomeFrom'))
    else:
        old_predictions.sort(lambda x, y: cmp((x.mQueryToken, x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
                                              (y.mQueryToken, y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

    nnew = 0
    last_prediction = None
    min_to = None

    for p in old_predictions:
        if not last_prediction:
            last_prediction = p
            min_to = p.mSbjctGenomeTo
            continue

        if last_prediction.mQueryToken == p.mQueryToken and \
            last_prediction.mSbjctToken == p.mSbjctToken and \
            last_prediction.mSbjctStrand == p.mSbjctStrand and \
            min_to > p.mSbjctGenomeFrom and \
            (min(last_prediction.mSbjctGenomeTo, p.mSbjctGenomeTo) -
                max(last_prediction.mSbjctGenomeFrom, p.mSbjctGenomeFrom)) > 0:

            min_to = min(min_to, p.mSbjctGenomeTo)
            if last_prediction.mQueryCoverage > p.mQueryCoverage:
                if param_loglevel >= 2:
                    print "# OVERLAP: kept %s(%s-%i, removed: %s" % (str(last_prediction.mPredictionId),
                                                                     last_prediction.mQueryToken,
                                                                     last_prediction.mSbjctGenomeFrom,
                                                                     str(p))
                if param_benchmarks:
                    CheckBenchmark(p)
                removed_predictions.append(p)
                continue
            else:
                if param_loglevel >= 2:
                    print "# OVERLAP: kept %s(%s-%i), removed: %s" % (str(p.mPredictionId),
                                                                      p.mQueryToken,
                                                                      p.mSbjctGenomeFrom,
                                                                      str(last_prediction))
                if param_benchmarks:
                    CheckBenchmark(last_prediction)
                removed_predictions.append(last_prediction)
                last_prediction = p
                continue

        new_predictions.append(last_prediction)
        min_to = p.mSbjctGenomeTo
        nnew += 1
        last_prediction = p

    new_predictions.append(last_prediction)
    nnew += 1

    return nnew

# ------------------------------------------------------------


def FilterConflicts(old_predictions, new_predictions, removed_predictions,
                    min_overlap, peptide_sequences):
    """remove conflicts.

    Remove overlapping entries between different queries.

    Only remove those sequences, which are alignable.

    If they are alignable, take the sequence with the highest score and highest coverage.
    (Take both, if score and coverage are not correlated.)
    """
    ##########################################################################
    # sort predictions by genomic region
    if isinstance(old_predictions, PredictionFile.PredictionFile):
        old_predictions.sort(
            ('mSbjctToken', 'mSbjctStrand', 'mSbjctGenomeFrom', 'mSbjctGenomeTo'))
    else:
        old_predictions.sort(lambda x, y: cmp((x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom, x.mSbjctGenomeTo),
                                              (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom, y.mSbjctGenomeTo)))

    ##########################################################################
    # filter predictions and resolve conflicts based on genomic overlap
    # deleted segments are put in a temporary storage space.
    alignator = alignlib_lite.makeAlignatorDPFull(
        alignlib_lite.ALIGNMENT_LOCAL, param_gop, param_gep)
    result = alignlib_lite.makeAlignmentVector()
    alignments = {}
    noverlaps = 0
    nredundants = 0

    nnew = 0
    last_prediction = None

    for this_prediction in old_predictions:
        try:
            (this_query_peptide, this_query_status, this_query_gene,
             this_query_transcript) = \
                re.split("\s+", this_prediction.mQueryToken)
        except ValueError:
            this_query_gene = None

        if not last_prediction:
            last_prediction = this_prediction
            last_query_gene = this_query_gene
            continue

        overlap = min(last_prediction.mSbjctGenomeTo,
                      this_prediction.mSbjctGenomeTo) -\
            max(last_prediction.mSbjctGenomeFrom,
                this_prediction.mSbjctGenomeFrom)
        union = max(last_prediction.mSbjctGenomeTo,
                    this_prediction.mSbjctGenomeTo) -\
            min(last_prediction.mSbjctGenomeFrom,
                this_prediction.mSbjctGenomeFrom)

        # resolve overlap between different genes
        if overlap > 0 and \
                (last_query_gene != this_query_gene or
                 last_query_gene is None):

            noverlaps += 1
            relative_overlap = 100 * overlap / union

            # Start conflict resolution, if overlap is above threshold.
            # Keep higher scoring segment.
            #
            # Check if queries are homologous.
            if relative_overlap >= param_max_percent_overlap:

                if peptide_sequences:
                    if last_prediction.mQueryToken < this_prediction.mQueryToken:
                        key = "%s-%s" % (last_prediction.mQueryToken,
                                         this_prediction.mQueryToken)
                    else:
                        key = "%s-%s" % (this_prediction.mQueryToken,
                                         last_prediction.mQueryToken)

                    if not alignments.has_key(key):
                        result.clear()
                        alignator.align(result,
                                        alignlib_lite.makeSequence(
                                            peptide_sequences[this_prediction.mQueryToken]),
                                        alignlib_lite.makeSequence(peptide_sequences[last_prediction.mQueryToken]))
                        alignments[key] = result.getScore()
                        if result.getScore() >= param_min_score_overlap:
                            nredundants += 1

                    if alignments[key] >= param_min_score_overlap:
                        is_overlap = 1
                    else:
                        is_overlap = 0
                else:
                    is_overlap = 1
            else:
                is_overlap = 0
        else:
            is_overlap = 0

        if is_overlap:
            # take best prediction. If difference is very small, set
            # difference to 0 (difference does not matter). In this case,
            # the first prediction is taken.
            d1 = last_prediction.mQueryCoverage - \
                this_prediction.mQueryCoverage
            if float(abs(d1)) / float(last_prediction.mQueryCoverage) < param_conflicts_min_difference:
                d1 = 0
            d2 = last_prediction.score - this_prediction.score
            if float(abs(d2)) / float(this_prediction.score) < param_conflicts_min_difference:
                d2 = 0
            if d1 >= 0 and d2 >= 0:
                if param_loglevel >= 2:
                    print "# CONFLICT: kept %i(%s-%i), overlap=%i(%5.2f), removed: %s" % (last_prediction.mPredictionId,
                                                                                          last_prediction.mQueryToken,
                                                                                          last_prediction.mSbjctGenomeFrom,
                                                                                          overlap, relative_overlap,
                                                                                          str(this_prediction))
                if param_benchmarks:
                    if CheckBenchmark(this_prediction, last_prediction):
                        print "# BENCHMARK KEPT with overlap=%i(%5.2f): %s" % (overlap, relative_overlap,
                                                                               str(last_prediction))

                removed_predictions.append(this_prediction)
                continue
            elif d1 <= 0 and d2 <= 0:
                if param_loglevel >= 2:
                    print "# CONFLICT: kept %i(%s-%i), overlap=%i(%5.2f), removed: %s" % (this_prediction.mPredictionId,
                                                                                          this_prediction.mQueryToken,
                                                                                          this_prediction.mSbjctGenomeFrom,
                                                                                          overlap, relative_overlap,
                                                                                          str(last_prediction))
                if param_benchmarks:
                    if CheckBenchmark(last_prediction, this_prediction):
                        print "# BENCHMARK KEPT with overlap=%i(%5.2f): %s" % (overlap, relative_overlap,
                                                                               str(this_prediction))
                removed_predictions.append(last_prediction)
                last_prediction = this_prediction
                last_query_gene = this_query_gene
                continue
            else:
                if param_loglevel >= 2:
                    print "# CONFLICT: non-correlated score/coverage. Keeping both %i(%s-%i) (%5.2f/%i/%i) and %i(%s-%i) (%5.2f/%i/%i)" % \
                          (this_prediction.mPredictionId,
                           this_prediction.mQueryToken, this_prediction.mSbjctGenomeFrom,
                           this_prediction.score, this_prediction.mQueryCoverage,
                           this_prediction.mPercentIdentity,
                           last_prediction.mPredictionId,
                           last_prediction.mQueryToken, last_prediction.mSbjctGenomeFrom,
                           last_prediction.score, last_prediction.mQueryCoverage,
                           last_prediction.mPercentIdentity)

        new_predictions.append(last_prediction)
        nnew += 1
        last_query_gene = this_query_gene
        last_prediction = this_prediction

    new_predictions.append(last_prediction)
    nnew += 1

    if param_loglevel >= 1:
        print "# calculated %i alignments for %i potential conflicts (%i above threshold)" % \
              (len(alignments), noverlaps, nredundants)

    return nnew

# ------------------------------------------------------------


def ReadBenchmarkingRegions(infile):
    """read benchmarking regions.
    """

    benchmarks = {}

    for line in infile:
        if line[0] == "#":
            continue

        (query_token, sbjct_token, sbjct_strand, sbjct_genome_from,
         sbjct_genome_to) = line[:-1].split("\t")
        sbjct_genome_from, sbjct_genome_to = int(
            sbjct_genome_from), int(sbjct_genome_to)

        key = "%s;%s" % (sbjct_token, sbjct_strand)
        if not benchmarks.has_key(key):
            benchmarks[key] = []
        benchmarks[key].append(
            (sbjct_genome_from, sbjct_genome_to, query_token))

    for key in benchmarks.keys():
        benchmarks[key].sort()

    return benchmarks


def ExchangeStreams(old_predictions, new_predictions):
    """close streams and exchange them, repopen them.
    """
    old_predictions.close()
    new_predictions.close()

    os.system("mv %s %s" %
              (new_predictions.GetFileName(), old_predictions.GetFileName()))

    old_predictions.open(mode="r")
    new_predictions.open(mode="w")

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
        elif o in ("-b", "--benchmark"):
            param_filename_benchmark = a
        elif o in ("-y", "--benchmark-synonyms"):
            param_filename_benchmark_synonyms = a
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-c", "--min-coverage-query"):
            param_min_coverage_query = float(a)
        elif o in ("-s", "--min-score"):
            param_min_total_score = float(a)
        elif o in ("-i", "--min-percent-identity"):
            param_min_percent_identity = float(a)
        elif o in ("-o", "--max-percent-overlap"):
            param_max_percent_overlap = float(a)
        elif o in ("-m", "--max-matches"):
            param_max_matches = int(a)
        elif o in ("-j", "--join-regions"):
            param_join_regions = int(a)
        elif o == "--disable-conflict":
            param_filter_conflicts = 0
        elif o == "--disable-suboptimal":
            param_filter_suboptimal = 0
        elif o == "--disable-overlap":
            param_filter_overlaps = 0
        elif o == "--disable-activation":
            param_reactivate_missed = 0
        elif o == "--conserve-memory":
            param_conserve_memory = 1
        elif o == "--min-interval-length":
            param_min_length = int(a)
        elif o == "--test":
            param_test = int(a)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()
    sys.stdout.flush()

    # read benchmarking regions
    if param_filename_benchmark:
        param_benchmarks = ReadBenchmarkingRegions(
            open(param_filename_benchmark, "r"))
        if param_loglevel >= 1:
            print "# read benchmarking regions for %i tokens" % len(param_benchmarks)
            sys.stdout.flush()
        if param_filename_benchmark_synonyms:
            infile = open(param_filename_benchmark_synonyms, "r")
            param_benchmark_synonyms = {}
            for line in infile:
                if line[0] == "#":
                    continue
                value, key = line[:-1].split("\t")
                param_benchmark_synonyms[key] = value
        else:
            param_benchmark_synonyms = {}
    else:
        param_benchmarks = {}
        param_benchmark_synonyms = {}

    # read peptide sequences
    if param_filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(param_filename_peptides, "r"))
    else:
        peptide_sequences = {}

    if param_conserve_memory:
        old_predictions, filename_old_predictions = tempfile.mkstemp()
        os.close(old_predictions)
        old_predictions = PredictionFile.PredictionFile()
        old_predictions.open(filename_old_predictions, "w")
    else:
        # array with final predictions
        old_predictions = []

    if param_loglevel >= 1:
        print "# reading predictions."
        sys.stdout.flush()

    nread = 0
    ninput = 0
    for line in sys.stdin:

        if line[0] == "#":
            continue

        entry = PredictionParser.PredictionParserEntry(expand=0)
        entry.Read(line)
        nread += 1

        # set prediction id
        if not entry.mPredictionId:
            entry.mPredictionId = nread

        # filter bad predictions right here in order to save memory:
        if entry.score < param_min_total_score:
            if param_loglevel >= 2:
                print "# PRUNING: reason: score below minimum: removing: %s" % str(entry)
            continue
        elif entry.mQueryCoverage < param_min_coverage_query:
            if param_loglevel >= 2:
                print "# PRUNING: reason: coverage below minimum: removing: %s" % str(entry)
            continue
        elif entry.mPercentIdentity < param_min_percent_identity:
            if param_loglevel >= 2:
                print "# PRUNING: reason: percent identity below minimum: removing: %s" % str(entry)
            continue
        elif entry.mSbjctTo - entry.mSbjctFrom + 1 < param_min_length:
            if param_loglevel >= 2:
                print "# PRUNING: reason: length of transcript below minimum: removing: %s" % str(entry)
            continue

        ninput += 1

        if param_test and ninput > param_test:
            break

        old_predictions.append(entry)

    if param_loglevel >= 1:
        print "# number of good predictions read: %i" % ninput
        sys.stdout.flush()

    if ninput == 0:
        print "# ERROR: no predictions"
        sys.exit(1)

    if param_conserve_memory:
        old_predictions.close()
        old_predictions.open(mode="r")
        removed_predictions, filename_removed_predictions = tempfile.mkstemp()
        os.close(removed_predictions)
        removed_predictions = PredictionFile.PredictionFile()
        removed_predictions.open(filename_removed_predictions, "w")

        new_predictions, filename_new_predictions = tempfile.mkstemp()
        os.close(new_predictions)
        new_predictions = PredictionFile.PredictionFile()
        new_predictions.open(filename_new_predictions, "w")
    else:
        removed_predictions = []
        new_predictions = []

    EvaluateBenchmark(old_predictions)

    ##########################################################################
    # remove bad predictions (low coverage, low score, etc).
    # predictions, x = FilterBadPredictions( predictions )

    # do not use removed predictions
    nbad = ninput

    ##########################################################################
    # remove overlapping entries
    if param_filter_overlaps:
        noverlaps = FilterOverlaps(
            old_predictions, new_predictions, removed_predictions)
        if param_conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []
        EvaluateBenchmark(old_predictions)
    else:
        noverlaps = nbad

    if param_loglevel >= 1:
        print "# number of predictions after overlap: %i" % noverlaps
        sys.stdout.flush()
    ##########################################################################
    # take only the best n predictions per query
    # sort by query_token, coverage_query and score
    if param_max_matches > 0:
        nmax_matches = FilterMaxMatches(
            old_predictions, new_predictions, param_max_matches)
        if param_conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []
        EvaluateBenchmark(old_predictions)
    else:
        nmax_matches = noverlaps

    if param_loglevel >= 1:
        print "# number of predictions after max matches: %i" % nmax_matches
        sys.stdout.flush()
    ##########################################################################
    # remove overlapping exons
    if param_filter_conflicts:
        nconflicts = FilterConflicts(old_predictions, new_predictions, removed_predictions,
                                     param_max_percent_overlap, peptide_sequences)
        if param_conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []
        EvaluateBenchmark(old_predictions)
    else:
        nconflicts = nmax_matches

    ##########################################################################
    # remove suboptimal matches
    if param_filter_suboptimal:
        nsuboptimal = FilterSuboptimal(old_predictions, new_predictions, removed_predictions,
                                       param_min_relative_coverage,
                                       param_min_relative_score,
                                       param_min_relative_percent_identity)
        if param_conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []
        EvaluateBenchmark(old_predictions)
    else:
        nsuboptimal = nconflicts

    ##########################################################################
    # get accepted transcripts and print out predictions
    if param_loglevel >= 1:
        print "# counting transcripts for %i predictions" % nsuboptimal
        sys.stdout.flush()

    accepted_transcripts = {}
    for p in old_predictions:
        accepted_transcripts[p.mQueryToken] = 1

    naccepted = len(accepted_transcripts)

    ##########################################################################
    # post-processing failed transcripts
    # print best predictions for transcripts that have no other match so far.
    ntactivated = 0
    nactivated = 0

    if param_reactivate_missed:

        if param_loglevel >= 1:
            print "# activating missed transcripts."
            sys.stdout.flush()

        # sort by peptide and score (descending)
        if param_conserve_memory:
            old_predictions.close()
            old_predictions.open(mode="a")
            removed_predictions.sort(('mQueryToken', 'score'))
        else:
            removed_predictions.sort(lambda x, y: cmp((x.mQueryToken, -x.score),
                                                      (y.mQueryToken, -y.score)))
        last_query_token = None
        for p in removed_predictions:

            # skip segments that have been printed out.
            if accepted_transcripts.has_key(p.mQueryToken):
                continue

            # skip lower scoring segments
            if last_query_token == p.mQueryToken:
                continue

            old_predictions.append(p)
            ntactivated += 1
            nactivated += 1

            last_query_token = p.mQueryToken

        if param_conserve_memory:
            old_predictions.close()
            old_predictions.open(mode="r")
        EvaluateBenchmark(old_predictions)

    nactivated += nsuboptimal

    ##########################################################################
    # join regions
    if param_join_regions:
        if param_loglevel >= 1:
            print "# joining regions" % nsuboptimal
            sys.stdout.flush()
        new_predictions = []
        njoined = JoinRegions(old_predictions, new_predictions)
        if param_conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []
    else:
        njoined = nactivated

    # finally: print everything
    for p in old_predictions:
        p.Write()

    if param_loglevel >= 1:
        print "# pairs: nread=%i, input=%i, bad=%i, overlaps=%i, max_matches=%i, conflicts=%i, suboptimal=%i, activated=%i, joined=%i" % \
              (nread, ninput, nbad, noverlaps, nmax_matches,
               nconflicts, nsuboptimal, nactivated, njoined)
        print "# transcripts: accepted=%i, activated=%i, total=%i" % (naccepted, ntactivated, naccepted + ntactivated)

    EvaluateBenchmark(old_predictions)

    if param_conserve_memory:
        os.remove(filename_old_predictions)
        os.remove(filename_new_predictions)
        os.remove(filename_removed_predictions)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
