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
optic/regions2graph.py - 
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

   python optic/regions2graph.py --help

Type::

   python optic/regions2graph.py --help

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
import tempfile
import gzip
import CGAT.Experiment as E
import CGAT.Exons as Exons
import CGAT.Genomics as Genomics
import CGAT.Intervalls as Intervalls
import CGAT.PredictionParser as PredictionParser
import CGAT.PredictionFile as PredictionFile
import alignlib_lite

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: optic/regions2graph.py 2754 2009-09-04 16:50:22Z andreas $

Read a graph of regions and dump out a bipartite graph.

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
--join-regions-max-coverage=    maximum coverage for regions not to be joined
-x, --conserve-memory           conserve memory.
--max-intron                    maximum intron length
--filter-regions                regions to remove
--filter-queries                ids not to test
Options for deciding when two queries are homologs:

--overlap-min-score=            minimum score for overlap
--overlap-max-coverage=         minimum maximum coverage for overlap
--overlap-min-coverage=         minimum mininum coverage for overlap
--overlap-min-identity=         minimum percent identity for overlap
""" % sys.argv[0]

global_alignments = {}
options = {}


class Region:
    mSbjctToken = ""
    mSbjctStrand = ""
    mSbjctGenomeFrom = 0
    mSbjctGenomeTo = 0

    def __init__(self):
        pass

    def __str__(self):
        return string.join((self.mSbjctToken, self.mSbjctStrand,
                            str(self.mSbjctGenomeFrom), str(self.mSbjctGenomeTo)), "\t")

# ------------------------------------------------------------


def CheckBenchmark(p, kept_p=None):

    key = "%s;%s" % (p.mSbjctToken, p.mSbjctStrand)

    if not options.benchmarks.has_key(key):
        return

    matches = options.benchmarks[key]
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
    if kept_p and options.benchmark_synonyms:
        if options.benchmark_synonyms.has_key(p.mQueryToken) and \
                options.benchmark_synonyms.has_key(kept_p.mQueryToken):
            if options.benchmark_synonyms[p.mQueryToken] == options.benchmark_synonyms[kept_p.mQueryToken]:
                status = "SYN "
            else:
                status = "DIF "
    if found:
        options.stdlog.write("# %s BENCHMARK REMOVED: %s (%i-%i): %s\n" %
                             (status, matches[x][2], matches[x][0], matches[x][1], str(p)))

    return found
# ------------------------------------------------------------


def EvaluateBenchmark(predictions, loglevel=2):
    """evaluate.
    """

    if not options.benchmarks:
        return

    # build hash of all benchmarking regions
    # hash is indexed by
    # sbjct_token;sbjct_strand;sbjct_genome_from;sbjct_genome_to;query_token
    benchmark_regions = {}
    for key in options.benchmarks.keys():
        for sbjct_genome_from, sbjct_genome_to, query_token in options.benchmarks[key]:
            benchmark_regions[
                "%s;%i;%i;%s" % (key, sbjct_genome_from, sbjct_genome_to, query_token)] = []

    found = 0
    nchecked = 0
    # count matches to benchmarking regions
    for p in predictions:

        nchecked += 1

        key = "%s;%s" % (p.mSbjctToken, p.mSbjctStrand)

        if not options.benchmarks.has_key(key):
            continue

        # iterate over matches to chromosome and strand
        # and record matches
        matches = options.benchmarks[key]

        x = 0
        while x < len(matches) and matches[x][1] < p.mSbjctGenomeFrom:
            x += 1

        while x < len(matches) and matches[x][0] < p.mSbjctGenomeTo:
            big_key = "%s;%i;%i;%s" % ((key,) + matches[x])
            benchmark_regions[big_key].append(p.mQueryToken)
            if loglevel >= 3:
                options.stdlog.write("# BENCHMARK FOUND: %s (%i-%i): %s\n" %
                                     (matches[x][2], matches[x][0], matches[x][1], str(p)))
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

        bench_synonym = options.benchmark_synonyms[bench_token]
        nidentical = 0
        nsyn = 0
        nnonsyn = 0
        for found_token in benchmark_regions[big_key]:

            if found_token == bench_token:
                nidentical += 1
            else:
                if bench_synonym == options.benchmark_synonyms[found_token]:
                    nsyn += 1
                else:
                    nnonsyn += 1

        nfound = nsyn + nidentical

        if loglevel >= 2:
            if n > 0:
                options.stdlog.write("# BENCHMARK %s: found=%i, identical=%i, synonymous=%i, mismatches=%i: %s\n" %
                                     (big_key, nfound, nidentical, nsyn, nnonsyn, string.join(benchmark_regions[big_key], ",")))

        if nfound:
            total_nfound += 1
        else:
            total_nmissed += 1

        if nfound > 1:
            total_nredundant += 1
        if nnonsyn > 0:
            total_nnonsyn += 1

    if loglevel >= 1:
        options.stdlog.write("# BENCHMARK SUMMARY: nchecked=%i, nloci=%i, nfound=%i, nmissed=%i, nredundant=%i, mismatches=%i\n" %
                             (nchecked, total_nloci, total_nfound, total_nmissed, total_nredundant, total_nnonsyn))

# ------------------------------------------------------------


def CombinePredictions(predictions):
    """combine predictions to a new prediction.
    """

    # combine segments
    s = predictions[0]
    parts = [(s.mQueryFrom, s.mQueryTo)]

    for x in predictions[1:]:
        parts.append((x.mQueryFrom, x.mQueryTo))
        s.mQueryFrom = min(x.mQueryFrom, s.mQueryFrom)
        s.mQueryTo = max(x.mQueryTo, s.mQueryTo)
        s.mSbjctGenomeFrom = min(x.mSbjctGenomeFrom, s.mSbjctGenomeFrom)
        s.mSbjctGenomeTo = max(x.mSbjctGenomeTo, s.mSbjctGenomeTo)
        s.mNIntrons += x.mNIntrons + 1
        s.mPercentIdentity += x.mPercentIdentity
        s.mPercentSimilarity += x.mPercentSimilarity
        s.score += x.score

    s.mSbjctFrom = 0
    s.mSbjctTo = 0
    s.mAlignmentString = ""
    s.mPercentIdentity /= len(predictions)
    s.mPercentSimilarity /= len(predictions)

    new_intervalls = Intervalls.CombineIntervallsLarge(parts)
    # calculate covered region
    covered = reduce(
        lambda x, y: x + y, map(lambda x: x[1] - x[0], new_intervalls))

    s.mQueryCoverage = 100 * covered / s.mQueryLength

    return s

# ------------------------------------------------------------


def JoinRegions(old_predictions, new_predictions):
    """join regions.

    Regions are joined, if they are within options.join_regions nucleotides within
    each other and the coverage of the query does not exceed a certain threshold.
    """

    # sort predictions by query and genomic region
    if isinstance(old_predictions, PredictionFile.PredictionFile):
        old_predictions.sort(
            ('mQueryToken', 'mSbjctToken', 'mSbjctStrand', 'mSbjctGenomeFrom'))
    else:
        old_predictions.sort(lambda x, y: cmp((x.mQueryToken, x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
                                              (y.mQueryToken, y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

    predictions = []
    last_p = None
    nnew = 0

    for p in old_predictions:

        is_new_chunk = False
        # check, if we are within the same "gene"
        # same gene is:
        # * same query, same chromosome, same strand
        # * gap not longer than options.join_regions
        # * do not join full lenght predictions
        if last_p:
            if last_p.mSbjctToken != p.mSbjctToken or \
               last_p.mSbjctStrand != p.mSbjctStrand or \
               last_p.mQueryToken != p.mQueryToken or \
               (p.mSbjctGenomeFrom - last_p.mSbjctGenomeTo) > options.join_regions or \
               p.mQueryCoverage > options.join_regions_max_coverage:
                is_new_chunk = True
        else:
            is_new_chunk = True

        if is_new_chunk:
            if predictions:
                new_predictions.append(CombinePredictions(predictions))
                nnew += 1
            predictions = []

        predictions.append(p)
        last_p = p

    new_predictions.append(CombinePredictions(predictions))
    nnew += 1

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

# ------------------------------------------------------------


def PrintEdges(region_id, region, predictions):
    """print a set of edges.

    This step prioritizes predictions:

        Edges are sorted by weight in descending manner, the 
                weight is the product of Coverage and Pide
        Doubles are removed
    """
    if len(predictions) == []:
        return

    predictions.sort(lambda x, y: cmp((x.mQueryCoverage * x.mPercentIdentity),
                                      (y.mQueryCoverage * y.mPercentIdentity)))
    predictions.reverse()

    queries = {}

    print_predictions = []
    # get number of matches:
    for p in predictions:
        if p.mQueryToken in queries:
            continue
        queries[p.mQueryToken] = 1
        print_predictions.append(p)

    edge_id = 0
    for p in print_predictions:
        weight = p.mQueryCoverage * p.mPercentIdentity
        edge_id += 1
        options.stdout.write("\t".join(map(str,
                                           (region_id,
                                            edge_id, len(print_predictions),
                                            region, p.mQueryToken, weight, p))) + "\n")
    return edge_id

# ------------------------------------------------------------


def ProcessRegion(predictions, region_id, region,
                  peptide_sequences=None,
                  filter_queries={}):
    """process a set of matches to a region.

    resolve region according to homology.
    """

    if options.loglevel >= 3:
        options.stdlog.write(
            "###################################################################\n")
        options.stdlog.write(
            "# resolving %i predictions in region %s\n" % (len(predictions), str(region)))
        sys.stdout.flush()

    predictions.sort(lambda x, y: cmp(x.score, y.score))
    predictions.reverse()

    alignator = alignlib_lite.makeAlignatorDPFull(
        alignlib_lite.ALIGNMENT_LOCAL, options.gop, options.gep)
    result = alignlib_lite.makeAlignmentVector()

    cluster = []

    map_sequence2cluster = range(0, len(predictions))
    edges = []

    noutput, nskipped = 0, 0

    if peptide_sequences:
        for x in range(len(predictions)):
            if options.loglevel >= 5:
                options.stdlog.write("# filtering from %i with prediction %i: %s\n" % (
                    x, predictions[x].mPredictionId, predictions[x].mQueryToken))
                sys.stdout.flush()

            if map_sequence2cluster[x] != x:
                continue

            region_id += 1
            edges = []

            if predictions[x].mQueryToken not in filter_queries:
                edges.append(predictions[x])
            else:
                nskipped += 1

            for y in range(x + 1, len(predictions)):

                if map_sequence2cluster[y] != y:
                    continue

                if predictions[x].mQueryToken < predictions[y].mQueryToken:
                    key = "%s-%s" % (predictions[x].mQueryToken,
                                     predictions[y].mQueryToken)
                else:
                    key = "%s-%s" % (predictions[y].mQueryToken,
                                     predictions[x].mQueryToken)

                # check if predictions are overlapping on the genomic sequence
                if min(predictions[x].mSbjctGenomeTo,   predictions[y].mSbjctGenomeTo) - \
                   max(predictions[x].mSbjctGenomeFrom, predictions[y].mSbjctGenomeFrom) < 0:
                    if options.loglevel >= 4:
                        options.stdlog.write("# alignment of predictions %i and %i: no overlap on genomic sequence, thus skipped\n" %
                                             (predictions[x].mPredictionId,
                                              predictions[y].mPredictionId))
                        sys.stdout.flush()
                    continue

                if not global_alignments.has_key(key):

                    seq1 = peptide_sequences[predictions[x].mQueryToken]
                    seq2 = peptide_sequences[predictions[y].mQueryToken]
                    result.clear()
                    s1 = alignlib_lite.makeSequence(seq1)
                    s2 = alignlib_lite.makeSequence(seq2)
                    alignator.align(result, s1, s2)

                    c1 = 100 * \
                        (result.getRowTo() - result.getRowFrom()) / len(seq1)
                    c2 = 100 * \
                        (result.getColTo() - result.getColFrom()) / len(seq2)
                    min_cov = min(c1, c2)
                    max_cov = max(c1, c2)

                    identity = alignlib_lite.calculatePercentIdentity(
                        result, s1, s2) * 100

                    # check if predictions overlap and they are homologous
                    if result.getScore() >= options.overlap_min_score and \
                       max_cov >= options.overlap_max_coverage and \
                       min_cov >= options.overlap_min_coverage and \
                       identity >= options.overlap_min_identity:
                        global_alignments[key] = True
                    else:
                        global_alignments[key] = False

                    if options.loglevel >= 4:
                        options.stdlog.write("# alignment=%s score=%i pid=%5.2f c1=%i c2=%i min_cov=%i max_cov=%i homolog=%s\n" %
                                             (key,
                                              result.getScore(),
                                              identity,
                                              c1, c2, min_cov, max_cov,
                                              global_alignments[key]))
                        sys.stdout.flush()

                if global_alignments[key]:
                    map_sequence2cluster[y] = x
                    if predictions[y].mQueryToken not in filter_queries:
                        edges.append(predictions[y])
                    else:
                        nskipped += 1

            noutput += PrintEdges(region_id, region, edges)

    return region_id, noutput, nskipped

# ------------------------------------------------------------


def CheckOverlap(regions, taboo_regions):
    """check whether exons lie in any of the taboo regions."""

    for first, last in regions:

        l = 0
        h = len(taboo_regions)
        while 1:
            m = (h + l) / 2

            if last < taboo_regions[m][0]:
                h = m
            elif first > taboo_regions[m][1]:
                l = m
            else:
                return True

            if h <= l + 1:
                break

    return False
# ------------------------------------------------------------


class MyBestEntry:
    mQueryCoverage = 0
    score = 0
    mPercentIdentity = 0

    def __init__(self):
        pass

# ------------------------------------------------------------


def GetBestPredictions(predictions):
    """return hash with best predictions for a query.
    """

    best_predictions = {}

    for p in predictions:
        if not best_predictions.has_key(p.mQueryToken):
            best_predictions[p.mQueryToken] = MyBestEntry()

        x = best_predictions[p.mQueryToken]
        x.mQueryCoverage = max(x.mQueryCoverage, p.mQueryCoverage)
        x.score = max(x.score, p.score)
        x.mPercentIdentity = max(x.mPercentIdentity, p.mPercentIdentity)

    return best_predictions

# ------------------------------------------------------------


def ExchangeStreams(old_predictions, new_predictions):
    """close streams and exchange them, reopen them.
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

    parser = E.OptionParser(
        version="%prog version: $Id: optic/regions2graph.py 2754 2009-09-04 16:50:22Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-b", "--benchmark", dest="filename_benchmark", type="string",
                      help="")

    parser.add_option("-y", "--benchmark-synonyms", dest="benchmark_synonyms", type="string",
                      help="")

    parser.add_option("-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="")

    parser.add_option("-c", "--min-coverage-query", dest="min_coverage_query", type="float",
                      help="")

    parser.add_option("-s", "--min-score", dest="min_total_score", type="float",
                      help="")

    parser.add_option("-i", "--min-percent-identity", dest="min_percent_identity", type="float",
                      help="")

    parser.add_option("-o", "--max-percent-overlap", dest="max_percent_overlap", type="float",
                      help="")

    parser.add_option("--overlap-min-score", dest="overlap_min_score", type="float",
                      help="")

    parser.add_option("--overlap-min-coverage", dest="overlap_min_coverage", type="float",
                      help="")

    parser.add_option("--overlap-min-identity", dest="overlap_min_identity", type="float",
                      help="")

    parser.add_option("--overlap-max-coverage", dest="overlap_max_coverage", type="float",
                      help="")

    parser.add_option("-m", "--max-matches", dest="max_matches", type="int",
                      help="")

    parser.add_option("-j", "--join-regions", dest="join_regions", type="int",
                      help="")

    parser.add_option("--join-regions-max-regions", dest="join_regions_max_regions", type="int",
                      help="")

    parser.add_option("--join-regions-max-coverage", dest="join_regions_max_coverage", type="float",
                      help="")

    parser.add_option("--min-interval-length", dest="min_length", type="int",
                      help="")

    parser.add_option("--test", dest="test", type="int",
                      help="")

    parser.add_option("--filter-queries", dest="filename_filter_queries", type="string",
                      help="")

    parser.add_option("--filter-regions", dest="filter_regions", type="string",
                      help="")

    parser.add_option("--conserve-memory", dest="conserve_memory", action="store_true",
                      help="")

    parser.add_option("--filter-suboptimal", dest="filter_suboptimal", action="store_true",
                      help="")

    parser.set_defaults(
        # overlap allowed for matches on genomic region
        max_percent_overlap=20,
        gop=-10.0,
        gep=-2.0,
        # thresholds for joining regions
        overlap_min_score=80,
        overlap_min_coverage=80,
        overlap_max_coverage=90,
        overlap_min_identity=50,
        # threshold for filtering bad predictions:
        # minimum score
        min_total_score=80,
        # joining regions
        join_regions=0,
        # maximum coverage of query for predictions to be joined
        # (This is to ensure not to join duplications. A range check
        # would be better, but runs into trouble with repeats).
        join_regions_max_coverage=90,
        # minimum coverage of query
        min_coverage_query=10,
        # conserve memory
        conserve_memory=0,
        # minimum percent identity
        min_percent_identity=0,
        # minimum length
        min_length=0,
        max_matches=0,
        filename_peptides=None,
        filename_filter_queries=None,
        # turn on/off various filters
        filter_suboptimal=False,
        filter_regions=False,
        # parameters for filter of suboptimal predictions
        min_relative_coverage=0.5,
        min_relative_score=0.5,
        min_relative_percent_identity=0.5,
        # minimum difference between non-correlated conflicts to keep them
        # both.
        conflicts_min_difference=0.1,
        # benchmarking data
        benchmarks=None,
        benchmark_synonyms=None,
        filename_benchmark=None,
        filename_benchmark_synonyms=None,
        test=None,
        max_intron=50000)

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    ##########################################################################
    # read filtering
    filter_queries = {}
    if options.filename_filter_queries:
        for line in open(options.filename_filter_queries, "r"):
            if line[0] == "#":
                continue
            query_token = line[:-1].split("\t")[0]
            filter_queries[query_token] = True

    if options.loglevel >= 1:
        options.stdlog.write(
            "# filtering for %i queries.\n" % len(filter_queries))

    ##########################################################################
    # read benchmarking regions
    if options.filename_benchmark:
        options.benchmarks = ReadBenchmarkingRegions(
            open(options.filename_benchmark, "r"))
        if options.loglevel >= 1:
            options.stdlog.write(
                "# read benchmarking regions for %i tokens\n" % len(options.benchmarks))
            sys.stdout.flush()
        if options.filename_benchmark_synonyms:
            infile = open(options.filename_benchmark_synonyms, "r")
            options.benchmark_synonyms = {}
            for line in infile:
                if line[0] == "#":
                    continue
                value, key = line[:-1].split("\t")
                options.benchmark_synonyms[key] = value
        else:
            options.benchmark_synonyms = {}
    else:
        options.benchmarks = {}
        options.benchmark_synonyms = {}

    ##########################################################################
    # read peptide sequences
    if options.filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(options.filename_peptides, "r"))
    else:
        peptide_sequences = {}

    if options.conserve_memory:
        old_predictions, filename_old_predictions = tempfile.mkstemp()
        os.close(old_predictions)
        old_predictions = PredictionFile.PredictionFile()
        old_predictions.open(filename_old_predictions, "w")
    else:
        # array with final predictions
        old_predictions = []

    if options.loglevel >= 1:
        options.stdlog.write("# reading predictions.\n")
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
        if entry.score < options.min_total_score:
            if options.loglevel >= 3:
                options.stdlog.write(
                    "# PRUNING: reason: score below minimum: removing: %s\n" % str(entry))
            continue
        elif entry.mQueryCoverage < options.min_coverage_query:
            if options.loglevel >= 3:
                options.stdlog.write(
                    "# PRUNING: reason: coverage below minimum: removing: %s\n" % str(entry))
            continue
        elif entry.mPercentIdentity < options.min_percent_identity:
            if options.loglevel >= 3:
                options.stdlog.write(
                    "# PRUNING: reason: percent identity below minimum: removing: %s\n" % str(entry))
            continue
        elif entry.mSbjctTo - entry.mSbjctFrom < options.min_length:
            if options.loglevel >= 3:
                options.stdlog.write(
                    "# PRUNING: reason: length of transcript below minimum: removing: %s\n" % str(entry))
            continue

        ninput += 1

        if options.test and ninput > options.test:
            break

        old_predictions.append(entry)

    if options.loglevel >= 1:
        options.stdlog.write("# predictions after input: %i\n" % ninput)
        sys.stdout.flush()

    if options.loglevel >= 10:

        options.stdlog.write(
            "############## start: predictions after input ###################################\n")
        for x in old_predictions:
            options.stdlog.write("# %s\n" % str(x))
        options.stdlog.write(
            "############## end: predictions after input #####################################\n")
        sys.stdout.flush()

    if ninput == 0:
        options.stdlog.write("# ERROR: no predictions\n")
        sys.exit(1)

    ##########################################################################
    # set up stacks of regions
    if options.conserve_memory:
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

    if options.benchmarks:
        EvaluateBenchmark(old_predictions)

    ##########################################################################
    # join regions
    if options.join_regions and options.join_regions_max_coverage:
        if options.loglevel >= 1:
            options.stdlog.write("# joining regions: maximum distance between segments = %i and maximum query coverage = %i\n" % (
                options.join_regions,
                options.join_regions_max_coverage))
            sys.stdout.flush()
        njoined = JoinRegions(old_predictions, new_predictions)
        if options.conserve_memory:
            ExchangeStreams(old_predictions, new_predictions)
        else:
            old_predictions = new_predictions
            new_predictions = []

        if options.loglevel >= 1:
            options.stdlog.write("# predictions after joining: %i\n" % njoined)
            sys.stdout.flush()

        if options.loglevel >= 10:
            options.stdlog.write(
                "############## start: predictions after joining ###################################\n")
            for x in old_predictions:
                options.stdlog.write("# %s" % str(x))
            options.stdlog.write(
                "############## end: predictions after joining #####################################\n")
            sys.stdout.flush()
    else:
        if options.loglevel >= 1:
            options.stdlog.write("# joining regions: skipped\n")
            sys.stdout.flush()

        njoined = ninput

    ##########################################################################
    # build map of best predictions
    if options.filter_suboptimal:
        if options.loglevel >= 1:
            options.stdlog.write("# calculating best predictions\n")
            sys.stdout.flush()
        best_predictions = GetBestPredictions(old_predictions)
    else:
        best_predictions = {}

    if options.loglevel >= 1:
        options.stdlog.write(
            "# calculated best predictions: %i\n" % len(best_predictions))
        sys.stdout.flush()

    ##########################################################################
    # get regions to eliminate
    filter_regions = {}
    if options.filter_regions:

        entry = PredictionParser.PredictionParserEntry(expand=0)

        filenames = options.filter_regions.split(",")

        for filename in filenames:
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# reading regions to filter from %s.\n" % (filename))
                sys.stdout.flush()

            if filename.endswith(".gz"):
                infile = gzip.open(filename, "r")
            else:
                infile = open(filename, "r")

            for line in infile:

                if line[0] == "#":
                    continue

                entry.Read(line)

                exons = Exons.Alignment2Exons(Genomics.String2Alignment(entry.mAlignmentString),
                                              entry.mQueryFrom, entry.mSbjctGenomeFrom)

                key = "%s-%s" % (entry.mSbjctToken, entry.mSbjctStrand)

                if key not in filter_regions:
                    filter_regions[key] = []

                for exon in exons:
                    filter_regions[key].append(
                        (exon.mGenomeFrom, exon.mGenomeTo))

            infile.close()

        for k in filter_regions.keys():
            filter_regions[k].sort()

    ##########################################################################
    # bipartite graph construction

    ##########################################################################
    # sort predictions by genomic region
    if options.conserve_memory:
        old_predictions.sort(
            ('mSbjctToken', 'mSbjctStrand', 'mSbjctGenomeFrom', 'mSbjctGenomeTo'))
    else:
        old_predictions.sort(lambda x, y: cmp((x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom, x.mSbjctGenomeTo),
                                              (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom, y.mSbjctGenomeTo)))

    ##########################################################################
    # filter predictions and resolve conflicts based on genomic overlap
    # deleted segments are put in a temporary storage space.
    min_from, max_from = None, None
    min_to, max_to = None, None
    region_id = 0
    noverlaps = 0
    last_prediction = None
    predictions = []
    region = Region()
    nclusters = 0
    neliminated_suboptimal = 0
    neliminated_overlap = 0

    noutput, nfiltered = 0, 0

    for this_prediction in old_predictions:

        # Filter 1: skip suboptimal predictions
        if this_prediction.mQueryToken in best_predictions:

            best_prediction = best_predictions[this_prediction.mQueryToken]

            neliminated_suboptimal += 1
            if float(this_prediction.mQueryCoverage) / best_prediction.mQueryCoverage < options.min_relative_coverage:
                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# PRUNING: reason: coverage below best: removing %s\n" % str(this_prediction))
                continue

            if float(this_prediction.score) / best_prediction.score < options.min_relative_score:
                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# PRUNING: reason: score below best: removing %s\n" % str(this_prediction))
                continue

            if float(this_prediction.mPercentIdentity) / best_prediction.mPercentIdentity < options.min_relative_percent_identity:
                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# PRUNING: reason: percent identity below best: removing %s\n" % str(this_prediction))
                continue

            neliminated_suboptimal -= 1

        # Filter 2: remove predictions overlapping with certain segments
        key = "%s-%s" % (this_prediction.mSbjctToken,
                         this_prediction.mSbjctStrand)

        if key in filter_regions:

            exons = Exons.Alignment2Exons(Genomics.String2Alignment(this_prediction.mAlignmentString),
                                          this_prediction.mQueryFrom, this_prediction.mSbjctGenomeFrom)

            if CheckOverlap(map(lambda x: (x.mGenomeFrom, x.mGenomeTo), exons), filter_regions[key]):
                if options.loglevel >= 2:
                    options.stdlog.write(
                        "# PRUNING: reason: overlapping with taboo region: removing %s\n" % str(this_prediction))
                neliminated_overlap += 1
                continue

        try:
            this_query_peptide, this_query_status, this_query_gene, this_query_transcript = \
                re.split("\s+", this_prediction.mQueryToken)
        except ValueError:
            this_query_gene = None

        # process first entry
        if min_from is None:
            min_from = this_prediction.mSbjctGenomeFrom
            max_from = this_prediction.mSbjctGenomeFrom
            max_to = this_prediction.mSbjctGenomeTo
            min_to = this_prediction.mSbjctGenomeTo
            predictions.append(this_prediction)
            last_prediction = this_prediction
            continue

        overlap = min_to > this_prediction.mSbjctGenomeFrom and \
            last_prediction.mSbjctToken == this_prediction.mSbjctToken and \
            last_prediction.mSbjctStrand == this_prediction.mSbjctStrand

        if options.loglevel >= 4:
            options.stdlog.write(
                "# from=%i, to=%i, working on: %s\n" % (min_from, max_to, str(this_prediction)))
            options.stdlog.flush()

        # resolve overlap between different genes
        if overlap:
            noverlaps += 1
        else:
            region.mSbjctToken = last_prediction.mSbjctToken
            region.mSbjctStrand = last_prediction.mSbjctStrand
            region.mSbjctGenomeFrom = min_from
            region.mSbjctGenomeTo = max_to

            region_id, nxoutput, nxfiltered = ProcessRegion(predictions,
                                                            region_id, region,
                                                            peptide_sequences,
                                                            filter_queries)

            noutput += nxoutput
            nfiltered += nxfiltered
            nclusters += 1
            predictions = []
            min_from = this_prediction.mSbjctGenomeFrom
            max_from = this_prediction.mSbjctGenomeFrom
            min_to = this_prediction.mSbjctGenomeTo
            max_to = this_prediction.mSbjctGenomeTo

        predictions.append(this_prediction)

        min_from = min(min_from, this_prediction.mSbjctGenomeFrom)
        max_from = max(max_from, this_prediction.mSbjctGenomeFrom)
        min_to = min(min_to, this_prediction.mSbjctGenomeTo)
        max_to = max(max_to, this_prediction.mSbjctGenomeTo)

        last_prediction = this_prediction

    if last_prediction:
        region.mSbjctToken = last_prediction.mSbjctToken
        region.mSbjctStrand = last_prediction.mSbjctStrand
        region.mSbjctGenomeFrom = min_from
        region.mSbjctGenomeTo = max_to

        region_id, nxoutput, nxfiltered = ProcessRegion(
            predictions, region_id, region, peptide_sequences, filter_queries)
        noutput += nxoutput
        nfiltered += nxfiltered

        nclusters += 1

    if options.conserve_memory:
        os.remove(filename_old_predictions)
        os.remove(filename_new_predictions)
        os.remove(filename_removed_predictions)

    if options.loglevel >= 1:
        options.stdlog.write("# pairs: nread=%i, input=%i, joined=%i, clusters=%i, regions=%i, eliminated_subopt=%i, eliminated_overlap=%i, noutput=%i, nfiltered=%i\n" %
                             (nread, ninput, njoined, nclusters, region_id, neliminated_suboptimal, neliminated_overlap, noutput, nfiltered))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
