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
gpipe/predictions2assembly.py - 
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

   python gpipe/predictions2assembly.py --help

Type::

   python gpipe/predictions2assembly.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import getopt
import tempfile
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.PredictionParser as PredictionParser
import CGAT.PredictionFile as PredictionFile


USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/predictions2assembly.py 698 2006-07-19 15:53:22Z andreas $

Compile collinear predictions that are spread accross several contigs into new predictions.
Old predictions are translated into the newly created contigs.

If --join-pattern is given, joined contigs are created..
If --joined-pattern contains a %%s, a new file is create for the assembled contig, otherwise
all are written to a single file.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-g, --genome-file=           pattern for filenames with the genomic DNA (FASTA).
-j, --join-pattern=             write joined contigs with pattern
-f, --format=                   input format [predictions]
-i, --max-intron                maximum intron length
-d, --max-difference            maximum difference between peptide and genomic gap
-c, --contigs-tsv-file=                  filename with contig sizes
-o, --max-overlap=              maximum overlap
-s, --filename-sizes=           filename with sizes
""" % sys.argv[0]

param_loglevel = 1

# maximum intron size
param_max_intron = 50000

param_format = "predictions"

param_long_options = ["verbose=", "help",
                      "format=", "max-intron=",
                      "join-pattern=", "genome-file=",
                      "max-overlap=",
                      "max-difference=", "contigs=",
                      "filename-sizes=", "version"]

param_short_options = "v:hf:i:d:c:o:jg:s:"

param_max_difference = 10

# relative permissive overlap
param_max_relative_overlap = 50

# absolute permissive overlap
param_max_overlap = 0

param_conserve_frame = 0

param_filename_contigs = "contig_sizes"

param_filename_join_pattern = None

# pattern for genomes, %s is substituted for the sbjct_token
param_genome_file = "genome_%s.fasta"

param_separator_contigs = "-"

param_filename_sizes = None

global_last_filename_genome = None
global_forward_sequences = None
global_reverse_sequences = None


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
                    max(segments[x].mSbjctGenomeFrom, segments[y].mSbjctGenomeFrom)) > 0:
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

                new_entries += ProcessSegments(segments)
                segments = []

            segments.append(entry)
            last_entry = entry

        new_entries += ProcessSegments(segments)

    if param_loglevel >= 2:
        print "# number of predictions: %i" % len(new_entries)

    return new_entries


class BoundaryPredictions:

    mPositiveMax = 0
    mPositiveMin = 1000000000
    mNegativeMax = 0
    mNegativeMin = 1000000000
    mPositiveMinPredictions = []
    mPositiveMaxPredictions = []
    mNegativeMinPredictions = []
    mNegativeMaxPredictions = []

    def __init__(self):
        pass

    def update(self, prediction):

        if prediction.mSbjctStrand == "+":
            if prediction.mSbjctGenomeFrom < self.mPositiveMin:
                self.mPositiveMinPredictions = filter(
                    lambda x: x.mSbjctGenomeFrom <= prediction.mSbjctGenomeFrom +
                    param_max_difference,
                    self.mPositiveMinPredictions)
                self.mPositiveMinPredictions.append(prediction)
                self.mPositiveMin = prediction.mSbjctGenomeFrom
            elif prediction.mSbjctGenomeFrom <= self.mPositiveMin + param_max_difference:
                self.mPositiveMinPredictions.append(prediction)

            if prediction.mSbjctGenomeTo > self.mPositiveMax:
                self.mPositiveMaxPredictions = filter(
                    lambda x: x.mSbjctGenomeTo >= prediction.mSbjctGenomeFrom -
                    param_max_difference,
                    self.mPositiveMaxPredictions)
                self.mPositiveMaxPredictions.append(prediction)
                self.mPositiveMax = prediction.mSbjctGenomeTo
            elif prediction.mSbjctGenomeTo >= self.mPositiveMax - param_max_difference:
                self.mPositiveMaxPredictions.append(prediction)

        else:
            if prediction.mSbjctGenomeFrom < self.mNegativeMin:
                self.mNegativeMinPredictions = filter(
                    lambda x: x.mSbjctGenomeFrom <= prediction.mSbjctGenomeFrom +
                    param_max_difference,
                    self.mNegativeMinPredictions)
                self.mNegativeMinPredictions.append(prediction)
                self.mNegativeMin = prediction.mSbjctGenomeFrom
            elif prediction.mSbjctGenomeFrom <= self.mNegativeMin + param_max_difference:
                self.mNegativeMinPredictions.append(prediction)

            if prediction.mSbjctGenomeTo > self.mNegativeMax:
                self.mNegativeMaxPredictions = filter(
                    lambda x: x.mSbjctGenomeTo >= prediction.mSbjctGenomeFrom -
                    param_max_difference,
                    self.mNegativeMaxPredictions)
                self.mNegativeMaxPredictions.append(prediction)
                self.mNegativeMax = prediction.mSbjctGenomeTo
            elif prediction.mSbjctGenomeTo >= self.mNegativeMax - param_max_difference:
                self.mNegativeMaxPredictions.append(prediction)

    def __str__(self):

        return "#" + string.join( map( str, (self.mPositiveMin, self.mPositiveMax, self.mNegativeMin, self.mNegativeMax))) + "\n" +\
               "# min-positive\n" + string.join( map( str, self.mPositiveMinPredictions), "\n") + "\n" +\
               "# max-positive\n" + string.join( map( str, self.mPositiveMaxPredictions), "\n") + "\n" +\
               "# min-negative\n" + string.join( map( str, self.mNegativeMinPredictions), "\n") + "\n" +\
               "# max-negative\n" + \
            string.join(map(str, self.mNegativeMaxPredictions), "\n") + "\n"


def CheckOverlap(l1, l2):
    """check if there are at least two predictions that are collinear."""

    results = []
    for p1 in l1:
        for p2 in l2:
            if p1.mQueryToken != p2.mQueryToken:
                continue
            overlap = min(p1.mQueryTo, p2.mQueryTo) - \
                max(p1.mQueryFrom, p2.mQueryFrom)
            if 100 * overlap / (p1.mQueryTo - p1.mQueryFrom + 1) >= param_max_relative_overlap or \
                    100 * overlap / (p2.mQueryTo - p2.mQueryFrom + 1) >= param_max_relative_overlap:
                continue

            if p1.mQueryTo < p2.mQueryFrom + param_max_overlap or \
               p1.mQueryFrom > p2.mQueryTo - param_max_overlap:
                results.append((p1, p2))

    return results


def CheckCollinearity(c1, c2):
    """check if there are at least two predictions that are collinear.

    Check in all 16 combinations.
    """

    results = []
    results += CheckOverlap(c1.mPositiveMinPredictions,
                            c2.mPositiveMinPredictions)
    results += CheckOverlap(c1.mPositiveMinPredictions,
                            c2.mPositiveMaxPredictions)
    results += CheckOverlap(c1.mPositiveMinPredictions,
                            c2.mNegativeMinPredictions)
    results += CheckOverlap(c1.mPositiveMinPredictions,
                            c2.mNegativeMaxPredictions)

    results += CheckOverlap(c1.mPositiveMaxPredictions,
                            c2.mPositiveMinPredictions)
    results += CheckOverlap(c1.mPositiveMaxPredictions,
                            c2.mPositiveMaxPredictions)
    results += CheckOverlap(c1.mPositiveMaxPredictions,
                            c2.mNegativeMinPredictions)
    results += CheckOverlap(c1.mPositiveMaxPredictions,
                            c2.mNegativeMaxPredictions)

    results += CheckOverlap(c1.mNegativeMinPredictions,
                            c2.mPositiveMinPredictions)
    results += CheckOverlap(c1.mNegativeMinPredictions,
                            c2.mPositiveMaxPredictions)
    results += CheckOverlap(c1.mNegativeMinPredictions,
                            c2.mNegativeMinPredictions)
    results += CheckOverlap(c1.mNegativeMinPredictions,
                            c2.mNegativeMaxPredictions)

    results += CheckOverlap(c1.mNegativeMaxPredictions,
                            c2.mPositiveMinPredictions)
    results += CheckOverlap(c1.mNegativeMaxPredictions,
                            c2.mPositiveMaxPredictions)
    results += CheckOverlap(c1.mNegativeMaxPredictions,
                            c2.mNegativeMinPredictions)
    results += CheckOverlap(c1.mNegativeMaxPredictions,
                            c2.mNegativeMaxPredictions)

    return results


def GetContig(prediction):
    """get contig sequence for prediction."""

    global global_last_filename_genome
    global global_forward_sequences
    global global_reverse_sequences

    if "%s" in param_genome_file:
        filename_genome = param_genome_file % prediction.mSbjctToken
    else:
        filename_genome = param_genome_file

    if global_last_filename_genome != filename_genome:
        if param_loglevel >= 2:
            print "# reading genome %s" % filename_genome

        try:
            global_forward_sequences, global_reverse_sequences = Genomics.ReadGenomicSequences(
                open(filename_genome, "r"))
        except IOError:
            raise "# WARNING: genome %s not found" % filename_genome

        global_last_filename_genome = filename_genome

    if prediction.mSbjctStrand == "+":
        return (prediction.mSbjctToken, global_forward_sequences[prediction.mSbjctToken], False)
    else:
        return (prediction.mSbjctToken, global_reverse_sequences[prediction.mSbjctToken], True)

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
        elif o in ("-f", "--format"):
            param_format = a
        elif o in ("-i", "--max-intron"):
            param_max_intron = int(a)
        elif o in ("-d", "--max-difference"):
            param_max_difference = int(a)
        elif o in ("-o", "--max-overlap"):
            param_max_overlap = int(a)
        elif o in ("-c", "--contigs-tsv-file"):
            param_filename_contigs = a
        elif o in ("-g", "--genome-file"):
            param_genome_file = a
        elif o in ("-j", "--join-pattern"):
            param_filename_join_pattern = a
        elif o in ("-s", "--filename-sizes"):
            param_filename_sizes = a

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()

    ninput = 0
    max_id = 0

    contig_sizes = Genomics.ReadContigSizes(open(param_filename_contigs, "r"))

    ##########################################################################
    # reading predictions
    contig = {}

    tmp_predictions, filename_tmp_predictions = tempfile.mkstemp()
    os.close(tmp_predictions)
    tmp_predictions = PredictionFile.PredictionFile()
    tmp_predictions.open(filename_tmp_predictions, "w")

    if param_format == "predictions":

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

            if entry.mSbjctToken not in contig:
                contig[entry.mSbjctToken] = BoundaryPredictions()

            contig[entry.mSbjctToken].update(entry)

            tmp_predictions.append(entry)

    if param_loglevel >= 4:
        for c in contig.keys():
            print "######start of %s #####################################################" % c
            print "#", str(contig[c])
            print "######end of %s #####################################################" % c

    tmp_predictions.close()

    max_id += 1
    first_pseudo_id = max_id

    cc = contig.keys()

    ##########################################################################
    # get pairs of colinear predictions on different contigs
    results = []

    if param_loglevel >= 1:
        print "# finished parsing %i contigs" % len(cc)
        sys.stdout.flush()

    for c1 in range(len(cc) - 1):

        if param_loglevel >= 1:
            print "# processing: %i/%i" % (c1 + 1, len(cc))
            sys.stdout.flush()

        for c2 in range(c1 + 1, len(cc)):

            r = CheckCollinearity(contig[cc[c1]], contig[cc[c2]])

            if r and param_loglevel >= 3:
                print "# --------------------------------------------------------"
                print "# %s and %s are collinear" % (cc[c1], cc[c2])
                for r1, r2 in r:
                    print "# ----------------------"
                    print "#", str(r1)
                    print "#", str(r2)

            results += r

    ##########################################################################
    # cluster co-linear predictions on different contigs by sbjct_token
    queries = {}

    for r1, r2 in results:
        if r1.mQueryToken not in queries:
            queries[r1.mQueryToken] = {}
        queries[r1.mQueryToken][r1.mPredictionId] = r1
        queries[r1.mQueryToken][r2.mPredictionId] = r2

    nnew = 0
    ncancelled = 0

    # set of contigs joined
    map_contig2new = {}

    # names of new contigs
    new_contigs = {}

    # remove old contig file, if it already exists.
    if param_filename_join_pattern and "%s" not in param_filename_join_pattern:
        if os.path.exists(param_filename_join_pattern):
            os.remove(param_filename_join_pattern)

    if param_filename_sizes:
        outfile_sizes = open(param_filename_sizes, "w")
    else:
        outfile_sizes = None

    ##########################################################################
    # join contigs
    for q in queries.keys():

        s = queries[q].values()

        s.sort(
            lambda x, y: cmp((x.mQueryFrom, x.mQueryTo), (y.mQueryFrom, y.mQueryTo)))

        if param_loglevel >= 2:
            print "# -----------------------------------------------"
            print "# predictions to be joined for query=", q
            for p in s:
                print "#", str(p)
            print "# -----------------------------------------------"

        new_prediction = s[0].GetCopy()
        last_contig_size = contig_sizes[new_prediction.mSbjctToken]
        do_cancel = False

        contigs = []

        if param_filename_join_pattern:
            contigs.append(GetContig(new_prediction))

        for p in s[1:]:
            overlap = new_prediction.mQueryTo - p.mQueryFrom + 1
            if overlap > 0:
                if overlap > param_max_overlap or \
                   100 * (p.mQueryTo - p.mQueryFrom + 1) / overlap > param_max_relative_overlap:
                    print "# dodgy prediction sequence (overlap = %i), joining of contigs cancelled." % overlap
                    sys.stdout.flush()
                    do_cancel = True
                    break

            if param_filename_join_pattern:
                contigs.append(GetContig(p))

            new_prediction.Add(p,
                               combine_contig=True,
                               allow_overlap=True,
                               contig_size=last_contig_size)

            last_contig_size += contig_sizes[p.mSbjctToken]

        if do_cancel:
            ncancelled += 1
            continue

        nnew += 1
        new_prediction.mPredictionId = max_id
        new_prediction.mSbjctStrand = "+"
        max_id += 1

        print "# joining\t" + string.join(map(lambda x: x.mSbjctToken + x.mSbjctStrand, s), "\t")

        if param_filename_join_pattern and len(contigs) > 0:

            new_contig = string.join(
                map(lambda x: x[0], contigs), param_separator_contigs)

            # do not write the same contig twice
            if new_contig not in new_contigs:

                new_contigs[new_contig] = 1
                lcontig = len(string.join(map(lambda x: x[1], contigs), ""))

                # check if contig already part of a different joined contig
                l = 0
                for id, sequence, switch in contigs:
                    if id in map_contig2new:
                        print "# WARNING: contig %s already joined" % id
                    map_contig2new[id] = (
                        new_contig, switch, l, lcontig - l - len(sequence))
                    l += len(sequence)

                # write new contig
                if "%s" in param_filename_join_pattern:
                    filename_genome = param_filename_join_pattern % new_contig
                    outfile = open(filename_genome, "w")
                else:
                    filename_genome = param_filename_join_pattern
                    outfile = open(filename_genome, "a")

                if outfile_sizes:
                    outfile_sizes.write("%s\t%i\t0\n" % (new_contig, lcontig))

                outfile.write(
                    ">" + new_contig + "\n" + string.join(map(lambda x: x[1], contigs), "") + "\n")
                outfile.close()

        print str(new_prediction)

    if outfile_sizes:
        outfile_sizes.close()

    ##########################################################################
    # move other predictions into the new contigs by translating their
    # coordinates
    tmp_predictions.open(mode="r")
    noutput = 0
    ntranslated = 0
    for p in tmp_predictions:

        if p.mSbjctToken in map_contig2new:

            p.mSbjctToken, switch, offset_pos, offset_neg = map_contig2new[
                p.mSbjctToken]

            if (switch and p.mSbjctStrand == "+") or \
               (not switch and p.mSbjctStrand == "-"):
                offset = offset_neg
            else:
                offset = offset_pos

            # change strand for inverted contigs
            if switch:
                if p.mSbjctStrand == "+":
                    p.mSbjctStrand = "-"
                else:
                    p.mSbjctStrand = "+"

            p.mSbjctGenomeFrom += offset
            p.mSbjctGenomeTo += offset
            ntranslated += 1

        noutput += 1

        print str(p)

    if param_loglevel >= 1:
        print "## nread=%i, nnew=%i, noutput=%i, ntranslated=%i, first_id=%i" %\
              (ninput, nnew, noutput, ntranslated, first_pseudo_id)
        print "# ncontigs=%i, npairs=%i, nqueries=%i, nnew=%i, njoined=%i, ncancelled=%i" %\
              (len(contig), len(results), len(queries),
               len(new_contigs), len(map_contig2new), ncancelled)

    os.remove(filename_tmp_predictions)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
