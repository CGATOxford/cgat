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
gpipe/exons2genes.py - 
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

   python gpipe/exons2genes.py --help

Type::

   python gpipe/exons2genes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import getopt
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import CGAT.PredictionParser as PredictionParser
import pgdb

USAGE = """python %s < predictions > genes

Version: $Id: gpipe/exons2genes.py 1799 2008-03-28 11:44:19Z andreas $

collect overlapping exons into genes.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-t, --filter-tokens             filter by query token
-p, --filter-predictions        filter by prediction_id
-P, --table-predictions         table with predictions
-Q, --table-quality             table with quality information
-G, --genes-tsv-file                     table with genic regions to be skipped
""" % sys.argv[0]

param_long_options = ["verbose=", "help", "overlap=", "filter-tokens=", "filter-predictions=",
                      "table-predictions=",
                      "table-quality=", "table-genes=", "quality=",
                      "version"]
param_short_options = "v:ho:t:p:P:T:G:Q:"

param_loglevel = 1

param_overlap_residues = 90
param_filename_filter_tokens = None
param_filename_filter_predictions = None
param_tablename_predictions = None
param_tablename_quality = None
param_tablename_genes = None
param_connection = "fgu202:andreas"

param_quality = None


def ResolveExonOverlaps(gene_id, predictions):
    """resolve overlaps between predictions based
    on exonic overlap."""

    all_exons = []
    n = 1
    for p in predictions:
        exons = Exons.Alignment2ExonBoundaries(Genomics.String2Alignment(p.mAlignmentString),
                                               query_from=0,
                                               sbjct_from=p.mSbjctGenomeFrom)

        for exon in exons:
            all_exons.append((exon.mGenomeFrom, exon.mGenomeTo, n))
        n += 1

    map_prediction2gene = range(0, len(predictions) + 1)
    map_gene2predictions = [None]
    for x in range(1, len(predictions) + 1):
        map_gene2predictions.append([x])

    all_exons.sort()
    # print all_exons

    # cluster exons by overlap
    last_exon_from, last_exon_to, last_p = all_exons[0]

    for exon_from, exon_to, p in all_exons[1:]:
        # if overlap
        if min(exon_to, last_exon_to) - max(exon_from, last_exon_from) > 0:
            # print "# overlap between %i and %i" % (p, last_p)
            # rewire pointers to point to gene of previous prediction
            # if they belong to different genes
            new_g = map_prediction2gene[last_p]
            old_g = map_prediction2gene[p]

            if new_g != old_g:
                for x in map_gene2predictions[old_g]:
                    map_gene2predictions[new_g].append(x)
                    map_prediction2gene[x] = new_g
                map_gene2predictions[old_g] = []

        # if no overlap: create new gene, if predictions has no gene
        # associated with it yet.
        else:
            # print "# no overlap between %i and %i" % (p, last_p)
            if not map_prediction2gene[p]:
                map_prediction2gene[p] = len(map_gene2predictions)
                map_gene2predictions.append([p])

        last_exon_to = max(last_exon_to, exon_to)
        last_p = p

    for x in range(1, len(map_gene2predictions)):
        if map_gene2predictions[x]:
            for p in map_gene2predictions[x]:
                print "%i\t%i" % (gene_id, predictions[p - 1].mPredictionId)
            gene_id += 1

    return gene_id

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
        elif o in ("-o", "--overlap"):
            param_overlap_residues = int(a)
        elif o in ("-t", "--filter-tokens"):
            param_filename_filter_tokens = a
        elif o in ("-p", "--filter-predictions"):
            param_filename_filter_predictions = a
        elif o in ("-P", "--table-predictions"):
            param_tablename_predictions = a
        elif o in ("-Q", "--table-quality"):
            param_tablename_quality = a
        elif o in ("-G", "--table-genes"):
            param_tablename_genes = a
        elif o in ("-q", "--quality="):
            param_quality = string.split(a, ",")

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()
    sys.stdout.flush()

    dbhandle = pgdb.connect(param_connection)

    filter_tokens = {}
    if param_filename_filter_tokens:
        infile = open(param_filename_filter_tokens, "r")
        for line in infile:
            if line[0] == "#":
                continue
            filter_tokens[line[:-1]] = 0
        if param_loglevel >= 1:
            print "# filtering for %i tokens" % len(filter_tokens)

    filter_predictions = {}
    if param_filename_filter_predictions:
        infile = open(param_filename_filter_predictions, "r")
        for line in infile:
            if line[0] == "#":
                continue
            filter_predictions[int(line[:-1])] = 0
        if param_loglevel >= 1:
            print "# filtering for %i predictions" % len(filter_predictions)

    # array with predictions
    predictions = []

    nfiltered_region = 0

    if param_tablename_genes:
        statement = """SELECT p.sbjct_token, p.sbjct_strand,
        MIN(p.sbjct_genome_from), MAX(p.sbjct_genome_to)
        FROM %s AS p, %s AS g
        WHERE p.prediction_id = g.prediction_id
        GROUP BY g.gene_id, p.sbjct_token, p.sbjct_strand
        """ % (param_tablename_predictions, param_tablename_genes)

        c = dbhandle.cursor()
        c.execute(statement)

        filter_genes = {}
        nregions = 0
        for line in c.fetchall():
            sbjct_token, sbjct_strand, sbjct_from, sbjct_to = line
            key = "%s%s" % (sbjct_token, sbjct_strand)
            if not filter_genes.has_key(key):
                filter_genes[key] = []
            filter_genes[key].append((sbjct_from, sbjct_to))
            nregions += 1

        for x in filter_genes.keys():
            filter_genes[x].sort()

        if param_loglevel >= 1:
            print "# filtering for %i regions" % nregions

    if param_tablename_predictions:

        if param_quality:
            statement = """SELECT p.* FROM %s AS p, %s AS q
            WHERE p.prediction_id = q.prediction_id AND
            q.class IN ('%s')""" % (param_tablename_predictions, param_tablename_quality, string.join(param_quality, "','"))
        else:
            statement = "SELECT p.* FROM %s AS p" % (
                param_tablename_predictions)

        c = dbhandle.cursor()
        c.execute(statement)

        for line in c.fetchall():
            entry = PredictionParser.PredictionParserEntry()

            entry.FillFromTable(line)

            if param_filename_filter_tokens:
                if filter_tokens.has_key(entry.mQueryToken):
                    filter_tokens[entry.mQueryToken] += 1
                else:
                    continue

            if param_filename_filter_predictions:
                if filter_predictions.has_key(entry.mPredictionId):
                    filter_predictions[entry.mPredictionId] += 1
                else:
                    continue

            if param_tablename_genes:
                key = "%s%s" % (entry.mSbjctToken, entry.mSbjctStrand)
                if filter_genes.has_key(key):
                    stop = 0
                    for sbjct_from, sbjct_to in filter_genes[key]:
                        if min(sbjct_to, entry.mSbjctGenomeTo) - max(sbjct_from, entry.mSbjctFrom) > 0:
                            stop = 1
                            nfiltered_region += 1
                            break
                    if stop:
                        continue

            predictions.append(entry)
    else:
        for line in sys.stdin:

            if line[0] == "#":
                continue

            entry = PredictionParser.PredictionParserEntry(expand=0)

            entry.Read(line)

            if param_filename_filter_tokens:
                if filter_tokens.has_key(entry.mQueryToken):
                    filter_tokens[entry.mQueryToken] += 1
                else:
                    continue

            if param_filename_filter_predictions:
                if filter_predictions.has_key(entry.mPredictionId):
                    filter_predictions[entry.mPredictionId] += 1
                else:
                    continue

            if param_tablename_genes:
                key = "%s%s" % (entry.mSbjctToken, entry.mSbjctStrand)
                if filter_genes.has_key(key):
                    stop = 0
                    for sbjct_from, sbjct_to in filter_genes[key]:
                        if min(sbjct_to, entry.mSbjctGenomeTo) - max(sbjct_from, entry.mSbjctFrom) > 0:
                            stop = 1
                            nfiltered_region += 1
                            break
                    if stop:
                        continue

            predictions.append(entry)

    if param_loglevel >= 1:
        print "# number of predictions read: %i" % len(predictions)

    ##########################################################################
    # sort predictions by genomic region
    predictions.sort(lambda x, y: cmp((x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom),
                                      (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeFrom)))

    sbjct_token = None
    sbjct_strand = None
    last_to = 0

    gene_id = 1
    collection = []
    nassigned = 0

    for p in predictions:

        if sbjct_token != p.mSbjctToken or\
           sbjct_strand != p.mSbjctStrand or\
           p.mSbjctGenomeFrom - (last_to - param_overlap_residues) > 0:
            if collection:
                gene_id = ResolveExonOverlaps(gene_id, collection)
            collection = []
            sbjct_token = p.mSbjctToken
            sbjct_strand = p.mSbjctStrand
            last_to = 0

        nassigned += 1
        collection.append(p)
        last_to = max(last_to, p.mSbjctGenomeTo)

    gene_id = ResolveExonOverlaps(gene_id, collection)

    total = 0

    if param_loglevel >= 1:
        print "# matched tokens:"
        for x in filter_tokens.keys():
            print "# %s\t%i" % (x, filter_tokens[x])
            total += filter_tokens[x]

    print "# genes=%i, assignments=%i, total=%i, filtered_region=%i" % (gene_id - 1, nassigned, total, nfiltered_region)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
