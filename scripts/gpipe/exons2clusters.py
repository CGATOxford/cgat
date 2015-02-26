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
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Exons as Exons
import CGAT.Genomics as Genomics
import alignlib_lite

USAGE = """python %s [OPTIONS] schema1 schema2 [...]

Version: $Id: gpipe/exonerate_combine_regions.py 2464 2009-02-02 10:28:09Z andreas $

Cluster exons into genes. Identical overlap is required.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --peptides-fasta-file=                 file with protein sequences (FASTA).
-o, --overlap=                  use overlap of # to collect clusters
-a, --alignment=                check alignment with minimum alignment score of #
-r, --regex-preferred=          regular expression for preferred gene identifiers.
-c, --contigs-tsv-file=                  filename with contig sizes
-g, --use-genome-length         use genome length for deciding which rep to use.
""" % sys.argv[0]

param_loglevel = 1

param_filename_peptides = None

param_long_options = ["verbose=", "help",
                      "peptides=", "overlap=", "alignment=", "regex-preferred=", "contigs=",
                      "use-genome-length", "genome-file=",
                      "version"]

param_short_options = "v:hp:o:a:r:c:gf:"

param_min_overlap = 0

# minimum minimum overlap
param_min_min_overlap = 0.1

# minimum maximum overlap
param_min_max_overlap = 0.5

param_min_alignment_score = 0.0
param_regex_preferred = None

param_min_terminal_exon_coverage = 0.8

# contig sizes
param_filename_contigs = None

# use genome length
param_use_genome_length = True

param_genome_file = None


def PrintCluster(cluster,
                 cluster_id,
                 lengths,
                 peptide_sequences=None,
                 regex_preferred=None):
    """print a cluster.

    Take longest sequence as representative. If preferred is given, only take
    genes matching preferred identifier.
    """

    if regex_preferred:
        rx = re.compile(regex_preferred)
    else:
        rx = None

    max_al = 0
    max_pl = 0
    rep_a = None
    rep_p = None
    for c in cluster:
        l = 0
        if c in lengths:
            l = lengths[c]

        if l > max_al:
            max_al = l
            rep_a = c

        if rx and rx.search(c) and l > max_pl:
            max_pl = l
            rep_p = c

    if max_pl > 0:
        max_l = max_pl
        rep = rep_p
    else:
        max_l = max_al
        rep = rep_a

    for mem in cluster:
        l = 0
        if mem in lengths:
            l = lengths[mem]
        if peptide_sequences:
            map_rep2mem = alignlib_lite.makeAlignmentVector()

            if rep == mem and rep in lengths:
                alignlib_lite.addDiagonal2Alignment(
                    map_rep2mem, 1, lengths[rep], 0)
            elif mem in peptide_sequences and \
                    rep in peptide_sequences:
                alignator = alignlib_lite.makeAlignatorDPFull(
                    alignlib_lite.ALIGNMENT_LOCAL, -10.0, -1.0)
                alignator.align(map_rep2mem,
                                alignlib_lite.makeSequence(
                                    peptide_sequences[rep]),
                                alignlib_lite.makeSequence(peptide_sequences[mem]))

            f = alignlib_lite.AlignmentFormatEmissions(map_rep2mem)
            print string.join(map(str, (rep, mem, l, f)), "\t")

        else:
            print string.join(map(str, (rep, mem, l)), "\t")

    sys.stdout.flush()

    return cluster_id


def CollectCluster(map_transcript2transcript, seed):

    if seed not in map_transcript2transcript:
        return []

    cluster = [seed]

    tt = map_transcript2transcript[seed]
    del map_transcript2transcript[seed]

    for t in tt:
        if t in map_transcript2transcript:
            cluster += CollectCluster(map_transcript2transcript, t)

    return cluster


def ClusterByExonCorrespondence(lengths={}, peptide_sequences=None):

    exons = Exons.ReadExonBoundaries(sys.stdin)
    if param_loglevel >= 1:
        print "# read exons for %i transcripts" % len(exons)

    if not lengths:
        for k in exons:
            lengths[k] = (exons[k][0].mPeptideTo / 3) + 1
            for e in exons[k][1:]:
                lengths[k] = max(lengths[k], (e.mPeptideTo / 3) + 1)

        if param_loglevel >= 1:
            print "# lengths for %i transcripts" % len(lengths)

    map_region2transcript = {}
    map_transcript2region = {}
    map_transcript2transcript = {}
    # build map of regions to transcripts
    for t in exons:
        map_transcript2region[t] = []
        for e in exons[t]:
            r = "%s-%s-%i-%i" % (e.mSbjctToken,
                                 e.mSbjctStrand, e.mGenomeFrom, e.mGenomeTo)
            if r not in map_region2transcript:
                map_region2transcript[r] = []
            map_region2transcript[r].append(t)
            map_transcript2region[t].append(r)

    # build map of transcript to transcript
    map_transcript2transcript = {}

    for t in map_transcript2region:
        map_transcript2transcript[t] = []
        for r in map_transcript2region[t]:
            for tt in map_region2transcript[r]:
                map_transcript2transcript[t].append(tt)

    for t in map_transcript2transcript:
        map_transcript2transcript[t].sort()
        l = None
        n = []
        for tt in map_transcript2transcript[t]:
            if t == tt:
                continue
            if l != tt:
                n.append(tt)
            l = tt
        map_transcript2transcript[t] = n

    # cluster greedily, take longest transcript
    cluster_id = 1
    for t in map_transcript2region:
        if t not in map_transcript2transcript:
            continue
        cluster = CollectCluster(map_transcript2transcript, t)
        PrintCluster(cluster, cluster_id, lengths, peptide_sequences,
                     param_regex_preferred)
        cluster_id += 1

    if param_loglevel >= 1:
        print "# RESULT: %i transcripts in %i genes" % (len(map_transcript2region), cluster_id - 1)


def CheckAlignments(peptide_sequences, query_token, other_tokens):
    """check wether query aligns to all others.
    """

    if param_loglevel >= 3:
        print "# checking query %s and sbjcts %s" % (query_token, str(other_tokens))
        sys.stdout.flush()

    if query_token not in peptide_sequences:
        return True

    result = alignlib_lite.makeAlignmentVector()
    alignator = alignlib_lite.makeAlignatorDPFull(alignlib_lite.ALIGNMENT_LOCAL,
                                                  -10.0, -1.0)
    row_seq = alignlib_lite.makeSequence(peptide_sequences[query_token])

    for x in other_tokens:
        if x not in peptide_sequences:
            continue
        col_seq = alignlib_lite.makeSequence(peptide_sequences[x])
        alignator.align(result, row_seq, col_seq)
        if param_loglevel >= 5:
            print "# %s - %s = %f" % (query_token, x, result.getScore())
        if result.getScore() > param_min_alignment_score:
            return True

    return False

# ------------------------------------------------------------------------


def ClusterByExonOverlap(exons, lengths={}, peptide_sequences={}, loglevel=1):
    """cluster transcripts by exon overlap.
    """

    map_transcript2cluster = {}
    map_cluster2transcripts = {}
    list_of_exons = []

    for k, ee in exons.items():
        if k not in map_transcript2cluster:
            map_transcript2cluster[k] = k
            map_cluster2transcripts[k] = [k, ]

        for e in ee:
            list_of_exons.append(e)

    list_of_exons.sort(lambda x, y: cmp((x.mSbjctToken, x.mSbjctStrand, x.mGenomeFrom),
                                        (y.mSbjctToken, y.mSbjctStrand, y.mGenomeFrom)))

    if loglevel >= 1:
        print "# sorted %i list_of_exons" % len(list_of_exons)

    last_id = list_of_exons[0].mQueryToken
    last_from = list_of_exons[0].mGenomeFrom
    last_to = list_of_exons[0].mGenomeTo
    last_token = list_of_exons[0].mSbjctToken
    last_strand = list_of_exons[0].mSbjctStrand

    for e in list_of_exons[1:]:

        if loglevel >= 3:
            print "# processing exon:", str(e)

        if last_token != e.mSbjctToken:

            if loglevel >= 2:
                print "# switching to strand %s" % e.mSbjctToken
                sys.stdout.flush()

        o = min(e.mGenomeTo, last_to) - max(e.mGenomeFrom, last_from)
        o1 = float(o) / (last_to - last_from)
        o2 = float(o) / (e.mGenomeTo - e.mGenomeFrom)

        if last_strand == e.mSbjctStrand and \
           last_token == e.mSbjctToken and \
           o >= param_min_overlap and \
           min(o1, o2) >= param_min_min_overlap and \
           max(o1, o2) >= param_min_max_overlap:

            last_cluster = map_transcript2cluster[last_id]

            if loglevel >= 4:
                print "# overlap %i(%5.2f/%5.2f) between %s and %s: %i-%i with %i-%i" % \
                      (o, o1, o2,
                       last_id, e.mQueryToken, last_from, last_to, e.mGenomeFrom, e.mGenomeTo)
                sys.stdout.flush()

            if e.mQueryToken not in map_cluster2transcripts[last_cluster]:
                join = True
                if param_min_alignment_score > 0 and \
                        not CheckAlignments(peptide_sequences, e.mQueryToken, map_cluster2transcripts[last_cluster]):
                    join = False

                if join:
                    map_cluster2transcripts[
                        last_cluster] += map_cluster2transcripts[e.mQueryToken]
                    for x in map_cluster2transcripts[e.mQueryToken]:
                        map_transcript2cluster[x] = last_cluster
                    map_cluster2transcripts[e.mQueryToken] = []

            last_from = min(e.mGenomeFrom, last_from)
            last_to = max(e.mGenomeTo, last_to)
        else:
            last_id = e.mQueryToken
            last_strand = e.mSbjctStrand
            last_token = e.mSbjctToken
            last_from = e.mGenomeFrom
            last_to = e.mGenomeTo

    return map_cluster2transcripts, map_transcript2cluster


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
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-o", "--overlap"):
            param_min_overlap = int(a)
        elif o in ("-a", "--alignment"):
            param_min_alignment_score = float(a)
        elif o in ("-r", "--regex-preferred"):
            param_regex_preferred = a
        elif o in ("-c", "--contigs-tsv-file"):
            param_filename_contigs = a
        elif o in ("-g", "--use-genome-length"):
            param_use_genome_length = True
        elif o in ("-f", "--genome-file"):
            param_genome_file = a

    print E.GetHeader()
    print E.GetParams()

    # read peptide sequences
    if param_filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(param_filename_peptides, "r"))
    else:
        peptide_sequences = {}

    if param_genome_file:
        # read from fasta file
        fasta = IndexedFasta.IndexedFasta(param_genome_file)
        contig_sizes = fasta.getContigSizes()
        delete_missing = True
    elif param_filename_contigs:
        # read contigs
        contig_sizes = Genomics.ReadContigSizes(
            open(param_filename_contigs, "r"))
        delete_missing = True
    else:
        contig_sizes = {"dummy": 1000000000}
        delete_missing = False

    if param_loglevel >= 1:
        print "# read %i peptide sequences" % len(peptide_sequences)
        sys.stdout.flush()

    exons = Exons.ReadExonBoundaries(sys.stdin,
                                     contig_sizes=contig_sizes,
                                     delete_missing=delete_missing,
                                     )

    if param_loglevel >= 1:
        print "# read exon information for %i transcripts" % len(exons)
        sys.stdout.flush()

    if len(exons) == 0:
        raise IOError("no exons in exon list.")

    Exons.SetRankToPositionFlag(exons)

    if param_use_genome_length:
        lengths = Exons.GetGenomeLengths(exons)
    else:
        lengths = Exons.GetPeptideLengths(exons)

    if param_min_overlap > 0:
        map_cluster2transcripts, map_transcript2cluster = ClusterByExonOverlap(exons,
                                                                               lengths,
                                                                               peptide_sequences,
                                                                               loglevel=param_loglevel)
    else:
        map_cluster2transcripts, map_transcript2cluster = \
            Exons.ClusterByExonIdentity(exons,
                                        max_terminal_num_exons=3,
                                        min_terminal_exon_coverage=param_min_terminal_exon_coverage,
                                        loglevel=param_loglevel)

    map_transcript2strand = {}
    for k, ee in exons.items():
        map_transcript2strand[k] = (ee[0].mSbjctStrand == "+")

    nnegatives, npositives = 0, 0
    # take longest transcript
    cluster_id = 1

    if peptide_sequences:
        print "rep\tmem\tlength\tquery_from\tquery_to\tquery_ali\tsbjct_from\tsbjct_to\tsbjct_ali"
    else:
        print "rep\tmem\tlength"

    for t in map_cluster2transcripts:
        if not map_cluster2transcripts[t]:
            continue
        PrintCluster(map_cluster2transcripts[t],
                     cluster_id,
                     lengths, peptide_sequences,
                     param_regex_preferred)
        cluster_id += 1

        if map_transcript2strand[t]:
            npositives += 1
        else:
            nnegatives += 1

    if param_loglevel >= 1:
        t = npositives + nnegatives
        print "# npositives=%i (%5.2f), nnegatives=%i (%5.2f)" % (npositives, npositives * 100.0 / t,
                                                                  nnegatives, nnegatives * 100.0 / t)
    if param_loglevel >= 1:
        print "# RESULT: %i transcripts in %i genes" % (len(map_transcript2cluster), cluster_id - 1)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
