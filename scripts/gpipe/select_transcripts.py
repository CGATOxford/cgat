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
gpipe/select_transcripts.py - 
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

   python gpipe/select_transcripts.py --help

Type::

   python gpipe/select_transcripts.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import time
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Exons as Exons
import CGAT.Intervals as Intervals
import alignlib_lite
import CGAT.GTF as GTF


USAGE = """python %s < predictions > genes

Version: $Id: gpipe/select_transcripts.py 2263 2008-11-17 16:36:29Z andreas $

select a set of transcripts from genes. Each transcript is
a representative of other predictions.

Transcripts are selected from predictions in the following order:

1. Remove redundant entries:
        Redundancy = same genomic region, same predicted peptide
        -> take the one with highest coverage of query
        
2. Resolve overlapping entries:
        If one entry is contained in another, take the longest one.

Algorithm: work with list elimination.

1. Build set of overlapping clusters (genes) for genes (gene_id/overlap_id).

2. for each class in (CG, ... )

        Get list of predictions sorted by length descendingly.
        
        Remove all predictions that are contained in it.
        
        Remove spanning predictions (=dubious predictions spanning good predictions
        on the same strand )


The codes are:

h: identity, same class
l: identity, lower class

e: exon swapping
g: gene spanners
m: master
i: identical sequence
p: partial identical sequence
h: eliminated due to high percent identity
l: eliminated due to lower percent identity
u: eliminated suboptimal spanning match
f: filter by list of transcripts

""" % sys.argv[0]


def getRangesFromExons(exons, both_strands=False, contig_sizes=None):
    """convert a list of exons into a dictionary of ranges.
    """

    if both_strands is True and contig_sizes is None:
        raise ValueError(
            "supply contig_sizes, if ranges are to computed on both strands.")

    filter_ranges = {}
    for key, ee in exons.items():
        min_exon_from = min(map(lambda x: x.mGenomeFrom, ee))
        max_exon_to = max(map(lambda x: x.mGenomeTo, ee))
        contig = ee[0].mSbjctToken
        strand = ee[0].mSbjctStrand

        k = "%s:%s" % (contig, strand)
        if k not in filter_ranges:
            filter_ranges[k] = []
        filter_ranges[k].append((min_exon_from, max_exon_to, key))
        if both_strands:
            # if filtering on both strands - add range on the other strand
            if strand == "+":
                strand == "-"
            else:
                strand = "+"
            l = contig_sizes[contig]
            min_exon_from, max_exon_to = l - max_exon_to, l - min_exon_from
            if k not in filter_ranges:
                filter_ranges[k] = []
            filter_ranges[k].append((min_exon_from, max_exon_to, key))

    return filter_ranges


def FilterEliminateOverlappingTranscripts(
        exons, filter_exons,
        eliminated_predictions, contig_sizes, options):
    """eliminate predictions that overlap or span a positive set of transcripts.
    """

    eliminated = []

    # convert list of filter exons into a list of ranges.
    filter_ranges = getRangesFromExons(
        filter_exons,
        both_strands=options.filter_remove_spanning_both_strands,
        contig_sizes=contig_sizes)

    for k, r in filter_ranges.items():
        filter_ranges[k] = Intervals.combineIntervals(map(lambda x: x[:2], r))

    exon_ranges = getRangesFromExons(exons,
                                     both_strands=False)

    # and now go through exons and delete transcripts whose
    # exons overlap one of the forbidden ranges
    for k, ee in exon_ranges.items():

        if k not in filter_ranges:
            continue

        ff = filter_ranges[k]
        ee.sort()

        # set exon index e and filter index f
        # (both are indices in sorted lists)
        e, f = 0, 0

        while e < len(ee):

            efrom, eto, id = ee[e]

            # increment filter, such that its extent
            # is larger than current range ee[e] to test.
            while f < len(ff) and ff[f][1] < efrom:
                f += 1
            if f == len(ff):
                break

            if eto < ff[f][0]:
                # no overlap
                pass
            else:
                options.stdout.write(
                    "%s\t%s\n" % (id, "eliminated: filtered by %s:%i:%i" % (k, ff[f][0], ff[f][1])))
                eliminated_predictions[id] = 0
                eliminated.append((id, "f"))

            e += 1

    return eliminated


def CheckSuboptimal(rep_id,
                    exons,
                    eliminated_predictions,
                    other_ids,
                    map_prediction2data,
                    options):

    overlaps = []

    # get predictions which overlap by exons (but not completely):
    for id in other_ids:
        if id == rep_id:
            continue
        if id in eliminated_predictions:
            continue
        if Exons.CheckOverlap(exons[rep_id], exons[id]) and \
            not Exons.CheckCoverage(exons[rep_id],
                                    exons[id],
                                    max_slippage=options.max_slippage):
            overlaps.append(id)

    rep = map_prediction2data[rep_id]
    identity = rep.mPid + options.suboptimal_min_identity_difference

    for x in range(0, len(overlaps) - 1):
        id1 = overlaps[x]
        d1 = map_prediction2data[id1]
        for y in range(x + 1, len(overlaps)):
            id2 = overlaps[y]
            d2 = map_prediction2data[id2]
            if options.loglevel >= 3:
                options.stdlog.write(
                    "# suboptimal: %s ? %s + %s: %s %s %s %s %i %i %i\n" %
                    (rep_id, id1, id2,
                     d1.mQuality in options.quality_remove_suboptimal,
                     d2.mQuality in options.quality_remove_suboptimal,
                     not Exons.CheckOverlap(
                         exons[id1], exons[id2]),
                     Exons.CheckCoverageAinB(exons[rep_id],
                                             exons[id1] + exons[id2],
                                             min_terminal_exon_coverage=0.0),
                     rep.mPid, d1.mPid, d2.mPid))

            if (d1.mQuality in options.quality_remove_suboptimal and
                d2.mQuality in options.quality_remove_suboptimal) and \
                not Exons.CheckOverlap(exons[id1], exons[id2]) and \
                Exons.CheckContainedAinB(
                    exons[rep_id], exons[id1] + exons[id2],
                    min_terminal_exon_coverage=0.0) and \
                    (identity < d1.mPid) and \
                    (identity < d2.mPid):
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# elimination: %s(%s) joins %s(%s) and %s(%s)\n" %
                        (rep_id, rep.mPid,
                         id1, d1.mPid,
                         id2, d2.mPid))
                return True

    return False


def EliminateSuboptimalPredictions(data,
                                   eliminated_predictions,
                                   map_transcript2cluster,
                                   map_cluster2transcripts,
                                   map_prediction2data,
                                   exons,
                                   options):
    """remove predictions that overlap other predictions, but have
    a much lower overall percent identity.
    """

    eliminated = []
    for d in data:

        rep_id, rep_quality = d.transcript_id, d.mQuality
        if rep_quality in options.quality_keep_suboptimal:
            continue

        if rep_id in eliminated_predictions:
            continue

        cluster = map_transcript2cluster[rep_id]

        if CheckSuboptimal(rep_id,
                           exons,
                           eliminated_predictions,
                           map_cluster2transcripts[cluster],
                           map_prediction2data,
                           options):

            options.stdout.write(
                "%s\t%s\n" % (rep_id, "eliminated: suboptimal"))

            eliminated_predictions[rep_id] = 0
            eliminated.append((rep_id, "u"))
            continue

    return eliminated


def CheckExonSwop(rep_id,
                  exons,
                  eliminated_predictions,
                  other_ids,
                  map_prediction2data,
                  options):
    """check for exon swop

    return true, if exon swop occurs.

    Exon swop occurs, if this prediction joins
    two predictions, one of which should be CG.

    None of the predictions should be fully contained
    in the master prediction.

    given:
        the rep_id to analyzse
        a map of rep_id to exons
        a list of rep_ids to check against

    -> is it an exon swopper?
      -> joining two CG predictions that do not overlap and
         contain no extra exons apart from the overlapping.
    -> is it large spanning prediction?
      -> spanning many predictions, including at least one CG?

    """
    overlaps = []
    # get predictions which overlap by exons (but not completely):

    for id in other_ids:
        if id == rep_id:
            continue
        if id in eliminated_predictions:
            continue
        if Exons.CheckOverlap(exons[rep_id], exons[id]) and \
            not Exons.CheckCoverage(exons[rep_id],
                                    exons[id],
                                    max_slippage=options.max_slippage):
            overlaps.append(id)

    if options.loglevel >= 3:
        options.stdlog.write(
            "# exon swop: %s overlaps with %i out of %i predictions\n" %
            (rep_id, len(overlaps), len(other_ids)))
        options.stdlog.flush()

    for x in range(0, len(overlaps) - 1):
        id1 = overlaps[x]
        for y in range(x + 1, len(overlaps)):
            id2 = overlaps[y]
            if options.loglevel >= 4:
                options.stdlog.write(
                    "# exon swop: %s ? %s + %s: %s %s %s %s\n" %
                    (rep_id, id1, id2,
                     map_prediction2data[id1].mQuality in options.quality_remove_exon_swopper,
                     map_prediction2data[id2].mQuality in options.quality_remove_exon_swopper,
                     not Exons.CheckOverlap(exons[id1], exons[id2]),
                     Exons.CheckCoverageAinB(
                         exons[rep_id],
                         exons[id1] + exons[id2],
                         min_terminal_num_exons=0,
                         min_terminal_exon_coverage=0.7,
                         max_slippage=options.max_slippage)))

            if (map_prediction2data[id1].mQuality in options.quality_remove_exon_swopper and
                map_prediction2data[id2].mQuality in options.quality_remove_exon_swopper) and \
                not Exons.CheckOverlap(exons[id1], exons[id2]) and \
                Exons.CheckCoverageAinB(exons[rep_id], exons[id1] + exons[id2],
                                        min_terminal_num_exons=0,
                                        min_terminal_exon_coverage=0.7,
                                        max_slippage=options.max_slippage):
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# elimination: %s(%s) joins %s(%s) and %s(%s)\n" %
                        (rep_id, map_prediction2data[rep_id].mQuality,
                         id1, map_prediction2data[id1].mQuality,
                         id2, map_prediction2data[id2].mQuality))
                return True

    return False


def EliminateExonSwoppers(data,
                          eliminated_predictions,
                          map_transcript2cluster,
                          map_cluster2transcripts,
                          map_prediction2data,
                          exons,
                          options):
    """remove predictions due to exon swop.
    """

    eliminated = []
    for d in data:

        rep_id, rep_quality = d.transcript_id, d.mQuality

        if rep_quality in options.quality_keep_exon_swopper:
            continue
        if rep_id in eliminated_predictions:
            continue

        cluster = map_transcript2cluster[rep_id]

        if CheckExonSwop(rep_id,
                         exons,
                         eliminated_predictions,
                         map_cluster2transcripts[cluster],
                         map_prediction2data,
                         options):

            options.stdout.write(
                "%s\t%s\n" % (rep_id, "eliminated: exon swopping"))

            eliminated_predictions[rep_id] = 0
            eliminated.append((rep_id, "e"))
            continue

    return eliminated


def EliminateSpanningPredictions(data,
                                 prediction_id,
                                 eliminated_predictions,
                                 options,
                                 filter_quality=None,
                                 this_quality=None):
    """eliminate dubious predictions that contain good predictions inside
    them without exon overlap.

    """

    eliminated = []

    for d in data:

        id = d.transcript_id
        if options.loglevel >= 4:
            options.stdlog.write(
                "# processing: id=%s class=%s\n" % (id, d.mQuality))

        if id in eliminated_predictions:
            continue

        if this_quality in options.quality_genes and \
                d.mQuality in options.quality_remove_dubious:
            eliminated_predictions[id] = prediction_id
            eliminated.append((id, "s"))

    return eliminated


def EliminateGeneSpanners(data, eliminated_predictions, exons, options):
    """eliminate predictions spanning full length other predictions
    without exon overlap.

    """

    def ProcessChunk(chunk, eliminated_predictions, exons):
        """process a cluster of overlapping predictions.

        Chunks are sorted by first position.
        Thus, only former can span later.
        """

        eliminated = []
        for x in range(0, len(chunk) - 1):
            xfrom, xto, xid, xquality = chunk[x]
            if xquality in options.quality_keep_gene_spanners:
                continue
            for y in range(x + 1, len(chunk)):
                yfrom, yto, yid, yquality = chunk[y]
                # print xid, yid, xfrom < yfrom, xto > yto,
                # Exons.CheckOverlap(exons[xid], exons[yid] ), xquality,
                # yquality
                if xfrom < yfrom and \
                        xto > yto and \
                        not Exons.CheckOverlap(
                            exons[str(xid)], exons[str(yid)]) and \
                        yquality in options.quality_remove_gene_spanners:
                    eliminated_predictions[xid] = 0
                    eliminated.append((xid, "g"))
                    if options.loglevel >= 1:
                        options.stdlog.write(
                            "# elimination: %s(%s) spans %s(%s)\n" %
                            (str(xid), xquality, str(yid), yquality))
                    break
        return eliminated

    # cluster by overlap
    eliminated = []
    chunk = []
    max_to = 0
    last_contig = None
    last_strand = None

    data.sort(key=lambda x: (x.contig, x.strand, x.start))

    for d in data:

        if last_contig != d.contig or \
                last_strand != d.strand or \
                max_to < d.end:
            if len(chunk) > 1:
                eliminated += ProcessChunk(chunk,
                                           eliminated_predictions, exons)
            max_to = 0
            chunk = []
            last_contig = d.contig
            last_strand = d.strand

        chunk.append((d.start, d.end, d.transcript_id, d.mQuality))
        max_to = max(max_to, d.end)

    eliminated += ProcessChunk(chunk, eliminated_predictions, exons)

    return eliminated


def EliminateRedundantEntries(rep,
                              data,
                              eliminated_predictions,
                              options,
                              peptides,
                              extended_peptides,
                              filter_quality=None,
                              this_quality=None):
    """eliminate redundant entries in a set."""

    eliminated = []

    rep_id = rep.transcript_id
    rep_coverage, rep_pid = rep.mQueryCoverage, rep.mPid

    alignator = alignlib_lite.makeAlignatorDPFull(
        alignlib_lite.ALIGNMENT_LOCAL, options.gop, options.gep)
    result = alignlib_lite.makeAlignmentVector()

    rep_seq = peptides[rep_id]
    rep_extended_seq = extended_peptides[rep_id]

    for entry in data:

        mem_id, mem_coverage, mem_pid, mem_quality = (entry.transcript_id,
                                                      entry.mQueryCoverage,
                                                      entry.mPid,
                                                      entry.mQuality)

        mem_seq = peptides[mem_id]
        mem_extended_seq = extended_peptides[mem_id]

        if options.loglevel >= 4:
            options.stdlog.write(
                "# processing: id=%s class=%s\n" % (mem_id, mem_quality))

        if mem_id in eliminated_predictions:
            continue

        if mem_extended_seq == rep_extended_seq:
            eliminated_predictions[mem_id] = rep_id
            eliminated.append((mem_id, "i"))

        elif mem_extended_seq in rep_extended_seq:
            eliminated_predictions[mem_id] = rep_id
            eliminated.append((mem_id, "p"))

        else:
            if mem_quality != this_quality or \
                    mem_quality in options.quality_exclude_same:

                seq1 = alignlib_lite.makeSequence(str(rep_seq))
                seq2 = alignlib_lite.makeSequence(str(mem_seq))

                alignator.align(result, seq1, seq2)

                if options.loglevel >= 5:
                    options.stdlog.write(
                        "# ali\n%s\n" %
                        alignlib_lite.AlignmentFormatExplicit(result, seq1, seq2))

                pidentity = 100 * \
                    alignlib_lite.calculatePercentIdentity(result, seq1, seq2)

                num_gaps = result.getNumGaps()

                if options.loglevel >= 4:
                    options.stdlog.write(
                        "# processing: id=%s class=%s pid=%5.2f rep_cov=%i mem_cov=%i\n" %
                        (mem_id, mem_quality, pidentity, rep_coverage, mem_coverage))

                if pidentity >= options.min_identity:

                    keep = False
                    if rep_coverage < mem_coverage - options.safety_coverage or \
                       rep_pid < mem_pid - options.safety_pide:
                        keep = True
                        reason = "covpid"
                    elif num_gaps >= options.max_gaps and \
                            mem_coverage > rep_coverage - options.safety_coverage:
                        keep = True
                        reason = "gaps"
                    elif mem_coverage >= rep_coverage - options.safety_coverage and \
                            100 * (result.getColTo() - result.getColFrom()) / len(mem_seq) < options.max_member_coverage:
                        keep = True
                        reason = "memcov"

                    if keep:
                        options.stdlog.write(
                            "# WARNING: not removing possibly good prediction: "
                            "%s: rep = %s, mem = %s, rep_cov=%i, rep_pid=%i, "
                            "mem_cov=%i, mem_pid=%i\n" %
                            (reason, rep_id, mem_id, rep_coverage, rep_pid,
                             mem_coverage, mem_pid))
                    else:
                        eliminated_predictions[mem_id] = rep_id
                        eliminated.append((mem_id, "h"))

                elif pidentity >= options.min_identity_non_genes and \
                        this_quality in options.quality_genes and \
                        mem_quality not in options.quality_genes:
                    if rep_coverage < mem_coverage - options.safety_coverage or \
                       rep_pid < mem_pid - options.safety_pide:
                        options.stdlog.write("# WARNING: not removing possibly good prediction: rep = %s, mem = %s, rep_cov=%i, rep_pid=%i, mem_cov=%i, mem_pid=%i\n" %
                                             (rep_id, mem_id, rep_coverage, rep_pid, mem_coverage, mem_pid))
                    else:
                        eliminated_predictions[mem_id] = rep_id
                        eliminated.append((mem_id, "l"))

    return eliminated


def EliminateRedundantEntriesByRange(rep, data,
                                     eliminated_predictions,
                                     options,
                                     peptides,
                                     extended_peptides,
                                     filter_quality=None,
                                     this_quality=None,
                                     exons={}):
    """eliminate redundant entries for prediction_id.

    Overlapping entries are defined by sequence ranges.
    """

    rep_id = rep.transcript_id

    if options.loglevel >= 2:
        options.stdlog.write("# eliminating transcripts for %s\n" % rep_id)
        options.stdlog.flush()

    rep_start, rep_end = rep.start, rep.end
    rep_contig, rep_strand = rep.contig, rep.strand

    # collect overlapping entries:
    entries = []
    for d in data:
        if d.transcript_id == rep_id:
            continue
        if d.transcript_id in eliminated_predictions:
            continue
        if d.contig != rep_contig:
            continue
        if d.strand != rep_strand:
            continue
        if filter_quality and d.mQuality in filter_quality:
            continue

        if max(d.mExtendedEnd, rep_end) - min(d.mExtendedStart, rep_start) > 0:
            entries.append(d)

    eliminated = EliminateRedundantEntries(rep, entries,
                                           eliminated_predictions,
                                           options,
                                           peptides,
                                           extended_peptides,
                                           filter_quality,
                                           this_quality)

    return eliminated


def EliminateRedundantEntriesByOverlap(rep,
                                       data,
                                       eliminated_predictions,
                                       options,
                                       peptides,
                                       extended_peptides,
                                       filter_quality=None,
                                       this_quality=None,
                                       exons={}):
    """eliminate redundant entries for prediction_id.

    Overlapping entries are defined overlap_id
    """

    if options.loglevel >= 3:
        options.stdlog.write(
            "# eliminating transcripts for %s out of %i\n" %
            (rep.transcript_id, len(data)))
        options.stdlog.flush()

    entries = []
    for mem in data[1:]:

        if mem.gene_id != rep.gene_id:
            continue
        if mem.mQuality in filter_quality:
            continue

        entries.append(mem)

    eliminated = EliminateRedundantEntries(rep, entries,
                                           eliminated_predictions,
                                           options,
                                           peptides,
                                           extended_peptides,
                                           filter_quality,
                                           this_quality)

    # get all overlaps no matter what strand (use export_sbjct_genome_from/to)
    entries = []
    if options.remove_spanning_predictions:

        for mem in data[1:]:

            if mem.contig != rep.contig:
                continue
            if min(mem.end, rep.end) - max(mem.start, rep.start) <= 0:
                continue
            if mem.mQuality in filter_quality:
                continue
            if mem.mQuality not in options.quality_remove_dubious:
                continue

            entries.append(mem)

        eliminated += EliminateSpanningPredictions(entries,
                                                   rep_id,
                                                   eliminated_predictions,
                                                   options,
                                                   filter_quality,
                                                   this_quality)

    if options.loglevel >= 3 and len(eliminated) > 0:
        options.stdlog.write(
            "# eliminated %i transcripts for %s\n" %
            (len(eliminated), rep_id))
        options.stdlog.flush()

    return eliminated


def PrintMembers(rep_id, outfile, eliminated, eliminated_by_method):
    """write members to outfile and keep counts.
    """

    nmembers = 0
    for mem_id, method in eliminated:
        nmembers += 1
        if method not in eliminated_by_method:
            eliminated_by_method[method] = 0
        eliminated_by_method[method] += 1
        outfile.write("%s\t%s\t%s\n" % (rep_id, mem_id, method))

    outfile.flush()

    return nmembers

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-o", "--overlap", dest="overlap_residues", type="int",
        help="overlap residues.")
    parser.add_option(
        "-t", "--filter-tokens", dest="filename_filter_tokens", type="string",
        help="filename to filter tokens.")
    parser.add_option(
        "-i", "--exon-identity", dest="exon_identity", action="store_true",
        help="exon identity.")
    parser.add_option(
        "--exons-file", dest="filename_exons", type="string",
        help="filename with exon information.")
    parser.add_option(
        "-m", "--output-members", dest="filename_members", type="string",
        help="output filename with members.")
    parser.add_option(
        "--overlap-id", dest="overlap_id", action="store_true",
        help="overlap id.")
    parser.add_option(
        "-s", "--remove-spanning", dest="remove_spanning_predictions", action="store_true",
        help="remove spanning predictions.")
    parser.add_option(
        "-c", "--remove-complement", dest="remove_complementary_predictions", action="store_true",
        help="remove complementary predictions.")
    parser.add_option(
        "--remove-exon-swoppers", dest="remove_exon_swoppers", action="store_true",
        help="remove exon swoppers.")
    parser.add_option(
        "--remove-gene-spanners", dest="remove_gene_spanners", action="store_true",
        help="remove gene spanners.")
    parser.add_option(
        "--remove-suboptimal", dest="remove_suboptimal", action="store_true",
        help="remove suboptimal predictions.")
    parser.add_option(
        "-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
        help="filename with peptide information.")
    parser.add_option(
        "--extended-peptides", dest="filename_extended_peptides", type="string",
        help="filename with peptide information - after extension.")
    parser.add_option("--test", dest="test_nids", type="string",
                      help="test nids.")
    # filter options
    parser.add_option(
        "--filter-transcripts", dest="filter_filename_transcripts", type="string",
        help="filename with transcripts that are used to filter.")
    parser.add_option(
        "--filter-remove-spanning", dest="filter_remove_spanning", action="store_true",
        help="remove all transcripts that span the filter set.")
    parser.add_option(
        "-g", "--genome-file", dest="genome_file", type="string",
        help="filename with genomic data (indexed).")
    parser.add_option(
        "--discard-large-clusters", dest="discard_large_clusters", type="int",
        help="if set discard clusters bigger than this size (patch) [default=%default].")

    parser.set_defaults(
        filename_members=None,
        filename_peptides=None,
        filename_extended_peptides=None,
        filename_exons=None,
        quality_hierarchy=("CG", "PG", "SG", "RG", "CP", "PP", "SP",
                           "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK"),
        # Classes, where redundancy is removed by similarity. When
        # exon structure is not conserved, I can't predict alternative
        # splice variants, so remove the redundancy.
        quality_exclude_same=("UG", "UP", "UF", "BF", "UK"),
        quality_genes=("CG", "SG", "PG", "RG", "UG"),
        # class that can be removed in spanning/complementary predictions
        quality_remove_dubious=("UG", "UP", "UF", "BF", "UK"),
        # class that is required for defining exon swopper event
        quality_remove_exon_swopper=("CG", "PG"),
        # class that will kept, in spite of being an exons swopper.
        quality_keep_exon_swopper=(),
        # class that is required for removing gene spanners
        quality_remove_gene_spanners=("CG"),
        # class that will kept, in spite of being a gene spanner
        quality_keep_gene_spanners=(),
        # class that is required for defining suboptimal matches
        quality_remove_suboptimal=("CG", "PG"),
        # class that will be kept, in spite of being a suboptimal match
        quality_keep_suboptimal=(),
        # gap penalties
        gop=-10.0,
        gep=-1.0,
        # maximum number of gaps to allow in alignment
        max_gaps=20,
        # threshold of percent identity that allows to remove a prediction
        # of a lower class.
        # This allows for insertions/deletions
        min_identity=98,
        # threshold of percent identity that allows to remove a prediction
        # of a non-gene by a gene
        min_identity_non_genes=80,
        # safety threshold: do not remove, if coverage of member is by x better
        # than representative
        safety_pide=10,
        safety_coverage=10,
        overlap_id=False,
        remove_spanning_predictions=False,
        remove_exon_swoppers=False,
        remove_gene_spanners=False,
        remove_suboptimal=False,
        # nids to use for testing
        test_nids=None,
        # remove members with less than maximum coverage
        max_member_coverage=90,
        # maximum allowable exon slippage
        max_slippage=9,
        # minimum difference in identity for suboptimal predictions to be
        # removed.
        suboptimal_min_identity_difference=10,
        # filter options
        filter_filename_transcripts=None,
        filter_remove_spanning=True,
        filter_remove_spanning_both_strands=True,
        genome_file=None,
        discard_large_clusters=None)

    (options, args) = E.Start(parser, add_database_options=True)

    if options.test_nids:
        options.test_nids = options.test_nids.split(",")

    # list of eliminated predictions
    eliminated_predictions = {}

    if options.filename_members:
        outfile_members = open(options.filename_members, "w")
    else:
        outfile_members = sys.stdout

    ######################################################
    ######################################################
    ######################################################
    # data
    ######################################################
    data = []

    class Entry:

        def __init__(self, gff):
            self.mPid = float(gff["pid"])
            self.mQueryCoverage = float(gff["qcov"])
            self.gene_id = gff['gene_id']
            self.transcript_id = gff['transcript_id']
            self.mExtendedStart = int(gff['xstart'])
            self.mExtendedEnd = int(gff['xend'])
            self.start = gff.start
            self.contig = gff.contig
            self.strand = gff.strand
            self.end = gff.end
            self.mQuality = gff['class']

    for gff in GTF.iterator(sys.stdin):
        data.append(Entry(gff))

    if options.loglevel >= 1:
        options.stdlog.write("# read %i transcripts.\n" % len(data))
        options.stdlog.flush()

    ######################################################
    ######################################################
    ######################################################
    # read peptide sequences
    ######################################################
    if options.loglevel >= 1:
        options.stdlog.write("# loading peptide databases ... ")
        options.stdlog.flush()

    if options.filename_peptides:
        peptides = IndexedFasta.IndexedFasta(options.filename_peptides)
        peptide_lengths = peptides.getContigSizes()
    else:
        peptide_lengths = {}
        peptides = {}

    ######################################################
    ######################################################
    ######################################################
    # read extended peptide sequences
    ######################################################
    if options.filename_extended_peptides:
        extended_peptides = IndexedFasta.IndexedFasta(
            options.filename_extended_peptides)
    else:
        extended_peptides = {}

    if options.loglevel >= 1:
        options.stdlog.write("finished\n")
        options.stdlog.flush()

    ######################################################
    ######################################################
    ######################################################
    # open genome file
    ######################################################
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
        contig_sizes = fasta.getContigSizes()
    else:
        contig_sizes = {}

    ######################################################
    ######################################################
    ######################################################
    # reading exons, clustering and formatting them.
    ######################################################
    if options.filename_exons:
        if options.loglevel >= 1:
            options.stdlog.write("# reading exon boundaries ... ")
            options.stdlog.flush()

        ids = [x.transcript_id for x in data]

        exons = Exons.ReadExonBoundaries(open(options.filename_exons, "r"),
                                         contig_sizes=contig_sizes,
                                         filter=set(ids))

        if options.loglevel >= 1:
            options.stdlog.write(
                "done - read exons for %i transcripts\n" % (len(exons)))

        if len(exons) == 0:
            raise ValueError("no exons found in table.")

        # flag terminal exons
        Exons.SetRankToPositionFlag(exons)

        identity_map_cluster2transcripts, identity_map_transcript2cluster =\
            Exons.ClusterByExonIdentity(exons,
                                        max_terminal_num_exons=3,
                                        max_slippage=options.max_slippage,
                                        loglevel=options.loglevel)

        overlap_map_cluster2transcripts, overlap_map_transcript2cluster =\
            Exons.ClusterByExonOverlap(exons,
                                       min_overlap=10,
                                       loglevel=options.loglevel)
    else:
        exons = {}

    ######################################################
    nrepresentatives, nmembers, neliminated = 0, 0, 0
    eliminated_by_method = {}

    ######################################################
    ######################################################
    ######################################################
    # read filter transcripts and apply filters
    ######################################################
    if options.filter_filename_transcripts:

        if options.loglevel >= 1:
            options.stdlog.write(
                "# reading exon boundaries for filter set ... ")
            options.stdlog.flush()

        filter_exons = Exons.ReadExonBoundaries(open(options.filter_filename_transcripts, "r"),
                                                delete_missing=True,
                                                contig_sizes=contig_sizes)

        if options.loglevel >= 1:
            options.stdlog.write(
                "done - read exons for %i transcripts\n" % (len(filter_exons)))

        t = time.time()
        eliminated = FilterEliminateOverlappingTranscripts(exons,
                                                           filter_exons,
                                                           eliminated_predictions,
                                                           contig_sizes,
                                                           options)

        n = PrintMembers(0, outfile_members, eliminated, eliminated_by_method)
        neliminated += n
        if options.loglevel >= 1:
            options.stdlog.write(
                "# removed %i transcripts overlapping or spanning transcripts in %i seconds.\n" % (n, time.time() - t))
            options.stdlog.flush()

    if options.remove_exon_swoppers and not exons:
        raise ValueError(
            "please specify exon table if using --remove-swoppers.")
    if options.remove_gene_spanners and not exons:
        raise ValueError(
            "please specify exon table if using --remove-gene-spanners.")

    ##########################################################################
    # remove predictions spanning other predictions but do not overlap with
    # them on an exon level.
    if options.remove_gene_spanners and exons:
        if options.loglevel >= 1:
            options.stdlog.write("# removing gene spanners\n")
            options.stdlog.flush()

        t = time.time()
        eliminated = EliminateGeneSpanners(data,
                                           eliminated_predictions,
                                           exons,
                                           options)

        n = PrintMembers(0, outfile_members, eliminated, eliminated_by_method)
        neliminated += n
        if options.loglevel >= 1:
            options.stdlog.write(
                "# removed %i gene spanners in %i seconds\n" % (n, time.time() - t))
            options.stdlog.flush()

    ##########################################################################
    # sort data by quality, length of prediction and coverage * pid

    if options.loglevel >= 1:
        options.stdlog.write("# sorting data\n")
        options.stdlog.flush()

    map2pos = {}
    for x in range(len(options.quality_hierarchy)):
        map2pos[options.quality_hierarchy[x]] = x

    data.sort(key=lambda x: (map2pos[x.mQuality], len(
        extended_peptides[x.transcript_id]), x.mQueryCoverage * x.mPid))

    # build map of prediction to quality
    map_prediction2data = {}
    for d in data:
        map_prediction2data[d.transcript_id] = d

    if options.loglevel >= 1:
        options.stdlog.write("# sorting data finished\n")
        options.stdlog.flush()

    ##########################################################################
    # remove predictions joining two other complete non-overlapping predictions
    if options.remove_exon_swoppers and exons:

        if options.loglevel >= 1:
            options.stdlog.write("# removing exon swoppers\n")
            options.stdlog.flush()

        eliminated = EliminateExonSwoppers(data,
                                           eliminated_predictions,
                                           identity_map_transcript2cluster,
                                           identity_map_cluster2transcripts,
                                           map_prediction2data,
                                           exons,
                                           options)

        n = PrintMembers(0, outfile_members, eliminated, eliminated_by_method)
        neliminated += n

        if options.loglevel >= 1:
            options.stdlog.write("# removed %i exon swoppers\n" % n)
            options.stdlog.flush()

    ##########################################################################
    # remove suboptimal predictions
    if options.remove_suboptimal and exons:

        if options.loglevel >= 1:
            options.stdlog.write("# removing suboptimal predictions\n")
            options.stdlog.flush()

        t = time.time()
        eliminated = EliminateSuboptimalPredictions(data,
                                                    eliminated_predictions,
                                                    overlap_map_transcript2cluster,
                                                    overlap_map_cluster2transcripts,
                                                    map_prediction2data,
                                                    exons,
                                                    options)

        n = PrintMembers(0, outfile_members, eliminated, eliminated_by_method)
        neliminated += n

        if options.loglevel >= 1:
            options.stdlog.write(
                "# removed %i suboptimal predictions in %i seconds\n" % (n, time.time() - t))
            options.stdlog.flush()

    ##########################################################################
    # remove redundant predictions
    l = len(data)

    options.report_step = max(1, int(l / 100))

    t2 = time.time()

    last_quality = None
    qualities = []

    options.stdout.write("%s\t%s\n" % ("rep", "comment"))

    for x in range(len(data)):

        if options.loglevel >= 1:
            if x % options.report_step == 0:
                options.stdlog.write("# process: %i/%i = %i %%, %i/%i = %i %% in %i seconds\n" %
                                     (x + 1, l,
                                      int(100 * (x + 1) / l),
                                      len(eliminated_predictions), l,
                                      100 * len(eliminated_predictions) / l,
                                      time.time() - t2))

                options.stdlog.flush()

        rep = data[x]

        rep_id, rep_quality = rep.transcript_id, rep.mQuality

        if rep_id in eliminated_predictions:
            continue

        if rep_quality != last_quality:
            if last_quality:
                qualities.append(last_quality)
            last_quality = rep_quality

        if options.loglevel >= 2:
            options.stdlog.write(
                "# processing prediction %s|%s\n" % (rep_id, rep_quality))
            options.stdlog.flush()

        eliminated = []

        if options.overlap_id:
            eliminated += EliminateRedundantEntriesByOverlap(rep,
                                                             data[x + 1:],
                                                             eliminated_predictions,
                                                             options,
                                                             peptides,
                                                             extended_peptides,
                                                             filter_quality=qualities,
                                                             this_quality=rep_quality)

        else:
            eliminated += EliminateRedundantEntriesByRange(rep,
                                                           data,
                                                           eliminated_predictions,
                                                           options,
                                                           peptides,
                                                           extended_peptides,
                                                           filter_quality=qualities,
                                                           this_quality=rep_quality)

        options.stdout.write("%s\t%i\n" % (rep_id, len(eliminated)))

        if outfile_members:
            outfile_members.write("%s\t%s\tm\n" % (str(rep_id), str(rep_id)))
            nrepresentatives += 1
            nmembers += PrintMembers(rep_id, outfile_members,
                                     eliminated, eliminated_by_method)

    if outfile_members != sys.stdout:
        outfile_members.close()

    options.stdlog.write("# representatives=%i, members=%i, eliminated=%i, total=%i\n" %
                         (nrepresentatives, nmembers, neliminated,
                          nrepresentatives + nmembers + neliminated))

    options.stdlog.write("# elimination by method:\n")

    for v, c in eliminated_by_method.items():
        options.stdlog.write("# method=%s, count=%i\n" % (v, c))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
