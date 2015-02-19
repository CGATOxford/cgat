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
gpipe/compare_predictions2exons.py - 
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

   python gpipe/compare_predictions2exons.py --help

Type::

   python gpipe/compare_predictions2exons.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import math
import alignlib_lite
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import CGAT.PredictionParser as PredictionParser
import CGAT.Exons as Exons

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/compare_predictions2exons.py 2011 2008-07-04 10:40:51Z andreas $

Evaluate genewise alignments.

Build a file with exon comparisions between genes.
""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: gpipe/compare_predictions2exons.py 2011 2008-07-04 10:40:51Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("-b", "--boundaries", dest="filename_boundaries", type="string",
                      help="filename with exon boundaries.")

    parser.add_option("-e", "--exons-file", dest="filename_exons", type="string",
                      help="filename with exons (output).")

    parser.add_option("-p", "--peptides-fasta-file", dest="filename_peptides", type="string",
                      help="filename with peptide sequences.")

    parser.add_option("-w", "--write-notfound", dest="write_notfound", action="store_true",
                      help="print exons for predictions not found in reference.")

    parser.add_option("-q", "--quality-pide", dest="quality_threshold_pide", type="int",
                      help="quality threshold (pide) for exons.")

    parser.set_defaults(
        genome_file="genome",
        filename_boundaries=None,
        filename_exons=None,
        filename_peptides=None,
        quality_threshold_pide=0,
        write_notfound=False,
        # allowed number of nucleotides for exon boundaries to
        # be considered equivalent.
        slipping_exon_boundary=9,
        # stop codons to search for
        stop_codons=("TAG", "TAA", "TGA"), )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    reference_exon_boundaries = {}
    if options.filename_boundaries:
        reference_exon_boundaries = Exons.ReadExonBoundaries(open(options.filename_boundaries, "r"),
                                                             do_invert=1,
                                                             remove_utr=1)
        E.info("read exon boundaries for %i queries" %
               len(reference_exon_boundaries))

    if options.filename_exons:
        outfile_exons = open(options.filename_exons, "w")
        outfile_exons.write("%s\n" % "\t".join((
            "prediction_id",
            "exon_id",
            "exon_from",
            "exon_to",
            "exon_frame",
            "reference_id",
            "reference_from",
            "reference_to",
            "reference_phase",
            "pidentity",
            "psimilarity",
            "nframeshifts",
            "ngaps",
            "nstopcodons",
            "is_ok",
            "genome_exon_from",
            "genome_exon_to")))

    else:
        outfile_exons = None

    if options.filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(options.filename_peptides, "r"))
        E.info("read peptide sequences for %i queries" %
               len(peptide_sequences))
    else:
        peptide_sequences = {}

    entry = PredictionParser.PredictionParserEntry()
    last_filename_genome = None

    nfound, nmissed_exons, nmissed_length = 0, 0, 0
    nempty_alignments = 0

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    options.stdout.write("%s\n" % "\t".join((
        "prediction_id",
        "number",
        "dubious_exons",
        "boundaries_sum",
        "boundaries_max",
        "identical_exons",
        "inserted_exons",
        "deleted_exons",
        "inserted_introns",
        "deleted_introns",
        "truncated_Nterminus",
        "truncated_Cterminus",
        "deleted_Nexons",
        "deleted_Cexons",
        "inserted_Nexons",
        "inserted_Cexons")))

    for line in sys.stdin:

        if line[0] == "#":
            continue

        try:
            entry.Read(line)
        except ValueError, msg:
            print "# parsing failed with msg %s in line %s" % (msg, line[:-1])
            sys.exit(1)

        exons = Genomics.Alignment2ExonBoundaries(entry.mMapPeptide2Genome,
                                                  query_from=entry.mQueryFrom,
                                                  sbjct_from=entry.mSbjctGenomeFrom,
                                                  add_stop_codon=0)

        if exons[-1][4] != entry.mSbjctGenomeTo:
            print "# WARNING: discrepancy in exon calculation!!!"
            for e in exons:
                print "#", str(e)
            print "#", str(entry)

        if options.loglevel >= 5:
            for e in exons:
                print "#", str(e)

        genomic_fragment = fasta.getSequence(entry.mSbjctToken, entry.mSbjctStrand,
                                             entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo)

        skip = False
        if peptide_sequences.has_key(entry.mQueryToken):

            query_sequence = alignlib_lite.makeSequence(
                peptide_sequences[entry.mQueryToken])
            sbjct_sequence = alignlib_lite.makeSequence(entry.mTranslation)

            percent_similarity, percent_identity = 0, 0
            if query_sequence.getLength() < entry.mMapPeptide2Translation.getRowTo():
                print "# WARNING: query sequence %s is too short: %i %i" % (entry.mQueryToken,
                                                                            query_sequence.getLength(
                                                                            ),
                                                                            entry.mMapPeptide2Translation.getRowTo())
                sys.stdout.flush()
                nmissed_length += 1
                skip = True

            elif sbjct_sequence.getLength() < entry.mMapPeptide2Translation.getColTo():
                print "# WARNING: sbjct sequence %s is too short: %i %i" % (entry.mSbjctToken,
                                                                            sbjct_sequence.getLength(
                                                                            ),
                                                                            entry.mMapPeptide2Translation.getColTo())
                sys.stdout.flush()
                nmissed_length += 1
                skip = True
            else:
                alignlib_lite.rescoreAlignment(entry.mMapPeptide2Translation,
                                               query_sequence,
                                               sbjct_sequence,
                                               alignlib_lite.makeScorer(query_sequence, sbjct_sequence))
                percent_identity = alignlib_lite.calculatePercentIdentity(entry.mMapPeptide2Translation,
                                                                          query_sequence,
                                                                          sbjct_sequence) * 100
                percent_similarity = alignlib_lite.calculatePercentSimilarity(
                    entry.mMapPeptide2Translation) * 100

            E.debug("prediction %s: percent identity/similarity: before=%5.2f/%5.2f, realigned=%5.2f/%5.2f" % (
                    str(entry.mPredictionId),
                    entry.mPercentSimilarity,
                    entry.mPercentIdentity,
                    percent_similarity,
                    percent_identity))

        else:
            query_sequence = None
            sbjct_sequence = None

        # default values
        exons_num_exons = "na"
        exons_boundaries_sum = "na"
        exons_boundaries_max = "na"
        dubious_exons = "na"

        ndeleted_exons, ninserted_exons, ndeleted_introns, ninserted_introns, nidentical_exons = 0, 0, 0, 0, 0
        truncated_Nterminal_exon, truncated_Cterminal_exon = 0, 0
        ndeleted_Nexons, ndeleted_Cexons = 0, 0
        ninserted_Nexons, ninserted_Cexons = 0, 0

        exons_offset = exons[0][3]

        if not reference_exon_boundaries.has_key(entry.mQueryToken):
            print "# WARNING: sequence %s has no exon boundaries" % (entry.mQueryToken)
            sys.stdout.flush()
            nmissed_exons += 1
            skip = True

        if not skip:

            nfound += 1

            ref_exons = reference_exon_boundaries[entry.mQueryToken]

            ref_exons_offset = ref_exons[0].mGenomeFrom

            exons_num_exons = len(ref_exons) - len(exons)
            exons_boundaries_sum = 0
            exons_phase = 0
            exons_boundaries_max = 0
            dubious_exons = 0

            inserted_exons = 0
            temp_inserted_exons = 0

            if options.loglevel >= 3:
                for e in exons:
                    options.stdlog.write("# %s\n" % str(e))
                for e in ref_exons:
                    options.stdlog.write("# %s\n" % str(e))

            min_pide = entry.mPercentIdentity * \
                options.quality_threshold_pide / 100

            in_sync = 0
            e, r = 0, 0

            while e < len(exons) and r < len(ref_exons):

                this_e, this_r = e + 1, r + 1
                percent_identity = 0
                percent_similarity = 0
                is_good_exon = 0

                if options.loglevel >= 4:
                    options.stdlog.write(
                        "# current exons: %i and %i\n" % (e, r))
                    sys.stdout.flush()

                exon_from, exon_to, exon_phase, exon_genome_from, exon_genome_to, exon_ali = exons[
                    e][0:6]
                ref_from, ref_to, ref_phase, ref_genome_from, ref_genome_to = (ref_exons[r].mPeptideFrom,
                                                                               ref_exons[
                                                                                   r].mPeptideTo,
                                                                               ref_exons[
                                                                                   r].frame,
                                                                               ref_exons[
                                                                                   r].mGenomeFrom,
                                                                               ref_exons[r].mGenomeTo)

                ref_genome_from -= ref_exons_offset
                ref_genome_to -= ref_exons_offset

                # get percent identity for exon
                exon_percent_identity = 0
                exon_percent_similarity = 0

                if query_sequence and sbjct_sequence:

                    tmp_ali = alignlib_lite.makeAlignmentVector()

                    xquery_from = exon_from / 3
                    xquery_to = exon_to / 3

                    alignlib_lite.copyAlignment(
                        tmp_ali, entry.mMapPeptide2Translation, xquery_from, xquery_to)

                    if tmp_ali.getLength() == 0:
                        options.stdlog.write("# WARNING: empty alignment %s\n" % str(
                            (ref_from, exon_from, ref_to, exon_to, xquery_from, xquery_to)))
                        nempty_alignments += 1
                    else:
                        if options.loglevel >= 5:
                            options.stdlog.write("# %s\n" % str(
                                alignlib_lite.AlignmentFormatExplicit(tmp_ali, query_sequence, sbjct_sequence)))

                        exon_percent_identity = alignlib_lite.calculatePercentIdentity(tmp_ali,
                                                                                       query_sequence,
                                                                                       sbjct_sequence) * 100
                        exon_percent_similarity = alignlib_lite.calculatePercentSimilarity(
                            tmp_ali) * 100

                if exon_percent_identity >= min_pide:
                    is_good_exon = 1
                else:
                    is_good_exon = 0

                if e < len(exons) - 1:
                    (next_exon_from, next_exon_to, next_exon_phase,
                     next_exon_genome_from, next_exon_genome_to, next_exon_ali) = exons[e + 1][0:6]
                else:
                    (next_exon_from, next_exon_to, next_exon_phase,
                     next_exon_genome_from, next_exon_genome_to, next_exon_ali) = 0, 0, 0, 0, 0, []

                if r < len(ref_exons) - 1:
                    next_ref_from, next_ref_to, next_ref_phase = (ref_exons[r + 1].mPeptideFrom,
                                                                  ref_exons[
                                                                      r + 1].mPeptideTo,
                                                                  ref_exons[r + 1].frame)
                else:
                    next_ref_from, next_ref_to, next_ref_phase = 0, 0, 0

                if options.loglevel >= 2:
                    options.stdlog.write("# %s\n" % "\t".join(map(str, (entry.mQueryToken,
                                                                        exon_from, exon_to, exon_phase,
                                                                        exon_genome_from, exon_genome_to,
                                                                        ref_from, ref_to, ref_phase))))
                    sys.stdout.flush()

                # beware of small exons.
                # if less than options.slipping_exon_boundary: boundary is 0
                # check if end is more than options.splipping_exon_boundary
                # apart as well.
                if exon_to - exon_from <= options.slipping_exon_boundary or \
                        ref_to - ref_from <= options.slipping_exon_boundary:
                    boundary = 0
                else:
                    boundary = options.slipping_exon_boundary

                if ref_to <= exon_from + boundary and \
                   ref_to <= exon_to - options.slipping_exon_boundary:
                    # no overlap
                    is_good_exon = 0
                    if e == 0:
                        ndeleted_Nexons += 1
                    else:
                        ndeleted_exons += 1
                    r += 1
                    exon_from, exon_to, exon_phase, exon_genome_from, exon_genome_to = 0, 0, 0, 0, 0
                    overlap = 0
                elif exon_to <= ref_from + boundary and \
                        exon_to <= ref_to - options.slipping_exon_boundary:
                    # no overlap
                    is_good_exon = 0
                    if r == 0:
                        ninserted_Nexons += 1
                    else:
                        ninserted_exons += 1
                    e += 1
                    ref_from, ref_to, ref_phase = 0, 0, 0
                    overlap = 0
                else:
                    # overlap
                    overlap = 1
                    dfrom = int(math.fabs(exon_from - ref_from))
                    dto = int(math.fabs(exon_to - ref_to))

                    # get percent identity for overlapping fragment
                    if query_sequence and sbjct_sequence:
                        # this the problem
                        tmp_ali = alignlib_lite.makeAlignmentVector()

                        xquery_from = max(ref_from / 3, exon_from / 3)
                        xquery_to = min(ref_to / 3, exon_to / 3)

                        alignlib_lite.copyAlignment(
                            tmp_ali, entry.mMapPeptide2Translation, xquery_from, xquery_to)

                        if tmp_ali.getLength() == 0:
                            options.stdlog.write("# warning: empty alignment %s\n" % str(
                                (ref_from, exon_from, ref_to, exon_to, xquery_from, xquery_to)))
                            percent_identity = 0
                            percent_similarity = 0
                        else:
                            if options.loglevel >= 5:
                                print str(alignlib_lite.AlignmentFormatExplicit(tmp_ali, query_sequence, sbjct_sequence))

                            percent_identity = alignlib_lite.calculatePercentIdentity(tmp_ali,
                                                                                      query_sequence,
                                                                                      sbjct_sequence) * 100
                            percent_similarity = alignlib_lite.calculatePercentSimilarity(
                                tmp_ali) * 100

                    if percent_identity >= min_pide:
                        is_good_exon = 1
                    else:
                        is_good_exon = 0
                        dubious_exons += 1

                    # adjust regions for terminal exons
                    if e == 0 and r == 0 and dfrom <= (entry.mQueryFrom - 1) * 3 and dfrom > 0:
                        if is_good_exon:
                            truncated_Nterminal_exon = dfrom
                        dfrom = 0

                    # truncated terminal exons
                    if e == len(exons) - 1 and r == len(ref_exons) - 1 and dto <= (entry.mQueryLength - entry.mQueryTo) * 3 and dto > 0:
                        if is_good_exon:
                            truncated_Cterminal_exon = dto
                        dto = 0

                    # do not count deviations for terminal query exons
                    if e == 0 and dfrom <= entry.mQueryFrom * 3 and dfrom > 0:
                        dfrom = 0

                    if e == len(exons) - 1 and dto <= (entry.mQueryLength - entry.mQueryTo) * 3 and dto > 0:
                        dto = 0

                    # permit difference of one codon (assumed to be stop)
                    if e == len(exons) - 1 and r == len(ref_exons) - 1 and dto == 3:
                        dto = 0

                    # deal with different boundary conditions:
                    if dfrom == 0 and dto == 0:
                        if is_good_exon:
                            nidentical_exons += 1
                        e += 1
                        r += 1
                    # next exon within this ref_exon
                    elif exon_to < ref_to and next_exon_to and next_exon_to <= ref_to + options.slipping_exon_boundary:
                        if is_good_exon:
                            ninserted_introns += 1
                        e += 1
                        in_sync = 1
                        dto = 0
                    # next ref_exon within this exon
                    elif ref_to < exon_to and next_ref_to and next_ref_to <= exon_to + options.slipping_exon_boundary:
                        if is_good_exon:
                            ndeleted_introns += 1
                        r += 1
                        in_sync = 1
                        dto = 0
                    else:
                        e += 1
                        r += 1
                        if in_sync:
                            dfrom = 0

                    if is_good_exon:
                        exons_boundaries_sum += dfrom + dto
                        exons_boundaries_max = max(dfrom, exons_boundaries_max)
                        exons_boundaries_max = max(dto, exons_boundaries_max)

                    ###########################################################
                    # count inserted/deleted introns and misplaced boundaries
                    ##
                    # if exon and next_exon in ref_exon: inserted intron
                    # if ref_exon and next_ref_exon in exon: deleted intron
                if outfile_exons:

                    if genomic_fragment and exon_genome_to:
                        nintrons, nframeshifts, ngaps, nsplits, nstopcodons, disruptions = Genomics.CountGeneFeatures(exon_genome_from - entry.mSbjctGenomeFrom,
                                                                                                                      exon_ali,
                                                                                                                      genomic_fragment,
                                                                                                                      border_stop_codon=0
                                                                                                                      )
                    else:
                        nintrons, nframeshifts, ngaps, nsplits, nstopcodons = 0, 0, 0, 0, 0

                    if exon_to == 0:
                        this_e = 0
                    if ref_to == 0:
                        this_r = 0
                    outfile_exons.write(string.join(map(str, (entry.mPredictionId,
                                                              this_e, exon_from, exon_to, exon_phase,
                                                              this_r, ref_from, ref_to, ref_phase,
                                                              percent_identity, percent_similarity,
                                                              nframeshifts, ngaps, nstopcodons,
                                                              is_good_exon,
                                                              exon_genome_from, exon_genome_to,
                                                              )), "\t") + "\n")

            while e < len(exons):
                exon_from, exon_to, exon_phase, exon_genome_from, exon_genome_to = exons[
                    e][0:5]
                e += 1
                ninserted_Cexons += 1

                if outfile_exons:
                    outfile_exons.write(string.join(map(str, (entry.mPredictionId,
                                                              e, exon_from, exon_to, exon_phase,
                                                              0, 0, 0, 0,
                                                              0, 0,
                                                              0, 0, 0,
                                                              1,
                                                              exon_genome_from, exon_genome_to,
                                                              )), "\t") + "\n")

            while r < len(ref_exons):
                ref_from, ref_to, ref_phase, ref_genome_from, ref_genome_to = (ref_exons[r].mPeptideFrom,
                                                                               ref_exons[
                                                                                   r].mPeptideTo,
                                                                               ref_exons[
                                                                                   r].frame,
                                                                               ref_exons[
                                                                                   r].mGenomeFrom,
                                                                               ref_exons[r].mGenomeTo)
                ndeleted_Cexons += 1
                ref_genome_from -= ref_exons_offset
                ref_genome_to -= ref_exons_offset
                r += 1
                if outfile_exons:
                    outfile_exons.write(string.join(map(str, (entry.mPredictionId,
                                                              0, 0, 0, 0,
                                                              r, ref_from, ref_to, ref_phase,
                                                              0, 0,
                                                              0, 0, 0,
                                                              0,
                                                              0, 0,
                                                              )), "\t") + "\n")
        else:
            if options.write_notfound:
                this_e = 0
                # use prediction's identity/similarity for exons.
                # This will still then flag stop-codons in later analysis
                percent_identity = entry.mPercentIdentity
                percent_similarity = entry.mPercentSimilarity

                for exon in exons:
                    this_e += 1
                    exon_from, exon_to, exon_phase, exon_genome_from, exon_genome_to, exon_ali = exon[
                        0:6]
                    if genomic_fragment:
                        nintrons, nframeshifts, ngaps, nsplits, nstopcodons, disruptions = Genomics.CountGeneFeatures(exon_genome_from - entry.mSbjctGenomeFrom,
                                                                                                                      exon_ali,
                                                                                                                      genomic_fragment)

                    outfile_exons.write(string.join(map(str, (entry.mPredictionId,
                                                              this_e, exon_from, exon_to, exon_phase,
                                                              0, 0, 0, 0,
                                                              percent_identity, percent_similarity,
                                                              nframeshifts, ngaps, nstopcodons,
                                                              1,
                                                              exon_genome_from, exon_genome_to,
                                                              )), "\t") + "\n")

        options.stdout.write("\t".join(map(str,
                                           (entry.mPredictionId,
                                            exons_num_exons,
                                            dubious_exons,
                                            exons_boundaries_sum,
                                            exons_boundaries_max,
                                            nidentical_exons,
                                            ninserted_exons, ndeleted_exons,
                                            ninserted_introns, ndeleted_introns,
                                            truncated_Nterminal_exon, truncated_Cterminal_exon,
                                            ndeleted_Nexons, ndeleted_Cexons,
                                            ninserted_Nexons, ninserted_Cexons))) + "\n")

    if outfile_exons:
        outfile_exons.close()

    E.info("found=%i, missed_exons=%i, missed_length=%i, empty_alis=%i" % (nfound,
                                                                           nmissed_exons,
                                                                           nmissed_length,
                                                                           nempty_alignments))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
