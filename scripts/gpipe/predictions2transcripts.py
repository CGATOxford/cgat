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
gpipe/predictions2transcripts.py - patch predictions into transcripts
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

attempt various ways to patch up predictions like

  * extending predictions to next start AUG/stopcodon.
  * filling introns (introns without frameshift and stop codon)

Usage
-----

Example::

   python gpipe/predictions2transcripts.py --help

Type::

   python gpipe/predictions2transcripts.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Exons as Exons
import CGAT.Prediction as Prediction
import CGAT.PredictionParser as PredictionParser
import CGAT.Stats as Stats


def findCodonReverse(sequence, start, found_codons, abort_codons=None):
    """find codon by tracking along sequence from start.

    This procedure will stop at completely masked codons and abort_codons

    """

    found = False

    if abort_codons:
        acodons = abort_codons + ("NNN", "XXX")
    else:
        acodons = ("NNN", "XXX")

    while start > 0 and genomic_sequence[start:start + 3] not in acodons:

        if genomic_sequence[start:start + 3] in found_codons:
            found = True
            break

        start -= 3

    return found, start


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: gpipe/predictions2transcripts.py 1841 2008-05-08 12:07:13Z andreas $",
                            usage=globals()["__doc__"])
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("-o", "--output-filename-summary", dest="output_filename_summary", type="string",
                      help="filename with summary information.")

    parser.add_option("--skip-header", dest="skip_header", action="store_true",
                      help="skip header.")

    parser.add_option("--start-codon-boundary", dest="start_codon_boundary", type="int",
                      help="maximum extension for start codon (make divisible by 3).")

    parser.add_option("--stop-codon-boundary", dest="stop_codon_boundary", type="int",
                      help="maximum extension for stop codon (make divisible by 3).")

    parser.add_option("--left-extension-mode", dest="left_extension_mode", type="choice",
                      choices=("first-start", "first-stop-backtrack"),
                      help="extension mode for 5' end.")

    parser.add_option("--fill-introns", dest="fill_introns", type="int",
                      help="fill intron if divisible by three and no stop codon up to a maximum length of #.")

    parser.add_option("--introns-max-stops", dest="introns_max_stops", type="int",
                      help="maximum number of stop codons to tolerate within an intron.")

    parser.add_option("--output-format", dest="output_format", type="choice",
                      choices=("predictions", "extensions", "filled-introns"),
                      help="output format.")

    parser.set_defaults(
        genome_file="genome",
        start_codons=("ATG"),
        stop_codons=("TAG", "TAA", "TGA"),
        start_codon_boundary=9999,
        stop_codon_boundary=9999,
        fill_introns=0,
        introns_max_stops=0,
        left_splice_signals=("GT",),
        right_splice_signals=("AG",),
        output_format="extensions",
        left_extension_mode="first-start",
        skip_header=False,
        output_filename_summary=None,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    options.start_codon_boundary = int(options.start_codon_boundary / 3)
    options.stop_codon_boundary = int(options.stop_codon_boundary / 3)

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    p = PredictionParser.PredictionParserEntry()

    ninput, noutput = 0, 0
    nfilled = 0
    nseqs_filled = 0
    nseqs_extended = 0
    left_extensions = []
    right_extensions = []
    filled_introns = []

    if not options.skip_header:
        if options.output_format == "predictions":
            options.stdout.write(Prediction.Prediction().getHeader() + "\n")
        elif options.output_format == "filled-introns":
            options.stdout.write("\t".join(("prediction_id",
                                            "intron",
                                            "peptide_sequence",
                                            "genomic_sequence")) + "\n")

    if options.output_filename_summary:
        outfile_summary = open(options.output_filename_summary, "w")
        outfile_summary.write("id\ttype\tnumber\tlength\tfrom\tto\tsequence\n")
    else:
        outfile_summary = None

    for line in options.stdin:

        if line[0] == "#":
            continue

        ninput += 1
        p.Read(line)

        lsequence = fasta.getLength(p.mSbjctToken)

        genome_from = max(0, p.mSbjctGenomeFrom - options.start_codon_boundary)
        genome_to = min(
            lsequence, p.mSbjctGenomeTo + options.stop_codon_boundary)

        genomic_sequence = fasta.getSequence(p.mSbjctToken, p.mSbjctStrand,
                                             genome_from,
                                             genome_to).upper()

        #######################################################################
        #######################################################################
        #######################################################################
        # Do extensions

        if options.start_codon_boundary or options.stop_codon_boundary:

            extension_start = p.mSbjctGenomeFrom - genome_from
            extension_stop = genome_to - p.mSbjctGenomeTo

            fragment_to = extension_start + \
                p.mSbjctGenomeTo - p.mSbjctGenomeFrom

            lfragment = len(genomic_sequence)

            ###################################################################
            ###################################################################
            ###################################################################
            # find start codon
            start = extension_start
            found_start = False
            if options.left_extension_mode == "first-start":

                found_start, start = findCodonReverse(genomic_sequence,
                                                      start,
                                                      options.start_codons,
                                                      options.stop_codons)

            elif options.left_extension_mode == "first-stop-backtrack":

                if genomic_sequence[start:start + 3] in options.start_codons:
                    found_start = True
                else:
                    found_start, start = findCodonReverse(genomic_sequence,
                                                          start,
                                                          options.stop_codons)

                    if found_start:
                        E.info("prediction %s: stop found at %i (%i) backtracking ..." % (
                            p.mPredictionId, start, extension_start - start))

                        # bracktrack to first start codon
                        found_start = False
                        while start < extension_start:
                            start += 3
                            if genomic_sequence[start:start + 3] in options.start_codons:
                                found_start = True
                                break
                        else:
                            start = extension_start

                        if found_start:
                            E.info("start codon found at %i (%i)." %
                                   (start, extension_start - start))
                        else:
                            E.info("no start codon found.")
                    else:
                        E.info("prediction %s: no stop found ... backtracking to start codon." % (
                            p.mPredictionId))

                        found_start, start = findCodonReverse(
                            genomic_sequence, start, options.start_codons)

                        E.info("prediction %s: no start codon found." %
                               (p.mPredictionId))

            if found_start:
                start += genome_from
            else:
                start = p.mSbjctGenomeFrom

            dstart = p.mSbjctGenomeFrom - start

            ###################################################################
            ###################################################################
            ###################################################################
            # find stop codon
            # stop points to the beginning of the codon, thus the stop codon will
            # not be part of the sequence.
            stop = fragment_to
            found_stop = 0
            while stop < lfragment and \
                    genomic_sequence[stop:stop + 3] not in ("NNN", "XXX"):
                if genomic_sequence[stop:stop + 3] in options.stop_codons:
                    found_stop = 1
                    break

                stop += 3

            if found_stop:
                stop += genome_from
            else:
                stop = p.mSbjctGenomeTo

            dstop = stop - p.mSbjctGenomeTo

            ###################################################################
            ###################################################################
            ###################################################################
            # build new prediction
            map_peptide2genome = []
            if dstart:
                map_peptide2genome.append(("G", 0, dstart))
            map_peptide2genome += p.mMapPeptide2Genome
            if dstop:
                map_peptide2genome.append(("G", 0, dstop))

            E.info("prediction %s: extension: found_start=%i, found_stop=%i, left=%i, right=%i" % (
                p.mPredictionId, found_start, found_stop, dstart, dstop))

            # save results
            p.mMapPeptide2Genome = map_peptide2genome
            p.mAlignmentString = Genomics.Alignment2String(map_peptide2genome)
            p.mSbjctGenomeFrom -= dstart
            p.mSbjctGenomeTo += dstop
            p.mSbjctFrom += dstart / 3
            p.mSbjctTo += dstart / 3 + dstop / 3

            if dstart or dstop:
                if dstart:
                    left_extensions.append(dstart)
                if dstop:
                    right_extensions.append(dstop)

                nseqs_extended += 1

        # update genomic sequence because borders might have changed.
        genomic_sequence = fasta.getSequence(p.mSbjctToken,
                                             p.mSbjctStrand,
                                             p.mSbjctGenomeFrom,
                                             p.mSbjctGenomeTo).upper()

        if options.fill_introns:

            has_filled = False

            exons = Exons.Alignment2Exons(p.mMapPeptide2Genome,
                                          query_from=0,
                                          sbjct_from=0)

            new_exons = []

            last_e = exons[0]

            nintron = 0

            for e in exons[1:]:

                nintron += 1
                lintron = e.mGenomeFrom - last_e.mGenomeTo

                if lintron > options.fill_introns or (lintron) % 3 != 0:
                    E.debug("prediction %s: intron %i of size %i discarded." %
                            (p.mPredictionId,
                             nintron, lintron))

                    new_exons.append(last_e)
                    last_e = e
                    continue

                # get sequence, include also residues from split codons
                # when checking for stop codons.
                if e.mAlignment[0][0] == "S":
                    offset_left = last_e.mAlignment[-1][2]
                    offset_right = e.mAlignment[0][2]
                else:
                    offset_left, offset_right = 0, 0

                sequence = genomic_sequence[
                    last_e.mGenomeTo - offset_left:e.mGenomeFrom + offset_right]

                # check for splice sites
                for signal in options.left_splice_signals:
                    if sequence[offset_left:offset_left + len(signal)] == signal:
                        left_signal = True
                        break
                else:
                    left_signal = False

                for signal in options.right_splice_signals:
                    if sequence[-(len(signal) + offset_right):-offset_right] == signal:
                        right_signal = True
                        break
                else:
                    right_signal = False

                nstops, ngaps = 0, 0
                for codon in [sequence[x:x + 3] for x in range(0, len(sequence), 3)]:
                    if codon in options.stop_codons:
                        nstops += 1
                    if "N" in codon.upper():
                        ngaps += 1

                    E.debug("prediction %s: intron %i of size %i (%i-%i) (%s:%s:%i:%i): stops=%i, gaps=%i, signals=%s,%s." %
                            (p.mPredictionId,
                             nintron, lintron,
                             offset_left, offset_right,
                             p.mSbjctToken, p.mSbjctStrand,
                             p.mSbjctGenomeFrom + last_e.mGenomeTo,
                             p.mSbjctGenomeFrom + e.mGenomeFrom,
                             nstops,
                             ngaps,
                             left_signal, right_signal))

                if nstops + ngaps > options.introns_max_stops:
                    new_exons.append(last_e)
                    last_e = e
                    continue

                E.info("prediction %s: filling intron %i of size %i: stops=%i, gaps=%i, signals=%s,%s" %
                       (p.mPredictionId,
                        nintron, lintron,
                        nstops,
                        ngaps,
                        left_signal, right_signal))

                e.Merge(last_e)
                has_filled = True
                nfilled += 1
                last_e = e

                if options.output_format == "filled-introns":
                    options.stdout.write("\t".join(map(str, (p.mPredictionId,
                                                             nintron,
                                                             Genomics.TranslateDNA2Protein(
                                                                 sequence),
                                                             sequence))) + "\n")

                filled_introns.append(lintron)
                p.mNIntrons -= 1

            new_exons.append(last_e)

            if has_filled:
                nseqs_filled += 1

            Exons.UpdatePeptideCoordinates(new_exons)

            p.mMapPeptide2Genome = Exons.Exons2Alignment(new_exons)
            p.mAlignmentString = Genomics.Alignment2String(
                p.mMapPeptide2Genome)

        # build translated sequence
        p.mMapPeptide2Translation, p.mTranslation = Genomics.Alignment2PeptideAlignment(
            p.mMapPeptide2Genome, p.mQueryFrom, 0, genomic_sequence)

        # output info
        if options.output_format == "predictions":
            options.stdout.write(str(p) + "\n")
        elif options.output_format == "extensions":
            if found_start:
                found_start = 1
            if found_stop:
                found_stop = 1
            options.stdout.write("\t".join(map(str, (p.mPredictionId,
                                                     found_start, found_stop,
                                                     dstart, dstop,
                                                     p.mTranslation,
                                                     p.mSbjctGenomeFrom, p.mSbjctGenomeTo,
                                                     p.mAlignmentString))) + "\n")

        noutput += 1
        options.stdout.flush()

    E.info("stats  : %s" %
           "\t".join(Stats.DistributionalParameters().getHeaders()))
    E.info("left   : %s" %
           str(Stats.DistributionalParameters(left_extensions)))
    E.info("right  : %s" %
           str(Stats.DistributionalParameters(right_extensions)))
    E.info("introns: %s" % str(Stats.DistributionalParameters(filled_introns)))
    E.info("ninput=%i, noutput=%i, nextended=%i, nfilled=%i, nexons_filled=%i" % (
        ninput, noutput, nseqs_extended, nseqs_filled, nfilled))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
