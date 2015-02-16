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
gpipe/gff2predictions.py - convert a gff or exons file to gpipe predictions
=====================================================================

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

   python gpipe/gff2predictions.py --help

Type::

   python gpipe/gff2predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Predictor2 as Predictor

USAGE = """python %s [OPTIONS] < psl > predictions

Convert GFF exon list to predictions format.

Version: $Id: gpipe/gff2predictions.py 2021 2008-07-10 16:00:48Z andreas $
"""


def checkIdentity(reference, translation, options):
    """check if the two sequences in reference and translation
    are (nearly) identical.

    Permitted deviations:
       * reference sequence is one residue longer
    """

    is_identical = reference == translation

    if not is_identical:
        # check for masked characters (stops)
        if abs(len(translation) == len(reference)) <= 1:
            nmismatch = 0
            # permit the first residue to be different, so start at residue 2
            for x in range(1, min(len(translation), len(reference))):
                if reference[x] != translation[x] and \
                   translation[x] != "X" and reference[x] != "U":
                    if options.loglevel >= 2:
                        options.stdlog.write(
                            "# residue mismatch at position %i: %s %s\n" %
                            (x, reference[x], translation[x]))
                    nmismatch += 1
                    if nmismatch > 10:
                        break

            return nmismatch == 0, nmismatch
        else:
            if options.loglevel >= 2:
                options.stdlog.write("# length mismatch: %i - %i\n" %
                                     (len(translation), len(reference)))
            return False, 0

    return is_identical, 0


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/gff2predictions.py 2021 2008-07-10 16:00:48Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-t", "--trans", dest="trans",
                      help="input is translated DNA.", action="store_true")

    parser.add_option("-f", "--format", dest="format",
                      help="input format.", type="choice",
                      choices=("exons", "psl", "gff"))

    parser.add_option("-o", "--output-format", dest="output_format",
                      help="output format", type="choice",
                      choices=('exontable', 'exons', 'predictions', 'cds', 'fasta'))

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genomic data (indexed).")

    parser.add_option("--predictions-file", dest="predictions_file", type="string",
                      help="filename with predictions. Use gene structures from this file if available.")

    parser.add_option("-i", "--gff-field-id", dest="gff_field_id", type="string",
                      help="field for the feature id in the gff info section.")

    parser.add_option("-p", "--filename-peptides", dest="filename_peptides", type="string",
                      help="Filename with peptide sequences. If given, it is used to check the predicted translated sequences.")

    parser.add_option("--no-realignment", dest="do_realignment", action="store_false",
                      help="do not re-align entries that do not parse correctly.")

    parser.add_option("--remove-unaligned", dest="remove_unaligned", action="store_true",
                      help="remove entries that have not been aligned correctly.")

    parser.add_option("--input-coordinates", dest="input_coordinates", type="string",
                      help="specify input format for input coordinates [forward|both-zero|one-closed|open].")

    parser.set_defaults(
        trans=False,
        output_format="predictions",
        format="psl",
        gff_field_id='id',
        input_coordinates="both-zero-open",
        filename_peptides=None,
        genome_file=None,
        do_realignment=True,
        predictions_file=None,
        remove_unaligned=False)

    (options, args) = E.Start(parser)

    if not options.genome_file:
        raise "please specify a genome file."

    fasta = IndexedFasta.IndexedFasta(options.genome_file)
    contig_sizes = fasta.getContigSizes()

    ninput, noutput, nskipped = 0, 0, 0
    nfound, nnotfound, nidentical, nmismatch, naligned, nunaligned = 0, 0, 0, 0, 0, 0

    if options.filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            IOTools.openFile(options.filename_peptides, "r"))
        predictor = Predictor.PredictorExonerate()
        predictor.mLogLevel = 0
    else:
        peptide_sequences = None
        predictor = None

    converter = IndexedFasta.getConverter(options.input_coordinates)

    predictions = {}
    if options.predictions_file:
        parser = PredictionParser.iterator_predictions(
            IOTools.openFile(options.predictions_file, "r"))
        for p in parser:
            predictions[p.mPredictionId] = p

    if options.output_format == "predictions":

        if options.format == "psl":

            if options.trans:
                parser = PredictionParser.PredictionParserBlatTrans()
            else:
                parser = PredictionParser.PredictionParserBlatCDNA()

            nmatches = 1
            for line in sys.stdin:
                if line[0] == "#":
                    continue
                if not re.match("^[0-9]", line):
                    continue

                try:
                    entries = parser.Parse((line,))
                except PredictionParser.AlignmentError, e:
                    print "# %s" % str(e)
                    print "#", line[:-1]
                    sys.exit(1)

                for entry in entries:
                    entry.mPredictionId = nmatches
                    nmatches += 1

                print str(entries)

        elif options.format == "exons":
            parser = PredictionParser.PredictionParserExons(
                contig_sizes=contig_sizes)
        else:
            raise"unknown format %s for output option %s" % (
                options.format, options.output_format)

        if options.loglevel >= 2:
            options.stdlog.write("# parsing.\n")
            options.stdlog.flush()

        results = parser.Parse(sys.stdin.readlines())

        if options.loglevel >= 2:
            options.stdlog.write("# parsing finished.\n")
            options.stdlog.flush()

        if options.loglevel >= 1:
            options.stdlog.write("# parsing: ninput=%i, noutput=%i, nerrors=%i\n" % (
                parser.GetNumInput(), parser.GetNumOutput(), parser.GetNumErrors()))

            for error, msg in parser.mErrors:
                options.stdlog.write("# %s : %s\n" % (str(error), msg))
                options.stdlog.flush()

        # if genomes are given: build translation
        if options.genome_file:

            results.Sort(lambda x, y: cmp(x.mSbjctToken, y.mSbjctToken))

            new_results = PredictionParser.Predictions()

            for entry in results:

                ninput += 1

                if options.loglevel >= 2:
                    options.stdlog.write("# processing entry %s:%s on %s:%s %i/%i.\n" % (entry.mPredictionId,
                                                                                         entry.mQueryToken,
                                                                                         entry.mSbjctToken,
                                                                                         entry.mSbjctStrand,
                                                                                         ninput, len(results)))
                    options.stdlog.flush()

                try:
                    lgenome = fasta.getLength(entry.mSbjctToken)
                    # added 3 residues - was a problem at split codons just before the stop.
                    # See for example the chicken sequence ENSGALP00000002741
                    genomic_sequence = fasta.getSequence(entry.mSbjctToken,
                                                         entry.mSbjctStrand,
                                                         entry.mSbjctGenomeFrom,
                                                         min(entry.mSbjctGenomeTo + 3, lgenome))

                except KeyError:
                    if options.loglevel >= 1:
                        options.stdlog.write("# did not find entry for %s on %s.\n" % (
                            entry.mPredictionId, entry.mSbjctToken))
                    nskipped += 1
                    continue

                if predictions and entry.mPredictionId in predictions:
                    if options.loglevel >= 2:
                        options.stdlog.write("# substituting entry %s on %s:%s.\n" % (entry.mPredictionId,
                                                                                      entry.mSbjctToken,
                                                                                      entry.mSbjctStrand))
                        options.stdlog.flush()
                    entry = predictions[entry.mPredictionId]

                exons = Exons.Alignment2Exons(
                    entry.mMapPeptide2Genome, 0, entry.mSbjctGenomeFrom)

                entry.mMapPeptide2Translation, entry.mTranslation = Genomics.Alignment2PeptideAlignment(
                    Genomics.String2Alignment(entry.mAlignmentString), entry.mQueryFrom, 0, genomic_sequence)

                entry.score = entry.mMapPeptide2Translation.getColTo(
                ) - entry.mMapPeptide2Translation.getColFrom() + 1

                (entry.mNIntrons, entry.mNFrameShifts, entry.mNGaps, entry.mNSplits, entry.mNStopCodons, entry.mNDisruptions ) = \
                    Genomics.CountGeneFeatures(0,
                                               entry.mMapPeptide2Genome,
                                               genomic_sequence)

                if peptide_sequences:

                    if str(entry.mPredictionId) in peptide_sequences:

                        reference = peptide_sequences[
                            str(entry.mPredictionId)].upper()

                        translation = entry.mTranslation
                        nfound += 1

                        is_identical, nmismatches = checkIdentity(
                            reference, translation, options)

                        if is_identical:
                            nidentical += 1
                        else:
                            nmismatch += 1

                            if options.do_realignment:
                                if options.loglevel >= 2:
                                    options.stdlog.write("# %s: mismatches..realigning in region %i:%i\n" % (
                                        entry.mPredictionId, entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo))
                                    options.stdlog.flush()

                                    result = predictor(entry.mPredictionId, reference,
                                                       entry.mSbjctToken, genomic_sequence,
                                                       "--subopt FALSE --score '%s'" % str(80))
                                    # "--exhaustive --subopt FALSE --score '%s'" % str(80) )

                                    if result:
                                        translation = result[0].mTranslation
                                        is_identical, nmismatches = checkIdentity(
                                            reference, translation, options)
                                    else:
                                        if options.loglevel >= 2:
                                            options.stdlog.write(
                                                "# %s: realignment returned empty result\n" % (entry.mPredictionId))
                                            options.stdlog.flush()
                                        is_identical = False

                                    if is_identical:
                                        naligned += 1
                                        prediction_id = entry.mPredictionId
                                        sbjct_genome_from = entry.mSbjctGenomeFrom
                                        entry = result[0]
                                        entry.mPredictionId = prediction_id
                                        entry.mSbjctGenomeFrom += sbjct_genome_from
                                    else:
                                        nunaligned += 1
                                        if options.loglevel >= 1:
                                            options.stdlog.write("# %s: mismatch on %s:%s:%i-%i after realignment\n# reference =%s\n# translated=%s\n# realigned =%s\n" %
                                                                 (entry.mPredictionId,
                                                                  entry.mSbjctToken,
                                                                  entry.mSbjctStrand,
                                                                  entry.mSbjctGenomeFrom,
                                                                  entry.mSbjctGenomeTo,
                                                                  reference,
                                                                  entry.mTranslation,
                                                                  translation))
                                            options.stdlog.flush()
                                        if options.remove_unaligned:
                                            nskipped += 1
                                            continue

                            else:
                                if options.loglevel >= 2:
                                    options.stdlog.write("# %s: mismatches on %s ... no realignment\n" %
                                                         (entry.mPredictionId,
                                                          entry.mSbjctToken,))
                                    if options.loglevel >= 3:
                                        options.stdlog.write("# %s: mismatch before realignment\n# reference =%s\n# translated=%s\n" %
                                                             (entry.mPredictionId,
                                                              reference,
                                                              translation))
                                    options.stdlog.flush()

                                if options.remove_unaligned:
                                    nskipped += 1
                                    continue

                    else:
                        nnotfound += 1

                new_results.append(entry)
                noutput += 1

            results = new_results
        if results:
            options.stdout.write(str(results) + "\n")

    elif options.output_format == "exontable":
        if options.format == "exons":
            exons = Exons.ReadExonBoundaries(sys.stdin, contig_sizes=contig_sizes,
                                             delete_missing=True)
        else:
            raise "unknown format."

        for k in exons.keys():
            ee = exons[k]

            id = 0
            for e in ee:
                id += 1
                print "\t".join(map(str, (e.mQueryToken,
                                          id,
                                          e.mPeptideFrom, e.mPeptideTo, e.frame,
                                          0, 0, 0, 0,
                                          0, 0,
                                          0, 0, 0,
                                          0,
                                          e.mGenomeFrom, e.mGenomeTo)))

    elif options.output_format == "exons":

        if options.format == "exons":
            parser = PredictionParser.PredictionParserExons(
                contig_sizes=contig_sizes)
        else:
            raise "unknown format %s." % options.format

        results = parser.Parse(sys.stdin.readlines())
        id = 0
        for entry in results:
            exons = Exons.Alignment2Exons(entry.mMapPeptide2Genome,
                                          entry.mQueryFrom,
                                          entry.mSbjctGenomeFrom, )

            for e in exons:
                id += 1
                print "\t".join(map(str, (entry.mQueryToken,
                                          entry.mSbjctToken, entry.mSbjctStrand,
                                          e.frame,
                                          e.mRank,
                                          e.mPeptideFrom, e.mPeptideTo,
                                          e.mGenomeFrom, e.mGenomeTo)))

    elif options.output_format == "cds":

        if options.format == "exons":

            exons = Exons.ReadExonBoundaries(sys.stdin,
                                             contig_sizes=contig_sizes,
                                             delete_missing=True)

            # sort by chromosome
            kk = exons.keys()
            kk.sort(
                lambda x, y: cmp(exons[x][0].mSbjctToken, exons[y][0].mSbjctToken))

            id = 0

            for k in kk:

                ee = exons[k]

                # attach sequence to cds
                for e in ee:
                    e.mSequence = fasta.getSequence(e.mSbjctToken, e.mSbjctStrand,
                                                    e.mGenomeFrom, e.mGenomeTo)

                exon_from = 0
                for e in ee:
                    id += 1
                    e.mPeptideFrom = exon_from
                    exon_from += e.mGenomeTo - e.mGenomeFrom
                    e.mPeptideTo = exon_from
                    print string.join(map(str, (id, e.mQueryToken,
                                                e.mPeptideFrom, e.mPeptideTo,
                                                e.frame,
                                                e.mGenomeFrom, e.mGenomeTo,
                                                e.mSequence
                                                )), "\t")

    elif options.output_format == "fasta":

        if options.format == "gff":
            gff_entries = GTF.readFromFile(sys.stdin)
            n = 0
            ninput = len(gff_entries)
            for e in gff_entries:

                sequence = fasta.getSequence(e.name, e.strand,
                                             e.start, e.end,
                                             converter)

                n += 1
                try:
                    id = e.fields[options.gff_field_id]
                except KeyError:
                    nskipped += 1
                    continue

                noutput += 1
                options.stdout.write(">%s %s:%s:%s:%s\n%s\n" % (e.fields[options.gff_field_id],
                                                                e.name, e.strand,
                                                                e.start, e.end, sequence))
        else:
            raise"unknown format %s for output option %s" % (
                options.format, options.output_format)

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))
        if peptide_sequences:
            options.stdlog.write("# nfound=%i, nnotfound=%i, nidentical=%i, nmismatch=%i, naligned=%i, nunaligned=%i\n" %
                                 (nfound, nnotfound, nidentical, nmismatch, naligned, nunaligned))
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
