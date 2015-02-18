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
gpipe/predictions2cds.py -
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

   python gpipe/predictions2cds.py --help

Type::

   python gpipe/predictions2cds.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import alignlib_lite
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Exons as Exons
import CGAT.PredictionParser as PredictionParser

USAGE = """python %s < predictions > genes

Version: $Id: gpipe/predictions2cds.py 1858 2008-05-13 15:07:05Z andreas $

extract cds or introns from predicions into various formats.

Output formats:

gff:
exons:        exon table output

""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-g", "--genome-file", dest="genome_file", type="string",
        help="filename with genome.")

    parser.add_option(
        "-o", "--is-forward-coordinates", dest="forward_coordinates",
        action="store_true",
        help="input uses forward coordinates.")

    parser.add_option(
        "-f", "--format", dest="format", type="choice",
        choices=(
            "default", "cds", "cdnas", "map", "gff", "intron-fasta", "exons"),
        help="output format.")

    parser.add_option(
        "-r", "--reset-to-start", dest="reset_to_start", action="store_true",
        help="move genomic coordinates to begin from 0.")

    parser.add_option("--reset-query", dest="reset_query", action="store_true",
                      help="move peptide coordinates to begin from 0.")

    parser.set_defaults(
        genome_file=None,
        forward_coordinates=False,
        format="default",
        reset_to_start=False,
        reset_query=False)

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    cds_id = 1

    entry = PredictionParser.PredictionParserEntry()

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    ninput, noutput, nskipped, nerrors = 0, 0, 0, 0

    for line in sys.stdin:

        if line[0] == "#":
            continue
        if line.startswith("id"):
            continue

        ninput += 1

        try:
            entry.Read(line)
        except ValueError, msg:
            options.stdlog.write(
                "# parsing failed with msg %s in line %s" % (msg, line))
            nerrors += 1
            continue

        cds = Exons.Alignment2Exons(entry.mMapPeptide2Genome,
                                    query_from=entry.mQueryFrom,
                                    sbjct_from=entry.mSbjctGenomeFrom,
                                    add_stop_codon=0)

        for cd in cds:
            cd.mSbjctToken = entry.mSbjctToken
            cd.mSbjctStrand = entry.mSbjctStrand

        if cds[-1].mGenomeTo != entry.mSbjctGenomeTo:
            options.stdlog.write(
                "# WARNING: discrepancy in exon calculation!!!\n")
            for cd in cds:
                options.stdlog.write("# %s\n" % str(cd))
            options.stdlog.write("# %s\n" % entry)

        lsequence = fasta.getLength(entry.mSbjctToken)
        genomic_sequence = fasta.getSequence(entry.mSbjctToken,
                                             entry.mSbjctStrand,
                                             entry.mSbjctGenomeFrom,
                                             entry.mSbjctGenomeTo)

        # deal with forward coordinates: convert them to negative strand
        # coordinates
        if options.forward_coordinates and \
                entry.mSbjctStrand == "-":
            entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo = lsequence - \
                entry.mSbjctGenomeTo, lsequence - entry.mSbjctGenomeFrom
            for cd in cds:
                cd.InvertGenomicCoordinates(lsequence)

        # attach sequence to cds
        for cd in cds:
            start = cd.mGenomeFrom - entry.mSbjctGenomeFrom
            end = cd.mGenomeTo - entry.mSbjctGenomeFrom
            cd.mSequence = genomic_sequence[start:end]

        # reset coordinates for query
        if options.reset_to_start:
            offset = entry.mPeptideFrom
            for cd in cds:
                cd.mPeptideFrom -= offset
                cd.mPeptideTo -= offset

        # play with coordinates
        if options.reset_to_start:
            offset = entry.mSbjctGenomeFrom
            for cd in cds:
                cd.mGenomeFrom -= offset
                cd.mGenomeTo -= offset
        else:
            offset = 0

        if options.format == "cds":
            rank = 0
            for cd in cds:
                rank += 1
                cd.mQueryToken = entry.mQueryToken
                cd.mSbjctToken = entry.mSbjctToken
                cd.mSbjctStrand = entry.mSbjctStrand
                cd.mRank = rank
                print str(cd)

        if options.format == "exons":
            rank = 0
            for cd in cds:
                rank += 1
                options.stdout.write("\t".join(map(str, (entry.mPredictionId,
                                                         cd.mSbjctToken,
                                                         cd.mSbjctStrand,
                                                         rank,
                                                         cd.frame,
                                                         cd.mPeptideFrom,
                                                         cd.mPeptideTo,
                                                         cd.mGenomeFrom,
                                                         cd.mGenomeTo))) + "\n")

        elif options.format == "cdnas":
            print string.join(map(str, (entry.mPredictionId,
                                        entry.mQueryToken,
                                        entry.mSbjctToken,
                                        entry.mSbjctStrand,
                                        entry.mSbjctGenomeFrom - offset,
                                        entry.mSbjctGenomeTo - offset,
                                        genomic_sequence)), "\t")

        elif options.format == "map":

            map_prediction2genome = alignlib_lite.makeAlignmentSet()

            for cd in cds:
                alignlib_lite.addDiagonal2Alignment(map_prediction2genome,
                                                    cd.mPeptideFrom + 1,
                                                    cd.mPeptideTo,
                                                    (cd.mGenomeFrom - offset) - cd.mPeptideFrom)

            print string.join(map(str, (entry.mPredictionId,
                                        entry.mSbjctToken,
                                        entry.mSbjctStrand,
                                        alignlib_lite.AlignmentFormatEmissions(map_prediction2genome))), "\t")

        elif options.format == "intron-fasta":
            rank = 0
            if len(cds) == 1:
                nskipped += 1
                continue

            last = cds[0].mGenomeTo
            for cd in cds[1:]:
                rank += 1
                key = "%s %i %s:%s:%i:%i" % (
                    entry.mPredictionId, rank, entry.mSbjctToken, entry.mSbjctStrand, last, entry.mSbjctGenomeFrom)
                sequence = genomic_sequence[
                    last - entry.mSbjctGenomeFrom:cd.mGenomeFrom - entry.mSbjctGenomeFrom]
                options.stdout.write(">%s\n%s\n" % (key, sequence))
                last = cd.mGenomeTo

        elif options.format == "gff-match":
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tTarget \"%s\" %i %i; Score %i; Introns %i; Frameshifts %i; Stops %i" % \
                  (entry.mSbjctToken,
                   "gpipe", "similarity",
                   entry.mSbjctGenomeFrom,
                   entry.mSbjctGenomeTo,
                   entry.mPercentIdentity,
                   entry.mSbjctStrand,
                   ".",
                   entry.mQueryToken,
                   entry.mQueryFrom,
                   entry.mQueryTo,
                   entry.score,
                   entry.mNIntrons,
                   entry.mNFrameShifts,
                   entry.mNStopCodons)

        elif options.format == "gff-exon":
            rank = 0
            for cd in cds:
                rank += 1
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tTarget \"%s\" %i %i; Score %i; Rank %i/%i; Prediction %i" % \
                      (entry.mSbjctToken,
                       "gpipe", "similarity",
                       cd.mGenomeFrom,
                       cd.mGenomeTo,
                       entry.mPercentIdentity,
                       entry.mSbjctStrand,
                       ".",
                       entry.mQueryToken,
                       cd.mPeptideFrom / 3 + 1,
                       cd.mPeptideTo / 3 + 1,
                       entry.score,
                       rank,
                       len(cds),
                       entry.mPredictionId)
        else:
            exon_from = 0
            for cd in cds:
                cd.mPeptideFrom = exon_from
                exon_from += cd.mGenomeTo - cd.mGenomeFrom
                cd.mPeptideTo = exon_from
                print string.join(map(str, (cds_id, entry.mPredictionId,
                                            cd.mPeptideFrom, cd.mPeptideTo,
                                            cd.frame,
                                            cd.mGenomeFrom, cd.mGenomeTo,
                                            cd.mSequence
                                            )), "\t")
                cds_id += 1

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i, nerrors=%i\n" % (
            ninput, noutput, nskipped, nerrors))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
