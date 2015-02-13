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
gpipe/predictions2introns.py -
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

   python gpipe/predictions2introns.py --help

Type::

   python gpipe/predictions2introns.py --help

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
import CGAT.PredictionParser as PredictionParser

USAGE = """python %s < predictions > introns

Version: $Id: gpipe/predictions2introns.py 2781 2009-09-10 11:33:14Z andreas $

Summarize information about introns in predictions.
""" % sys.argv[0]


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

    parser = E.OptionParser(version="%prog version: $Id: gpipe/predictions2introns.py 2781 2009-09-10 11:33:14Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("-o", "--output-filename-summary", dest="output_filename_summary", type="string",
                      help="filename with summary information.")

    parser.add_option("--skip-header", dest="skip_header", action="store_true",
                      help="skip header.")

    parser.add_option(
        "--fill-introns", dest="fill_introns", type="int",
        help="fill intron if divisible by three and no stop codon up to a "
        "maximum length of #.")

    parser.add_option(
        "--introns-max-stops", dest="introns_max_stops", type="int",
        help="maximum number of stop codons to tolerate within an intron.")

    parser.add_option("--output-format", dest="output_format", type="choice",
                      choices=("predictions", "extensions", "filled-introns"),
                      help="output format.")

    parser.set_defaults(
        genome_file="genome",
        start_codons=("ATG"),
        stop_codons=("TAG", "TAA", "TGA"),
        skip_header=False,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

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
        options.stdout.write("\t".join(("prediction_id",
                                        "intron",
                                        "contig",
                                        "strand",
                                        "start",
                                        "end",
                                        "length",
                                        "nstops",
                                        "type",
                                        "prime5",
                                        "prime3",
                                        )) + "\n")

    for line in sys.stdin:

        if line[0] == "#":
            continue

        ninput += 1
        p.Read(line)

        lsequence = fasta.getLength(p.mSbjctToken)

        genomic_sequence = fasta.getSequence(p.mSbjctToken, p.mSbjctStrand,
                                             p.mSbjctGenomeFrom,
                                             p.mSbjctGenomeTo).upper()

        exons = Exons.Alignment2Exons(p.mMapPeptide2Genome,
                                      query_from=0,
                                      sbjct_from=0)

        new_exons = []

        last_e = exons[0]

        nintron = 0

        for e in exons[1:]:

            nintron += 1
            lintron = e.mGenomeFrom - last_e.mGenomeTo

            intron_is_l3 = lintron % 3 != 0

            if intron_is_l3:
                # get sequence, include also residues from split codons
                # when checking for stop codons.
                # note that e.mAlignment can sometimes be empty. This might
                # be an exonerate bug. In the alignment string there are two
                # consecutive exons.
                if e.mAlignment and last_e.mAlignment and e.mAlignment[0][0] == "S":
                    offset_left = last_e.mAlignment[-1][2]
                    offset_right = e.mAlignment[0][2]
                else:
                    offset_left, offset_right = 0, 0

                sequence = genomic_sequence[
                    last_e.mGenomeTo - offset_left:e.mGenomeFrom + offset_right]

                intron_nstops = 0
                for codon in [sequence[x:x + 3] for x in range(0, len(sequence), 3)]:
                    if codon in options.stop_codons:
                        intron_nstops += 1
            else:
                intron_nstops = 0

            # check for splice signals
            sequence = genomic_sequence[last_e.mGenomeTo:e.mGenomeFrom]

            intron_type, prime5, prime3 = Genomics.GetIntronType(sequence)

            if options.loglevel >= 2:
                options.stdlog.write("\t".join(map(str, (p.mPredictionId,
                                                         nintron,
                                                         lintron,
                                                         intron_nstops,
                                                         intron_type,
                                                         genomic_sequence[last_e.mGenomeTo - 6:last_e.mGenomeTo].lower() + "|" + sequence[:5] + "..." +
                                                         sequence[-5:] + "|" + genomic_sequence[e.mGenomeFrom:e.mGenomeFrom + 6].lower()))) + "\n")

            options.stdout.write("\t".join(map(str, (p.mPredictionId,
                                                     nintron,
                                                     p.mSbjctToken,
                                                     p.mSbjctStrand,
                                                     last_e.mGenomeTo +
                                                     p.mSbjctGenomeFrom,
                                                     e.mGenomeFrom +
                                                     p.mSbjctGenomeFrom,
                                                     lintron,
                                                     intron_nstops,
                                                     intron_type, prime5, prime3))) + "\n")

            last_e = e

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i.\n" % (
            ninput, noutput))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
