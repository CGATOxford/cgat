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
gpipe/prediction2pairs.py - 
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

   python gpipe/prediction2pairs.py --help

Type::

   python gpipe/prediction2pairs.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.PredictionParser as PredictionParser
import CGAT.IndexedFasta as IndexedFasta
import alignlib_lite

USAGE = """python %s [OPTIONS] < assignments > pairs

Version: $Id: gpipe/prediction2pairs.py 2031 2008-07-15 09:19:05Z andreas $

Take genewise predictions and write aligned pairs of genomic dnas. This
step assumes that there are no frameshifts in the cds sequences. Frameshifts
in the predictions are removed.
""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/prediction2pairs.py 2031 2008-07-15 09:19:05Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genomic data (indexed).")

    parser.add_option("-c", "--cds-gtf-file", dest="filename_cds", type="string",
                      help="filename with cds seguences.")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("paired_fasta", ),
                      help="output format, valid options are: paired_fasta: concatenated pairwise alignments in FASTA format")

    parser.set_defaults(
        genome_file="genome",
        filename_cds="cds.fasta",
        format="paired_fasta",
        filename_suffix=".fasta",
        filename_prefix="",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(1)

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    # reading CDS sequences
    if options.filename_cds:
        cds_sequences = Genomics.ReadPeptideSequences(
            open(options.filename_cds, "r"))
    else:
        cds_sequences = {}

    if options.loglevel >= 1:
        options.stdlog.write("# read %i CDS sequences\n" % len(cds_sequences))

    last_filename_genome = None

    p = PredictionParser.PredictionParserEntry()

    ninput, noutput, nsanity, n3, nlength = 0, 0, 0, 0, 0

    for line in options.stdin:

        if line[0] == "#":
            continue
        if line[0] == '"':
            continue

        p.Read(line)

        ninput += 1

        genomic_fragment = fasta.getSequence(p.mSbjctToken, p.mSbjctStrand,
                                             p.mSbjctGenomeFrom, p.mSbjctGenomeTo)

        if len(genomic_fragment) == 0:
            raise "ERROR: empty fragment %s:%s for line" % (
                p.mSbjctGenomeFrom, p.mSbjctGenomeTo), line

        try:
            cds_fragment = cds_sequences[p.mQueryToken]
        except KeyError:
            options.stdlog.write(
                "# ERROR: cds not found: query %s.\n" % p.mQueryToken)
            continue

        map_query2sbjct, genomic_fragment = Genomics.Alignment2CDNA(p.mMapPeptide2Genome,
                                                                    query_from=p.mQueryFrom,
                                                                    sbjct_from=0,
                                                                    genome=genomic_fragment)

        # check for errors:
        if map_query2sbjct.getRowTo() != p.mQueryTo * 3:
            options.stdlog.write("# ERROR: boundary shift in query at line %s\n# %i %i\n" % (
                line, map_query2sbjct.getRowTo(), p.mQueryTo * 3))

        if map_query2sbjct.getColTo() > len(genomic_fragment):
            options.stdlog.write("# ERROR: length mismatch in line %s\n# genomic fragment (%i) shorter than last aligned residue (%i)\n" %
                                 (line, len(genomic_fragment), map_query2sbjct.getColTo()))
            options.stdlog.write(
                "# cds     %s\n# genomic %s\n" % (str(cds_fragment), genomic_fragment))
            nlength += 1
            continue

        if map_query2sbjct.getRowTo() > len(cds_fragment):
            options.stdlog.write("# ERROR: length mismatch in line %s\n# cds fragment (%i) shorter than last aligned residue (%i)\n" %
                                 (line, len(cds_fragment), map_query2sbjct.getRowTo()))
            options.stdlog.write(
                "# cds     %s\n# genomic %s\n" % (str(cds_fragment), genomic_fragment))
            nlength += 1
            continue

        cds_seq = alignlib_lite.makeSequence(cds_fragment)
        genomic_seq = alignlib_lite.makeSequence(genomic_fragment)

        f = alignlib_lite.AlignmentFormatExplicit(
            map_query2sbjct, cds_seq, genomic_seq)
        row_ali = f.mRowAlignment
        col_ali = f.mColAlignment

        row_ali, col_ali = Genomics.RemoveFrameShiftsFromAlignment(
            row_ali, col_ali)

        row_ali = Genomics.MaskStopCodons(row_ali)
        col_ali = Genomics.MaskStopCodons(col_ali)

        if len(row_ali) != len(col_ali):
            options.stdlog.write("# ERROR: wrong alignment lengths.\n")
            sys.exit(1)

        if len(row_ali) % 3 or len(col_ali) % 3:
            options.stdlog.write(
                "# ERROR: sequences are not a multiple of 3 in line: %s\n" % line)
            options.stdlog.write("# %6i %s\n# %6i %s\n" % (
                len(row_ali), str(row_ali), len(col_ali), str(col_ali)))
            n3 += 1

        input = re.sub("[-X]", "", p.mTranslation)
        ref = re.sub("[-X]", "", Genomics.TranslateDNA2Protein(col_ali))
        if input != ref:
            if options.loglevel >= 1:
                options.stdlog.write("# sanity check failed for %s - %s\n# %6i %s\n# %6i %s\n" % (p.mPredictionId, p.mQueryToken,
                                                                                                  len(input), input,
                                                                                                  len(ref), ref))
            nsanity += 1
            continue

        options.stdout.write(">%s\n%s\n" % (p.mPredictionId, row_ali))
        options.stdout.write(">%s_vs_%s_%s_%i_%i\n%s\n" %
                             (p.mQueryToken, p.mSbjctToken, p.mSbjctStrand, p.mSbjctGenomeFrom, p.mSbjctGenomeTo, col_ali))
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nsanity=%i, nlength=%i, n3=%i\n" % (
            ninput, noutput, nsanity, nlength, n3))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
