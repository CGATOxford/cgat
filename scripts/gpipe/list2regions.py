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
gpipe/list2regions.py - predict genes from a list of associations
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Options::

    -h, --help                      print this message.
    -v, --verbose=                  loglevel.
    -p, --peptides-fasta-file=                 file with peptide sequences (FASTA).
    -g, --genome-file=              pattern for filenames with the genomic DNA (FASTA).
    -m, --map=                      map of peptide identifiers
    --disable-conflict              turn of resolution of conflicts
    --disable-overlap               turn of resolution of overlaps
    --disable-suboptimal            turn of elimination of suboptimal predictions
    --disable-activation            turn of reactivation of eliminated queries

Usage
-----

Example::

   python gpipe/list2regions.py 

Type::

   python gpipe/list2regions.py --help

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
import CGAT.Genomics as Genomics
from CGAT.Predictor2 import PredictorExonerate

param_long_options = ["verbose=", "help", "max-percent-overlap=",
                      "min-coverage-query=", "min-score=", "min-percent-identity=",
                      "max-matches=", "peptides=", "genome-file=",
                      "map=", "version"]

param_short_options = "v:ho:c:s:i:m:p:"


# pattern for genomes, %s is substituted for the sbjct_token
param_genome_file = "genome_%s.fasta"
param_filename_peptides = None
param_filename_map = None
param_options = "--subopt FALSE"
param_loglevel = 2

param_border = 100


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ("-p", "--peptides-fasta-file"):
            param_filename_peptides = a
        elif o in ("-m", "--map"):
            param_filename_map = a
        elif o in ("-g", "--genome_file"):
            param_genome_file = a

    last_filename_genome = None
    # read peptide sequences
    if param_filename_peptides:
        peptide_sequences = Genomics.ReadPeptideSequences(
            open(param_filename_peptides, "r"))
    else:
        peptide_sequences = {}

    map_a2b = {}

    if param_filename_map:
        infile = open(param_filename_map, "r")
        for line in infile:
            a, b = string.split(line[:-1], "\t")
            map_a2b[a] = b

    predictor = PredictorExonerate()

    nmissed, nfound, nfailed = 0, 0, 0

    for line in sys.stdin:

        gene, sbjct_genome_from, sbjct_genome_to, sbjct_strand, chromosome, query_token = re.split(
            "\s+", line[:-1])

        sbjct_token = "chr" + chromosome
        sbjct_genome_from, sbjct_genome_to = map(
            int, (sbjct_genome_from, sbjct_genome_to))
        if sbjct_strand == "1":
            sbjct_strand = "+"
        else:
            sbjct_strand = "-"

        filename_genome = param_genome_file % sbjct_token

        if last_filename_genome != filename_genome:
            E.debug("reading genome %s" % filename_genome)

            forward_sequences, reverse_sequences = Genomics.ReadGenomicSequences(
                open(filename_genome, "r"))
            last_filename_genome = filename_genome

            lgenome = len(forward_sequences[sbjct_token])

        if sbjct_strand == "+":
            genomic_sequence = forward_sequences[sbjct_token]
        else:
            genomic_sequence = reverse_sequences[sbjct_token]
            x = sbjct_genome_to
            sbjct_genome_to = lgenome - sbjct_genome_from
            sbjct_genome_from = lgenome - x

        sbjct_genome_from -= 100
        sbjct_genome_to += 100

        if map_a2b.has_key(query_token):
            query_token = map_a2b[query_token]

        E.debug("aligning: %s to %s:%s:%i-%i" %
                (query_token, sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to))

        if not peptide_sequences.has_key(query_token):
            E.warn("peptides sequence not found for %s" % query_token)
            nmissed += 1
            continue

        result = predictor(query_token,
                           peptide_sequences[query_token],
                           sbjct_token,
                           genomic_sequence[sbjct_genome_from:sbjct_genome_to],
                           param_options,
                           0, sbjct_genome_to - sbjct_genome_from)

        if result:
            result.ShiftGenomicRegion(sbjct_genome_from)
            result.SetStrand(sbjct_strand)
            print str(result)
            nfound += 1
        else:
            nfailed += 1

    if param_loglevel >= 1:
        print "# found=%i, missed=%i, failed=%i" % (nfound, nmissed, nfailed)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
