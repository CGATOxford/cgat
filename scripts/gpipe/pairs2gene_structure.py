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
gpipe/pairs2gene_structure.py - 
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

   python gpipe/pairs2gene_structure.py --help

Type::

   python gpipe/pairs2gene_structure.py --help

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
import CGAT.PredictionParser as PredictionParser
import alignlib_lite

USAGE = """python %s [OPTIONS] < assignments > pairs

Version: $Id: gpipe/pairs2gene_structure.py 1799 2008-03-28 11:44:19Z andreas $

Take a list of orthologous transcripts and write out a list
of orthologous transcripts.

Options:
-h, --help                      print this message.
-v, --verbose                   loglevel.
-g, --genome-file=           pattern for filenames with the genomic DNA (FASTA).
-c, --cds-gtf-file=                      filename with coding sequences
-f, --format=                   output format, valid options are:
                                paired_fasta: concatenated pairwise alignments in FASTA format
                                
""" % sys.argv[0]

param_long_options = ["verbose=", "help", "genome-file=", "format=",
                      "cds=", "version"]
param_short_options = "v:hg:f:c:"

param_loglevel = 0

# pattern for genomes, %s is substituted for the sbjct_token
param_genome_file = "genome_%s.fasta"

# filename with cdss
param_filename_cdss = "cds.fasta"

# output format
param_format = "paired_fasta"

# prefix/suffix for output files
param_filename_suffix = ".fasta"
param_filename_prefix = ""

# ------------------------------------------------------------


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
        elif o in ("-g", "--genome-file"):
            param_genome_file = a
        elif o in ("-c", "--cds-gtf-file"):
            param_filename_cds = a

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(1)

    print E.GetHeader()
    print E.GetParams()

    # reading CDS sequences
    if param_filename_cds:
        cds_sequences = Genomics.ReadPeptideSequences(
            open(param_filename_cds, "r"))
    else:
        cds_sequences = {}

    if param_loglevel >= 1:
        print "# read %i CDS sequences" % len(cds_sequences)

    last_filename_genome = None

    p = PredictionParser.PredictionParserEntry()
    for line in sys.stdin:

        if line[0] == "#":
            continue
        if line[0] == '"':
            continue

        p.Read(line)

        # read genomic sequence
        if "%s" in param_genome_file:
            filename_genome = param_genome_file % p.mSbjctToken
        else:
            filename_genome = param_genome_file

        if last_filename_genome != filename_genome:
            if param_loglevel >= 2:
                print "# reading genome %s" % filename_genome

            forward_sequences, reverse_sequences = Genomics.ReadGenomicSequences(
                open(filename_genome, "r"))
            last_filename_genome = filename_genome

        if p.mSbjctStrand == "+":
            genomic_sequence = forward_sequences[p.mSbjctToken]
        else:
            genomic_sequence = reverse_sequences[p.mSbjctToken]

        try:
            cds_fragment = cds_sequences[p.mQueryToken]
        except KeyError:
            print "# ERROR: cds not found: query %s." % p.mQueryToken
            continue

        genomic_fragment = genomic_sequence[
            p.mSbjctGenomeFrom:p.mSbjctGenomeTo]

        if len(genomic_fragment) == 0:
            raise "ERROR: empty fragment %s:%s for line" % (
                p.mSbjctGenomeFrom, p.mSbjctGenomeTo), line

        map_query2sbjct, genomic_fragment = Genomics.Alignment2CDNA(p.mMapPeptide2Genome,
                                                                    query_from=p.mQueryFrom -
                                                                    1,
                                                                    sbjct_from=0,
                                                                    genome=genomic_fragment)

        # check for errors:
        if map_query2sbjct.getRowTo() != p.mQueryTo * 3:
            print str(p)
            raise "# ERROR: boundary shift in query: %i %i" % (
                map_query2sbjct.getRowTo(), p.mQueryTo * 3)

        if map_query2sbjct.getColTo() > len(genomic_fragment):
            print "# ERROR: length mismatch: genomic fragment (%i) shorter than last aligned residue (%i)" %\
                (len(genomic_fragment), map_query2sbjct.getColTo())
            print "#", line
            print "# cds"
            print "#", cds_fragment
            print "# genomic"
            print "#", genomic_fragment
            continue

        if map_query2sbjct.getRowTo() > len(cds_fragment):
            print "# ERROR: length mismatch: cds fragment (%i) shorter than last aligned residue (%i)" %\
                (len(cds_fragment), map_query2sbjct.getRowTo())
            print "#", line
            print "# cds"
            print "#", cds_fragment
            print "# genomic"
            print "#", genomic_fragment
            continue

        cds_seq = alignlib_lite.makeSequence(cds_fragment)
        genomic_seq = alignlib_lite.makeSequence(genomic_fragment)

        data = map(lambda x: string.split(x, "\t"),
                   string.split(alignlib_lite.writePairAlignment(cds_seq,
                                                                 genomic_seq,
                                                                 map_query2sbjct), "\n"))

        row_ali, col_ali = Genomics.RemoveFrameShiftsFromAlignment(
            data[0][1], data[1][1])

        row_ali = Genomics.MaskStopCodons(row_ali)
        col_ali = Genomics.MaskStopCodons(col_ali)

        if len(row_ali) != len(col_ali):
            print "# ERROR: wrong alignment lengths."
            sys.exit(1)

        if len(row_ali) % 3 or len(col_ali) % 3:
            print line
            print row_ali
            print col_ali
            print len(row_ali), len(col_ali)
            print " ERROR: non-codons in alignment."
            sys.exit(1)

        print ">%i\n%s" % (p.mPredictionId, row_ali)
        print ">%s_vs_%s_%s_%i_%i\n%s" % \
              (p.mQueryToken, p.mSbjctToken, p.mSbjctStrand,
               p.mSbjctGenomeFrom, p.mSbjctGenomeTo, col_ali)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
