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
gpipe/exons2map.py - 
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

   python gpipe/exons2map.py --help

Type::

   python gpipe/exons2map.py --help

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
import CGAT.Exons as Exons
import alignlib_lite


USAGE = """python %s [OPTIONS] < psl > predictions

Convert exon list to a map of prediction to genome.

Note: This file takes in forward strand coordinates, but
returns forward/backward coordinates.

Version: $Id: gpipe/exons2map.py 1799 2008-03-28 11:44:19Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-c, --contigs-tsv-file=                  filename with contig lengths
""" % sys.argv[0]


param_long_options = ["verbose=", "help", "contigs=", "version"]
param_short_options = "v:hc:"

param_trans = None

param_filename_contigs = None


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
        elif o in ("-c", "--contigs-tsv-file"):
            param_filename_contigs = a

    print E.GetHeader()
    print E.GetParams()

    last_exon = Exons.Exon()

    contig_sizes = {}
    if param_filename_contigs:

        infile = open(param_filename_contigs, "r")
        for line in infile:
            if line[0] == "#":
                continue

            sbjct_token, size = line[:-1].split("\t")[:2]
            contig_sizes[sbjct_token] = int(size)

    map_prediction2genome = alignlib_lite.makeAlignmentSet()
    nexons, npairs = 0, 0

    for line in sys.stdin:

        if line[0] == "#":
            continue

        this_exon = Exons.Exon()
        this_exon.Read(line)

        if this_exon.mSbjctStrand == "-":
            this_exon.InvertGenomicCoordinates(
                contig_sizes[this_exon.mSbjctToken])

        nexons += 1

        if last_exon.mQueryToken != this_exon.mQueryToken:

            if last_exon.mQueryToken:
                f = alignlib_lite.AlignmentFormatEmissions(
                    map_prediction2genome)
                print string.join(map(str, (last_exon.mQueryToken,
                                            last_exon.mSbjctToken,
                                            last_exon.mSbjctStrand,
                                            f)), "\t")

                npairs += 1
            map_prediction2genome.clear()

        alignlib_lite.addDiagonal2Alignment(map_prediction2genome,
                                            this_exon.mPeptideFrom + 1,
                                            this_exon.mPeptideTo + 1,
                                            this_exon.mGenomeFrom - this_exon.mPeptideFrom)

        last_exon = this_exon

    f = alignlib_lite.AlignmentFormatEmissions(map_prediction2genome)
    print string.join(map(str, (last_exon.mQueryToken,
                                last_exon.mSbjctToken,
                                last_exon.mSbjctStrand,
                                f)), "\t")
    npairs += 1

    print "# nexons=%i, npairs=%i" % (nexons, npairs)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
