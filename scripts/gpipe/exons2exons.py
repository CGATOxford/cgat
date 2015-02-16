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
gpipe/exons2exons.py - 
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

   python gpipe/exons2exons.py --help

Type::

   python gpipe/exons2exons.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import CGAT.IndexedFasta as IndexedFasta


USAGE = """python %s [OPTIONS] < exons.in > exons.out

Modify an exon list. Methods to apply are:

remove-stop: removes a traling stop codon from the last exon.

"""


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/exons2exons.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method",
                      help="method to apply.", type="choice",
                      choices=("remove-stop",))

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genomic data (indexed).")

    parser.add_option("--is-forward-coordinates", dest="forward_coordinates", action="store_true",
                      help="work in forward coordinates.")

    parser.set_defaults(
        method=None,
        forward_coordinates=False,
        genome_file=None)

    (options, args) = E.Start(parser)

    if options.method == "remove-stop" and not options.genome_file:
        raise "please supply genome file for method %s" % options.method

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
        contig_sizes = fasta.getContigSizes()
        exons = Exons.ReadExonBoundaries(sys.stdin,
                                         contig_sizes=contig_sizes)
    else:
        exons = Exons.ReadExonBoundaries(sys.stdin)

    ninput, noutput, nremoved_stops, nremoved_exons = 0, 0, 0, 0
    for id, ee in exons.items():

        if options.loglevel >= 3:
            for e in ee:
                options.stdlog.write("# %s\n" % str(e))

        if options.method == "remove-stop":
            e = ee[-1]
            d = min(3, e.mPeptideTo - e.mPeptideFrom)
            if d < 3:
                codon2 = fasta.getSequence(
                    e.mSbjctToken, e.mSbjctStrand, e.mGenomeTo - d, e.mGenomeTo)
                prev_e = ee[-2]
                codon1 = fasta.getSequence(
                    prev_e.mSbjctToken, prev_e.mSbjctStrand, prev_e.mGenomeTo - (3 - d), prev_e.mGenomeTo)
                codon = codon1 + codon2
            else:
                codon = fasta.getSequence(
                    e.mSbjctToken, e.mSbjctStrand, e.mGenomeTo - d, e.mGenomeTo)

            if codon.upper() in Genomics.StopCodons:

                if d < 3:
                    nremoved_exons += 1
                    d = 3 - d
                    del ee[-1]
                    e = ee[-1]

                e.mGenomeTo -= d
                e.mPeptideTo -= d
                nremoved_stops += 1

                if e.mGenomeTo == e.mGenomeFrom:
                    nremoved_exons += 1
                    del ee[-1]
                    e = ee[-1]

            assert(e.mGenomeTo > e.mGenomeFrom)
            assert(e.mPeptideTo > e.mPeptideFrom)

        if options.forward_coordinates:

            l = contig_sizes[ee[0].mSbjctToken]
            for e in ee:
                e.InvertGenomicCoordinates(l)

        for e in ee:
            options.stdout.write(str(e) + "\n")

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nremoved_stops=%i, nremoved_exons=%i\n" % (
            ninput, noutput, nremoved_stops, nremoved_exons))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
