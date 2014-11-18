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
gtf2exons.py - convert gtf to exons file format
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a gtf formatted file to exons format used
by gpipe.

The gtf file should contain only the query token as
its last option.

Usage
-----

Example::

   python gtf2exons.py --help

Type::

   python gtf2exons.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse

import CGAT.Experiment as E
import CGAT.Exons as Exons
import CGAT.IndexedFasta as IndexedFasta


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gtf2exons.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genomic data (indexed).")

    parser.add_option("--coordinate-format", dest="coordinate_format", type="string",
                      help="input type of coordinates.")

    parser.add_option("--is-forward-coordinates", dest="forward_coordinates", action="store_true",
                      help="output forward coordinates.")

    parser.add_option("-e", "--extract-id", dest="extract_id", type="string",
                      help="""regular expression to extract id from id column, e.g. 'transcript_id "(\S+)"'.""" )

    parser.set_defaults(
        coordinate_format="zero-forward",
        forward_coordinates=False,
        genome_file=None,
        extract_id=None)

    (options, args) = E.Start(parser)

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
        contig_sizes = fasta.getContigSizes()
    else:
        contig_sizes = {}

    if options.extract_id:
        extract_id = re.compile(options.extract_id)
    else:
        extract_id = None

    converter = IndexedFasta.getConverter(options.coordinate_format)

    exons = Exons.ReadExonBoundaries(sys.stdin,
                                     contig_sizes=contig_sizes,
                                     converter=converter,
                                     do_invert=True,
                                     format="gtf",
                                     gtf_extract_id=extract_id)

    ntranscripts, nexons, nerrors = 0, 0, 0
    for id, ee in exons.items():
        ntranscripts += 1
        has_error = False
        for e in ee:
            if options.forward_coordinates and e.mSbjctToken in contig_sizes and \
                    e.mSbjctStrand == "-":
                l = contig_sizes[e.mSbjctToken]
                e.mGenomeFrom, e.mGenomeTo = l - e.mGenomeTo, l - e.mGenomeFrom

            if e.mGenomeFrom < 0:
                has_error = True
                if options.loglevel >= 1:
                    options.stderr.write("# Error: %s\n" % str(e))
                break

            options.stdout.write(str(e) + "\n")
            nexons += 1

        if has_error:
            nerrors += 1
            continue

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ntranscripts=%i, nexons=%i, nerrors=%i\n" % (ntranscripts, nexons, nerrors))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
