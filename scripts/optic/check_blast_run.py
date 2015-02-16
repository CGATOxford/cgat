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
optic/check_blast_run.py -
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

check blast run.

Input: one or more sets of vertices.

Usage
-----

Example::

   python optic/check_blast_run.py --help

Type::

   python optic/check_blast_run.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments

parser = E.OptionParser(
    version="%prog version: $Id: optic/check_blast_run.py 2781 2009-09-10 11:33:14Z andreas $")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-n", "--vertices", dest="vertices", action="append",
                      help="filename with vertices.")
    parser.add_option("-e", "--extra", dest="filename_extra", type="string",
                      help="filename to store extra vertices in.")
    parser.add_option("-m", "--missed", dest="filename_missed", type="string",
                      help="filename to store missed vertices in.")

    parser.set_defaults(
        vertices=[],
        filename_extra=None,
        filename_missed=None,
    )

    (options, args) = E.Start(parser)

    if len(options.vertices) == "":
        raise "please specify one set of vertices."

    vertices = {}
    index = 0
    missed_queries = []
    nvertices = [0] * len(options.vertices)
    for x in range(len(options.vertices)):
        f = options.vertices[x]
        vv = map(lambda x: x[:-1].split("\t")[0],
                 filter(lambda x: x[0] != "#", open(f, "r").readlines()))
        nvertices[x] = len(vv)
        missed_queries.append([])
        for v in vv:
            vertices[v] = [x, 0, 0]
        if options.loglevel >= 1:
            print "# read %i vertices from %s" % (len(vv), f)
            sys.stdout.flush()

    l = BlastAlignments.Link()
    extra_vertices = {}
    for line in sys.stdin:

        if line[0] == "#":
            continue

        l.Read(line)

        if l.mQueryToken in vertices:
            vertices[l.mQueryToken][1] += 1
        else:
            extra_vertices[l.mQueryToken] = 1

        if l.mSbjctToken in vertices:
            vertices[l.mSbjctToken][2] += 1
        else:
            extra_vertices[l.mSbjctToken] = 1

    found_queries = [0] * len(options.vertices)
    found_sbjcts = [0] * len(options.vertices)

    for v, vv in vertices.items():
        index, nquery, nsbjct = vv
        if nquery:
            found_queries[index] += 1
        else:
            missed_queries[index].append(v)

        if nsbjct:
            found_sbjcts[index] += 1

    headers = ("set", "name", "tvertex", "nmissed", "pmissed",
               "nquery", "pquery", "nsbjct", "psbjct")

    print "\t".join(headers)

    for x in range(len(options.vertices)):
        print "%i\t%s\t%i\t%i\t%5.2f\t%i\t%5.2f\t%i\t%5.2f" % (x, options.vertices[x],
                                                               nvertices[x],
                                                               len(missed_queries[
                                                                   x]),
                                                               100 *
                                                               float(
                                                                   len(missed_queries[x])) / nvertices[x],
                                                               found_queries[
                                                                   x],
                                                               100 *
                                                               float(
                                                                   found_queries[x]) / nvertices[x],
                                                               found_sbjcts[x],
                                                               100 * float(found_sbjcts[x]) / nvertices[x])

    print "//"
    print "%i vertices not in set" % len(extra_vertices)

    if options.filename_extra and len(extra_vertices) > 0:
        outfile = open(options.filename_extra, "w")
        for x in extra_vertices.keys():
            outfile.write("%s\n" % x)
        outfile.close()

    if options.filename_missed:
        outfile = open(options.filename_missed, "w")
        for x in range(len(options.vertices)):
            for y in missed_queries[x]:
                outfile.write("%i\t%s\t%s\n" % (x, options.vertices[x], y))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
