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
optic/transcripts2links.py - 
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

   python optic/transcripts2links.py --help

Type::

   python optic/transcripts2links.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import warnings
import CGAT.Experiment as E
import CGAT.Genomics as Genomics


def WriteLinks(outfile, chunk, weight=0.0, separator="|"):

    for x in range(len(chunk) - 1):
        for y in range(x + 1, len(chunk)):
            outfile.write(chunk[x] + "\t" +
                          chunk[y] + "\t" +
                          str(weight) + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/transcripts2links.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-w", "--weight", dest="weight", type="float",
                      help="weight to assign each pair.")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method to use [genes|sequences|togenes].")

    parser.add_option("-f", "--filter-tsv-file", dest="filter", type="string",
                      help="filename with filter (take only those ids in filter).")

    parser.set_defaults(
        method=None,
        weight=0.0,
        separator="|",
        filter=None,
    )

    (options, args) = E.Start(parser)

    subset = {}
    if options.filter:
        x = map(lambda x: x[:-1],
                filter(lambda x: x[0] != "#", open(options.filter, "r").readlines()))
        for xx in x:
            subset[xx] = 1

    if options.method == "genes":

        transcripts = map(
            lambda x: x[:-1].split("\t")[0], filter(lambda x: x[0] != "#", sys.stdin.readlines()))

        if subset:
            transcripts = filter(lambda x: x in subset, transcripts)

        transcripts = map(
            lambda x: x[:-1].split(options.separator), transcripts)
        transcripts.sort(lambda x, y: cmp(x[2], y[2]))

        last_gene = 0
        chunk = []
        for species, id, gene, status in transcripts:
            if last_gene != gene:
                WriteLinks(sys.stdout, chunk, options.weight)
                chunk = []
                last_gene = gene
            chunk.append(options.separator.join((species, id, gene, status)))

        WriteLinks(sys.stdout, chunk, options.weight)

    elif options.method == "sequences":

        sequences = Genomics.ReadPeptideSequences(sys.stdin, filter=subset)
        ids = sequences.keys()

        hids = {}

        for id, seq in sequences.items():
            hid = Genomics.GetHID(seq)
            if hid not in hids:
                hids[hid] = []

            hids[hid].append(id)

        for c in hids.values():
            WriteLinks(sys.stdout, c, options.weight, options.separator)

    elif options.method == "subsequences":

        # ignore warnings from networkx/matplotlib that a display
        # can not be found
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import networkx

        sequences = Genomics.ReadPeptideSequences(sys.stdin)
        ids = sequences.keys()

        # build links between all genes
        graph = networkx.Graph()
        graph.add_nodes_from(ids)
        for x in range(len(ids) - 1):
            idx = ids[x]
            for y in range(x + 1, len(ids)):
                idy = ids[y]
                if sequences[idx] == sequences[idy] or \
                        sequences[idx] in sequences[idy] or \
                        sequences[idy] in sequences[idx]:
                    graph.add_edge(idx, idy)

        components = networkx.connected_components(graph)
        for c in components:
            WriteLinks(sys.stdout, c, options.weight, options.separator)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
