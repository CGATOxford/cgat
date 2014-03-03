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
optic/orthologs2list.py - 
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

   python optic/orthologs2list.py --help

Type::

   python optic/orthologs2list.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import time
import warnings

import CGAT.Experiment as E
import CGAT.Orthologs as Orthologs
import CGAT.GraphTools as GraphTools
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/orthologs2list.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-s", "--species-regex", dest="species_regex", type="string",
                      help="regular expression to extract species from identifier.")

    parser.add_option("-g", "--gene-regex", dest="gene_regex", type="string",
                      help="regular expression to extract gene from identifier.")

    parser.add_option("-b", "--only-best", dest="only_best", action="store_true",
                      help="write only the best pair for a pairing.")

    parser.add_option("-w", "--no-within", dest="within", action="store_false",
                      help="do not write within species pairs.")

    parser.add_option("-d", "--distances", dest="filename_distances", type="string",
                      help="filename with distances between transcripts.")

    parser.add_option("-c", "--no-combine-genes", dest="combine_genes", action="store_false",
                      help="do not combine orthologous clusters which contain the same gene.")

    parser.add_option("--filename-restrict-filter1", dest="filename_restrict_filter1", type="string",
                      help="filename with ids to filter out.")
    parser.add_option("--filename-restrict-filter2", dest="filename_restrict_filter2", type="string",
                      help="filename with ids to filter out.")

    parser.add_option("-f", "--format", dest="format", type="choice", choices=("graph", "components"),
                      help="output format.")

    parser.add_option("-m", "--mode", dest="mode", type="choice", choices=("orthologs", "orphans"),
                      help="analyze either 'orthologs' or 'orphans'.")

    parser.add_option("--genome1", dest="genome1", type="string",
                      help="first genome.")
    parser.add_option("--genome2", dest="genome2", type="string",
                      help="second genome.")

    parser.set_defaults(
        species_regex="^([^|]+)\|",
        gene_regex="^[^|]+\|[^|]+\|([^|]+)\|",
        only_best=None,
        filename_distances=None,
        within=True,
        combine_genes=True,
        report_step=100000,
        use_networkx=False,
        separator="|",
        genome1=None,
        genome2=None,
        mode="orthologs",
        filename_restrict_filter1=None,
        filename_restrict_filter2=None,
        format="graph",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    rs = re.compile(options.species_regex)
    rg = re.compile(options.gene_regex)

    t0 = time.time()
    # retrieve matches between pairs:
    pairs = {}
    max_dist = 0

    if options.filename_distances and options.only_best:
        infile = open(options.filename_distances, "r")
        for line in infile:
            if line[0] == "#":
                continue
            a, b, d = line[:-1].split("\t")[:3]
            d = float(d)
            if a < b:
                key = "%s-%s" % (a, b)
            else:
                key = "%s-%s" % (b, a)

            max_dist = max(d, max_dist)
            pairs[key] = d

        infile.close()

    cluster_id = 0
    ninput, noutput, nmissed, nskipped, nsingletons = 0, 0, 0, 0, 0

    # Read positive filter information:
    filter_restrict1 = {}
    if options.filename_restrict_filter1:
        xx, e = IOTools.ReadList(open(options.filename_restrict_filter1, "r"))
        for x in xx:
            filter_restrict1[Orthologs.Transcript(x).mTranscript] = True

    filter_restrict2 = {}
    if options.filename_restrict_filter2:
        xx, e = IOTools.ReadList(open(options.filename_restrict_filter2, "r"))
        for x in xx:
            filter_restrict2[Orthologs.Transcript(x).mTranscript] = True

    if options.loglevel >= 1:
        options.stdlog.write("# read filtering information: %i/%i\n" %
                             (len(filter_restrict1), len(filter_restrict2)))

    t1 = time.time()

    if options.loglevel >= 1:
        options.stdlog.write("# finished input in %i seconds.\n" % (t1 - t0))

    orthologs = []

    if options.mode == "orthologs":
        orthologs = Orthologs.ReadInterpretation(sys.stdin, options.separator,
                                                 genome1=options.genome1,
                                                 genome2=options.genome2,
                                                 filter_restrict_transcripts1=filter_restrict1,
                                                 filter_restrict_transcripts2=filter_restrict2)
    else:
        orthologs = Orthologs.ReadOrphans(sys.stdin, options.separator,
                                          genome1=options.genome1,
                                          genome2=options.genome2,
                                          filter_restrict_transcripts1=filter_restrict1,
                                          filter_restrict_transcripts2=filter_restrict2)

    ninput = len(orthologs)

    max_dist = map(lambda x: x[4], orthologs)

    t2 = time.time()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# reading %i groups in %i seconds.\n" % (ninput, t2 - t1))

    if options.combine_genes:

        if options.use_networkx:

            nclusters = len(orthologs)
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# before combining genes: %i clusters\n" % len(orthologs))
                options.stdlog.flush()

            # build links between all genes
            # ignore warnings from networkx/matplotlib that a display
            # can not be found
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                import networkx

            graph = networkx.Graph()

            # This procedure skips genes with "0". This is a patch, because
            # these genes should not be there in the first place.
            iteration = 0
            for transcripts1, transcripts2, genes1, genes2, weight in orthologs:
                iteration += 1

                if options.loglevel >= 1:
                    if (iteration % options.report_step == 0):
                        options.stdlog.write("# iteration: %i/%i (%i%%) in %i seconds.\n" %
                                             (iteration, nclusters,
                                              100 * iteration / nclusters, time.time() - t2))
                        options.stdlog.flush()

                for g in genes1.keys():
                    graph.add_node((1, g))
                for g in genes2.keys():
                    graph.add_node((2, g))
                for g1 in genes1.keys():
                    if g1 == "0":
                        continue
                    for g2 in genes1.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge((1, g1), (2, g2))
                    for g2 in genes2.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge((1, g1), (2, g2))
                for g1 in genes2.keys():
                    if g1 == "0":
                        continue
                    for g2 in genes2.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge((2, g1), (2, g2))

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# created graph in %i seconds.\n" % (time.time() - t2))
                options.stdlog.flush()

            tt2 = time.time()

            components = networkx.connected_components(graph)

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# calculated connected components in %i seconds\n" % (time.time() - tt2))
                options.stdlog.flush()

        else:

            graph = GraphTools.ExternalGraph()

            iteration = 0
            nclusters = len(orthologs)

            for transcripts1, transcripts2, genes1, genes2, weight in orthologs:

                iteration += 1

                if options.loglevel >= 1:
                    if (iteration % options.report_step == 0):
                        options.stdlog.write("# iteration: %i/%i (%i%%) in %i seconds.\n" %
                                             (iteration, nclusters,
                                              100 * iteration / nclusters, time.time() - t1))
                        options.stdlog.flush()

                f = "%s;%s"

                for g1 in genes1.keys():
                    if g1 == "0":
                        continue
                    for g2 in genes1.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge(f % (1, g1), f % (2, g2))
                    for g2 in genes2.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge(f % (1, g1), f % (2, g2))
                for g1 in genes2.keys():
                    if g1 == "0":
                        continue
                    for g2 in genes2.keys():
                        if g2 == "0":
                            continue
                        graph.add_edge(f % (2, g1), f % (2, g2))

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# created graph in %i seconds\n" % (time.time() - t2))
                options.stdlog.flush()

            tt2 = time.time()

            graph.finalize()
            components = graph.connected_components()

            if options.loglevel >= 1:
                options.stdlog.write("# retrieved %i connected components in %i seconds\n" % (
                    len(components), time.time() - tt2))
                options.stdlog.flush()

            for x in range(len(components)):
                components[x] = map(lambda y: y.split(";"), components[x])

        tt2 = time.time()

        map_gene2cluster = {}
        for x in range(len(components)):
            for a, b in components[x]:
                map_gene2cluster[b] = x

        new_orthologs = [[[], [], 0] for x in range(len(components))]

        singletons = []

        for transcripts1, transcripts2, genes1, genes2, weight in orthologs:
            if genes1:
                try:
                    cluster_id = map_gene2cluster[genes1.keys()[0]]
                except KeyError:
                    singletons.append(genes1)
            elif genes2:
                try:
                    cluster_id = map_gene2cluster[genes2.keys()[0]]
                except KeyError:
                    singletons.append(genes2)
            else:
                raise "Error, both genes1 and genes2 are emtpy."

            new_orthologs[cluster_id][0] += transcripts1
            new_orthologs[cluster_id][1] += transcripts2
            new_orthologs[cluster_id][2] = weight

        nsingletons = len(singletons)

        orthologs = map(lambda x: (x[0], x[1],
                                   Orthologs.GetGenes(x[0]),
                                   Orthologs.GetGenes(x[1]),
                                   weight), new_orthologs)

        if options.loglevel >= 1:
            options.stdlog.write(
                "# combining genes in %i seconds\n" % (time.time() - tt2))
            options.stdlog.flush()

        if options.loglevel >= 1:
            options.stdlog.write("# after combining genes: %i clusters, %i singletons\n" % (
                len(orthologs), nsingletons))

    t3 = time.time()

    if options.loglevel >= 1:
        options.stdlog.write("# gene clustering in %i seconds.\n" % (t3 - t2))

    cluster_id = 0

    def getCode(s):
        if len(s) == 1:
            return "1"
        elif len(s) == 0:
            return "0"
        else:
            return "m"

    for transcripts1, transcripts2, genes1, genes2, weight in orthologs:

        cluster_id += 1

        g1 = getCode(genes1)
        g2 = getCode(genes2)
        t1 = getCode(transcripts1)
        t2 = getCode(transcripts2)

        if options.format == "graph":

            # find best transcripts
            best_transcripts = {}
            if options.only_best:
                # print only best match between each possible set of genes in
                # ortholog pair
                for gg1, tt1 in genes1.items():
                    for gg2, tt2 in genes2.items():
                        best = max_dist
                        best_pair = None
                        for x in tt1:
                            for y in tt2:
                                if x < y:
                                    key = "%s-%s" % (x, y)
                                else:
                                    key = "%s-%s" % (y, x)

                            if key in pairs:
                                if best > pairs[key]:
                                    best = pairs[key]
                                    best_pair = (x, y)
                        if best_pair:
                            best_transcripts[x] = 1
                            best_transcripts[y] = 1
                            options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                                best_pair[0], best_pair[1], weight, g1, g2, str(t1), str(t2), cluster_id))
                            noutput += 1
                        else:
                            options.stdlog.write(
                                "# missed link between: %s %s\n" % (str(genes1), str(genes2)))
                            nmissed += 1
            else:
                for x in transcripts1:
                    for y in transcripts2:
                        options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                            x, y, weight, g1, g2, str(t1), str(t2), cluster_id))
                        noutput += 1

            if options.within:

                # add self links for first species.
                for x in range(len(transcripts1) - 1):
                    for y in range(x + 1, len(transcripts1)):
                        if not best_transcripts or \
                           (transcripts1[x] in best_transcripts and
                                transcripts1[y] in best_transcripts):
                            options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                                str(transcripts1[x]), str(transcripts1[y]), weight, g1, g2, str(t1), str(t2), cluster_id))
                            noutput += 1

                # add self links for second species
                for x in range(len(transcripts2) - 1):
                    for y in range(x + 1, len(transcripts2)):
                        if not best_transcripts or \
                           (transcripts2[x] in best_transcripts and
                                transcripts2[y] in best_transcripts):
                            options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                                str(transcripts2[x]), str(transcripts2[y]), weight, g1, g2, str(t1), str(t2), cluster_id))
                            noutput += 1

                # If orphans, also add links for genes with a single
                # transripts.
                if options.mode == "orphans":
                    if len(transcripts1) == 1:
                        x, y = 0, 0
                        options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                            str(transcripts1[x]), str(transcripts1[y]), weight, g1, g2, str(t1), str(t2), cluster_id))
                    elif len(transcripts2) == 1:
                        x, y = 0, 0
                        options.stdout.write("%s\t%s\t%6.4f\t%s%s\t%s%s\t%i\n" % (
                            str(transcripts2[x]), str(transcripts2[y]), weight, g1, g2, str(t1), str(t2), cluster_id))

        elif options.format == "components":

            for gg1, tt1 in genes1.items():
                for t in tt1:
                    options.stdout.write("%s\t%i\n" % (str(t), cluster_id))

            for gg2, tt2 in genes2.items():
                for t in tt2:
                    options.stdout.write("%s\t%i" % (str(t), cluster_id))

    if options.loglevel >= 1:
        options.stdout.write("# ninput=%i, noutput=%i, nmissed=%i, skipped=%i, nsingletons=%i\n" % (
            ninput, noutput, nmissed, nskipped, nsingletons))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
