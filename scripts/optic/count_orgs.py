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
optic/count_orgs.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

count presence or absence of species/transcripts/genes in clusters.

Clusters can be given as a map, links or as trees. If a reference tree is given,
absences/presences are checked, if they make sense phylogenetically.

Usage
-----

Example::

   python optic/count_orgs.py --help

Type::

   python optic/count_orgs.py --help

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
import getopt
import optparse

import CGAT.Orthologs as Orthologs
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools


class SpeciesCounts:

    def __init__(self):
        self.mGenes = set()
        self.mTranscripts = set()
        self.mTrees = set()


def GetPattern(data, l):

    pattern = ["0"] * l
    for x in range(len(data)):
        if data[x] > 0:
            pattern[x] = "1"
    pattern = string.join(pattern, "")
    return pattern


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: optic/count_orgs.py 1706 2007-12-11 16:46:11Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--reference-tree", dest="reference_tree", type="string",
                      help="reference tree to read.")

    parser.add_option("-p", "--filename-patterns", dest="filename_patterns", type="string",
                      help="filename with patterns to output.")

    parser.add_option("-u", "--filename-summary", dest="filename_summary", type="string",
                      help="filename with summary to output.")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("map", "links", "trees"),
                      help="output format.")

    parser.add_option("-o", "--organisms", dest="column2org", type="string",
                      help="sorted list of organisms.")

    parser.add_option("-s", "--species-regex", dest="species_regex", type="string",
                      help="regular expression to extract species from identifier.")

    parser.add_option("-g", "--gene-regex", dest="gene_regex", type="string",
                      help="regular expression to extract gene from identifier.")

    parser.set_defaults(
        reference_tree=None,
        format="map",
        filename_patterns=None,
        column2org=None,
        species_regex="^([^|]+)\|",
        gene_regex="^[^|]+\|[^|]+\|([^|]+)\|",
        separator="|",
        filename_summary=None,
    )

    (options, args) = E.Start(parser)

    if options.reference_tree:
        if options.reference_tree[0] == "(":
            nexus = TreeTools.Newick2Nexus(options.reference_tree)
        else:
            nexus = TreeTools.Newick2Nexus(open(options.reference_tree, "r"))
        reference_tree = nexus.trees[0]

        if options.loglevel >= 3:
            print "# reference tree:"
            print reference_tree.display()

    else:
        reference_tree = None

    clusters = {}
    if options.format == "map":

        for line in sys.stdin:
            if line[0] == "#":
                continue
            id, r = line[:-1].split("\t")
            if r not in clusters:
                clusters[r] = []
            clusters[r].append(id)

    elif options.format == "trees":

        nexus = TreeTools.Newick2Nexus(sys.stdin)

        for tree in nexus.trees:
            clusters[tree.name] = tree.get_taxa()

    elif options.format == "links":
        members = set()
        id = None
        for line in sys.stdin:
            if line[0] == "#":
                continue

            if line[0] == ">":
                if id:
                    clusters[id] = members
                x = re.match(">cluster #(\d+)", line[:-1])
                if x:
                    id = x.groups()[0]
                else:
                    id = line[1:-1]
                members = set()
                continue

            data = line[:-1].split("\t")[:2]
            members.add(data[0])
            members.add(data[1])

        if id:
            clusters[id] = members

    if len(clusters) == 0:
        raise "empty input."

    ########################################################################
    ########################################################################
    ########################################################################
    # sort out reference tree
    ########################################################################
    rs = re.compile(options.species_regex)
    rg = re.compile(options.gene_regex)
    extract_species = lambda x: rs.search(x).groups()[0]

    # prune tree to species present
    species_set = set()
    for cluster, members in clusters.items():
        species_set = species_set.union(set(map(extract_species, members)))

    if reference_tree:

        TreeTools.PruneTree(reference_tree, species_set)

        if options.loglevel >= 1:
            options.stdlog.write(
                "# Tree after pruning: %i taxa.\n" % len(reference_tree.get_taxa()))

    if options.column2org:
        options.column2org = options.column2org.split(",")
    elif reference_tree:
        options.column2org = []
        for nx in reference_tree.get_terminals():
            options.column2org.append(reference_tree.node(nx).get_data().taxon)
    else:
        options.column2org = []
        for x in species_set:
            options.column2org.append(x)

    options.org2column = {}
    for x in range(len(options.column2org)):
        options.org2column[options.column2org[x]] = x

    if reference_tree:
        reference_patterns = TreeTools.calculatePatternsFromTree(
            reference_tree, options.column2org)

        if options.loglevel >= 3:
            print "# reference patterns:"
            print reference_patterns

    ##########################################################################
    notus = len(options.column2org)
    patterns = {}
    species_counts = [SpeciesCounts() for x in options.column2org]

    # first genes, then transcripts
    options.stdout.write("mali\tpattern\tpresent\tngenes\t%s\tntranscripts\t%s\n" % (
        "\t".join(options.column2org), "\t".join(options.column2org)))

    keys = clusters.keys()
    keys.sort()
    for cluster in keys:
        members = clusters[cluster]

        count_genes = [{} for x in range(len(options.org2column))]
        count_transcripts = [0] * len(options.org2column)

        for m in members:
            data = m.split(options.separator)

            if len(data) == 4:
                s, t, g, q = data
            elif len(data) == 2:
                s, g = data
                t = g

            if s not in options.org2column:
                raise "unknown species %s" % s

            col = options.org2column[s]

            count_transcripts[col] += 1
            if g not in count_genes[col]:
                count_genes[col][g] = 0

            count_genes[col][g] += 1

            species_counts[col].mGenes.add(g)
            species_counts[col].mTranscripts.add(t)
            species_counts[col].mTrees.add(cluster)

        ntotal_transcripts = reduce(lambda x, y: x + y, count_transcripts)
        npresent_transcripts = len(filter(lambda x: x > 0, count_transcripts))
        ntotal_genes = reduce(lambda x, y: x + y, map(len, count_genes))
        npresent_genes = len(filter(lambda x: x > 0, map(len, count_genes)))

        pattern = GetPattern(count_transcripts, notus)
        if pattern not in patterns:
            patterns[pattern] = 0
        patterns[pattern] += 1
        options.stdout.write(string.join((cluster, pattern, str(npresent_genes),
                                          str(ntotal_genes),
                                          string.join(
                                              map(str, map(len, count_genes)), "\t"),
                                          str(ntotal_transcripts),
                                          string.join(map(str, count_transcripts), "\t")), "\t") + "\n")

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # write pattern summary
    ##########################################################################
    xx = patterns.keys()
    xx.sort()
    if options.filename_patterns:
        outfile = open(options.filename_patterns, "w")
    else:
        outfile = sys.stdout

    for x in range(len(options.column2org)):
        outfile.write("# %i = %s\n" % (x, options.column2org[x]))

    if reference_tree:
        outfile.write("pattern\tcounts\tisok\n")
    else:
        outfile.write("pattern\tcounts\n")

    for x in xx:
        if reference_tree:
            if x in reference_patterns:
                is_ok = "1"
            else:
                is_ok = "0"
            outfile.write("%s\t%s\t%s\n" % (x, patterns[x], is_ok))
        else:
            outfile.write("%s\t%s\n" % (x, patterns[x]))

    if outfile != sys.stdout:
        outfile.close()

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # write summary counts per species
    ##########################################################################
    if options.filename_summary:
        outfile = open(options.filename_summary, "w")
    else:
        outfile = sys.stdout

    outfile.write("species\tntranscripts\tngenes\tntrees\n")

    for species, col in options.org2column.items():
        outfile.write("%s\t%i\t%i\t%i\n" % (species,
                                            len(species_counts[
                                                col].mTranscripts),
                                            len(species_counts[col].mGenes),
                                            len(species_counts[col].mTrees)))

    if outfile != sys.stdout:
        outfile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
