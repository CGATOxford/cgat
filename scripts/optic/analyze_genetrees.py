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
optic/analyze_genetrees.py - 
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

   python optic/analyze_genetrees.py --help

Type::

   python optic/analyze_genetrees.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import re

from types import *
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.Tree as Tree
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats
import CGAT.TreeReconciliation as TreeReconciliation


def printCounts(heights_per_species, relheights_per_species,
                heights_per_tree, relheights_per_tree,
                options,
                prefix_header, prefix_row):
    """print counts for each section."""

    def printHeightsPerTree(values, section, options, prefix_header,
                            prefix_row):

        if not values:
            return

        outfile, is_new = TreeReconciliation.getFile(options, section)
        if is_new:
            outfile.write("%s%s\theights\n" % (
                prefix_header, "\t".join(Stats.DistributionalParameters().getHeaders())))

        s = Stats.DistributionalParameters(values)
        s.setFormat(options.format_branch_length)
        outfile.write("%s%s\t%s\n" % (prefix_row,
                                      str(s),
                                      ",".join(map(lambda x: options.format_branch_length % x, values))))

    def printHeightsPerSpecies(values, section, options, prefix_header,
                               prefix_row):

        if not values:
            return

        # distributions of distance to node
        outfile, is_new = TreeReconciliation.getFile(options, section)
        if is_new:
            outfile.write("%sspecies\t%s\theights\n" % (
                prefix_header, "\t".join(Stats.DistributionalParameters().getHeaders())))

        for species in sorted(values.keys()):
            s = Stats.DistributionalParameters(values[species])
            s.setFormat(options.format_branch_length)
            outfile.write("%s%s\t%s\t%s\n" % (prefix_row,
                                              species,
                                              str(s),
                                              ",".join(map(lambda x: options.format_branch_length % x, values[species]))))

    printHeightsPerSpecies(
        heights_per_species, "species_heights", options, prefix_header, prefix_row)
    printHeightsPerSpecies(
        relheights_per_species, "species_relheights", options, prefix_header, prefix_row)

    printHeightsPerTree(
        heights_per_tree, "tree_heights", options, prefix_header, prefix_row)
    printHeightsPerTree(
        relheights_per_tree, "tree_relheights", options, prefix_header, prefix_row)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_genetrees.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-r", "--species-regex", dest="species_regex", type="string",
                      help="regular expression to extractspecies from identifier.")

    parser.add_option("--gene-regex", dest="gene_regex", type="string",
                      help="regular expression to extract gene from identifier.")

    parser.add_option("--filename-filter-positives", dest="filename_filter_positives", type="string",
                      help="filename with positive list of trees to analyze.")

    parser.add_option("-s", "--filename-species-tree", dest="filename_species_tree", type="string",
                      help="filename with species tree.")

    parser.add_option("--filename-species2colour", dest="filename_species2colour", type="string",
                      help="filename with map of species to colours. If not given, random colours are assigned to species.")

    parser.add_option("-t", "--species-tree", dest="species_tree", type="string",
                      help="species tree.")

    parser.add_option("-e", "--locations-tsv-file", dest="filename_locations", type="string",
                      help="filename with map of transcript information to location information.")

    parser.add_option("--no-create", dest="create", action="store_false",
                      help="do not create files, but append to them.")

    parser.add_option("--max-separation", dest="max_separation", type="int",
                      help="maximum allowable separation between syntenic segments for border plot (set to 0, if syntey is enough).")

    parser.add_option("--filename-species2url", dest="filename_species2url", type="string",
                      help="filename with mapping information of species to URL.")

    parser.add_option("--column-prefix", dest="prefix", type="string",
                      help="prefix to add as first column.")

    parser.add_option("--outgroup-species", dest="outgroup_species", type="string",
                      help="species to used as outgroups. Separate multiple species by ','.")

    parser.add_option("--subtrees-trees", dest="subtrees_trees", action="store_true",
                      help="write trees for subtrees.")

    parser.add_option("--subtrees-identifiers", dest="subtrees_identifiers", action="store_true",
                      help="write identifiers of subtrees.")

    parser.add_option("--svg-add-ids", dest="svg_add_ids", action="store_true",
                      help="add node ids to svg plot.")

    parser.add_option("--svg-otus", dest="svg_otus", type="string",
                      help="otus to output in svg species tree.")

    parser.add_option("--svg-branch-lenghts", dest="svg_branch_lengths", type="choice",
                      choices=("contemporary", "uniform", "median"),
                      help="branch lengths in species tree.")

    parser.add_option("--print-totals", dest="print_totals", action="store_true",
                      help="output totals sections.")

    parser.add_option("--print-subtotals", dest="print_subtotals", action="store_true",
                      help="output subtotals sections.")

    parser.add_option("--print-best", dest="print_best", action="store_true",
                      help="output best node assignment for each node in gene tree.")

    parser.add_option("--print-svg", dest="print_svg", action="store_true",
                      help="output svg files.")

    parser.add_option("--print-species-svg", dest="print_species_svg", action="store_true",
                      help="output species svg files.")

    parser.add_option("--output-filename-pattern", dest="output_pattern", type="string",
                      help="""output pattern for separate output of sections [default: %default].
                      Set to None, if output to stdout. Can contain one %s to be substituted with section.""")

    parser.add_option("--output-pattern-svg", dest="output_pattern_svg", type="string",
                      help="filename for svg output. If it contains %s, this is replaced by gene_tree name.")

    parser.add_option("--filename-node-types", dest="filename_node_types", type="string",
                      help="filename with node type information from a previous run.")

    parser.add_option("--analyze-resolution-data", dest="analyze_resolution_data", type="choice", action="append",
                      choices=("stats", "histograms"),
                      help="stdin is resolution data.")

    parser.add_option("--filter-quality", dest="filter_quality", type="choice",
                      choices=("all", "genes", "pseudogenes"),
                      help="filter predictions by gene type.")

    parser.add_option("--filter-location", dest="filter_location", type="choice",
                      choices=("all", "local", "non-local", "cis", "unplaced"),
                      help="filter predictions by location.")

    parser.add_option("--remove-unplaced", dest="remove_unplaced", action="store_true",
                      help="remove predictions on unplaced contigs.")

    parser.add_option("--skip-without-outgroups", dest="skip_without_outgroups", action="store_true",
                      help="skip clusters without outgroups.")

    parser.set_defaults(
        filter_quality="all",
        filter_location="all",
        remove_unplaced=False,
        species_regex="^([^|]+)\|",
        gene_regex="^[^|]+\|[^|]+\|([^|]+)\|",
        filename_species_tree=None,
        priority={"Speciation": 0,
                   "SpeciationDeletion": 1,
                  "Transcripts": 2,
                   "DuplicationLineage": 3,
                  "Duplication": 4,
                  "DuplicationDeletion": 5,
                  "DuplicationInconsistency": 6,
                  "Outparalogs": 7,
                  "InconsistentTranscripts": 8,
                  "Inconsistency": 9,
                  "Masked": 10},
        species_tree=None,
        filename_species2colour=None,
        filename_locations=None,
        max_separation=0,
        filename_species2url=None,
        separator="|",
        prefix=None,
        output_pattern=None,
        output_pattern_svg=None,
        outgroup_species=None,
        svg_add_ids=False,
        svg_branch_lengths="median",
        svg_otus=None,
        subtrees=False,
        print_svg=False,
        print_subtotals=False,
        print_totals=False,
        print_best=False,
        subtrees_identifiers=False,
        create=True,
        min_branch_length=0.00,
        filename_node_types=None,
        format_branch_length="%6.4f",
        nodetypes_inconsistency=("InconsistentTranscripts", "Inconsistency"),
        analyze_resolution_data=None,
        warning_small_branch_length=0.01,
        filename_filter_positives=None,
        skip_without_outgroups=False,
    )

    (options, args) = E.Start(
        parser, add_database_options=True, add_csv_options=True)

    if options.outgroup_species:
        options.outgroup_species = set(options.outgroup_species.split(","))

    if options.svg_otus:
        options.svg_otus = set(options.svg_otus.split(","))

    rx_species = re.compile(options.species_regex)
    extract_species = lambda x: rx_species.match(x).groups()[0]
    if options.gene_regex:
        rx_gene = re.compile(options.gene_regex)
        extract_gene = lambda x: rx_gene.match(x).groups()[0]
    else:
        extract_gene = None

    extract_quality = lambda x: x.split(options.separator)[3]

    #########################################################################
    #########################################################################
    #########################################################################
    # read positive list of malis
    #########################################################################
    if options.filename_filter_positives:
        filter_positives, nerrors = IOTools.ReadList(
            open(options.filename_filter_positives, "r"))
        filter_positives = set(filter_positives)
    else:
        filter_positives = None

    #########################################################################
    #########################################################################
    #########################################################################
    # read location info
    #########################################################################
    if options.filename_locations:
        map_id2location = TreeReconciliation.readLocations(open(options.filename_locations, "r"),
                                                           extract_species)
    else:
        map_id2location = {}

    if (options.remove_unplaced or options.filter_location != "all") and not options.filename_locations:
        raise "please supply a file with location information."

    #########################################################################
    #########################################################################
    #########################################################################
    # delete output files
    #########################################################################
    if options.create and options.output_pattern:
        for section in ("details", "subtrees", "subids", "details", "trees", "nodes", "categories"):
            fn = options.output_pattern % section
            if os.path.exists(fn):
                if options.loglevel >= 1:
                    options.stdlog.write("# deleting file %s.\n" % fn)
                os.remove(fn)

    if options.loglevel >= 1:
        options.stdlog.write("# reading gene trees.\n")
        options.stdlog.flush()

    gene_nexus = TreeTools.Newick2Nexus(sys.stdin)

    Tree.updateNexus(gene_nexus)

    if options.loglevel >= 1:
        options.stdlog.write(
            "# read %i gene trees from stdin.\n" % len(gene_nexus.trees))
        options.stdlog.flush()

    #########################################################################
    #########################################################################
    #########################################################################
    # main loop over gene trees
    #########################################################################
    ninput, nfiltered, nskipped, noutput = 0, 0, 0, 0
    nskipped_filter, nskipped_outgroups = 0, 0

    # total counts
    total_heights_per_species = {}
    total_relheights_per_species = {}
    total_heights_per_tree = []
    total_relheights_per_tree = []

    for gene_tree in gene_nexus.trees:

        ninput += 1

        xname = re.sub("_tree.*", "", gene_tree.name)
        xname = re.sub("subtree_", "", xname)

        if filter_positives and xname not in filter_positives:
            nskipped_filter += 1
            continue

        if options.loglevel >= 6:
            gene_tree.display()

        #######################################################################
        #######################################################################
        #######################################################################
        # get identifier for this tree and update prefixes accordingly
        #######################################################################
        if options.prefix:
            if len(gene_nexus.trees) > 0:
                prefix_header = "prefix1\tprefix2\t"
                prefix_row = options.prefix + "\t" + gene_tree.name + "\t"
                prefix_prefix = options.prefix + "_" + gene_tree.name + "_"
                prefix_name = options.prefix + "_" + gene_tree.name
            else:
                prefix_header = "prefix\t"
                prefix_row = options.prefix + "\t"
                prefix_prefix = options.prefix + "_"
                prefix_name = options.prefix
        else:
            if len(gene_nexus.trees) > 0:
                prefix_header = "prefix\t"
                prefix_row = gene_tree.name + "\t"
                prefix_prefix = gene_tree.name + "\t"
                prefix_name = gene_tree.name
            else:
                prefix_header, prefix_row, prefix_prefix, prefix_name = "", "", "", ""

        #######################################################################
        #######################################################################
        #######################################################################
        # apply filters to gene tree
        #######################################################################
        TreeReconciliation.filterTree(gene_tree, options, map_id2location)

        otus = TreeTools.GetTaxa(gene_tree)

        if len(otus) <= 1:
            nfiltered += 1
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# tree %s: empty after filtering - skipped.\n" % gene_tree.name)
            continue

        this_species_list = map(extract_species, otus)
        # check, if only outgroups
        if options.outgroup_species:
            if not set(this_species_list).difference(options.outgroup_species):
                nfiltered += 1
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# tree %s: only outgroups after filtering - skipped.\n" % gene_tree.name)
                continue

            if options.skip_without_outgroups and not set(this_species_list).intersection(options.outgroup_species):
                nskipped_outgroups += 1
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# tree %s: no outgroups - skipped.\n" % gene_tree.name)
                continue

        #######################################################################
        #######################################################################
        #######################################################################
        # reroot gene tree, if outgroups have been given.
        #######################################################################
        if options.outgroup_species:
            TreeReconciliation.rerootTree(gene_tree, extract_species, options)

        #######################################################################
        #######################################################################
        #######################################################################
        # compute distance to root for each node
        #######################################################################
        distance_to_root = TreeTools.GetDistanceToRoot(gene_tree)

        #######################################################################
        #######################################################################
        #######################################################################
        # compute counts
        #######################################################################
        # heights per tree
        heights_per_tree = []
        # relative heights per tree
        relheights_per_tree = []
        # distance to root
        heights_per_species = {}
        # distance to root (relative to maximum distance to root)
        relheights_per_species = {}

        analysis_set, gene_set, pseudogene_set, other_set = TreeReconciliation.getAnalysisSets(
            gene_tree, extract_quality, options)

        if len(analysis_set) == 0:
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# tree %s: empty analysis set - skipped.\n" % gene_tree.name)
            nskipped += 1
            continue

        reference_height = TreeReconciliation.getReferenceHeight(distance_to_root,
                                                                 gene_tree,
                                                                 gene_set,
                                                                 options,
                                                                 extract_species,
                                                                 method="median")

        if reference_height is None:
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# tree %s: reference height not computable or 0 - skipped.\n" % gene_tree.name)
            nskipped += 1
            continue

        for node_id in analysis_set:

            node = gene_tree.node(node_id)
            species = extract_species(node.data.taxon)
            height = distance_to_root[node_id]

            if height < options.warning_small_branch_length:
                options.stdlog.write("# tree %s: small distance %s to root at node %i: %s\n" %
                                     (gene_tree.name, options.format_branch_length % height,
                                      node_id, node.data.taxon))

            relheight = height / reference_height
            try:
                heights_per_species[species].append(height)
            except KeyError:
                heights_per_species[species] = [height]
                relheights_per_species[species] = []

            relheights_per_species[species].append(relheight)

            # do not use outgroup species
            if options.outgroup_species and species in options.outgroup_species:
                continue

            heights_per_tree.append(height)
            relheights_per_tree.append(relheight)

        if options.loglevel >= 1:
            options.stdlog.write("# tree %s: reference_height=%s\n" % (
                gene_tree.name, options.format_branch_length % reference_height))
            options.stdlog.flush()

        if options.print_subtotals:
            printCounts(heights_per_species, relheights_per_species,
                        heights_per_tree, relheights_per_tree,
                        options, prefix_header, prefix_row)

        #######################################################################
        #######################################################################
        #######################################################################
        # update total counts
        #######################################################################
        TreeReconciliation.appendCounts(
            total_heights_per_species, heights_per_species)
        TreeReconciliation.appendCounts(
            total_relheights_per_species, relheights_per_species)

        TreeReconciliation.appendCounts(
            total_heights_per_tree, heights_per_tree)
        TreeReconciliation.appendCounts(
            total_relheights_per_tree, relheights_per_tree)

        noutput += 1

    if options.print_totals:

        if options.prefix:
            prefix_header = "prefix1\tprefix2\t"
            prefix_row = options.prefix + "\t" + "total" + "\t"
            prefix_prefix = options.prefix + "_" + "total" + "_"
            prefix_name = options.prefix + "_" + "total"
        else:
            prefix_header = "prefix\t"
            prefix_row = "total" + "\t"
            prefix_prefix = "total" + "_"
            prefix_name = "total"

        printCounts(total_heights_per_species, total_relheights_per_species,
                    total_heights_per_tree, total_relheights_per_tree,
                    options, prefix_header, prefix_row)

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, nfiltered=%i, nskipped=%i, nskipped_filter=%i, nskipped_outgroups=%i, noutput=%i\n" % (
            ninput, nfiltered, nskipped, nskipped_filter, nskipped_outgroups, noutput))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
