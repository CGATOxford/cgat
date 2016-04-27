'''
lca2table.py
====================================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Summarise results of LCA analysis - from mtools lcamapper.sh

Usage
-----

Example::

   python lca2table.py < infile > outfile

Type::

   python lca2table.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.LCA as LCA
import CGAT.Experiment as E
import collections


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--summarise", dest="summarise", type="choice",
                      choices=("level-counts", "taxa-counts", "individual"),
                      help="summarise the taxa counts - no. phyla etc")

    parser.add_option("--output-map", dest="output_map", action="store_true",
                      help="ouput map of taxonomy")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.output_map:
        found = []
        options.stdout.write("""kingdom\t
        phylum\t
        class\t
        order\t
        family\t
        genus\t
        species\n""")
        # only output the mapping file - do not continue
        # summarise regardless of the specified options
        for lca in LCA.iterate(options.stdin):
            hierarchy = [lca.kingdom,
                         lca.phylum,
                         lca._class,
                         lca.order,
                         lca.family,
                         lca.genus,
                         lca.species]
            if hierarchy in found:
                continue
            else:
                found.append(hierarchy)
                options.stdout.write("\t".join(hierarchy) + "\n")
        return

    if options.summarise == "level-counts":
        level_counts = collections.defaultdict(set)
        total = 0
        nreads_kingdom = 0
        nreads_kingdom_plus = 0
        nreads_phylum = 0
        nreads_phylum_plus = 0
        nreads_class = 0
        nreads_class_plus = 0
        nreads_order = 0
        nreads_order_plus = 0
        nreads_family = 0
        nreads_family_plus = 0
        nreads_genus = 0
        nreads_genus_plus = 0
        nreads_species = 0
        nreads_species_plus = 0
        nreads_subspecies = 0
        nreads_subspecies_plus = 0

        c = E.Counter()
        for lca in LCA.iterate(options.stdin):
            total += 1
            if lca.kingdom != "NA":
                nreads_kingdom += 1
                level_counts["kingdom"].add(lca.kingdom)
            else:
                c.kingdom_unmapped += 1

            if lca.kingdom_plus != "NA":
                nreads_kingdom_plus += 1
                level_counts["kingdom+"].add(lca.kingdom_plus)
            else:
                c.kingdom_plus_unmapped += 1

            if lca.phylum != "NA":
                nreads_phylum += 1
                level_counts["phylum"].add(lca.phylum)
            else:
                c.phylum_unmapped += 1

            if lca.phylum_plus != "NA":
                nreads_phylum_plus += 1
                level_counts["phylum+"].add(lca.phylum_plus)
            else:
                c.phylum_plus_unmapped += 1

            if lca._class != "NA":
                nreads_class += 1
                level_counts["class"].add(lca._class)
            else:
                c.class_unmapped += 1

            if lca._class_plus != "NA":
                nreads_class_plus += 1
                level_counts["class+"].add(lca._class_plus)
            else:
                c.class_plus_unmapped += 1

            if lca.order != "NA":
                nreads_order += 1
                level_counts["order"].add(lca.order)
            else:
                c.order_unmapped += 1

            if lca.order_plus != "NA":
                nreads_order_plus += 1
                level_counts["order+"].add(lca.order_plus)
            else:
                c.order_plus_unmapped += 1

            if lca.family != "NA":
                nreads_family += 1
                level_counts["family"].add(lca.family)
            else:
                c.family_unmapped += 1

            if lca.family != "NA":
                nreads_family_plus == 1
                level_counts["family+"].add(lca.family_plus)
            else:
                c.family_plus_unmapped += 1

            if lca.genus != "NA":
                nreads_genus += 1
                level_counts["genus"].add(lca.genus)
            else:
                c.genus_unmapped += 1

            if lca.genus_plus != "NA":
                nreads_genus_plus == 1
                level_counts["genus+"].add(lca.genus_plus)
            else:
                c.genus_plus_unmapped += 1

            if lca.species != "NA":
                nreads_species += 1
                level_counts["species"].add(lca.species)
            else:
                c.species_unmapped += 1

            if lca.species_plus != "NA":
                nreads_species_plus += 1
                level_counts["species+"].add(lca.species_plus)
            else:
                c.species_plus_unmapped += 1

            if lca.subspecies != "NA":
                nreads_subspecies += 1
                level_counts["subspecies"].add(lca.subspecies)
            else:
                c.subspecies_unmapped += 1

            if lca.subspecies_plus != "NA":
                nreads_subspecies_plus += 1
                level_counts["subspecies+"].add(lca.subspecies_plus)
            else:
                c.subspecies_plus_unmapped += 1

        options.stdout.write("\t".join(["nkingdom",
                                        "nkingdom+",
                                        "nphylum",
                                        "nphylum+",
                                        "nclass",
                                        "nclass+",
                                        "norder",
                                        "norder+",
                                        "nfamily",
                                        "nfamily+",
                                        "ngenus",
                                        "ngenus+",
                                        "nspecies",
                                        "nspecies+",
                                        "nsubspecies",
                                        "nsubspecies+",
                                        "nseqkingdom",
                                        "nseqkingdom+",
                                        "nseqphylum",
                                        "nseqphylum+",
                                        "nseqclass",
                                        "nseqclass+",
                                        "nseqorder",
                                        "nseqorder+",
                                        "nseqfamily",
                                        "nseqfamily+",
                                        "nseqgenus",
                                        "nseqgenus+",
                                        "nseqspecies",
                                        "nseqspecies+",
                                        "nseqsubspecies",
                                        "nseqsubspecies+"]) + "\n")

        options.stdout.write("\t".join(map(
            str, [len(level_counts["kingdom"]),
                  len(level_counts["kingdom+"]),
                  len(level_counts["phylum"]),
                  len(level_counts["phylum+"]),
                  len(level_counts["class"]),
                  len(level_counts["class+"]),
                  len(level_counts["order"]),
                  len(level_counts["order+"]),
                  len(level_counts["family"]),
                  len(level_counts["family+"]),
                  len(level_counts["genus"]),
                  len(level_counts["genus+"]),
                  len(level_counts["species"]),
                  len(level_counts["species+"]),
                  len(level_counts["subspecies"]),
                  len(level_counts["subspecies+"]),
                  nreads_kingdom,
                  nreads_phylum,
                  nreads_phylum_plus,
                  nreads_class,
                  nreads_class_plus,
                  nreads_order,
                  nreads_order_plus,
                  nreads_family,
                  nreads_family_plus,
                  nreads_genus,
                  nreads_genus_plus,
                  nreads_species,
                  nreads_species_plus,
                  nreads_subspecies,
                  nreads_subspecies_plus])) + "\n")
    elif options.summarise == "taxa-counts":
        unmapped = collections.defaultdict(int)
        total = 0
        taxa_counts = {"kingdom": collections.defaultdict(int),
                       "kingdom+": collections.defaultdict(int),
                       "phylum": collections.defaultdict(int),
                       "phylum+": collections.defaultdict(int),
                       "class": collections.defaultdict(int),
                       "class+": collections.defaultdict(int),
                       "order": collections.defaultdict(int),
                       "order+": collections.defaultdict(int),
                       "family": collections.defaultdict(int),
                       "family+": collections.defaultdict(int),
                       "genus": collections.defaultdict(int),
                       "genus+": collections.defaultdict(int),
                       "species": collections.defaultdict(int),
                       "species+": collections.defaultdict(int),
                       "subspecies": collections.defaultdict(int),
                       "subspecies+": collections.defaultdict(int)}

        c = E.Counter()
        for lca in LCA.iterate(options.stdin):
            total += 1
            if lca.kingdom != "NA":
                taxa_counts["kingdom"][lca.kingdom] += 1
            else:
                c.kingdom_unmapped += 1
                unmapped["kingdom"] += 1
            if lca.kingdom_plus != "NA":
                taxa_counts["kingdom+"][lca.kingdom_plus] += 1
            else:
                c.kingdom_plus_unmapped += 1
                unmapped["kingdom+"] += 1
            if lca.phylum != "NA":
                taxa_counts["phylum"][lca.phylum] += 1
            else:
                c.phylum_unmapped += 1
                unmapped["phylum"] += 1
            if lca.phylum_plus != "NA":
                taxa_counts["phylum+"][lca.phylum_plus] += 1
            else:
                c.phylum_plus_unmapped += 1
                unmapped["phylum+"] += 1
            if lca._class != "NA":
                taxa_counts["class"][lca._class] += 1
            else:
                c.class_unmapped += 1
                unmapped["class"] += 1
            if lca._class_plus != "NA":
                taxa_counts["class+"][lca._class_plus] += 1
            else:
                c.class_plus_unmapped += 1
                unmapped["class+"] += 1
            if lca.order != "NA":
                taxa_counts["order"][lca.order] += 1
            else:
                c.order_unmapped += 1
                unmapped["order"] += 1
            if lca.order_plus != "NA":
                taxa_counts["order+"][lca.order_plus] += 1
            else:
                c.order_plus_unmapped += 1
                unmapped["order+"] += 1
            if lca.family != "NA":
                taxa_counts["family"][lca.family] += 1
            else:
                c.family_unmapped += 1
                unmapped["family"] += 1
            if lca.family_plus != "NA":
                taxa_counts["family+"][lca.family_plus] += 1
            else:
                c.family_plus_unmapped += 1
                unmapped["family+"] += 1
            if lca.genus != "NA":
                taxa_counts["genus"][lca.genus] += 1
            else:
                c.genus_unmapped += 1
                unmapped["genus"] += 1
            if lca.genus_plus != "NA":
                taxa_counts["genus+"][lca.genus_plus] += 1
            else:
                c.genus_plus_unmapped += 1
                unmapped["genus+"] += 1
            if lca.species != "NA":
                taxa_counts["species"][lca.species] += 1
            else:
                c.species_unmapped += 1
                unmapped["species"] += 1
            if lca.species_plus != "NA":
                taxa_counts["species+"][lca.species_plus] += 1
            else:
                c.species_plus_unmapped += 1
                unmapped["species+"] += 1
            if lca.subspecies != "NA":
                taxa_counts["subspecies"][lca.subspecies] += 1
            else:
                c.subspecies_unmapped += 1
                unmapped["subspecies"] += 1
            if lca.subspecies_plus != "NA":
                taxa_counts["subspecies+"][lca.subspecies_plus] += 1
            else:
                c.subspecies_plus_unmapped += 1
                unmapped["subspecies+"] += 1

        options.stdout.write("level\ttaxa\tcount\tproportion\trpm\n")
        for level, taxa_count in taxa_counts.iteritems():
            total_level = total - unmapped[level]
            for taxa, count in taxa_count.iteritems():
                options.stdout.write("\t".join(
                    [level,
                     taxa,
                     str(count),
                     str(float(count)/total_level),
                     str(float(count)/(float(total_level)/1000000))]) + "\n")

        E.info(c)

    elif options.summarise == "individual":
        # each read is output with its respective
        # taxon assignments
        options.stdout.write("\t".join(["id",
                                        "kingdom",
                                        "kingdom+",
                                        "phylum",
                                        "phylum+",
                                        "class",
                                        "class+",
                                        "order",
                                        "order+",
                                        "family",
                                        "family+",
                                        "genus",
                                        "genus+",
                                        "species",
                                        "species+",
                                        "subspecies",
                                        "subspecies+"]) + "\n")
        for lca in LCA.iterate(options.stdin):
            options.stdout.write("\t".join([lca.identifier,
                                            lca.kingdom,
                                            lca.kingdom_plus,
                                            lca.phylum,
                                            lca.phylum_plus,
                                            lca._class,
                                            lca._class_plus,
                                            lca.order,
                                            lca.order_plus,
                                            lca.family,
                                            lca.family_plus,
                                            lca.genus,
                                            lca.genus_plus,
                                            lca.species,
                                            lca.species_plus,
                                            lca.subspecies,
                                            lca.subspecies_plus]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
