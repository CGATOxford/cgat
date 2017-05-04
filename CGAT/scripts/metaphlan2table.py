'''
script_template.py
=============================================

:Tags: Python

Purpose
-------

Convert the output of a metaphlan analysis to a preferred table format


Usage
-----

Example::

   python metaphlan2table.py --help

Type::

   python metaphlan2table.py --help

for command line help.

Documentation
-------------

Code
----

'''

import sys
import optparse
import CGAT.Metaphlan as Metaphlan
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(
        version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
        usage=globals()["__doc__"])

    parser.add_option("-t", "--sequence-type", dest="type", type="choice",
                      choices=("read_map", "rel_ab"), help="type of file to be parsed to a table")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    assert options.type, "must specify infile type"
    if options.type == "read_map":
        options.stdout.write(
            "seq_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
        for entry in Metaphlan.read_map_iterator(sys.stdin):
            options.stdout.write("\t".join(
                [entry.seq_id, entry.kingdom, entry.phylum, entry.c_lass, entry.order, entry.family, entry.genus, entry.species]) + "\n")

    elif options.type == "rel_ab":
        options.stdout.write("taxon_level\ttaxon\trel_abundance\n")
        for entry in Metaphlan.relative_abundance_iterator(sys.stdin):
            options.stdout.write(
                "\t".join([entry.taxon_level, entry.taxon, entry.abundance]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
