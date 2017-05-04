'''
Metagenemark.py
=============================================

:Tags: Python

Purpose
-------

Parses the 'gff' output that is not strictly in gff format
from MetaGeneMark analysis into gff, fasta and amino acid
sequence files

Usage
-----

Example::

   python formatMetagenemark.py --help

Type::

   python formatMetagenemark.py --help

for command line help.

Documentation
-------------

Takes input from stdin and outputs to stdout

Command line options
--------------------

'''

import os
import sys
import optparse
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                                   usage=globals()["__doc__"])

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("gff", "fasta", "aa"), help="supply help")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.format == "gff":
        os.system("grep -v '#' | grep -v '^$'")

    elif options.format == "aa":
        pattern = "Protein"

    elif options.format == "fasta":
        pattern = "DNA"

    # assume that each sequence line does not just contain the
    # amino acids A, C, T or G - THIS IS NOT OPTIMAL
    # This is a list of the amino acids that are not in either "DNA" or "A", "C", "T", "G"
    # This assumes that each sequence (line) will contain at least one of these
    amino_acids = ["M", "R", "Q", "E", "H", "I",
                   "L", "K", "F", "P", "S", "W", "Y", "V"]
    result = []
    name = None
    for line in options.stdin.readlines():
        if not line.startswith("##") or line.find("date") != -1 or line.find("gff") != -1 or line.find("source") != -1:
            continue
        data = line[2:-1]
        if data.startswith("%s" % pattern):
            name = data
            if result:
                options.stdout.write(
                    ">%s\n%s\n" % (prot_name, "".join(result)))
                result = []
        else:
            if pattern == "Protein":
                if "".join(map(str, [data.find(x) != -1 for x in amino_acids])).find("True") != -1 and data.find("end") == -1 and data.find("%s" % pattern):
                    result.append(data)
                    prot_name = name
            elif pattern == "DNA":
                if name:
                    if "".join(map(str, [data.find(x) != -1 for x in amino_acids])).find("True") == -1 and data.find("end") == -1 and data.find("%s" % pattern):
                        result.append(data)
                        prot_name = name
    if result:
        options.stdout.write(">%s\n%s\n" % (prot_name, "".join(result)))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
