'''
tree_species2genes.py - collapse leaves with the same species.
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

-p, --pattern-species=          regex pattern to extract species from identifier
-g, --genes-tsv-file=                    filename with list of genes per species

Usage
-----

Example::

   python tree_species2genes.py --help

Type::

   python tree_species2genes.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import getopt
import tempfile
import time
import popen2

from Bio.Nexus import Nexus
from Bio.Nexus.Nodes import Node

import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments
import CGAT.Genomics as Genomics

param_loglevel = 1

param_long_options = ["verbose=", "help",
                      "pattern-species=",
                      "genes=", "version"]

param_short_options = "v:hp:g:"

param_pattern_species = "^([^@:]+)[@:]"

param_filename_genes = None


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ("-p", "--pattern-species"):
            param_pattern_species = a
        elif o in ("-g", "--genes-tsv-file"):
            param_filename_genes = a

    if not param_filename_genes:
        print "please specify a filename with genes"
        print globals()["__doc__"]
        sys.exit(1)

    print E.GetHeader()
    print E.GetParams()

    rx = re.compile(param_pattern_species)
    infile = open(param_filename_genes, "r")
    map_species2genes = {}
    for line in infile:
        if line[0] == "#":
            continue
        gene = line[:-1].split("\t")[0]
        species = rx.search(gene).groups()[0]
        if species not in map_species2genes:
            map_species2genes[species] = []
        map_species2genes[species].append(gene)

    nexus = TreeTools.Newick2Nexus(sys.stdin)

    TreeTools.Species2Genes(nexus, map_species2genes)

    print TreeTools.Nexus2Newick(nexus)

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
