'''
tree_map_leaves.py - 
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

   python tree_map_leaves.py --help

Type::

   python tree_map_leaves.py --help

for command line help.

Command line options
--------------------

'''
import sys
import getopt
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

USAGE = """python %s [OPTIONS] < tree.in > tree.out

Version: $Id: tree_map_leaves.py 2782 2009-09-10 11:40:29Z andreas $

collapse leafs with the same species.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-a, --map-tsv-file=                    apply mapping in file
-c, --create=                   create mapping in file
-i, --invert=                   invert map
-b, --strip-branches            remove branch lengths
""" % sys.argv[0]

param_loglevel = 1

param_long_options = ["verbose=", "help", "pattern-species=",
                      "invert", "create=", "apply=", "strip-brances",
                      "version"]

param_short_options = "v:hp:a:c:i"

param_apply = None
param_create = None
param_invert = False
param_remove_branch_lengths = False


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
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-a", "--map-tsv-file"):
            param_apply = a
        elif o in ("-c", "--create"):
            param_create = a
        elif o in ("-i", "--invert"):
            param_invert = True
        elif o in ("-b", "--strip-branches"):
            param_remove_branch_lengths = True

    print E.GetHeader()
    print E.GetParams()

    keys = {}
    if param_apply:
        infile = open(param_apply, "r")
        for line in infile:
            if line[0] == "#":
                continue
            a, b = line[:-1].split("\t")[:2]
            if param_invert:
                a, b = b, a
            keys[a] = b

    nexus = TreeTools.Newick2Nexus(sys.stdin)

    notu = 0

    for tree in nexus.trees:
        if param_loglevel >= 2:
            tree.display()

        for nx in tree.get_terminals():
            t1 = tree.node(nx).get_data().taxon

            if param_create:
                if t1 not in keys:
                    keys[t1] = "otu%i" % notu
                    notu += 1

            if t1 in keys:
                tree.node(nx).get_data().taxon = keys[t1]

    print TreeTools.Nexus2Newick(nexus)

    if param_create:
        outfile = open(param_create, "w")
        for key in keys:
            outfile.write("%s\t%s\n" % (key, keys[key]))
        outfile.close()

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
