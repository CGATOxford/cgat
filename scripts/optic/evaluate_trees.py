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
optic/evaluate_trees.py -
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

   python optic/evaluate_trees.py --help

Type::

   python optic/evaluate_trees.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

USAGE = """python %s [OPTIONS] < tree.in > tree.out

Version: $Id: optic/evaluate_trees.py 2781 2009-09-10 11:33:14Z andreas $

compare a set of trees to reference trees.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --pattern-species=          regex pattern to extract species from identifier
-g, --genes-tsv-file=                    filename with list of genes per species
""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/evaluate_trees.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-r", "--reference=", dest="filename_reference_tree",
                      help="filename with reference tree.", type="string")

    parser.set_defaults(
        filename_reference_tree=None
    )

    (options, args) = E.Start(parser)

    if not options.filename_reference_tree:
        print "please supply reference tree."

    if options.loglevel >= 1:
        print "# reading reference tree."

    nexus = TreeTools.Newick2Nexus(open(options.filename_reference_tree, "r"))
    reference_tree = nexus.trees[0]

    if options.loglevel >= 1:
        print "# reading sample trees."

    nexus2 = TreeTools.Newick2Nexus(sys.stdin)

    ntotal, nok, nfailed = 0, 0, 0
    ntopology, ntaxa, nleaves = 0, 0, 0
    for t in nexus2.trees:
        ntotal += 1
        is_ok, reason = TreeTools.IsCompatible(reference_tree, t)
        if is_ok:
            nok += 1
        else:
            nfailed += 1
            if reason == "topology":
                ntopology += 1
            elif reason == "taxa":
                ntaxa += 1
            elif reason == "leaves":
                nleaves += 1

    print "# total=%i, compatible=%i, failed=%i, topology=%i, taxa=%i, leaves=%i" %\
          (ntotal, nok, nfailed, ntopology, ntaxa, nleaves)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
