################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree2patterns.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
#################################################################################
"""
tree2patterns.py - convert trees to taxon patterns
==================================================

:Author: Andreas Heger
:Release: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collections of trees from stdin and - based
on a reference tree - outputs for each tree a pattern indicating
if a certain taxon has been found or net.

``1`` indicates that a pattern has been found, ``0`` otherwise.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

""" 

import os
import sys
import string
import re
import getopt
import time
import sets
import optparse
import math
import tempfile


from Bio.Nexus import Nexus
from Bio.Nexus.Nodes import Node

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree2patterns.py 2781 2009-09-10 11:33:14Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-t", "--reference-tree", dest="reference_tree", type="string",
                      help="reference tree to read."  )

    parser.add_option("-s", "--sort-order", dest="sort_order", type="string",
                      help="output order of OTU."  )

    parser.set_defaults(
        reference_tree = None,
        sort_order = [],
        )

    (options, args) = E.Start( parser )

    if not options.sort_order:
        for nx in reference_tree.get_terminals():
            options.sort_order.append( reference_tree.node(nx).get_data().taxon )
    else:
        options.sort_order = options.sort_order.split(",")

    if not options.reference_tree:
        raise "no reference tree defined."

    nexus = TreeTools.Newick2Nexus( options.reference_tree )
    reference_tree = nexus.trees[0]
    
    if options.loglevel >= 3:
        print "# reference tree:"
        print reference_tree.display()

    patterns = TreeTools.calculatePatternsFromTree( tree, options.sort_order )

    for p in patterns:
        print p

    E.Stop()



