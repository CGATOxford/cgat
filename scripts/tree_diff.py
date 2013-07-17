################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree_diff.py 2782 2009-09-10 11:40:29Z andreas $
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
tree_diff.py - compare two sets of trees
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts reads two collections of trees and checks if the
trees in the two collections are identical.

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
import optparse
import math
import tempfile
import subprocess

from types import *

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree_diff.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-1", "--filename-tree1", dest="filename_tree1", type="string",
                      help="filename with first tree(s)."  )
    parser.add_option("-2", "--filename-tree2", dest="filename_tree2", type="string",
                      help="filename with second tree(s)."  )
    parser.add_option("-o", "--outgroup", dest="outgroup", type="string",
                      help="reroot with outgroup before processing.")

    parser.set_defaults(
        filename_tree1 = None,
        filename_tree2 = None,
        outgroup = None
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if (len(args) == 2):
        options.filename_tree1, options.filename_tree2 = args

    if not options.filename_tree1 or not options.filename_tree2:
        raise ValueError( "please specify two trees." )
    
    ## take first trees
    nexus = TreeTools.Newick2Nexus( open(options.filename_tree1, "r") )
    trees1 = nexus.trees
    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees from %s.\n" % (len(trees1), options.filename_tree1))

    ## take first trees
    nexus = TreeTools.Newick2Nexus( open(options.filename_tree2, "r") )
    trees2 = nexus.trees
    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees from %s.\n" % (len(trees2), options.filename_tree2))
        
    ntotal, nsame, ndiff = 0,0,0

    if options.outgroup:
        for tree in trees1:
            tree.root_with_outgroup( options.outgroup )
        for tree in trees2:
            tree.root_with_outgroup( options.outgroup )
    
    for x in range(len(trees1)):
        for y in range(len(trees2)):
            if options.loglevel >= 2:
                print trees1[x]
                print trees2[y]
            if trees1[x].is_identical(trees2[y]):
                code = "="
                nsame += 1
            else:
                code = "<>"
                ndiff += 1
            options.stdout.write("%s\t%i\t%i\n" % (code, x, y ))
            ntotal += 1
            
    options.stdlog.write( "# n1=%i, n2=%i, ntotal=%i, nsame=%i, ndiff=%i\n" % (len(trees1), len(trees2), ntotal, nsame, ndiff))
        
    E.Stop()
