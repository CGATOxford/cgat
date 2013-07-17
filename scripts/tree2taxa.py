################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree2taxa.py 2782 2009-09-10 11:40:29Z andreas $
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
tree2taxa.py - extract taxa in a tree
=====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collections of trees from stdin
and for each tree outputs the taxa found within the
tree.

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
import optparse

from types import *
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree2taxa.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option( "--skip-trees", dest="skip_trees", action="store_true",
                       help="do not output tree names in third field [default=%default]."  )

    parser.set_defaults(
        skip_trees = False
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    nexus = TreeTools.Newick2Nexus( sys.stdin )
    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees from stdin.\n" % len(nexus.trees))

    ntree = 0
    ntotal = len(nexus.trees)

    if ntotal == 1:
        options.stdout.write("taxon\n")
    else:
        if options.skip_trees:
            options.stdout.write("taxon\ttree\n")
        else:
            options.stdout.write("taxon\ttree\tname\n")
            
    for tree in nexus.trees:
        ntree += 1
        taxa = TreeTools.GetTaxa( tree )

        if ntotal == 1:
            for t in taxa:
                options.stdout.write("%s\n" % (t))
        elif options.skip_trees:
            for t in taxa:
                options.stdout.write("%s\t%i\n" % (t, ntree))
        else:
            for t in taxa:
                options.stdout.write("%s\t%i\t%s\n" % (t, ntree, tree.name))

    if options.loglevel >= 1:
        options.stdlog.write( "# ntotal=%i\n" % (ntotal))
        
    E.Stop()
