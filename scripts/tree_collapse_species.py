################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree_collapse_species.py 2782 2009-09-10 11:40:29Z andreas $
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
tree_collapse_species.py - collapse single species nodes 
========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of trees from stdin and 
collapses nodes that have children of the same species.

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
import tempfile
import time
import popen2

from Bio.Nexus import Nexus

USAGE="""python %s [OPTIONS] < tree.in > tree.out

Version: $Id: tree_collapse_species.py 2782 2009-09-10 11:40:29Z andreas $



Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --pattern-species=          regex pattern to extract species from identifier
""" % sys.argv[0]

import CGAT.Experiment as E

param_loglevel = 1

param_long_options=["verbose=", "help",
                    "pattern-species=",
                    "version"]

param_short_options="v:hp:"

param_pattern_species = "^([^@:]+)[@:]"


def PruneTree( tree, id ):

    if id not in tree.get_terminals():
        raise "Not a terminal taxon: %i" % id
    
    prev=tree.unlink(id)
    tree.kill(id)
    if not prev==tree.root and len(tree.node(prev).succ)==1:
        
        succ=tree.node(prev).get_succ()[0]
        new_bl=tree.node(prev).data.branchlength+tree.node(succ).data.branchlength
        tree.collapse(prev)
        tree.node(succ).data.branchlength=new_bl
        
    return prev

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-p", "--pattern-species"):
            param_pattern_species = a

    print E.GetHeader()
    print E.GetParams()

    lines = ["#NEXUS\nBegin trees;\ntree tree = "] + sys.stdin.readlines() + ["End;"]
    nexus = Nexus.Nexus( string.join(lines, "") )

    if len(nexus.trees) != 1:
        raise "no tree found in file."
    
    tree = nexus.trees[0]
    
    if param_loglevel >= 2:
        tree.display()

    rx = re.compile( param_pattern_species )
    changed = True
    
    while changed:
        changed = False
        leaves = tree.get_terminals()
        for x in range(0, len(leaves) -1 ):
            nx = leaves[x]
            t1 = tree.node(nx).get_data().taxon
            s1 = rx.search( t1 ).groups()[0]
            p1 = tree.node(nx).get_prev()            
            for y in range( x + 1, len(leaves)):
                ny = leaves[y]
                t2 = tree.node(ny).get_data().taxon                
                s2 = rx.search( t2 ).groups()[0]
                p2 = tree.node(ny).get_prev()
                if s1 == s2 and tree.is_monophyletic( (nx, ny) ) != -1:
                    print "collapsing nodes", t1, t2, nx, ny, p1, p2                    
                    PruneTree( tree, nx )
                    if param_loglevel >= 2:
                        tree.display()
                    changed = True
        
    print E.GetFooter()
