################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree2matrix.py 2782 2009-09-10 11:40:29Z andreas $
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
tree2matrix.py - convert a tree to a distance matrix
====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert one or more trees into distance matrix. The nodes
in the tree between distances are output is defined by
the --pairs option. Possible values are:

+--------------------+----------------------------------------------------------------------+
|                    |compute distances between all nodes in the tree                       |
|all                 |                                                                      |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+
|leaves              |compute distance between all leaves in the tree                       |
|                    |                                                                      |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+
|branches            |computed distances between each branches                              |
|                    |                                                                      |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+
|terminals           |get distance of each leave to its parent                              |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+
|lineage             |get distance of each leaf to the next speciation node                 |
|                    |                                                                      |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+
|between-species     |get distance between different species, excluding same species nodes  |
|                    |                                                                      |
|                    |                                                                      |
|                    |                                                                      |
+--------------------+----------------------------------------------------------------------+

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
""" program $Id: tree2matrix.py 2782 2009-09-10 11:40:29Z andreas $

"""
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools

def TranslateNode( node, tree, terminals, options ):
    if options.do_translate:
        return options.separator.join( sorted(TreeTools.GetLeaves( tree, node )) )
    elif node in terminals:
        return tree.node( node ).data.taxon
    else:
        return str(node)

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree2matrix.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="number format to use."  )
    parser.add_option("-g", "--graph", dest="to_graph", action="store_true",
                      help="convert tree(s) to graph(s)."  )
    parser.add_option("-a", "--table", dest="to_table", action="store_true",
                      help="convert tree(s) to table."  )
    parser.add_option("-t", "--translate", dest="do_translate", action="store_true",
                      help="translate internal nodes to clades."  )
    parser.add_option( "--output-pattern", dest="output_filename_pattern", type="string",
                       help="pattern for output file if there are multiple trees in the file.""")
    parser.add_option( "--pairs", dest="pairs", type="choice",
                       choices=("all", "leaves", "branches", "terminals", "lineage", "between-species" ),
                       help="choose pairs of nodes to output.""")
    parser.add_option( "--species", dest="species", type="string",
                       help="comma separated list of species that are considered. All others are ignored." )


    parser.set_defaults(
        format = "%6.4f",
        to_graph = False,
        to_table = False,
        do_translation = False,
        separator = ":",
        do_all_on_all = False,
        do_branches = False,
        do_terminals = False,
        output_filename_pattern = None,
        pairs = "branches",
        species = None,
        regex_species = ("^([^|]+)" ),
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if options.species: options.species = options.species.split(",")
    nexus = TreeTools.Newick2Nexus( sys.stdin )

    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees from stdin.\n" % len(nexus.trees))

    ntree = 0
    outfile = None

    ## The table is a hash of lists
    table = {}

    extract_species = lambda x: re.search( options.regex_species, x).groups()[0]

    for tree in nexus.trees:

        if len(nexus.trees) == 1:
            outfile = options.stdout
        elif options.output_filename_pattern:
            ntree += 1
            if outfile != None: outfile.close()            
            outfile = open( options.output_filename_pattern % ntree, "w" )
        else:
            outfile = options.stdout
            
        ## prune tree, if an explicit species list is given
        if options.species:
            species = set(options.species)
            terminals = tree.get_terminals()
            for x in terminals:
                taxon = tree.node(x).data.taxon
                if extract_species(taxon) not in species:
                    tree.prune( taxon )

        ## define node list
        terminals = tree.get_terminals()
        set_terminals = set(terminals)
        node_list = []
        
        if options.pairs == "all":
            nodes = TreeTools.GetAllNodes(tree)
            for x in range(len(nodes)):
                for y in range(0,x):
                    node_list.append( (nodes[x], nodes[y]) )
        elif options.pairs == "terminals":
            for x in terminals:
                node_list.append( (x, tree.node(x).prev) )
        elif options.pairs == "leaves":
            nodes = terminals
            for x in range(len(nodes)):
                for y in range(0,x):
                    node_list.append( (nodes[x], nodes[y]) )
        elif options.pairs == "branches":
            nodes = TreeTools.GetAllNodes(tree)
            for x in range(len(nodes)):
                if tree.node(x).prev:
                    node_list.append( (x, tree.node(x).prev) )
        elif options.pairs == "between-species":
            nodes = terminals
            for x in range(len(nodes)):
                for y in range(0,x):
                    s1 = extract_species( tree.node(nodes[x]).data.taxon) 
                    s2 = extract_species( tree.node(nodes[y]).data.taxon) 
                    if s1 != s2:
                        node_list.append( (nodes[x], nodes[y]) )

        elif options.pairs == "lineage":
            raise "not implemented."
        
        if options.to_graph:
            outfile.write( "node1\tnode2\tdistance\n" )
            
            links = TreeTools.Tree2Graph( tree )
            for n1, n2, weight in links:

                node1 = TranslateNode( n1, tree, set_terminals, options )
                node2 = TranslateNode( n2, tree, set_terminals, options )

                if node1 > node2: node1, node2 = node2, node1
                outfile.write( "%s\t%s\t%s\n" % (
                    node1, node2,
                    options.format % weight ) )

        elif options.to_table:
            
            if options.do_all_on_all:
                nodes = TreeTools.GetAllNodes(tree)
            else:
                nodes = terminals

            for n1, n2 in node_list:
                
                node1 = TranslateNode( n1, tree, set_terminals, options )
                node2 = TranslateNode( n2, tree, set_terminals, options )

                if node1 > node2: node1, node2 = node2, node1

                if options.do_terminals:
                    key = "%s" % node2
                else:
                    key = "%s-%s" % (node1, node2)
                    
                if key not in table: table[key] = []

                table[key].append( options.format % tree.distance( n1, n2 ) )
        else:
            outfile.write( "node1\tnode2\tdistance\n" )
            
            for n1, n2 in node_list:
                node1 = TranslateNode( n1, tree, set_terminals, options )
                node2 = TranslateNode( n2, tree, set_terminals, options )
                
                if node1 > node2: node1, node2 = node2, node1
                
                outfile.write( "%s\t%s\t%s\n" % ( \
                        node1, node2, 
                        options.format % tree.distance( n1, n2 )))

    if options.to_table:
        outfile = sys.stdout
        outfile.write("branch\t%s\n" % ("\t".join(map(str, range(0,len(nexus.trees))))))
                                                      
        for key, values in table.items():
            outfile.write( "%s\t%s\n" % (key, "\t".join(values)))
        
    if outfile != sys.stdout:
        outfile.close()
    
    E.Stop()
