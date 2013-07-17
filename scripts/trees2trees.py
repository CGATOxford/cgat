################################################################################
#   Gene prediction pipeline 
#
#   $Id: trees2trees.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2006 Andreas Heger 
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
'''
trees2trees.py - manipulate a collection of trees
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of trees from stdin and
outputs them on stdout. The trees can be filtered by
branch length, monophyly, number of taxa and many more
criteria.

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
'''


import os
import sys
import string
import re
import optparse
import math
import time
import random

USAGE="""python %s [OPTIONS]


""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.TreeTools as TreeTools

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: trees2trees.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-c", "--output-filename-map", dest="output_filename_map", type="string",
                      help="filename of map to output."  )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("filter", "split"),
                      help="method to use: filter removed trees, while split writes them to individual files. DEFAULT=%default"  )

    parser.add_option("-d", "--output-pattern", dest="output_pattern", type="string",
                      help="filename pattern for output multiple alignment files."  )

    parser.add_option("--filter-terminal-max-length", dest="filter_max_length", type="float",
                      help="remove terminal branches with a branch length larger than this."  )

    parser.add_option("--filter-terminal-min-length", dest="filter_min_length", type="float",
                      help="remove any branches with a branch length smaller than this."  )

    parser.add_option("--filter-min-length", dest="filter_min_length", type="float",
                      help="remove terminal branches with a branch length smaller than this."  )

    parser.add_option("--filter-max-length", dest="filter_min_length", type="float",
                      help="remove any branches with a branch length smaller than this."  )

    parser.add_option("--filter-by-trees", dest="filter_by_trees", type="string", action="append",
                      help="mask branches according to trees. Give filenames with mask trees. These trees need to have the same names and structure as the input trees, but can be in any order."  )

    parser.add_option("--filter-by-monophyly", dest="filter_by_monophyly", type="string",
                      help="only retain trees where the given taxa are monphyletic. Supply taxa as a comma-separated list."  )

    parser.add_option("--min-support", dest="min_support", type="float",
                      help="for monophyly filtering, only accept trees with minimum support."  )

    parser.add_option("--filter-ntaxa", dest="filter_ntaxa", type="int", 
                      help="filter by number of taxa."  )

    parser.add_option("--filter-simple-orthologs", dest="filter_simple_orthologs", action="store_true", 
                      help="filter for trees for simple orhtologs. This works by counting the number of taxa."  )

    parser.add_option("--filter", dest="filter", type="choice",
                      choices=("taxa", "trees"),
                      help="filter removes taxa or whole trees." )

    parser.set_defaults(
        output_pattern="%s.tree",
        output_filename_map = None,
        filter_terminal_max_length = None,
        filter_terminal_min_length = None,
        filter_max_length = None,
        filter_min_length = None,
        method ="split",
        filter = "taxa",
        filtered_branch_length = -999,
        filter_by_trees = [],
        filter_by_monophyly = None,
        filter_ntaxa = None,
        filter_simple_orthologs = None,
        min_support = 0.0,
        regex_species = ("^([^|]+)" ),
        )

    (options, args) = E.Start( parser )

    nexus = TreeTools.Newick2Nexus( sys.stdin )
    
    if options.loglevel >= 1:
        options.stdlog.write("# read %i trees from stdin.\n" % len(nexus.trees))

    ninput, noutput, nskipped = 0, 0, 0
    ndiscarded = 0
    ndiscarded_taxa = 0
    ndiscarded_branches = 0

    extract_species = lambda x: re.search( options.regex_species, x).groups()[0]
    
    if options.filter_by_trees:
        nexus_filter = []
        nexus_maps = []
        for filename in options.filter_by_trees:
            nexus_filter.append( TreeTools.Newick2Nexus( open( filename, "r") ) )
            trees = nexus_filter[-1].trees
            if options.loglevel >=1 :
                options.stdlog.write("# read %i trees for filtering from %s\n" % (len(trees), filename))

            nexus_map = {}
            for x in range( len(trees)):
                nexus_map[trees[x].name] = x
            nexus_maps.append( nexus_map )

    if options.filter_by_monophyly:
        monophyly_taxa = options.filter_by_monophyly.split(",")
        if len(monophyly_taxa) == 0:
            raise "please supply at least two taxa for the monophyly test."
            
    if options.output_filename_map:
        outfile_map = open(options.output_filename_map, "a" )
    else:
        outfile_map = None

    for tree in nexus.trees:

        ninput += 1
        id = tree.name
        has_discarded = False

        if options.filter_ntaxa != None:

            ntaxa = len(tree.get_terminals())
            if ntaxa != options.filter_ntaxa:
                if options.loglevel >= 2:
                    options.stdlog.write("# tree %s: removed because number of taxa (%i) different\n" % \
                                         (id, ntaxa ) )
                has_discarded = True
                
        if options.filter_simple_orthologs:
            ntaxa = len(tree.get_terminals())
            nspecies = len(set(map( lambda x: extract_species(tree.node(x).data.taxon), tree.get_terminals() )))
            if nspecies != ntaxa:
                if options.loglevel >= 2:
                    options.stdlog.write("# tree %s: removed because not a simple ortholog cluster: ntaxa!=nspecies (%i!=%i)\n" % \
                                             (id, ntaxa, nspecies ) )

                has_discarded = True

        if options.filter_terminal_max_length != None:
            for x in tree.get_terminals():
                node = tree.node(x)
                if node.data.branchlength >= options.filter_terminal_max_length:
                    has_discarded = True
                    ndiscarded_taxa += 1                    
                    tree.prune( node.data.taxon )
                    if options.loglevel >= 2:
                        options.stdlog.write("# tree %s: removed taxon %s because terminal branchlength to large: %s\n" % \
                                             (id, node.data.taxon, str(node.data.branchlength)) )

        if options.filter_terminal_min_length != None:
            for x in tree.get_terminals():
                node = tree.node(x)
                if node.data.branchlength <= options.filter_terminal_min_length:
                    has_discarded = True
                    ndiscarded_taxa += 1                    
                    tree.prune( node.data.taxon )
                    if options.loglevel >= 2:
                        options.stdlog.write("# tree %s: removed taxon %s because terminal branchlength to small: %s\n" % \
                                             (id, node.data.taxon, str(node.data.branchlength)) )
                    
        if options.filter_max_length != None:
            for x in tree.get_nodes(tree.root):
                if x == tree.root: continue                
                node = tree.node(x)
                if node.data.branchlength >= options.filter_max_length:
                    has_discarded = True
                    ndiscarded_branches += 1                    
                    if options.loglevel >= 2:
                        options.stdlog.write("# tree %s: removed branch %i because branchlength to large: %s\n" % \
                                             (id, x, tree.name, str(node.data.branchlength)) )
                    node.data.branchlength = options.filtered_branch_length
                    
        if options.filter_min_length != None:
            for x in tree.get_nodes(tree.root):
                if x == tree.root: continue
                node = tree.node(x)
                if node.data.branchlength <= options.filter_min_length:
                    has_discarded = True
                    ndiscarded_branches += 1
                    if options.loglevel >= 2:
                        options.stdlog.write("# tree %s: removed branch %i because internal branchlength too small: %s\n" % \
                                             (id, x, str(node.data.branchlength)) )
                    node.data.branchlength = options.filtered_branch_length
                    
        if options.filter_by_trees:
            found = []
            for y in range(len(nexus_maps)):
                if id in nexus_maps[y]:
                    found.append( (y, nexus_filter[y].trees[nexus_maps[y][id]]) )

            if not found:
                ndiscarded += 1
                continue

            for x in tree.get_nodes(tree.root):
                if x == tree.root: continue
                for y, other_tree in found:
                    other_node = other_tree.node( x )
                    if other_node.data.branchlength == options.filtered_branch_length:
                        node = tree.node(x)
                        if options.loglevel >= 2:
                            options.stdlog.write("# tree %s: removed branch %i because internal branchlength masked by tree %i:%s.\n" % \
                                                 (id, x, y, other_tree.name) )
                        
                        node.data.branchlength = options.filtered_branch_length
                        has_discarded = True
                        ndiscarded_branches += 1
                        break

        if options.filter_by_monophyly:

            terminals = set(map( lambda x: tree.node(x).data.taxon, tree.get_terminals()))
            
            for t in monophyly_taxa:
                if t not in terminals:
                    if options.loglevel >= 2:
                        options.stdlog.write( "taxon %s not in tree %s\n" % (t, tree.name))
                    nskipped += 1
            succ = tree.node(tree.root).succ
            ## use minimum support at root, if it is not the same (if trees
            ## are rooted)
            if len(succ) == 2:
                m = min( map( lambda x: tree.node(x).data.support, succ) )
                for x in succ:
                    tree.node(x).data.support = m
                
            if not TreeTools.IsMonophyleticForTaxa( tree, monophyly_taxa, support=options.min_support ):
                ndiscarded += 1
                continue
            
        if has_discarded:
            ndiscarded += 1
            if options.filter=="trees" or options.filter_ntaxa:
                continue

        if options.method == "split":

            output_filename = re.sub( "%s", id, options.output_pattern )

            dirname = os.path.dirname(output_filename)

            if dirname and not os.path.exists( dirname ):
                os.makedirs( dirname )

            if not os.path.exists( output_filename ):
                outfile = open(output_filename, "w" )
                outfile.write( TreeTools.Tree2Newick( tree ) + "\n" )
                noutput += 1
            else:
                if options.loglevel >= 1:
                    options.stdlog.write("# skipping because output for tree %s already exists: %s\n" % (id, output_filename))                        
                nskipped += 1
                continue

        elif options.method == "filter":
            options.stdout.write( ">%s\n%s\n" % (tree.name, TreeTools.Tree2Newick( tree )) )
            noutput += 1
            
        if outfile_map:
            for t in TreeTools.GetTaxa( tree ):
                outfile_map.write( "%s\t%s\n" % (t, id) )

    if outfile_map:
        outfile_map.close()

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i, with_discarded=%i, discarded_taxa=%i, discarded_branches=%i.\n" %\
                             (ninput, noutput, nskipped,
                              ndiscarded, ndiscarded_taxa, ndiscarded_branches))
        
    E.Stop()
    
