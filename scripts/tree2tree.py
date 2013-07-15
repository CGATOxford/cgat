################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree2tree.py 2782 2009-09-10 11:40:29Z andreas $
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
tree2tree.py - manipulate trees
===============================

:Author: Andreas Heger
:Release: $Id$

:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of trees from stdin and outputs the again on stdout
after manipulating them. Manipulations include

   * renaming taxa

   * normalizing branch lengths

The complete list of methods is:

normalize
+++++++++

divide-by-tree
++++++++++++++

rename
++++++

set-uniform-branch-length
+++++++++++++++++++++++++

extract-with-pattern
++++++++++++++++++++

build-map
++++++++++++++++++++

remove-pattern
++++++++++++++++++++

unroot
++++++++++++++++++++

midpoint-root
++++++++++++++++++++

balanced-root
++++++++++++++++++++

add-node-names
++++++++++++++++++++

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

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools
import CGAT.WrapperPhylip as WrapperPhylip

def Process( lines, other_trees, options, map_old2new, ntree):

    nexus = TreeTools.Newick2Nexus( map(lambda x: x[:-1], lines) )

    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees.\n" % len(nexus.trees))

    nskipped = 0
    ntotal = len(nexus.trees)
    extract_pattern = None
    species2remove = None
    write_map = False

    phylip_executable = None
    phylip_options = None

    index = 0

    ## default: do not output internal node names
    write_all_taxa = False
    
    for tree in nexus.trees:

        if options.outgroup:
            tree.root_with_outgroup( options.outgroup )

        for method in options.methods:

            if options.loglevel >=3:
                options.stdlog.write("# applying method %s to tree %i.\n" % (method, index))
                
            if method == "midpoint-root":
                tree.root_midpoint()

            elif method == "balanced-root":
                tree.root_balanced()

            elif method == "unroot":
                TreeTools.Unroot(tree)
                
            elif method=="phylip":
                if not phylip_executable:
                    phylip_executable=options.parameters[0]
                    del options.parameters[0]
                    phylip_options = re.split("@", options.parameters[0] )
                    del options.parameters[0]
                    
                    phylip = WrapperPhylip.Phylip()
                    phylip.setProgram( phylip_executable )
                    phylip.setOptions( phylip_options )

                phylip.setTree( tree )

                result = phylip.run()

                nexus.trees[index] = result.mNexus.trees[0]
                
            elif method == "normalize":
                if options.value == 0:
                    v = 0
                    for n in tree.chain.keys():
                        v = max( v, tree.node(n).data.branchlength )
                else:
                    v = options.value

                for n in tree.chain.keys():
                    tree.node(n).data.branchlength /= float( options.value )
                    
            elif method == "divide-by-tree":

                if len(other_trees) > 1:
                    other_tree = other_trees[ntree]
                else:
                    other_tree = other_trees[0]
                    
                ## the trees have to be exactly the same!!
                if options.loglevel >= 2:
                    print tree.display()
                    print other_tree.display()

                if not tree.is_identical( other_tree ):
                    nskipped += 1
                    continue

                ## even if the trees are the same (in topology), the node numbering might not be
                ## the same. Thus build a map of node ids.
                map_a2b = TreeTools.GetNodeMap( tree, other_tree )

                for n in tree.chain.keys():
                    try:
                        tree.node(n).data.branchlength /= float( other_tree.node(map_a2b[n]).data.branchlength )
                    except ZeroDivisionError:
                        options.stdlog.write("# Warning: branch for nodes %i and %i in tree-pair %i: divide by zero\n" %\
                                             ( n, map_a2b[n], ntree ))
                        continue

            elif method == "rename":
                if not map_old2new:
                    
                    map_old2new = IOTools.ReadMap( open(options.parameters[0], "r"), columns=(0,1) )

                    if options.invert_map:
                        map_old2new = IOTools.getInvertedDictionary( map_old2new, make_unique = True)
                        
                    del options.parameters[0]

                unknown = []
                for n, node in tree.chain.items():
                    if node.data.taxon:
                        try:
                            node.data.taxon = map_old2new[node.data.taxon]
                        except KeyError:
                            unknown.append( node.data.taxon )

                for taxon in unknown:
                    tree.prune( taxon )

            ## reformat terminals
            elif method == "extract-with-pattern":
                
                if not extract_pattern:
                    extract_pattern = re.compile(options.parameters[0])
                    del options.parameters[0]
                
                for n in tree.get_terminals():
                    node = tree.node( n )
                    node.data.taxon = extract_pattern.search( node.data.taxon ).groups()[0]

            elif method == "set-uniform-branchlength":
                for n in tree.chain.keys():
                    tree.node(n).data.branchlength = options.value

            elif method == "build-map":
                # build a map of identifiers
                options.write_map = True
                for n in tree.get_terminals():
                    node = tree.node( n )
                    if node.data.taxon not in map_old2new:
                        new = options.template_identifier % (len(map_old2new) + 1)
                        map_old2new[node.data.taxon] = new                        
                    node.data.taxon = map_old2new[node.data.taxon]

            elif method == "remove-pattern":
                if species2remove == None:
                    species2remove = re.compile(options.parameters[0])
                    del options.parameters
                taxa = []
                for n in tree.get_terminals():
                    t = tree.node(n).data.taxon
                    skip = False
                    if species2remove.search( t ):
                        continue
                    if not skip:
                        taxa.append(t)
                TreeTools.PruneTree( tree, taxa )

            elif method == "add-node-names":

                inode = 0
                write_all_taxa = True
                for n, node in tree.chain.items():
                    if not node.data.taxon:
                        node.data.taxon = "inode%i" % inode
                        inode += 1

            elif method == "newick2nhx":
                ## convert names to species names
                for n in tree.get_terminals():
                    t = tree.node(n).data.taxon
                    d = t.split("|")
                    if len(d) >= 2:
                        tree.node(n).data.species = d[0]
                        
        index += 1
        ntree += 1

    if options.output_format == "nh":
        options.stdout.write(TreeTools.Nexus2Newick( nexus, write_all_taxa = True,
                                                     with_branchlengths = options.with_branchlengths) + "\n" )
    else:
        for tree in nexus.trees:
            tree.writeToFile( options.stdout, format = options.output_format )

    return ntotal, nskipped, ntree

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree2tree.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-d", "--value", dest="value", type="float",
                      help="normalizing value."  )
    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="""methods to apply [normalize|divide-by-tree|divide-by-tree|rename|set-uniform-branch-length|extract-with-pattern|build-map|remove-pattern|unroot|midpoint-root|balanced-root|add-node-names"""  )
    parser.add_option("-2", "--filename-tree2", dest="filename_tree2", type="string",
                      help="filename with second tree."  )
    parser.add_option("-o", "--outgroup", dest="outgroup", type="string",
                      help="reroot with outgroup before processing.")
    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameters for methods.")
    parser.add_option("-e", "--template-identifier", dest="template_identifier", type="string",
                      help="""template identifier [%default]. A %i is replaced by the position
                      of the sequence in the file."""  )
    parser.add_option("-i", "--invert-map", dest="invert_map", action="store_true",
                      help="""invert map.""")
    parser.add_option("-f", "--filter", dest="filter", type="choice",
                      choices=("max-branch-length",),
                      help="filter trees")
    parser.add_option("--output-format", dest="output_format", type="choice",
                      choices=( "nh", "nhx" ),
                      help=("output format for trees."))
    parser.add_option("-b", "--no-branch-lengths", dest="with_branchlengths", action="store_false",
                      help="""do not write branchlengths. Per default, 0 branch lengths are added.""")

    parser.set_defaults(
        value = 0,
        methods = "",
        filename_tree2 = None,
        outgroup = None,
        parameters = "",
        template_identifier="ID%06i",
        write_map = False,
        invert_map = False,
        filter = None,
        output_format = "nh",
        with_branchlengths = True,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    options.methods = options.methods.split(",")
    options.parameters = options.parameters.split(",")    

    other_trees = []
    ## read other trees
    if options.filename_tree2:
        other_nexus = TreeTools.Newick2Nexus( open(options.filename_tree2, "r") )
        if len(other_nexus.trees) > 0:
            other_trees = other_nexus.trees
        else:
            other_tree = other_nexus.trees[0]
            other_trees = [other_tree]

    lines = sys.stdin.readlines()

    ntotal, nskipped, ntree = 0, 0, 0

    if options.filter:

        nexus = TreeTools.Newick2Nexus( lines )
        
        new_trees = []

        value = float(options.parameters[0])
        del options.parameters[0]

        ## decision functions: return true, if tree
        ## is to be skipped
        if options.filter == "max-branch-length":
            f = lambda x: x >= value

        for tree in nexus.trees:
            ntotal += 1

            for id, node in tree.chain.items():
                if f(node.data.branchlength):
                    nskipped += 1
                    break
            else:
                new_trees.append(tree)
                ntree += 1
                
        nexus.trees = new_trees

        options.stdout.write(TreeTools.Nexus2Newick( nexus, with_names = True ) + "\n" )
        
    else:

        ## iterate over chunks
        chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)))

        map_old2new = {}

        if chunks:
            for c in range(len(chunks)-1):
                a, b = chunks[c], chunks[c+1]
                options.stdout.write( lines[a] )
                a += 1
                Process( lines[a:b], other_trees, options, map_old2new,ntree )

            options.stdout.write( lines[chunks[-1]] )
            t, s, ntree = Process( lines[chunks[-1]+1:], other_trees, options, map_old2new, ntree )
            ntotal += t
            nskipped += s
        else:
            ntotal, nskipped, ntree = Process( lines, other_trees, options, map_old2new,ntree )

        if options.write_map:
            p = options.parameters[0]
            if p:
                outfile = open(p, "w")
            else:
                outfile = options.stdout

            outfile.write("old\tnew\n")            
            for old_id, new_id in map_old2new.items():
                outfile.write("%s\t%s\n" % (old_id, new_id) )
            if p:
                outfile.close()

    if options.loglevel >= 1:
        options.stdlog.write( "# ntotal=%i, nskipped=%i\n" % (ntotal, nskipped))
        
    E.Stop()
