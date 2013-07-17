################################################################################
#   Gene prediction pipeline 
#
#   $Id: trees2tree.py 2782 2009-09-10 11:40:29Z andreas $
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
trees2tree.py - aggregate trees into a single tree
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of trees from stdin and aggregates
them into fewer or a single tree.

Possible aggregation functions are:
   * min, max, sum, mean: min, max, sum, mean of branch-lengths
   * nr: write a set of non-redundant trees.

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
import math

import scipy
import scipy.stats
import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools
import CGAT.WrapperPhylip as WrapperPhylip

def getBestTree( trees, method = "select-largest"):
    """select best tree out of a set of trees."""

    if method == "select-largest":
        sizes = zip( map( lambda x: len(x.get_taxa()), trees ), range(len(trees)))
        sizes.sort()

        best_tree = sizes[-1][1]
        if options.loglevel >= 3:
            for x in range(len(trees)):
                if x == best_tree: continue
                options.stdlog.write("# skipped tree: %s: %s\n" % (trees[x].name, TreeTools.Tree2Newick( trees[x])))
        
        return trees[best_tree]
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: trees2tree.py 2782 2009-09-10 11:40:29Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("counts", "min", "max", "sum", "mean", "median", "stddev", "non-redundant", "consensus",
                               "select-largest"), 
                      help="aggregation function."  )

    parser.add_option("-r", "--regex-id", dest="regex_id", type="string",
                      help="regex pattern to extract identifier from tree name for the selection functions." )

    parser.add_option("-w", "--write-values", dest="write_values", type="string",
                      help="if processing multiple trees, write values to file." )

    parser.add_option("-e", "--error-branchlength", dest="error_branchlength", type="float",
                      help="set branch length without counts to this value." )

    parser.set_defaults(
        method = "mean",
        regex_id = None,
        filtered_branch_lengths = (-999.0,999.0),
        write_values = None,
        error_branchlength = None,
        separator=":",
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if options.loglevel >= 2:
        options.stdlog.write( "# reading trees from stdin.\n" )
        options.stdlog.flush()

    nexus = TreeTools.Newick2Nexus( sys.stdin )
    if options.loglevel >= 1:
        options.stdlog.write( "# read %i trees from stdin.\n" % len(nexus.trees))

    nskipped = 0
    ninput = len(nexus.trees)
    noutput = 0
    nerrors = 0

    if options.method == "non-redundant":
        ## compute non-redudant trees
        template_trees = []
        template_counts = []
        ntree = 0
        for tree in nexus.trees:

            for x in range(0, len(template_trees)):
                is_compatible, reason = TreeTools.IsCompatible( tree, template_trees[x] )
                if is_compatible:
                    template_counts[x] += 1
                    break
            else:
                template_counts.append( 1 )
                template_trees.append( tree )

            if options.loglevel >= 2:
                options.stdlog.write("# tree=%i, ntemplates=%i\n" % (ntree, len(template_trees)))
                
            ntree += 1
            
        for x in range(0, len(template_trees)):
            if options.loglevel >= 1:
                options.stdlog.write("# tree: %i, counts: %i, percent=%5.2f\n" %\
                                     (x, template_counts[x], template_counts[x] * 100.0 / ntotal))
            options.stdout.write(TreeTools.Tree2Newick( template_trees[x] ) + "\n")

    elif options.method in ("select-largest",):
        ## select one of the trees with the same name.
        clusters = {}
        for x in range(0, len(nexus.trees)):
            n = nexus.trees[x].name

            if options.regex_id:
                n = re.search(options.regex_id, n).groups()[0]
                
            if n not in clusters: clusters[n] = []
            clusters[n].append(x)

        new_trees = []
        
        for name, cluster in clusters.items():
            new_trees.append( getBestTree( [ nexus.trees[x] for x in cluster], options.method ) )
            
        for x in range(0, len(new_trees)):
            options.stdout.write(">%s\n" % new_trees[x].name )
            options.stdout.write(TreeTools.Tree2Newick( new_trees[x],  ) + "\n")
            noutput += 1

        nskipped = ntotal - noutput
            
    elif options.method == "consensus":
        
        phylip = WrapperPhylip.Phylip()
        phylip.setLogLevel( options.loglevel - 2 )
        phylip.setProgram( "consense" )
        phylip_options = []
        phylip_options.append("Y")
        
        phylip.setOptions( phylip_options )
        phylip.setTrees( nexus.trees )

        result = phylip.run()

        options.stdout.write("# consensus tree built from %i trees\n" % (phylip.mNInputTrees) )
        options.stdout.write(TreeTools.Tree2Newick( result.mNexus.trees[0] ) + "\n")
        noutput = 1

    else:
        if options.method in ("min", "max", "sum", "mean", "counts"):

            xtree = nexus.trees[0]
            for n in xtree.chain.keys():
                if xtree.node(n).data.branchlength in options.filtered_branch_lengths:
                    xtree.node(n).data.branchlength = 0
                ntotals = [1] * len(xtree.chain.keys())
                    
            if options.method == "min":
                f = min
            elif options.method == "max":
                f = max
            elif options.method == "sum":
                f = lambda x, y: x + y
            elif options.method == "mean":
                f = lambda x, y: x + y
            elif options.method == "counts":
                f = lambda x, y: x + 1
                for n in xtree.chain.keys():
                    if xtree.node(n).data.branchlength not in options.filtered_branch_lengths:
                        xtree.node(n).data.branchlength = 1
                    else:
                        xtree.node(n).data.branchlength = 0                        
            else:
                raise "unknown option %s" % options.method
            
            for tree in nexus.trees[1:]:
                
                for n in tree.chain.keys():
                    if tree.node(n).data.branchlength not in options.filtered_branch_lengths:
                        xtree.node(n).data.branchlength = f(xtree.node(n).data.branchlength, tree.node(n).data.branchlength)
                        ntotals[n] += 1

            if options.method == "mean":
                for n in xtree.chain.keys():
                    if ntotals[n] > 0:
                        xtree.node(n).data.branchlength = float(xtree.node(n).data.branchlength) / ntotals[n]
                    else:
                        if options.error_branchlength != None:
                            xtree.node(n).data.branchlength = options.error_branchlength
                            if options.loglevel >= 1:
                                options.stdlog.write("# no counts for node %i - set to %f\n" % (n, options.error_branchlength) )
                                nerrors += 1
                        else:
                            raise "no counts for node %i" % n
                        
        else:
            ## collect all values for trees
            values = [ [] for x in range(TreeTools.GetSize( nexus.trees[0] )) ]
            
            for tree in nexus.trees:
                for n, node in tree.chain.items():
                    if node.data.branchlength not in options.filtered_branch_lengths:
                        values[n].append( node.data.branchlength )

            tree = nexus.trees[0]
            for n, node in tree.chain.items():
                if len(values[n]) > 0:
                    if options.method == "stddev":
                        node.data.branchlength = scipy.std( values[n] )
                    elif options.method == "median":
                        node.data.branchlength = scipy.median( values[n] )
                else:
                    if options.error_branchlength != None:
                        node.data.branchlength = options.error_branchlength
                        if options.loglevel >= 1:
                            options.stdlog.write("# no counts for node %i - set to %f\n" % (n, options.error_branchlength) )
                            nerrors += 1
                    else:
                        raise "no counts for node %i" % n
                    
            if options.write_values:
                outfile = open(options.write_values, "w" )
                for n, node in tree.chain.items():
                    values[n].sort()
                    id = options.separator.join( sorted(TreeTools.GetLeaves( tree, n )) )
                    outfile.write( "%s\t%s\n" % (id, ";".join( map(str, values[n] ) ) ))
                outfile.close()
                
        del nexus.trees[1:]
        options.stdout.write(TreeTools.Nexus2Newick( nexus ) + "\n")
        noutput = 1
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ntotal=%i, nskipped=%i, noutput=%i, nerrors=%i\n" % (ninput, nskipped, noutput, nerrors))
        
    E.Stop()
