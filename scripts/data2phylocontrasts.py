####
####
##
## Project PairsDBTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: data2phylocontrasts.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####
"""
data2phylocontrasts.py - compute phylogenetic independent contrasts
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a tab-separated table and computes phylogenetic independent 
contrasts according to Felsenstein (1985).

In addition to the table, this script requires a tree with branch lengths.
For each taxon in the tree, the tab-separated table should contain a
column.

This script uses :file:`CONTRAST` of the :file:`PHYLIP` package as a back-end.

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

import sys
import re
import string
import os
import optparse
import math

import CGAT.Experiment as E
import scipy
import scipy.stats
import CGAT.WrapperPhylip as WrapperPhylip
import CGAT.TreeTools as TreeTools

def calculateCorrelationCoefficient( a, b):
    """calculate correlation coefficient for regression through the origin."""

    s1 = 0
    s2 = 0
    s3 = 0
    for i in range(len(a)):
        s1 += a[i] * b[i]
        s2 += a[i] * a[i]
        s3 += b[i] * b[i]
    return s1 / math.sqrt( abs(s2) * abs(s3) )

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: data2phylocontrasts.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )
    parser.add_option("-t", "--filename-tree", dest="filename_tree", type="string",
                      help="filename with tree(s)." )
    parser.add_option("--skip-header", dest="add_header", action="store_false",
                      help="do not add header to flat format."  )
    parser.add_option("--write-header", dest="write_header", action="store_true",
                      help="write header and exit."  )
    parser.add_option("--debug", dest="debug", action="store_true",
                      help="debug mode")
    parser.add_option("--display-tree", dest="display_tree", action="store_true",
                      help="display the tree")
    
    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("contrasts", "spearman", "pearson", "compute" ),
                      help="methods to perform on contrasts." )
    
    parser.set_defaults(
        columns = "all",
        filename_tree = None,
        add_header = True,
        write_header = False,
        debug = False,
        methods = [],
        value_format = "%6.4f",
        pvalue_format = "%e",
        display_tree = False,
        )

    (options, args) = E.Start( parser, quiet = True )

    if options.columns not in ( "all", "all-but-first"):
        options.columns = map(lambda x: int(x) -1 , options.columns.split(","))

    phylip = WrapperPhylip.Phylip()
    
    if options.debug:
        phylip.setLogLevel( options.loglevel )

    phylip.setProgram( "contrast" )

    ##########################################################
    ##########################################################
    ##########################################################        
    # retrieve data and give to phylip
    data = []    
    headers = []
    first = True
    for line in sys.stdin:
        if line[0] == "#": continue
        d = line[:-1].strip().split("\t") 
        if first:
            first = False
            headers = d[1:]
            continue
        data.append(d)
        
    phylip.setData( data )
    ncolumns = len(headers)
    nrows = len(data)
    
    ##########################################################
    ##########################################################
    ##########################################################        
    # read trees
    nexus = None
    if options.filename_tree:
        nexus = TreeTools.Newick2Nexus( open(options.filename_tree, "r") )

    if not nexus:
        raise ValueError("please provide trees with branchlenghts")
    
    ##########################################################
    ##########################################################
    ##########################################################
    # set up phylip
    phylip_options = []
    ## print out contrasts
    phylip_options.append("C")        
    phylip_options.append("Y")        
    phylip.setOptions( phylip_options )
    
    ##########################################################
    ##########################################################
    ##########################################################
    ## main loop
    ##########################################################    
    for tree in nexus.trees:

        if options.display_tree:
            tree.display()

        ## compute this before giving the tree to the phylip module,
        ## as it remaps taxon names.
        map_node2data = {}
        for x in range(nrows):
            taxon = data[x][0]
            map_node2data[tree.search_taxon( taxon )] = x
        
        phylip.setTree( tree )
        
        result = phylip.run()

        for method in options.methods:

            if method in ("pearson", "spearman"):

                options.stdout.write( "header1\theader2\tr\tp\tcode\n" )
                
                n = len(result.mContrasts)
                columns = []
                for c in range( ncolumns ):
                    columns.append( map( lambda x: x[c], result.mContrasts ) )

                for x in range( 0, ncolumns - 1):
                    for y in range( x+1, ncolumns) :
                        
                        # phylip value
                        phy_r = result.mCorrelations[x][y]

                        import rpy
                        from rpy import r as R
                        
                        ## Various ways to calculate r. It is not possible to use
                        ## cor.test or lsfit directly, as you have to perform a
                        ## regression through the origin.
                        
                        ## uncomment to check pearson r against phylip's value
                        ## r = calculateCorrelationCoefficient( columns[x], columns[y] )

                        ## for significance, use linear regression models in R
                        rpy.set_default_mode(rpy.NO_CONVERSION)
                        linear_model = R.lm(R("y ~ x - 1"), data = R.data_frame(x=columns[x], y=columns[y]))
                        rpy.set_default_mode(rpy.BASIC_CONVERSION)

                        ss = R.summary(linear_model)

                        ## extract the p-value
                        p = ss['coefficients'][-1][-1]

                        if p < 0.001:
                            code = "***"
                        elif p < 0.01:
                            code = "**"
                        elif p < 0.05:
                            code = "*"
                        else:
                            code = ""

                        options.stdout.write( "\t".join( (headers[x], headers[y],
                                                          options.value_format % phy_r,
                                                          options.pvalue_format % p,
                                                          code)) + "\n" )

            elif method == "contrasts":
                
                options.stdout.write( "\t".join(headers) + "\n" )
                for d in result.mContrasts:
                    options.stdout.write( "\t".join(map( lambda x: options.value_format % x, d )) + "\n ")

            elif method == "compute":

                ## make room for all internal nodes and one dummy node
                ## for unrooted trees.
                max_index = TreeTools.GetMaxIndex( tree ) + 2
                variances = [None] * max_index
                values = [ [None] * nrows for x in range(max_index) ]
                contrasts = []
                for x in range(max_index):
                    contrasts.append( [None] * ncolumns )
                branchlengths = [None] * max_index
                
                def update_data( node_id, bl, c1, c2, ):

                    b1, b2 = branchlengths[c1], branchlengths[c2]
                    rb1 = 1.0/b1
                    rb2 = 1.0/b2
                    # compute variance
                    variance = math.sqrt(b1+b2)

                    # extend branch length of this node to create correct variance for parent
                    branchlengths[node_id] = bl + (b1 * b2) / (b1 + b2)
                    variances[node_id] = variance

                    for c in range(ncolumns):
                        v1, v2 = values[c1][c], values[c2][c]
                        # save ancestral value as weighted mean
                        values[node_id][c] = (( rb1*v1 + rb2*v2) ) / (rb1+rb2)
                        # compute normalized contrast
                        contrasts[node_id][c] = (v1-v2) / variance
                
                def update_contrasts( node_id ):
                    """update contrasts for a node."""
                    node = tree.node(node_id)
                    if node.succ:
                        if len(node.succ) == 2:
                            c1, c2 = node.succ
                            update_data( node_id, node.data.branchlength, c1, c2 )
                        else:
                            assert( node_id == tree.root )
                            assert( len(node.succ) == 3 )
                            update_data( node_id, node.data.branchlength, node.succ[0], node.succ[1] )
                            update_data( max_index - 1, node.data.branchlength, node_id, node.succ[2] )
                    else:
                        for c in range(ncolumns):
                            values[node_id][c] = float(data[map_node2data[node_id]][c+1])
                            
                        branchlengths[node_id] = node.data.branchlength
                    
                tree.dfs( tree.root, post_function = update_contrasts )

                options.stdout.write( "node_id\tvariance\t%s\n" % "\t".join(headers) )
                for node_id in range(max_index):
                    if variances[node_id] == None: continue
                    options.stdout.write( "%s\t%s\t%s\n" % ( node_id,
                                                             options.value_format % variances[node_id],                                                             
                                                             "\t".join( map( lambda x: options.value_format % x, contrasts[node_id])),
                                                             ))
                                                                                 
                                                                                      
                
    E.Stop()
    
    










