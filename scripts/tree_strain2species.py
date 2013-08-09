################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree_strain2species.py 2782 2009-09-10 11:40:29Z andreas $
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
trees_train2species.py - convert trees containing strains to species trees
===========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of trees from stdin and collapses all nodes that
contain only children of the same species but different strains.

This script was used for the yeast project.

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
import subprocess

import scipy.stats

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools

def parseIdentifier( id, options ):

    data = id.split(options.separator)
    if len(data) == 4:
        return data
    elif len(data) == 2:
        s,g = data
        t = g
        q = "UK"
        return s,t,g,q

def buildIdentifier( schema, transcript, gene, quality, options ):

    if transcript == None:
        return options.separator.join( (schema, gene ) )
    else:
        return options.separator.join( (schema, transcript, gene, quality ) )

def getMergers( tree, map_strain2species, options):
    """merge strains to species.

    returns the new tree with species merged and
    a dictionary of genes including the genes that have been merged.

    Currently, only binary merges are supported.
    """
    
    n = TreeTools.GetSize( tree ) + 1
    all_strains = map_strain2species.keys()
    all_species = map_strain2species.values()
    
    genes = []
    for x in range(n):
        g = {}
        for s in all_strains: g[s] = set()
        genes.append( g )
        
    ## build list of species pairs that can be joined.
    map_species2strain = IOTools.getInvertedDictionary( map_strain2species )
    pairs = []

    for species, strains in map_species2strain.items():
        for x in range(len(strains)):
            for y in range(0,x):
                pairs.append( (strains[x], strains[y]) )

    ## map of genes to new genes
    ## each entry in the list is a pair of genes of the same species
    ## but different strains to be joined.
    map_genes2new_genes = []

    ## dictionary of merged genes. This is to ensure that no gene
    ## is merged twice
    merged_genes = {}

    def count_genes( node_id ):
        """record number of genes per species for each node

        This is done separately for each strain. The counts are aggregated for each species
        over strains by taking the maximum gene count per strain. This ignores any finer
        tree structure below a species node.
        """
        node = tree.node(node_id)
        
        if node.succ:
            this_node_set = genes[node_id]
            ## process non-leaf node
            for s in node.succ:

                ## propagate: terminated nodes force upper nodes to terminate (assigned to None).
                if not genes[s]:
                    this_node_set = None
                    break

                ## check if node merges genes that are not part of the positive set
                for strain in all_strains:
                    if strain in map_strain2species:
                        # merge genes from all children
                        this_node_set[strain] = this_node_set[strain].union(genes[s][strain])
                        
                        if len(this_node_set[strain]) > 1:
                            # more than two genes for a single species, so no join
                            this_node_set = None
                            break
                        
                    elif strain not in map_strain2species and \
                           this_node_set[strain] > 0:
                        this_node_set = None
                        break

            if this_node_set == None:
                genes[node_id] = None
                return

            for strain_x, strain_y in pairs:
                if len(this_node_set[strain_x]) == 1 and len(this_node_set[strain_y]) == 1:
                    species = map_strain2species[strain_x]
                    gene_x, gene_y = tuple(this_node_set[strain_x])[0], tuple(this_node_set[strain_y])[0]

                    ## check if these to genes have already been merged or are
                    ## merged with other partners already
                    ## The merged genes are assigned the same node_id, if they have
                    ## been already merged.
                    key1 = strain_x + gene_x
                    key2 = strain_y + gene_y
                    if key1 > key2: key1, key2 = key2, key1
                    
                    merge = False
                    if key1 in merged_genes and key2 in merged_genes:
                        if merged_genes[key1] == merged_genes[key2]:
                            merge = True
                    elif key1 not in merged_genes and key2 not in merged_genes:
                        merge = True
                        merged_genes[key1] = node_id
                        merged_genes[key2] = node_id
                        
                    if merge:
                        map_genes2new_genes.append( (node_id, species, strain_x, gene_x, strain_y, gene_y ) )
                        
                    ## once two genes have been joined, they can not be remapped further
                    genes[node_id] = None
                    return
        else:
            ## process leaf
            strain, t, g, q = parseIdentifier( node.data.taxon, options )

            if strain in map_strain2species:
                genes[node_id][strain].add( g )
            else:
                ## do not process nodes that do not need to be mapped
                genes[node_id] = None
    
    tree.dfs( tree.root, post_function = count_genes )

    return map_genes2new_genes

def applyMergers( tree, mergers, counters, map_strain2species, options ):
    """apply mergers to a tree."""
    
    new_genes = {}
    for node_id, species, strain_x, gene_x, strain_y, gene_y in mergers:
        if species not in counters:
            counters[species] = 0
        else:
            counters[species] += 1
            
        new_name = buildIdentifier( species, None, options.pattern_gene % counters[species], None, options )
        tree.truncate( node_id, new_name )
        new_genes[new_name] = [ (strain_x, gene_x), (strain_y, gene_y) ]

    ## rename all remaining taxa
    for n in tree.get_terminals():
        strain, t, g, q = parseIdentifier( tree.node(n).data.taxon, options )

        if strain in map_strain2species:
            species = map_strain2species[strain]
            if options.keep_old_names:
                new_name =  buildIdentifier( species, t, g, q)
            else:
                if species not in counters:
                    counters[species] = 0
                else:
                    counters[species] += 1
                new_name = buildIdentifier( species, None, options.pattern_gene % counters[species], None, options )

            tree.node(n).data.taxon = new_name
            new_genes[new_name] = [(strain, g),]
            
    return new_genes

def processGeneTrees( chunks, lines, map_strain2species, options ):
    """process gene trees."""

    if options.output_filename_genes:
        output_genes = open(options.output_filename_genes, "w" )
    else:
        output_genes = options.stdout

    ## for counting genes
    counters = {}
    ## dictionary of merged genes, used to test if some genes appear more than once
    merged = {}
    
    def processChunk( lines, map_strain2species, options ):

        nexus = TreeTools.Newick2Nexus( lines )
        global ninput, noutput, nskipped, nmerged
        
        for tree in nexus.trees:
            ninput += 1

            if options.loglevel >= 3: tree.display()

            mergers = getMergers( tree, map_strain2species, options )

            if options.loglevel >= 3:
                options.stdlog.write("# found %i pairs of genes that will be merged.\n" % (len(mergers) ) )
            
            if len(mergers) > 0:
                nmerged += 1

            n = applyMergers( tree, mergers, counters, map_strain2species, options )

            if len(tree.get_terminals()) <= 1:
                nskipped += 1
                continue

            for new_name, values in n.items():
                for strain, gene in values:
                    if (strain,gene) in merged:
                        options.stdlog.write( "# warning: strain %s and gene %s already appeared in tree %s" % (merged[(strain,gene)]))
                        nwarnings += 1
                    merged[(strain,gene)] = None
                    output_genes.write( "%s\t%s\n" % (options.separator.join((strain,gene)), new_name ) )

            tree.writeToFile( options.stdout, format = options.output_format )
            noutput += 1
            
    if chunks:
        for c in range(len(chunks)-1):
            a, b = chunks[c], chunks[c+1]
            processChunk( lines[a:b], map_strain2species, options )
    else:
        processChunk( lines, map_strain2species, options )

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i, nmerged=%i\n" % (ninput, noutput, nskipped, nmerged))

def getSpeciesTreeMergers( tree, full_map_strain2species, options ):
    """merge strains to species.

    Simply rename all taxa of strains to the species.
    """
    
    nnodes = TreeTools.GetSize( tree ) + 1

    map_strain2species = {}
    for n in tree.get_terminals():
        node = tree.node(n)
        taxon = node.data.taxon
        if taxon in full_map_strain2species:
            map_strain2species[taxon] = full_map_strain2species[taxon]
            node.data.taxon = map_strain2species[taxon]
            
    if len(map_strain2species) == 0: return []

    all_species = tree.get_taxa()
    mapped_species = set(map_strain2species.values())

    species_at_node = []
    for x in range(nnodes):
        g = {}
        for s in all_species: g[s] = 0
        species_at_node.append( g )
        
    def count_species( node_id ):
        """record species for each node
        """
        node = tree.node(node_id)
        if node.succ:
            ## process non-leaf node
            for s in node.succ:
                for species in all_species:
                    species_at_node[node_id][species] += species_at_node[s][species]
        else:
            ## process leaf
            species = node.data.taxon
            species_at_node[node_id][species] = 1
    
    tree.dfs( tree.root, post_function = count_species )
    
    # now merge all those that contain only a single species
    # proceed top-down
    
    nodes_to_skip = set()
    mergers = []

    def merge_species( node_id ):
        
        if node_id in nodes_to_skip: return

        total = sum( species_at_node[node_id].values() )
        for species in mapped_species:
            if species_at_node[node_id][species] <= 1 or \
                    species_at_node[node_id][species] != total: continue
            
            ## merge species
            children = tree.get_leaves( node_id )
            for child in children:
                nodes_to_skip.add( child)

            mergers.append( (node_id, children) )
                
    tree.dfs( tree.root, pre_function = merge_species )

    return mergers

def applySpeciesTreeMergers( tree, mergers, map_strain2species, options ):
    """apply mergers to a tree."""

    for node_id, children in mergers:

        node = tree.node(node_id)

        branch_lengths = [ tree.node(c).data.branchlength for c in children ]

        # copy taxon name from first child
        node.data.taxon = tree.node(children[0]).data.taxon

        # set new branch length
        if options.merge_mode == "ignore":
            pass
        elif options.merge_mode == "add-mean":
            node.data.branchlength += scipy.mean( branch_lengths )
        elif options.merge_mode == "add-max":
            node.data.branchlength += max( branch_lengths )
        elif options.merge_mode == "add-min":
            node.data.branchlength += min( branch_lengths )
            
        # remove all children
        for child in children:
            tree.truncate( child, keep_node = node_id )

def processSpeciesTrees( chunks, lines, map_strain2species, options ):
    """process gene trees."""

    def processChunk( lines, map_strain2species, options ):

        nexus = TreeTools.Newick2Nexus( lines )
        global ninput, noutput, nskipped, nmerged
        
        for tree in nexus.trees:
            ninput += 1

            if options.loglevel >= 3: tree.display()

            mergers = getSpeciesTreeMergers( tree, map_strain2species, options )

            if options.loglevel >= 3:
                options.stdlog.write("# found %i nodes in the tree that will be merged.\n" % (len(mergers) ) )
            
            if len(mergers) > 0:
                nmerged += 1

            n = applySpeciesTreeMergers( tree, mergers, map_strain2species, options )

            if len(tree.get_terminals()) <= 1:
                nskipped += 1
                continue

            tree.writeToFile( options.stdout, format = options.output_format )
            noutput += 1
            
    if chunks:
        for c in range(len(chunks)-1):
            a, b = chunks[c], chunks[c+1]
            processChunk( lines[a:b], map_strain2species, options )
    else:
        processChunk( lines, map_strain2species, options )

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i, nmerged=%i\n" % (ninput, noutput, nskipped, nmerged))

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree_strain2species.py 2782 2009-09-10 11:40:29Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option( "--filename-synonyms", dest="filename_synonyms", type="string" ,
                       help="filename with synonyms. Use this to aggregate several strains for a species.")

    parser.add_option( "--filename-genes", dest="output_filename_genes", type="string" ,
                       help="output filename with new gene names.")

    parser.add_option( "--species-tree", dest="species_tree", action="store_true",
                       help="input tree are species trees. If not given, the trees are assumed to be gene trees." )
                       
    parser.add_option( "--merge-mode", dest="merge_mode", type="choice",
                       choices=("ignore", "add-mean", "add-max", "add-min" ),
                       help="how to deal with branch lengths of merged nodes." )
    
    parser.set_defaults(
        filename_synonyms = "map_strain2species",
        pattern_gene = "J%06i",
        output_format = "nh", 
        separator = "|",
        output_filename_genes = None,
        keep_old_names = False,
        species_tree = False,
        merge_mode = "ignore",
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    ########################################################################
    ########################################################################
    ########################################################################
    ## read synonyms
    if options.filename_synonyms:
        infile = open(options.filename_synonyms, "r" )
        map_strain2species = IOTools.ReadMap( infile )
        infile.close()
    else:
        map_strain2species = {}

    lines = map(lambda x: x[:-1], sys.stdin.readlines())

    ninput, noutput, nskipped, nmerged = 0, 0, 0, 0

    ## iterate over chunks
    chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)))
    if len(chunks) == 0: chunks=[0]
    chunks.append(len(lines))

    if options.species_tree:
        processSpeciesTrees( chunks, lines, map_strain2species, options )
    else:
        processGeneTrees( chunks, lines, map_strain2species, options )

    E.Stop()
