################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: trees2sets.py 2782 2009-09-10 11:40:29Z andreas $
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
trees2sets.py - compute ortholog sets from trees
================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

extract sets from a list of trees. The trees need to be rooted already.

Given a reference tree and a collection of trees, this script outputs
(depending on the option --method)

   * strict: 1:1 ortholog clusters for a selected set of species
     combinations. All species need to be present for a species combination
   
   * degenerate: 1:m ortholog cluster for a selected set of species
     combinations. All species need to be present for a species combination.

   * lineage: lineage specific duplications.

   * any: 1:1 and 1:m ortholog clusters

The default set of species combinations consists of all monophyletic clades.
Different sets can be specified with the option --enumeration:

   * exhaustive: all possible species combinations
   
   * monophyletic: all monophyletic species combinations (provide a species tree)
   
   * full: one set of all species

   * explicit: species combination given on command line

   * lineage: lineage specific duplications

   * pairwise: all pairwise species combinations

The script can map strains to species.

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
import optparse

import CGAT.Orthologs as Orthologs

USAGE="""python trees2sets.py < stdin 


"""

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.SetTools as SetTools
import CGAT.IOTools as IOTools

def buildPattern( nspecies, s, members = None):
    """build a present/absent pattern for species
    in set s.
    """
    pattern = []
    for x in range(nspecies):
        if x in s:
            if members:
                pattern.append( str(len(members[x])) )
            else:
                pattern.append( "1" )                
        else:
            pattern.append( "0" )
    return pattern

def parseIdentifier( id, options ):

    data = id.split(options.separator)
    if len(data) == 4:
        return data
    elif len(data) == 2:
        s,g = data
        t = g
        q = "UK"
        return s,t,g,q

##--------------------------------------------------------------------------------
def getOrthologNodes( tree, positive_set, options, selector = "strict", outgroups = None ):
    """get all ortholog nodes in tree for species in positive_set.

    Depending on the selector function, different sets are returned:
    
    If selector is "strict", only strict orthologs are returned. These
    contain exactly one gene per species for all species in the positive_set.
    
    If selector is "degenerate", only degenerate orthologs are returned.
    These contain at least gene per species for species in the positive_set.
    
    Collect genes in tree for each species.

    Returns the node_id for which a set fulfills the criteria
    and the set for which it fulfills it.

    Avoid double counting: if you are interested in species A and B,
    any branches involving others species should be ignored. Make sure to only
    count once and not every time a discarded branch is removed.

    Thus, as soon as A and B merge, any node further up the tree have to be
    ignored.

    total_genes_function:   if true, node is recorded
    total_species_function: if true, iteration stops

    """

    nspecies = len(options.org2column)
    
    if selector == "strict":
        ## strict orthologs: at most one gene per species
        exit_function = lambda num_genes_for_species: num_genes_for_species > 1
        keep_function = lambda num_genes_for_species: num_genes_for_species == 1
        total_genes_function = lambda num_genes_at_node,num_species_in_pattern: num_genes_at_node == num_species_in_pattern
        total_species_function = lambda num_species_at_node,num_species_in_pattern: num_species_at_node == num_species_in_pattern
        check_outgroup_function = lambda x: False        
        negative_set = set()
    elif selector == "degenerate":
        ## degenerate orthologs: any number of genes per species,
        exit_function = lambda num_genes_for_species: False
        keep_function = lambda num_genes_for_species: num_genes_for_species > 0
        total_genes_function = lambda num_genes_at_node,num_species_in_pattern: num_genes_at_node > num_species_in_pattern
        total_species_function = lambda num_species_at_node,num_species_in_pattern: num_species_at_node == num_species_in_pattern        
        check_outgroup_function = lambda x: False
        negative_set = set()
    elif selector == "lineage":
        ## lineage specific duplications: at least 1 gene
        exit_function = lambda num_genes_for_species: False
        keep_function = lambda num_genes_for_species: num_genes_for_species > 1
        total_genes_function = lambda num_genes_at_node,num_species_in_pattern: num_genes_at_node >= num_species_in_pattern
        total_species_function = lambda num_species_at_node,num_species_in_pattern: False
        check_outgroup_function = lambda x: False
        negative_set = set( range(nspecies) ).difference( positive_set )
    elif selector == "any":
        ## any number of orthologs, including
        ## orphans
        exit_function = lambda num_genes_for_species: False
        keep_function = lambda num_genes_for_species: True
        total_genes_function = lambda num_genes_at_node,num_species_in_pattern: True
        total_species_function = lambda num_species_at_node,num_species_in_pattern: num_species_at_node == num_species_in_pattern
        check_outgroup_function = lambda x: False
        negative_set = set()
    elif selector == "outgroup":
        ## group selector
        exit_function = lambda num_genes_for_species: False
        keep_function = lambda num_genes_for_species: num_genes_for_species > 0
        total_genes_function = lambda num_genes_at_node,num_species_in_pattern: False        
        total_species_function = lambda num_species_at_node,num_species_in_pattern: False
        ## check for outgrup: needs to have outgroup and at least one other species
        ## ie.: sum of all genes in outgroups larger than sum of all genes in all species
        if not outgroups:
            raise "usage error: please supply outgroups if 'outgroup'-selector is chosen."
        check_outgroup_function = lambda genes: 0 < sum( [ len(genes[x]) for x in outgroups ] ) < sum( map( lambda x: len(x), genes) )
        negative_set = set()
    else:
        raise "unknown selector %s" % selector

    ## work here: set genes[node_id] to None,
    ## 1. if the gene count for a species of interest is > 1
    ## 2. if the gene count for all species of interest is 1 in
    ## the child node.

    if options.loglevel >= 5:
        options.stdlog.write( "# gene tree\n" )
        tree.display()

    n = TreeTools.GetSize( tree ) + 1
    genes = []
    for x in range(n):
        genes.append( [ set() for x in range(nspecies)] )

    ortholog_nodes = []

    def count_genes( node_id ):
        """record number of genes per species for each node
        """
        node = tree.node(node_id)
        
        if options.loglevel >= 6:
            options.stdlog.write("# node_id=%i\n" % node_id)
            if options.loglevel >= 10:
                options.stdlog.write("# sets=%s\n" % (str(genes)))

        # species in pattern
        num_species_in_pattern = len(positive_set)
        
        if node.succ:
            ## process non-leaf node
            for s in node.succ:

                ## propagate: terminated nodes force upper nodes to terminate (assigned to None).
                if not genes[s]:
                    genes[node_id] = None
                    return
                
                ## total number of genes at node
                num_genes_at_node = 0
                ## total number of species at node
                num_species_at_node = 0

                ## compute new gene set for each species at node
                for x in positive_set:
                    genes[node_id][x] = genes[node_id][x].union(genes[s][x])
                    
                    num_genes_for_species = len(genes[node_id][x])
                    if exit_function(num_genes_for_species):
                        genes[node_id] = None
                        return
                    num_genes_at_node += num_genes_for_species
                    if num_genes_for_species: num_species_at_node += 1
                    
            if options.loglevel >= 6:
                print "node=", node_id, "species_at_node", num_species_at_node, "genes_at_node=", num_genes_at_node, \
                    "num_genes_for_species=", num_genes_for_species, "ngenes=", sum( map( lambda x: len(x), genes[node_id]) )
                options.stdlog.write("# genes at node %i\t%s\n" % (node_id, genes[node_id]))
                if outgroups:
                    print sum( [ len(genes[node_id][x]) for x in outgroups ] )
                    print check_outgroup_function(genes[node_id])

            # check stop criterion
            if total_species_function( num_species_at_node, num_species_in_pattern ):
                ## check if positive requirements are fulfilled
                for x in positive_set:
                    if not keep_function(len(genes[node_id][x])):
                        if options.loglevel >= 6:
                            options.stdlog.write("# keep function false for species %i\n" % x)
                        break
                else:
                    if total_genes_function( num_genes_at_node, num_species_in_pattern):
                        if options.loglevel >= 6:
                            options.stdlog.write("# recording node %i\n" % x)
                        ortholog_nodes.append( (node_id, genes[node_id]) )
                genes[node_id] = None
                return
            elif check_outgroup_function( genes[node_id] ):
                ortholog_nodes.append( (node_id, genes[node_id]) )
                genes[node_id] = None
                return
            elif negative_set:
                if total_genes_function( num_genes_at_node, num_species_in_pattern):
                    if options.loglevel >= 6:
                        options.stdlog.write("# recording node %i\n" % node_id)
                    ortholog_nodes.append( (node_id, genes[node_id]) )

        else:
            ## process leaf
            s, t, g, q = parseIdentifier( node.data.taxon, options )
            c = options.org2column[s]
            if c in positive_set:
                genes[node_id][c].add( g )
            elif c in negative_set:
                genes[node_id] = None

    tree.dfs( tree.root, post_function = count_genes )
    
    return ortholog_nodes

##------------------------------------------------------------------------------------------                
##------------------------------------------------------------------------------------------                
##------------------------------------------------------------------------------------------
def rerootTree( gene_tree, extract_species, options ):
    """reroot tree."""

    otus = TreeTools.GetTaxa( gene_tree )
    
    ## find monophyletic trees of outgroup_species
    try:
        outgroup_taxa = filter( lambda x: extract_species(x) in options.outgroups, otus)
    except AttributeError:
        raise "error while rerooting tree in tree %s with %s" % (gene_tree.name, str(otus))

    if gene_tree.is_monophyletic( outgroup_taxa ):
        r = outgroup_taxa
    else:
        r = [outgroup_taxa[0],]

    if r:
        if options.loglevel >= 1:
            options.stdlog.write("# tree %s: rerooting with %i outgroups:  %s.\n" % (gene_tree.name, len(r), ",".join(r)))
            options.stdlog.flush()            
    else:
        if options.loglevel >= 1:
            options.stdlog.write("# tree %s: no outgroup found, tree will not be rerooted.\n" % gene_tree.name)
            options.stdlog.flush()
            
    gene_tree.root_with_outgroup( r )

    if options.loglevel >= 5:
        gene_tree.display()

##--------------------------------------------------------------------------------
def writeOrthologSets( outfile, nexus,
                       extract_species,
                       extract_gene,
                       options,
                       reference_tree = None,
                       method = "strict",
                       outgroups = None ):
    """output ortholog sets.

    A "strict" ortholog set contains exactly one gene for each species,
    while a "degenerate" ortholog set contains at least one gene for each
    species. 
    """

    ######################################################################
    # build species set to compare
    sets = []
    species = options.column2org
    nspecies = len(species)
    
    if options.enumeration == "monophyletic":
        if reference_tree:
            for members, h1, h2 in TreeTools.GetSubsets( reference_tree ):
                if len(members) > 1:
                    sets.append( members )
        else:
            raise "please specify a species tree for monophyletic enumeration"
        
    elif options.enumeration == "exhaustive":
        for x in range(2, len(species) ):
            sets += list(SetTools.xuniqueCombinations( species, x))
        sets.append(species)
                
    elif options.enumeration == "pairwise":
        
        for x in range(len(species)-1):
            for y in range(x+1, len(species)):
                sets.append( (species[x], species[y]) )
                
    elif options.enumeration == "full":                
        sets.append( species )

    elif options.enumeration == "lineage":
        for s in species: 
            sets.append( (s,) )

    elif options.enumeration == "explicit":
        for x in range(2, len(options.species_set)):
            sets += list(SetTools.xuniqueCombinations( options.species_set, x ))
        sets.append(options.species_set)

    ######################################################################
    # build sets with positional information
    xsets = []
    map_frozenset2set = {}
    for x in range(len(sets)):
        ss = frozenset(map( lambda x: options.org2column[x], sets[x]) )
        xsets.append( ss )
        map_frozenset2set[ss] = x

    ######################################################################
    ## collect outgroups
    if outgroups:
        noutgroups = set()
        for x in outgroups:
            noutgroups.add( options.org2column[x] )
    else:
        noutgroups = None
        
    ######################################################################
    # loop over each tree and set
    # I did not see a way to loop a tree once for all sets without doing
    # complicated counting. The problem is that counting has to be stopped
    # at different tree heights for different sets.
    ninput, noutput, nempty, nskipped = 0, 0, 0, 0

    counts = [0] * len(sets)

    options.stdout.write( "nspecies\tname\tid\tcluster\tpattern\t%s\tnode_id\tmembers\n" % "\t".join(species))

    cluster_id = 0
    nerrors = 0
    
    for tree in nexus.trees:
        
        ninput += 1
        ntotal_tree = 0

        if options.loglevel >= 3:
            options.stdlog.write("# processing tree %s\n" % tree.name )

        if options.reroot:
            rerootTree( tree, extract_species, options )
        
        for c in range(len(xsets)):
            ## numbered species set: 0,1,...
            sn = xsets[c]
            ## literal species set: species1, species2, ...
            sl = sets[c]
            
            ortholog_nodes = getOrthologNodes( tree, sn, options, selector = method,
                                               outgroups = noutgroups )
            ntotal_tree += len(ortholog_nodes)

            n = 0

            pattern =  buildPattern( nspecies, sn )

            ## check for inconsistent partitions (the same gene in different
            ## ortholog clusters) within the current tree
            found_genes = set()
            ortho_sets = set()
            
            ## reverse ortholog_node - work in top-down manner.
            ortholog_nodes.reverse()

            for node_id, members in ortholog_nodes:
                n += 1
                cluster_id += 1
                
                otus = filter( lambda x: extract_species(x) in sl, tree.get_taxa( node_id ) )
                genes = set( map( extract_gene, otus ) )

                if found_genes.intersection( genes ):

                    ## only take largest cluster for lineage specific duplications
                    if method == "lineage":
                        continue

                    if frozenset(genes) in ortho_sets:
                        nskipped += 1
                        if options.loglevel >= 1:
                            options.stdlog.write("# %s: cluster %i: redundant node: %i - skipped because already present: %s\n" %
                                                 (tree.name, n, node_id, str(found_genes.intersection( genes ))))
                    else:
                        nerrors += 1
                        if options.loglevel >= 1:
                            options.stdlog.write("# %s: cluster %i: inconsistent node: %i - the same gene in different clusters: %s\n" %
                                                 (tree.name, n, node_id, str(found_genes.intersection( genes ))))
                    
                found_genes = found_genes.union( genes )
                ortho_sets.add( frozenset(genes) )
                
                xpattern = buildPattern( nspecies, sn, members )
                
                options.stdout.write( "%i\t%s\t%i\t%i\t%s\t%s\t%i\t%s\n" % (len(sl),
                                                                            tree.name,
                                                                            n,
                                                                            cluster_id, 
                                                                            "".join(pattern),
                                                                            "\t".join(xpattern),
                                                                            node_id,
                                                                            ";".join(otus) ))
                
            counts[c] += n
            
        if ntotal_tree == 0:
            nempty += 1
        else:
            noutput += 1

    if options.loglevel >= 1:
        options.stdout.write("# ninput=%i, nempty=%i, noutput=%i, nskipped=%i, nerrors=%i\n" % (ninput, nempty, noutput, nskipped, nerrors) )

    ## write summary information
        
    if options.filename_summary:
        outfile = open(options.filename_summary, "w")
    else:
        outfile = options.stdout
        outfile.write("//\n")
        
    outfile.write( "cluster\tpattern\tcounts\t%s\n" % ("\t".join(species)) )
    
    for c in range(len(xsets)):
        pattern = buildPattern( nspecies, xsets[c])
        outfile.write("%i\t%s\t%i\t%s\n" % ( c,
                                             "".join(pattern),
                                             counts[c],
                                             "\t".join(pattern) ))
                                             
    if outfile != options.stdout:
        outfile.close()


if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: trees2sets.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])
    
    parser.add_option("-t", "--reference-tree", dest="reference_tree", type="string",
                      help="reference tree to read."  )

    parser.add_option("-e", "--enumeration", dest="enumeration", type="choice",
                      choices=("monophyletic", "full", "pairwise", "exhaustive", "explicit", "lineage" ),
                      help="enumeration of ortholog groups."  )

    parser.add_option( "-o", "--organisms", dest="column2org", type="string" ,
                       help="sorted list of organisms.")

    parser.add_option("-p", "--filename-patterns", dest="filename_patterns", type="string",
                      help="filename with patterns to output."  )

    parser.add_option("-u", "--filename-summary", dest="filename_summary", type="string",
                      help="filename with summary to output."  )

    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("strict", "degenerate", "any", "outgroup", "lineage" ),
                      help="sets to extract."  )

    parser.add_option("-s", "--species-set", dest="species_set", type="string",
                      help="comma separated list of species." )

    parser.add_option("-g", "--outgroups", dest="outgroups", type="string",
                      help="comma separated list of outgroup species." )

    parser.add_option( "--species-regex", dest="species_regex", type="string" ,
                       help="regular expression to extract species from identifier.")

    parser.add_option( "--gene-regex", dest="gene_regex", type="string" ,
                       help="regular expression to extract gene from identifier.")

    parser.add_option( "--reroot", dest="reroot", type="choice" ,
                       choices=("outgroup", "midpoint"),
                       help="reroot trees before computing sets.")

    parser.set_defaults(
        reference_tree = None,
        enumeration = "full",
        column2org=None,
        separator = "|",
        species_regex ="^([^|]+)\|",
        gene_regex = "^[^|]+\|[^|]+\|([^|]+)\|",
        filename_summary = None,
        methods = [],
        species_set = None,
        outgroups = None,
        reroot = None,
        )

    (options, args) = E.Start( parser )

    if len(options.methods) == 0:
        options.methods.append("strict")

    if options.species_set:
        options.species_set = options.species_set.split(",")
        options.enumeration = "explicit"

    #######################################################################
    ## warning: outgroup method is useless, as it requires
    ## only a single outgroup per tree and the tree rooted
    ## with the outgroup.
    if "outgroup" in options.methods and not options.outgroups:
        raise "please supply --outgroups if method 'outgroup' is chosen."

    if options.outgroups:
        options.outgroups = options.outgroups.split(",")
        
    ########################################################################
    ########################################################################
    ########################################################################
    if options.reference_tree:
        if options.reference_tree[0] == "(":
            nexus = TreeTools.Newick2Nexus( options.reference_tree )
        else:
            nexus = TreeTools.Newick2Nexus( open(options.reference_tree, "r") )
        reference_tree = nexus.trees[0]
    
        if options.loglevel >= 3:
            options.stdlog.write("# reference tree:\n%s\n" % reference_tree.display())
    else:
        reference_tree = None
        raise ValueError( "please supply a reference tree" )

    ########################################################################
    ########################################################################
    ########################################################################
    ## read all trees
    ########################################################################    
    nexus = TreeTools.Newick2Nexus( sys.stdin )

    ########################################################################
    ########################################################################
    ########################################################################
    ## sort out reference tree
    ########################################################################
    rs = re.compile( options.species_regex )
    rg = re.compile( options.gene_regex )
    extract_species = lambda x: parseIdentifier(x, options)[0]
    extract_gene = lambda x:parseIdentifier(x, options)[2]
    
    ## prune reference tree to species present
    species_set = set()
    for tree in nexus.trees:
        try:
            species_set = species_set.union( set(map( extract_species, tree.get_taxa()) ) )
        except AttributeError:
            raise "parsing error while extracting species from %s" % str(tree.get_taxa())

    TreeTools.PruneTree( reference_tree, species_set )
    
    if options.loglevel >= 1:
        options.stdlog.write("# reference tree after pruning has %i taxa.\n" % len(reference_tree.get_taxa()))

    if options.column2org:
        options.column2org = options.column2org.split(",")
    elif reference_tree:
        options.column2org = []
        for nx in reference_tree.get_terminals():
            options.column2org.append( reference_tree.node(nx).get_data().taxon )
    
    options.org2column = {}
    for x in range(len(options.column2org)):
        options.org2column[options.column2org[x]] = x

    for method in options.methods:

            #######################################################################################
            #######################################################################################
            #######################################################################################
            ## print out a list of ortholog clusters
            #######################################################################################
            writeOrthologSets( options.stdout,
                               nexus,
                               extract_species,
                               extract_gene,
                               options = options,
                               reference_tree = reference_tree,
                               method = method,
                               outgroups = options.outgroups )
    
    E.Stop()
