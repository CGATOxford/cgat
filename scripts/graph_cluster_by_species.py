################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
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
graph_cluster_by_species.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python graph_cluster_by_species.py --help

Type::

   python graph_cluster_by_species.py --help

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
import getopt
import tempfile
import time
import popen2
import optparse

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: graph_cluster_by_species.py 2782 2009-09-10 11:40:29Z andreas $

Cluster transcripts based on a distance measure. The transcripts can be clustered by

   * species: Find best matching clusters including as many species as possible. Used to
        select a subset of 1:1 orthologous genes. Better though is to use a tree.
        
   * genes: build clusters out of best matching transcripts. Used to build clusters of orthologous
        transcripts from a cluster of orthologs. Clustering is done by reciprocal best hit.
        Identical transcripts have a distance of 0 and are grouped together.
        
""" % sys.argv[0]

import CGAT.Experiment as E
import networkx

###########################################################
def AddSynonym( rep_transcript, mem_transcript, map_synonyms, map_to_synonyms):

    if rep_transcript not in map_synonyms:
        map_synonyms[rep_transcript] = []

    if mem_transcript in map_synonyms:            

        for s in map_synonyms[mem_transcript]:
            map_to_synonyms[s] = rep_transcript

        map_synonyms[rep_transcript] += map_synonyms[mem_transcript]                    
        del map_synonyms[mem_transcript]

    map_synonyms[rep_transcript].append( mem_transcript )
    map_to_synonyms[mem_transcript] = rep_transcript

###########################################################
def clusterBySpecies( links, options ):
    """cluster links in infile by species."""
    
    map_priority = {}
    for x in range( len(options.quality_priority) ):
        map_priority[options.quality_priority[x]] = x
        
    transcripts = {}
    genomes = {}
    graph = {}
    reciprocal_best = {}
    map_synonyms = {}
    map_to_synonyms = {}
    

    ## build synonym list. Use priority list
    for transcript1, transcript2, weight in links:
    
        weight = float(weight)

        species1, t1, gene1, quality1 = transcript1.split(options.separator)
        species2, t2, gene2, quality2 = transcript2.split(options.separator)

        if weight <= options.synonym_weight and species1 == species2:

            ## use priority to decide about synonym
            if map_priority[quality1] > map_priority[quality2]:
                transcript1, transcript2 = transcript2, transcript1

            AddSynonym( transcript1, transcript2, map_synonyms, map_to_synonyms )
            
    ## second parse for actual processing. Synonyms are mapped
    ninput, noutput, nskipped = 0, 0, 0
    
    for transcript1, transcript2, weight in links:

        weight = float(weight)
        ninput += 1
        
        if options.min_weight != None and weight < options.min_weight:
            nskipped += 1
            continue

        if options.max_weight != None and weight > options.max_weight:
            nskipped += 1
            continue

        if transcript1 in map_to_synonyms: transcript1 = map_to_synonyms[transcript1]
        if transcript2 in map_to_synonyms: transcript2 = map_to_synonyms[transcript2]        

        if transcript1 > transcript2: transcript1, transcript2 = transcript2, transcript1
        
        species1 = transcript1.split(options.separator)[0]
        species2 = transcript2.split(options.separator)[0]

        if species1 not in transcripts: transcripts[species1] = {}
        if species2 not in transcripts: transcripts[species2] = {}
        
        transcripts[species1][transcript1] = 1
        transcripts[species2][transcript2] = 1        

        genomes[species1] = 1
        genomes[species2] = 1

        key = "%s-%s" % (species1, species2)
        if key not in graph: graph[key] = []
        graph[key].append( (weight, transcript1, transcript2) )
        
        key = "%s-%s" % (species2, species1)
        if key not in graph: graph[key] = []
        graph[key].append( (weight, transcript1, transcript2) )

        ## compute reciprocal best matches
        if transcript1 not in reciprocal_best: reciprocal_best[transcript1] = {}
        if species2 not in reciprocal_best[transcript1]:
            reciprocal_best[transcript1][species2] = (weight, transcript2)
        else:
            if weight < reciprocal_best[transcript1][species2][0]:
                reciprocal_best[transcript1][species2] = (weight, transcript2)

        if transcript2 not in reciprocal_best: reciprocal_best[transcript2] = {}
        if species1 not in reciprocal_best[transcript2]:
            reciprocal_best[transcript2][species1] = (weight, transcript1)
        else:
            if weight < reciprocal_best[transcript2][species1][0]:
                reciprocal_best[transcript2][species1] = (weight, transcript1)

        noutput += 1
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i, ngenomes=%i\n" % (ninput, noutput, nskipped, len(genomes) ))

    ## Compute reciprocal best match graph
    subgraph = networkx.Graph()        
    for k, v in reciprocal_best.items():
        subgraph.add_node( k )
        
        this_genome = k.split(options.separator)[0]
        
        for other_genome, best in v.items():
            
            if this_genome == other_genome: continue
            
            best_weight, best_transcript = best
            if reciprocal_best[best_transcript][this_genome][1] == k:

                subgraph.add_edge( k, best_transcript )

    ## compute components
    components = networkx.connected_components( subgraph )

    if options.loglevel >= 1:
        options.stdlog.write( "# component sizes before pruning: %s\n" % ",".join(map( lambda x: str(len(x)), components)))

    ## prune components: eliminate sequences from the same species in each component
    new_components = []

    ## sort by species and quality, take first entry, the rest are synonyms
    for component in components:
        component = map( lambda x: (x.split(options.separator)[0], x.split(options.separator)[1], x), component)
        component.sort()
        new_component = []
        last_species = None
        rep_transcript = None
        for species, quality, transcript in component:
            
            if last_species == species:
                AddSynonym( rep_transcript, transcript, map_synonyms, map_to_synonyms )                
            else:
                new_component.append( transcript )
                last_species = species
                rep_transcript = transcript
                
        new_components.append( new_component )
        
    components = new_components

    return components, map_synonyms

###########################################################
def updateGraphs( transcript1, transcript2,
                  weight,
                  map_transcript2id,                  
                  map_id2transcripts,
                  ids, graph,
                  reciprocal_best):

    id1 = map_transcript2id( transcript1 )
    id2 = map_transcript2id( transcript2 )    

    if id1 not in map_id2transcripts: map_id2transcripts[id1] = {}
    if id2 not in map_id2transcripts: map_id2transcripts[id2] = {}

    map_id2transcripts[id1][transcript1] = 1
    map_id2transcripts[id2][transcript2] = 1        

    ids[id1] = 1
    ids[id2] = 1

    key = "%s-%s" % (id1, id2)
    if key not in graph: graph[key] = []
    graph[key].append( (weight, transcript1, transcript2) )

    key = "%s-%s" % (id2, id1)
    if key not in graph: graph[key] = []
    graph[key].append( (weight, transcript1, transcript2) )

    ## compute reciprocal best matches
    if transcript1 not in reciprocal_best: reciprocal_best[transcript1] = {}
    if id2 not in reciprocal_best[transcript1]:
        reciprocal_best[transcript1][id2] = (weight, transcript2)
    else:
        if weight < reciprocal_best[transcript1][id2][0]:
            reciprocal_best[transcript1][id2] = (weight, transcript2)

    if transcript2 not in reciprocal_best: reciprocal_best[transcript2] = {}
    if id1 not in reciprocal_best[transcript2]:
        reciprocal_best[transcript2][id1] = (weight, transcript1)
    else:
        if weight < reciprocal_best[transcript2][id1][0]:
            reciprocal_best[transcript2][id1] = (weight, transcript1)

###########################################################
def getReciprocalBestComponents( reciprocal_best, map_transcript2id ):
                
    ## compute reciprocal best match graph
    subgraph = networkx.Graph()        
    for transcript, ids in reciprocal_best.items():
        
        subgraph.add_node( transcript )

        this_id = map_transcript2id( transcript )
        
        for other_id, best in ids.items():

            if this_id == other_id: continue

            best_weight, best_transcript = best
            if reciprocal_best[best_transcript][this_id][1] == transcript:
                subgraph.add_edge( transcript, best_transcript )
                
    ## compute components
    return networkx.connected_components( subgraph )

#############################################################################################
def getPrunedComponents( components, map_transcript2id,
                         map_priority,
                         map_synonyms, map_to_synonyms, options ):
    """eliminate sequences from the same species in each component.

    Take the highest quality in case of conflict.
    """
    new_components = []

    ## sort by transcript and quality
    for component in components:

        cc = []
        for transcript in component:
            s, t, g, q = transcript.split(options.separator)
            id = map_transcript2id( transcript )
            cc.append( (id, map_priority[q], transcript ) )
        
        cc.sort()
        new_component = []
        last_id = None
        rep_transcript = None
        for transcript in component:

            id = map_transcript2id(transcript)
            
            if last_id == id:
                AddSynonym( rep_transcript, transcript, map_synonyms, map_to_synonyms )                
            else:
                new_component.append( transcript )
                last_id = id
                rep_transcript = transcript
                
        new_components.append( new_component )
        
    return new_components
    
###########################################################
def clusterByReciprocity( links, map_transcript2id, options ):
    """given a graph of links between transcripts, split into
    clusters by transcripts."""

    ## second parse for actual processing. Synonyms are mapped

    map_id2transcripts = {}
    ids = {}
    graph = {}
    reciprocal_best = {}
    map_synonyms = {}
    map_to_synonyms = {}

    map_priority = {}
    for x in range( len(options.quality_priority) ):
        map_priority[options.quality_priority[x]] = x

    ninput, noutput, nskipped = 0, 0, 0
    for transcript1, transcript2, weight in links:

        ninput += 1
        weight = float(weight)

        if options.min_weight != None and weight < options.min_weight:
            nskipped += 1
            continue

        if options.max_weight != None and weight > options.max_weight:
            nskipped += 1
            continue

        if transcript1 > transcript2: transcript1, transcript2 = transcript2, transcript1
        noutput += 1
        updateGraphs( transcript1, transcript2, weight,
                      map_transcript2id,
                      map_id2transcripts, ids, graph, reciprocal_best )
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i, nids=%i\n" % (ninput, noutput, nskipped, len(ids) ))

    components = getReciprocalBestComponents( reciprocal_best, map_transcript2id )

    if options.loglevel >= 1:
        options.stdlog.write( "# component sizes before pruning: %s\n" % ",".join(map( lambda x: str(len(x)), components)))

    components = getPrunedComponents( components, map_transcript2id,
                                      map_priority,
                                      map_synonyms, map_to_synonyms,
                                      options )
    
    if options.loglevel >= 1:
        options.stdlog.write( "# component sizes after pruning: %s\n" % ",".join(map( lambda x: str(len(x)), components)))

    return components, map_synonyms

##-------------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph_cluster_by_species.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-m", "--genome-master", dest="genome_master", type="string",
                      help="genome to use as master." )

    parser.add_option("-s", "--filename-synonyms", dest="filename_synonyms", type="string",
                      help="output filename for synonyms." )

    parser.add_option("-u", "--filename-summary", dest="filename_summary", type="string",
                      help="output filename of component summary." )

    parser.add_option("-w", "--synonym-weight", dest="synonym_weight", type="float",
                      help="weight for grouping together synonyms." )

    parser.add_option("-q", "--quality-priority", dest="quality_priority", type="string",
                      help="priority list of status codes for sorting.")

    parser.add_option("--min-weight", dest="min_weight", type="float",
                      help="minimum weight to use.")

    parser.add_option("--max-weight", dest="max_weight", type="float",
                      help="maximum weight to use.")

    parser.add_option("-o", "--mode", dest="mode", type="choice",
                      choices=("species", "genes"),
                      help="clustering mode - species or genes.")
    
    parser.set_defaults( \
        genome_master = None,
        synonym_weight = 0,
        filename_synonyms = None,
        filename_summary = None,
        quality_priority="CG,PG,SG,RG,CP,PP,SP,RP,CF,PF,SF,UG,UP,UF,BF,UK",
        mode="species",
        min_weight = None,
        max_weight = None,
        separator = "|")
        
    (options, args) = E.Start( parser )

    options.quality_priority = options.quality_priority.split(",")

    ###############################################################
    ###############################################################
    ###############################################################        
    ## input
    ###############################################################                

    links = map(lambda x: x[:-1].split("\t"), filter( lambda x: x[0] != "#", sys.stdin.readlines()))

    ###############################################################
    ###############################################################
    ###############################################################        
    ## clustering
    ###############################################################                

    if options.mode == "species":

        def map_transcript2id( transcript, separator= "|" ):
            s, t, g, q = transcript.split(separator)
            return s
        
        components, map_synonyms = clusterBySpecies( links, options )

    elif options.mode == "genes":
        
        def map_transcript2id( transcript, separator= "|" ):
            s, t, g, q = transcript.split(separator)
            return "%s%s%s" % (s, separator, g )
        
        components, map_synonyms = clusterByReciprocity( links, map_transcript2id, options )

    ###############################################################
    ###############################################################
    ###############################################################        
    ## output
    ###############################################################                
    n = 1
    options.stdout.write( "id\tcomponent\n" )    
    for component in components:
        for c in component:
            options.stdout.write( "%s\t%i\n" % (c, n))
        n += 1

    n = 1
    if options.filename_summary:
        outfile = open(options.filename_summary, "w")
        outfile.write("component\tsize\tnspecies\tnmaster\tnsynonyms\n" )
        for component in components:
            species = map( lambda x: x.split(options.separator)[0], component)
            nsynonyms = 0
            for c in component:
                if c in map_synonyms: nsynonyms += len(map_synonyms[c])
            outfile.write( "%i\t%i\t%i\t%i\t%i\n" % ( n,
                                                      len(component),
                                                      len(species),
                                                      len(filter( lambda x: x == options.genome_master, species)),
                                                      nsynonyms))

            n += 1
        
    if options.filename_synonyms and len(map_synonyms):
        outfile = open(options.filename_synonyms, "w")
        outfile.write("rep\tmem\n")
        for k,v in map_synonyms.items():
            for m in v:
                outfile.write( "%s\t%s\n" % (k, m))
        outfile.close()
        
    if options.loglevel >= 1:
        
        options.stdlog.write( "# ncomponents=%i, nsynonyms=%i\n" % (len(components), len(map_synonyms)))
        options.stdlog.write( "# final component sizes: %s\n" % ",".join(map( lambda x: str(len(x)), components)))

    E.Stop()
