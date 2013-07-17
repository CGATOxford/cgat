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
compare_clusters.py - compare two partitions of the same data
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Input: two maps.

Method: build a graph between clusters in the two partitions.
Clusters are connected, if they share common identifiers.

Usage
-----

Example::

   python compare_clusters.py --help

Type::

   python compare_clusters.py --help

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
import time
import random
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import networkx

def getFile( section, options ):
    if options.output_pattern:
        if options.loglevel >= 1:
            options.stdlog.write( "# output for section %s goes to %s\n" % (section, options.output_pattern % section))
        return open(options.output_pattern % section, "w" )
    else:
        return options.stdout

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: compare_clusters.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-o", "--output-pattern", dest="output_pattern", type="string" ,
                       help="output pattern for filenames.")

    parser.set_defaults(
        output_pattern = None,
        format = "%5.2f",
        )
    
    (options, args) = E.Start( parser,
                                        add_pipe_options = True )

    if len(args) != 2:
        raise "please supply to filenames with the clusters."

    map_id2cluster1, map_cluster2ids1 = IOTools.ReadMap( open( args[0] ), both_directions = True )
    map_id2cluster2, map_cluster2ids2 = IOTools.ReadMap( open( args[1] ), both_directions = True )

    graph = networkx.Graph()
    
    for a in map_cluster2ids1.keys(): graph.add_node( (1,a) )
    for b in map_cluster2ids2.keys(): graph.add_node( (2,b) )

    ## build graph between clusters
    for cluster1, ids1 in map_cluster2ids1.items():
        for id1 in ids1:
            if id1 in map_id2cluster2:
                graph.add_edge( (1,cluster1), (2,map_id2cluster2[id1]) )

    components = networkx.connected_components( graph )

    #######################################################
    #######################################################
    #######################################################    
    ## write components and compute counts
    #######################################################    
    outfile = getFile( "components", options )
    outfile.write("id\ttotal\tn1\tn2\tmembers1\tmembers2\n" )
    n = 0
    counts = {}
    subsets = []
    for component in components:

        m1, m2 = [], []
        
        for x in component:
            if x[0] == 1:
                m1.append(x[1])
            else:
                m2.append(x[1])
                
        t = len(component)
        n1 = len(m1)
        n2 = len(m2)
        cc = ( n1, n2 )
        if cc not in counts: counts[cc] = 0
        counts[cc] += 1
        
        if cc == (1,1): subsets.append( n )
        
        n += 1
        outfile.write( "%i\t%i\t%i\t%i\t%s\t%s\n" % (n, t, n1, n2, ",".join(m1), ",".join(m2) ))

    if outfile != options.stdout:
        outfile.close()
    else:
        outfile.write ("//\n")

    #######################################################
    #######################################################
    #######################################################    
    ## write counts
    #######################################################    
    outfile = getFile( "counts", options )
    outfile.write("n1\tn2\tcounts\tpcounts1\tpcounts2\n" )
    for cc, c in counts.items():
        outfile.write("%i\t%i\t%i\t%s\t%s\n" % ( cc[0], cc[1], c,
                                                 options.format % (100.0 * float(c) / len(map_cluster2ids1)),
                                                 options.format % (100.0 * float(c) / len(map_cluster2ids2)) ))

    if outfile != options.stdout:
        outfile.close()
    else:
        outfile.write ("//\n")

    #######################################################
    #######################################################
    #######################################################    
    ## analyze subsets - how many of the 1:1 clusters
    ## contain the exact members?
    #######################################################    
    outfile = getFile( "subsets", options )
    outfile.write("id\tn1\tn2\tunion\tinter\tunique1\tunique2\n" )    

    ntrue = 0
    nrest1 = 0
    nrest2 = 0
    nother = 0
    
    for component_id in subsets:
        component = components[component_id]
        if component[0][0] == 1:
            id1, id2 = component[0][1], component[1][1]
        else:
            id1, id2 = component[1][1], component[0][1]

        members1 = set(map_cluster2ids1[id1])
        members2 = set(map_cluster2ids2[id2])

        union = len(members1.union(members2))
        intersection = len(members1.intersection(members2))
        rest1 = len(members1.difference(members2))
        rest2 = len(members2.difference(members1))

        if rest1 == 0 and rest2 == 0:
            ntrue += 1
        elif rest1 == 0:
            nrest1 += 1
        elif rest2 == 0:
            nrest2 += 1
        else:
            nother += 1
        
        outfile.write( "%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % ( component_id,
                                                          len(members1),
                                                          len(members2),
                                                          union,
                                                          intersection,
                                                          rest1, rest2))

    if outfile != options.stdout:
        outfile.close()
    else:
        outfile.write ("//\n")

    ## write subset statistics
    ntotal = len(subsets)
    options.stdout.write("# subset statistics of 1:1 corresponding clusters\n" )
    options.stdout.write( "class\tcounts\ttotal\n")
    options.stdout.write( "%s\t%i\t%s\n" % ("total", ntotal, options.format % 100 ))
    options.stdout.write( "%s\t%i\t%s\n" % ("true", ntrue, options.format % (100.0 * ntrue / ntotal )))
    options.stdout.write( "%s\t%i\t%s\n" % ("unique1", nrest1, options.format % (100.0 * nrest1 / ntotal )))
    options.stdout.write( "%s\t%i\t%s\n" % ("unique2", nrest2, options.format % (100.0 * nrest2 / ntotal )))
    options.stdout.write( "%s\t%i\t%s\n" % ("other", nother, options.format % (100.0 * nother / ntotal )))            

    E.Stop()
