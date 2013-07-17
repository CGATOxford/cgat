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
graph2stats.py - calculate statistics for a (redundant) graph 
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a :term:`graph` in :term:`edge list` format and 
computes stats on each edges. This only make sense, if the graph
contains multiple edges between a pair of vertices.


Usage
-----

Example::

   python graph2stats.py --help

Type::

   python graph2stats.py --help

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
import time
import optparse
import math
import tempfile

""" program $Id: graph2stats.py 2782 2009-09-10 11:40:29Z andreas $

python graph2stats.py < graph.in

calculate statistics for a (redundant) graph.

Possibilities are: 
"""
import CGAT.Experiment as E
import CGAT.Histogram as Histogram

def PrintHistograms( outfile, titles, histograms, options ):

    combined_histogram = Histogram.Combine( hists )

    outfile.write( "\t".join( ("bin",) + titles ) )
    Histogram.Print( combined_histogram, nonull = options.nonull )        

    
    

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph2stats.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-r", "--range", dest="range", type="string",
                      help="range to calculate histogram for."  )
    parser.add_option("-b", "--bin-size", dest="bin_size", type="string",
                      help="bin size."  )
    parser.add_option("-i", "--titles", dest ="titles", action="store_true",
                      help="use supplied column titles." )
    parser.add_option("-s", "--make-symmetric", dest ="make_symmetric", action="store_true",
                      help="symmetrize graph." )
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )
    parser.add_option("-p", "--output-pattern", dest="output_pattern", type="string",
                      help="pattern for output files." )
    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method." )
    parser.add_option("-o", "--output-format", dest="output_format", type="string",
                      help="output format." )
    
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    
    parser.set_defaults(
        bin_size = None,
        range = None,
        titles = False,
        columns = "all",
        append = (),
        empty_bins = False,
        min_value = None,
        max_value = None,
        normalize = False,
        cumulative = False,
        reverse_cumulative = False,
        nonull = None,
        make_symmetric = False,
        output_pattern = "%s.hist",
        method = "histograms",
        output_format = "semi"
        )

    (options, args) = E.Start( parser )

    if options.columns != "all":
        options.columns = map(lambda x: int(x) -1 , options.columns.split(","))

    if options.range:
        options.min_value, options.max_value = map(float, options.range(split(",")))
        
    # retrieve data
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    vals = {}

    if options.method == "histograms":

        ## read data
        for line in lines:

            v1, v2, w = line[:-1].split("\t")[:3]

            try:
                w[3] = float(w[3])
            except ValueError:
                nerrors += 1
                continnue

            if v1 not in vals: vals[v1] = {}
            if v2 not in vals[v1]: vals[v1][v2] = []
            vals[v1][v2].append( w )
            if options.make_symmetric:
                if v2 not in vals: vals[v2] = {}
                if v1 not in vals[v2]: vals[v2][v1] = []
                vals[v2][v1].append( w )

        ## convert to histograms
        for k1, vv in vals.items():
            for k2 in vv.keys():
                if len(vv[k2]) == 0: continue

                h = Histogram.Calculate( vv[k2],
                                         no_empty_bins = options.empty_bins,
                                         increment = options.bin_size,
                                         min_value = options.min_value,
                                         max_value = options.max_value)
                
                if options.normalize: h = Histogram.Normalize( h )
                if options.cumulative: h = Histogram.Cumulate( h )
                if options.reverse_cumulative: h = Histogram.Cumulate( h, direction = 0 )

                vv[k2] = h
        
        ## write output
        if options.output == "semi":
            for k1, vv in vals.items():

                outfile = open( options.output_pattern % k1 )

                kk2 = vv.keys()
                kk2.sort()

                hists = []
                for k2 in kk2: hists.append( vv[k2] )

                PrintHistograms( outfile, kk2, hists, options )

                outfile.close()
                
    elif options.method == "counts":

        ## read data
        for line in lines:

            v1, v2 = line[:-1].split("\t")[:2]

            if v1 not in vals: vals[v1] = {}
            if v2 not in vals[v1]: vals[v1][v2] = 0
            vals[v1][v2] += 1
            if options.make_symmetric:
                if v2 not in vals: vals[v2] = {}
                if v1 not in vals[v2]: vals[v2][v1] = 0
                vals[v2][v1] += 1

        ## convert to histograms
        for k1, vv in vals.items():
            for k2 in vv.keys():
                options.stdout.write( "%s\t%s\t%i\n" % (k1, k2, vv[k2] ) )
                
    E.Stop()

