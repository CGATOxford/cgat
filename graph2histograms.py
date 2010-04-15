################################################################################
#   Gene prediction pipeline 
#
#   $Id: graph2histograms.py 2782 2009-09-10 11:40:29Z andreas $
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
import os, sys, string, re, getopt, time, optparse, math, tempfile

""" program $Id: graph2histograms.py 2782 2009-09-10 11:40:29Z andreas $

python graph2stats.py < graph.in

calculate statistics for a (redundant) graph.

Possibilities are: 
"""
import Experiment, Histogram

import scipy, scipy.stats

def PrintHistograms( outfile, histograms, options ):

    hists = []
    if options.sort:
        titles = []
        for k2 in options.sort:
            if k2 in histograms:
                hists.append( vv[k2] )
                titles.append( k2 )
    else:
        titles = histograms.keys()
        titles.sort()
        for k2 in titles:
            hists.append( vv[k2] )

    combined_histogram = Histogram.Combine( hists )

    outfile.write( "\t".join( ("bin",) + tuple(titles) ) + "\n" )

    Histogram.Write( outfile, combined_histogram, nonull = options.nonull )        


if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: graph2histograms.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-r", "--range", dest="range", type="string",
                      help="range to calculate histogram for."  )
    parser.add_option("-b", "--bin-size", dest="bin_size", type="string",
                      help="bin size."  )
    parser.add_option("-i", "--titles", dest ="titles", action="store_true",
                      help="use supplied column titles." )
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )
    parser.add_option("-p", "--output-pattern", dest="output_pattern", type="string",
                      help="pattern for output files." )
    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="method." )
    parser.add_option("-o", "--output-format", dest="output_format", type="string",
                      help="output format." )
    parser.add_option("-s", "--no-self", dest="no_self", action="store_true",
                      help="do not allow self links." )
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    parser.add_option("--sort", dest="sort", type="string",
                      help="sort order for categories")
    
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
        make_symmetric = True,
        output_pattern = "%s.hist",
        methods = "histograms",
        output_format = "semi",
        output_number = "%6.4f",
        no_self = False,
        sort = "",
        )

    (options, args) = Experiment.Start( parser )

    if options.columns != "all":
        options.columns = map(lambda x: int(x) -1 , options.columns.split(","))

    if options.sort: options.sort=options.sort.split(",")
    
    if options.range:
        options.min_value, options.max_value = map(float, options.range(split(",")))

    if options.methods: options.methods = options.methods.split(",")
    # retrieve data
    lines = filter( lambda x: x[0] not in ( "#", ">"), sys.stdin.readlines())

    vals = {}
    nerrors = 0
    
    ## read data
    for line in lines:

        v1, v2, w = line[:-1].split("\t")[:3]

        if options.no_self and v1 == v2:
            continue

        try:
            w = float(w)
        except ValueError:
            nerrors += 1
            continnue

        if v1 not in vals: vals[v1] = {}
        if v2 not in vals[v1]: vals[v1][v2] = []
        vals[v1][v2].append( w )
        if options.make_symmetric:
            if v1 == v2: continue

            if v2 not in vals: vals[v2] = {}
            if v1 not in vals[v2]: vals[v2][v1] = []
            vals[v2][v1].append( w )

    for method in options.methods:
        
        if method == "histograms":

            ## convert to histograms
            for k1, vv in vals.items():
                for k2 in vv.keys():
                    if len(vv[k2]) == 0: continue

                    if options.loglevel >= 1:
                        print "# calculating histogram for %s %s: %i values" % (k1, k2, len(vv[k2]))

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
            if options.output_format == "semi":
                for k1, vv in vals.items():

                    outfile = open( options.output_pattern % k1, "w" )
                    PrintHistograms( outfile, vv, options )
                    outfile.close()

        elif method == "summary":

            outfile = sys.stdout

            ## get summary statistics
            outfile.write( "token1\ttoken2\tcount\tmin\tmax\tmean\tmedian\tstd\tsum\n" )

            if options.sort:
                vv1 = options.sort
            else:
                vv1 = vals.keys()
                vv1.sort()

            for v1 in vv1:
                if v1 not in vals: continue
                vv = vals[v1]

                if options.sort:
                    vv2 = options.sort
                else:
                    vv2 = vv.keys()
                    vv2.sort()

                for v2 in vv2:
                    if v2 not in vv: continue
                    myvals = vv[v2]
                    outfile.write( "\t".join( (
                        v1, v2,
                        "%i" % len(myvals),
                        options.output_number % min(myvals),
                        options.output_number % max(myvals),                    
                        options.output_number % scipy.mean(myvals),
                        options.output_number % scipy.median(myvals),                    
                        options.output_number % scipy.std(myvals),
                        options.output_number % reduce( lambda x,y: x+y, myvals) ) ) + "\n" )

            if outfile != sys.stdout:
                outfile.close()
        
    Experiment.Stop()

