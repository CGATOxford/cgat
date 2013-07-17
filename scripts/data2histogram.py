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
data2histogram.py - histogram data in a table
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes histograms over one or more
columns of a table.

Usage
-----

Example::

   python data2histogram.py --help

Type::

   python data2histogram.py --help

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
import types

import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import numpy

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: data2histogram.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-r", "--range", dest="range", type="string",
                      help="range to calculate histogram for."  )
    parser.add_option("-b", "--bin-size", dest="bin_size", type="string",
                      help="bin size."  )
    parser.add_option("-i", "--titles", dest ="titles", action="store_true",
                      help="use supplied column titles." )
    parser.add_option( "--no-titles", dest ="titles", action="store_false",
                      help="no column titles given." )
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )
    parser.add_option("--min-data", dest="min_data", type="int",
                      help="minimum amount of data required, if less data, then the histogram will be empty [default=%default].")
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    parser.add_option("--no-empty-bins", dest="no_empty_bins", action="store_true",
                      help="do not display empty bins.")
    parser.add_option("--with-empty-bins", dest="no_empty_bins", action="store_false",
                      help="display empty bins.")
    parser.add_option("--normalize", dest="normalize", action="store_true",
                      help="normalize histogram.")
    parser.add_option("--cumulative", dest="cumulative", action="store_true",
                      help="calculate cumulative histogram.")
    parser.add_option("--reverse-cumulative", dest="reverse_cumulative", action="store_true",
                      help="calculate reverse cumulative histogram.")
    parser.add_option("--headers", dest="headers", type="string",
                      help="use the following headers.")
    parser.add_option("--ignore-out-of-range", dest="ignore_out_of_range", action="store_true",
                      help="ignore values that are out of range (as opposed to truncating them to range border.")
    parser.add_option("--missing", dest="missing_value", type="string",
                      help="entry for missing values [%default]." )
    parser.add_option("--dynamic-bins", dest="dynamic_bins", action="store_true",
                      help="each value constitutes its own bin." )
    parser.add_option("--on-the-fly", dest="on_the_fly", action="store_true",
                      help="on the fly computation of histograms. Requires setting of min-value, max-value and bin_size." )
    
    parser.set_defaults(
        bin_size = None,
        range = None,
        titles = True,
        columns = "all",
        append = (),
        no_empty_bins = True,
        min_value = None,
        max_value = None,
        normalize = False,
        cumulative = False,
        reverse_cumulative = False,
        nonull = None,
        ignore_out_of_range = False,
        min_data = 1,
        headers = None,
        missing_value = "na",
        dynamic_bins = False,
        on_the_fly = False,
        bin_format = "%.2f",
        value_format = "%6.4f",
        )

    (options, args) = E.Start( parser )

    if options.columns != "all":
        options.columns = map(lambda x: int(x) -1 , options.columns.split(","))

    if options.range:
        options.min_value, options.max_value = map(float, options.range.split(","))

    if options.headers:
        options.headers = options.headers.split(",")
        
    if options.on_the_fly:
        if options.min_value == None or options.max_value == None or options.bin_size == None:
            print USAGE
            raise "please supply columns, min-value, max-value and bin-size for on-the-fly computation."

        # try to glean titles from table:
        if options.titles:
            while 1:
                line = sys.stdin.readline()
                if not line: break
                if line[0] == "#": continue
                data = line[:-1].split("\t")
                break

            if options.columns == "all":
                options.titles = data
                options.columns = range( len(data) )
            else:
                options.titles = [ data[x] for x in options.columns ]
                
        bins = numpy.arange( options.min_value, options.max_value, float(options.bin_size))
        hh = Histogram.fillHistograms( sys.stdin, options.columns, [ bins for x in range(len(options.columns) ) ] )
        n = len(hh)

        titles = ['bin']

        if options.headers:
            titles.append( options.headers[x] )
        elif options.titles:
            titles.append( options.titles[x] )
        else:
            for x in options.columns:
                titles.append( "col%i" % (x+1) )

        if len(titles) > 1:
            options.stdout.write( "\t".join(titles) + "\n" )
        
        for x in range(len(bins)):
            v = []
            v.append( options.bin_format % bins[x] )
            for c in range(n):
                v.append( options.value_format % hh[c][x] )

            options.stdout.write( "\t".join( v ) + "\n" )

    else:
        ## in-situ computation of histograms
        # retrieve data
        first = True
        vals = []

        # parse data, convert to floats
        for l in options.stdin:

            if l[0] == "#": continue

            data = string.split(l[:-1], "\t")

            if first:
                first = False
                ncols = len(data)
                if options.columns == "all":
                    options.columns = range(ncols)

                vals = [ [] for x in options.columns ]

                if options.titles:
                    try:
                        options.titles = [ data[x] for x in options.columns ]
                    except IndexError:
                        raise IndexError, "not all columns %s found in data %s" % (str(options.columns), str(data))
                    continue

            for x in range(len(options.columns)):

                try:
                    v = string.atof(data[options.columns[x]])
                except IndexError:
                    print "# IndexError in line:", l[:-1]
                    continue
                except ValueError:
                    continue

                vals[x].append( v )

        lines = None

        hists = []
        titles = []
        
        if not vals:
            if options.loglevel >= 1:
                options.stdlog.write( "# no data\n" )
            E.Stop()
            sys.exit(0)

        for x in range(len(options.columns)):

            if options.loglevel >= 1:
                options.stdlog.write( "# column=%i, num_values=%i\n" % (options.columns[x], len(vals[x])) )

            if len(vals[x]) < options.min_data: continue

            h = Histogram.Calculate( vals[x],
                                     no_empty_bins = options.no_empty_bins,
                                     increment = options.bin_size,
                                     min_value = options.min_value,
                                     max_value = options.max_value,
                                     dynamic_bins = options.dynamic_bins,
                                     ignore_out_of_range = options.ignore_out_of_range )

            if options.normalize: h = Histogram.Normalize( h )
            if options.cumulative: h = Histogram.Cumulate( h )
            if options.reverse_cumulative: h = Histogram.Cumulate( h, direction = 0 )

            hists.append(h)

            for m in options.append:
                if m == "normalize":
                    hists.append( Histogram.Normalize( h ) )

            if options.headers:
                titles.append( options.headers[x] )
            elif options.titles:
                titles.append( options.titles[x] )
            else:
                titles.append( "col%i" % options.columns[x] )

        if titles:
            options.stdout.write( "bin\t" + "\t".join(titles) + "\n" )

        if len(hists) == 1:
            Histogram.Print( hists[0], nonull = options.nonull )
        else:
            combined_histogram = Histogram.Combine( hists, missing_value = options.missing_value )
            Histogram.Print( combined_histogram, nonull = options.nonull )        

    E.Stop()






