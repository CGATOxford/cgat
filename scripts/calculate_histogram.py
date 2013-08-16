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
calculate_histogram.py - calculate histogram from data
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script calculates histograms from data in a
tab-separated table.

read in data and build histogram of column.

Options::

    -c, --column            column to take [default = 0]
    -a, --append=           append columns [normalize]
    -n, --normalize         normalize column
    --cumulative            cumulative histogram
    --reverse-cumulative    reverse cumulative histogram
    -i, --titles            use supplied titles
    # at start of line is a comment

Usage
-----

Example::

   python calculate_histogram.py --help

Type::

   python calculate_histogram.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import string
import os
import getopt
import time

import CGAT.Experiment as E
import CGAT.Histogram as Histogram

param_loglevel = 1
param_separator = "//"
param_take = None
param_fill = 0
param_nonull = None
param_columns = [0,]
param_empty_bins = 1
param_titles = False

param_lower_limit = None
param_upper_limit = None
param_bin_size = None
param_scale = None
param_normalize = None
param_append = []

param_cumulative = False
param_reverse_cumulative = False

param_long_options = ["Verbose=", "nonull", "fill","take=", "column=", "show_empty",
                      "upper=", "lower=", "bin-size=", "scale=", "normalize", "append=", "titles",
                      "cumulative", "reverse-cumulative", "help",
                      "version"]

param_short_options = "v:nft:c:eu:l:b:a:ioh"

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)

    except getopt.error, msg:
        print globals()["__doc__"]
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-t", "--take"):
            param_take = map(string.atoi, string.split(a, ","))
        elif o in ("-f", "--fill"):
            param_fill = 1
        elif o in ("-n", "--nonull"):
            param_nonull = ""
        elif o in ("-c", "--column"):
            if a == "all":
                param_columns = "all"
            else:
                param_columns = map(lambda x: int(x) - 1, string.split( a, ","))
        elif o in ("-e", "--show_empty"):
            param_empty_bins = 0
        elif o in ("-u", "--upper"):
            param_upper_limit = string.atof(a)
        elif o in ("-l", "--lower"):
            param_lower_limit = string.atof(a)
        elif o in ("-b", "--bin-size"):
            param_bin_size = string.atof(a)
        elif o in ("-s", "--scale"):
            param_scale = float(a)
        elif o in ("-o", "--normalize"):
            param_normalize = True
        elif o in ("-a", "--append"):
            param_append = string.split( a, ",")
        elif o in ("-i", "--titles"):
            param_titles = True
        elif o == "--cumulative":
            param_cumulative = True
        elif o == "--reverse-cumulative":
            param_reverse_cumulative = True
        
    histograms = []

    if param_loglevel >= 1:
        print E.GetHeader()
        print E.GetParams()    
    
    vals = []
    
    for x in param_columns: vals.append( [] )
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    ncols = len(string.split(lines[0][:-1], "\t"))
    if param_columns == "all":
        param_columns = range(ncols)
        for x in param_columns: vals.append( [] )

    if param_titles:
        data = lines[0][:-1].split("\t")
        del lines[0]
        param_titles = map( lambda x: data[x], param_columns)
        
    for l in lines:
        data = string.split(l[:-1], "\t")
            
        for x in range(len(param_columns)):
            try:
                v = string.atof(data[param_columns[x]])
            except IndexError:
                print "# IndexError in line:", l[:-1]
                continue
            except ValueError:
                continue

            if param_scale:
                v *= param_scale

            if param_upper_limit != None and v > param_upper_limit:
                v = param_upper_limit

            if param_lower_limit != None and v < param_lower_limit:
                v = param_lower_limit

            vals[x].append( v )

    lines = None

    hists = []
    titles = []
    
    for x in range(len(param_columns)):
        if param_loglevel >= 1:
            print "# column=%i, num_values=%i" % (param_columns[x], len(vals[x]))

        if len(vals[x]) == 0: continue
        
        h = Histogram.Calculate( vals[x], no_empty_bins = param_empty_bins, increment = param_bin_size)
        if param_scale: h = Histogram.Scale( h, 1.0 / param_scale )

        if param_normalize: h = Histogram.Normalize( h )
        if param_cumulative: h = Histogram.Cumulate( h )
        if param_reverse_cumulative: h = Histogram.Cumulate( h, direction = 0 )
        
        hists.append(h)

        for m in param_append:
            if m == "normalize":
                hists.append( Histogram.Normalize( h ) )

        if param_titles:
            titles.append( param_titles[x] )

    if titles:
        print "bin\t" + "\t".join(titles)

    if len(hists) == 1:
        Histogram.Print( hists[0], nonull = param_nonull )
    else:
        combined_histogram = Histogram.Combine( hists )
        Histogram.Print( combined_histogram, nonull = param_nonull )        

    if param_loglevel >= 1:        
        print E.GetFooter()
    










