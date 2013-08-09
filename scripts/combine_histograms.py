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
combine_histograms.py - combine several histograms
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

OPTIONS:
-n, --normalize         normalize each histogram
-m, --missing=          set missing values to #
-t, --titles            use supplied titles
-h, --headers=          set header to #, "auto" = filenames
-f, --format=           format to use when combining histograms
--format-value=         format for numbers
--format-bin=           format for bin values

-> relative frequencies
-> cumulative counts and frequencies in both directions

Usage
-----

Example::

   python combine_histograms.py --help

Type::

   python combine_histograms.py --help

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

param_long_options = ["missing=", "headers=", "titles", "normalize", 
                      "format=", "format-bin=", "format-value=", "sort=", "help",
                      "version"]
param_short_options = "v:ht:m:h:s:f:"

param_headers = None
param_normalize = 0

param_missing_value = 0

param_titles = True

param_format = None
param_format_bin = None
param_format_value = None

param_sort = None

import CGAT.Experiment as E
import CGAT.Histogram as Histogram

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)
                                      

    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(1)

    for o,a in optlist:
        if o in ( "--help",):
            print globals()["__doc__"]
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--headers"):
            param_headers = string.split( a, "," )
        elif o in ("-n", "--normalize"):
            param_normalize = 1
        elif o in ("-m", "--missing"):
            param_missing_value = a
        elif o == "--no-titles":
            param_titles = False
        elif o == "--no-titles":
            param_titles = False
        elif o in ("-f", "--format"):
            param_format = a
        elif o == "--format-value":
            param_format_value = a
        elif o == "--format-bin":
            param_format_bin = a
        elif o in ("-s", "--sort"):
            if a in ("numerical", "alphabetic"):
                param_sort = a
            else:
                param_sort = string.split( a, "," )

    if len(args) < 1:
        print globals()["__doc__"], "please specify at one histogram."
        sys.exit(1)

    param_filenames = args

    print E.GetHeader()
    print E.GetParams()    

    histograms = []

    ## first
    headers = ['bin',]
    if param_headers and headers != "auto":
        headers = [param_headers[0],]
        del param_headers[0]
        
    for x in range(len(param_filenames)):
        
        filename = param_filenames[x]
        if not os.path.exists(filename):
            print "# skipped because file not present: %s" % filename
            continue

        file = open(filename, "r")

        lines = filter( lambda x: x[0] <> "#", file)
        
        if len(lines) == 0:
            continue

        if param_titles:
            h = lines[0][:-1].split("\t")[1:]
            del lines[0]
                
        if param_headers == "auto":
            headers.append( os.path.basename( filename ) )
        elif param_headers:
            headers.append( param_headers[x] )
        elif param_titles:
            headers += h

        data = map( lambda x: map(string.atof, string.split(x[:-1], "\t")), lines)

        ## add empty data point for empty histograms
        if len(data) == 0: data = [ (0,0) ] 
        
        histograms.append( data )

    ## sort the whole thing:
    if param_sort:
        sort_order = []

        if param_sort == "numerical":
            t = zip( map(int, headers[1:]), range( 1, len(headers) + 1) )
            t.sort()

            for tt in t:
                sort_order.append( headers[tt[1]] )

        elif param_sort == "alphabetical":
            t = zip( headers[1:], range( 1, len(headers) + 1) )
            t.sort()

            for tt in t:
                sort_order.append( headers[tt[1]] )
        else:
            sort_order = param_sort

        ## map header to old position
        map_header2pos = {}
        for x in range(1, len(headers)):
            map_header2pos[headers[x]] = x

        order = []
        for x in sort_order:
            if x in map_header2pos:                
                order.append( map_header2pos[x] )

        new_headers = [headers[0]]
        new_histograms = []

        for x in order:
            new_headers.append( headers[x] )
            new_histograms.append( histograms[x-1] )

        histograms = new_histograms
        headers = new_headers
        
    combined_histogram = Histogram.Combine( histograms, param_missing_value )

    if headers:
        print "\t".join(headers)

    if param_normalize:
        combined_histogram = Histogram.Normalize( combined_histogram )

    Histogram.Print( combined_histogram,
                     format_bin=param_format_bin,
                     format_value=param_format_value,
                     )

    print E.GetFooter()


