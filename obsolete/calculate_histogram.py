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

Usage
-----

Example::

   python calculate_histogram.py < in.data > out.tsv

Type::

   python calculate_histogram.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import string
import os
import getopt
import time

import CGAT.Experiment as E
import CGAT.Histogram as Histogram

##---------------------------------------------------------------------------------------------------------        
def main( argv = None ):
    
    if argv == None: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                             usage = globals()["__doc__"] )

    parser.add_option("-n", "--nonull", dest="nonull", action = "store_true",
                      help="no null [default=%default]"  )

    parser.add_option("-e", "--show-empty", dest="empty_bins", action = "store_true",
                      help="show empty bins [default=%default]"  )

    parser.add_option("-o", "--normalize", dest="normalize", action = "store_true",
                      help="normalize histogram [default=%default]"  )

    parser.add_option("-i", "--titles", dest="titles", action = "store_true",
                      help="use titles supplied in ... [default=%default]"  )

    parser.add_option( "--cumulative", dest="cumulative", action = "store_true",
                      help="compute cumulative histogram [default=%default]"  )

    parser.add_option( "--reverse-cumulative", dest="reverse_cumulative", action = "store_true",
                      help="compute reverse cumulative histogram [default=%default]"  )

    parser.add_option( "-c", "--column", dest="column", type = "int",
                      help="columns to take [default=%default]"  )
    
    parser.add_option( "-b", "--bin-size", dest="bin_size", type = "float",
                      help="bin size to use [default=%default]"  )

    parser.add_option( "-u", "--upper", dest="upper_limit", type = "float",
                      help="upper limit to use [default=%default]"  )

    parser.add_option( "-l", "--lower", dest="lower_limit", type = "float",
                      help="lower limit to use [default=%default]"  )

    parser.add_option( "-s", "--scale", dest="scale", type = "float",
                      help="scale to use [default=%default]"  )

    parser.add_option( "-a", "--append", dest="append", type = "choice", action="append",
                       choices = ("normalize", ),
                       help="append columns [default=%default]"  )

    parser.set_defaults(
        nonull = None,
        columns = [0,],
        empty_bins = True,
        titles = False,
        lower_limit = None,
        upper_limit = None,
        bin_size = None,
        scale = None,
        normalize = None,
        append = [],
        cumulative = False,
        reverse_cumulative = False )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if options.columns:
        if options.columns != "all":
            options.columns = [ int(x) - 1 for x in options.columns.split( ",") ]
    else:
        options.columns.append( 0 )

    histograms = []
    
    vals = []
    
    for x in options.columns: vals.append( [] )
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    ncols = len(string.split(lines[0][:-1], "\t"))
    if options.columns == "all":
        options.columns = range(ncols)
        for x in options.columns: vals.append( [] )

    if options.titles:
        data = lines[0][:-1].split("\t")
        del lines[0]
        options.titles = map( lambda x: data[x], options.columns)
        
    for l in lines:
        data = string.split(l[:-1], "\t")
            
        for x in range(len(options.columns)):
            try:
                v = string.atof(data[options.columns[x]])
            except IndexError:
                print "# IndexError in line:", l[:-1]
                continue
            except ValueError:
                continue

            if options.scale:
                v *= options.scale

            if options.upper_limit != None and v > options.upper_limit:
                v = options.upper_limit

            if options.lower_limit != None and v < options.lower_limit:
                v = options.lower_limit

            vals[x].append( v )

    lines = None

    hists = []
    titles = []
    
    for x in range(len(options.columns)):
        E.info( "column=%i, num_values=%i" % (options.columns[x], len(vals[x])) )

        if len(vals[x]) == 0: continue
        
        h = Histogram.Calculate( vals[x], no_empty_bins = options.empty_bins, increment = options.bin_size)
        if options.scale: h = Histogram.Scale( h, 1.0 / options.scale )

        if options.normalize: h = Histogram.Normalize( h )
        if options.cumulative: h = Histogram.Cumulate( h )
        if options.reverse_cumulative: h = Histogram.Cumulate( h, direction = 0 )
        
        hists.append(h)

        for m in options.append:
            if m == "normalize":
                hists.append( Histogram.Normalize( h ) )

        if options.titles:
            titles.append( options.titles[x] )

    if titles:
        options.stdout.write( "bin\t" + "\t".join(titles) + "\n" )

    if len(hists) == 1:
        Histogram.Print( hists[0], nonull = options.nonull )
    else:
        combined_histogram = Histogram.Combine( hists )
        Histogram.Print( combined_histogram, nonull = options.nonull )        

    E.Stop()

if __name__ == '__main__':
    sys.exit(main(sys.argv))








