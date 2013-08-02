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
data2stats.py - summary statistics on table columns/rows
========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes summary statistics (min, mean, median)
on rows/columns of a table.

Usage
-----

Example::

   python data2stats.py < table.in > stats.out

Type::

   python data2stats.py --help

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
import optparse

import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import scipy

def PrintValues( outfile, values,  options, prefix = "",titles = None):

    if options.flat or options.aggregate_column:

        if options.add_header:
            if prefix: outfile.write( "prefix\t" )
            
            if titles: outfile.write( "column\t" )
                
            print "\t".join( ("nval", "min", "max", "mean", "median", "stddev", "sum", "q1", "q3" ) )
        
        for x in range(len(values)):

            vals = values[x]

            if len(vals) == 0:

                if options.output_empty:
                    if titles: outfile.write( titles[x] + "\t" )
                    if prefix: outfile.write( prefix + "\t" )

                    outfile.write( "0" + "\tna" * 8  + "\n" )

                continue

            if titles: outfile.write( titles[x] + "\t" )
            if prefix: outfile.write( prefix + "\t" )

            vals.sort()
            if len(vals) > 4:
                q1 = options.value_format % vals[len(vals) // 4]
                q3 = options.value_format % vals[len(vals) * 3 // 4]
            else:
                q1 = options.value_format % vals[0]
                q3 = options.value_format % vals[-1]

            outfile.write( "\t".join( ( "%i" % len(vals),
                                        options.value_format % float(min(vals)),
                                        options.value_format % float(max(vals)),
                                        options.value_format % scipy.mean(vals),
                                        options.value_format % scipy.median(vals),
                                        options.value_format % scipy.std(vals),                                      
                                        options.value_format % reduce( lambda x, y: x+y, vals),
                                        q1, q3,
                                        )) + "\n")
            
    else:

        if titles:
            print "category\t%s" % string.join(titles,"\t")

        print "count\t%s"  % (string.join( map(lambda v: "%i" % len(v), values), "\t"))
        print "min\t%s"    % (string.join( map(lambda v: options.value_format % min(v), values), "\t"))
        print "max\t%s"    % (string.join( map(lambda v: options.value_format % max(v), values), "\t"))
        print "mean\t%s"   % (string.join( map(lambda v: options.value_format % scipy.mean(v), values), "\t"))
        print "median\t%s" % (string.join( map(lambda v: options.value_format % scipy.median(v), values), "\t"))
        print "stddev\t%s" % (string.join( map(lambda v: options.value_format % scipy.std(v), values), "\t"))
        print "sum\t%s"    % (string.join( map(lambda v: options.value_format % reduce( lambda x,y: x+y, v), values), "\t"))
        print "q1\t%s"     % (string.join( map(lambda v: options.value_format % scipy.stats.scoreatpercentile(v,per=25), values), "\t"))
        print "q3\t%s"     % (string.join( map(lambda v: options.value_format % scipy.stats.scoreatpercentile(v,per=75), values), "\t"))

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: data2stats.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    parser.add_option("--scale", dest="scale", type="float",
                      help="scale values.")
    parser.add_option("-a", "--aggregate-column", dest="aggregate_column", type="int",
                      help="use column to aggregate."  )
    parser.add_option("-i", "--no-title", dest ="titles", action="store_false",
                      help="do not use supplied column titles." )
    parser.add_option("-e", "--headers", dest ="headers", type="string",
                      help="headers." )
    parser.add_option("-r", "--rows", dest="rows", action="store_true",
                      help="data is in rows."  )
    parser.add_option("--ignore-zeros", dest="ignore_zeros", action="store_true",
                      help="ignore zero values."  )
    parser.add_option("-f", "--format", dest="value_format", type="string",
                      help="number format."  )
    parser.add_option("-x", "--flat", dest="flat", action="store_true",
                      help="flat format."  )
    parser.add_option("--skip-header", dest="add_header", action="store_false",
                      help="do not add header to flat format."  )
    parser.add_option("--write-header", dest="write_header", action="store_true",
                      help="write header and exit."  )
    parser.add_option("--skip-empty", dest="output_empty", action="store_false",
                      help="do not output empty columns."  )
    parser.add_option("--output-empty", dest="output_empty", action="store_true",
                      help="output empty columns."  )

    parser.set_defaults(
        columns = "all",
        min_value = None,
        max_value = None,
        scale = None,
        aggregate_column = None,
        titles = True,
        headers = None,
        rows = False,
        value_format = "%6.4f",
        flat = False,
        test = "%5.2f",
        add_header = True,
        write_header = False,
        ignore_zeros = False,
        output_empty = False,
        separator = "\t"
        )
    


    (options, args) = E.Start( parser, quiet = True )

    if options.columns not in ( "all", "all-but-first", "variable" ):
        options.columns = map(lambda x: int(x) -1 , options.columns.split(","))

    if options.headers: options.headers = options.headers.split(",")
    
    ## write header for flat output
    if options.write_header:
        print "\t".join( ("nval", "min", "max", "mean", "median", "stddev", "sum", "q1", "q3") )
        sys.exit(0)

    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    outfile = sys.stdout

    if len(lines) > 0:
        
        ncols = len(string.split(lines[0][:-1], "\t"))
        
        if options.columns == "all":
            options.columns = range(0,ncols)
        elif options.columns == "all-but-first":
            options.columns = range(1,ncols)
        elif options.columns == "variable":
            pass
            
        if options.rows:

            ## ignore first value: is row title            
            if options.columns != "variable":
                del options.columns[0]

            if options.titles: del lines[0]
                
            ## write header for flat output
            if options.flat:
                if options.headers:
                    x = options.headers[0]
                else:
                    x = "row"
                print "\t".join( (x, "nval", "min", "max", "mean", "median", "stddev", "sum", "q1", "q3") )
                options.add_header = False
                
            for l in lines:
                data = l[:-1].split(options.separator)

                if options.columns == "variable":
                    vals = data[1:]
                else:
                    vals = map(lambda x: data[x], options.columns)

                ## remove unknown values
                vals= [ float(x) for x in vals if x and x.lower() not in ("na", "nan" ) ]

                if options.ignore_zeros:
                    vals = [ x for x in vals if x != 0.0 ]

                ## now convert to float
                vals = map(float, vals)
                
                PrintValues( outfile, [vals], options, data[0] )
                
        else:

            last_aggregate = None

            if options.titles:
                data = lines[0][:-1].split("\t")
                
                if not options.headers:
                    options.headers = map( lambda x: data[x], options.columns )
                    del lines[0]

                if options.aggregate_column != None:
                    print "category\t%s" % string.join(options.headers,"\t")
                
            vals = [ [] for x in range(len(options.columns)) ]

            for l in lines:

                data = l[:-1].split("\t")

                for c in range(len(options.columns)):

                    try:
                        val = string.atof(data[options.columns[c]])

                    except IndexError:
                        print "# IndexError in line:", l[:-1]
                        continue
                    except ValueError:
                        continue

                    if options.aggregate_column != None:

                        if last_aggregate != data[options.aggregate_column]:

                            if last_aggregate:
                                PrintValues( outfile, vals, options, last_aggregate )

                            vals = [ [] for x in range(len(options.columns)) ]
                            last_aggregate = data[options.aggregate_column]

                    if options.scale:
                        val *= options.scale

                    if options.max_value != None and val > options.max_value:
                        val = options.max_value

                    if options.min_value != None and val < options.min_value:
                        val = options.min_value

                    vals[c].append( val )

            lines = None

            ## remove empty columns
            nvals = []
            titles = []
            for c in range(len(options.columns)):
                if vals[c] or options.output_empty:
                    nvals.append( vals[c] )

                    if options.headers:
                        titles.append( options.headers[c] )

            vals = nvals
            
            PrintValues( outfile, vals, options, last_aggregate, titles )
    
    else:
        if options.titles:
            titles = ["missing",]
        else:
            titles = []

        if options.output_empty:
            PrintValues( outfile, [ [], ], options, None, titles )

    if options.loglevel >= 1:
        E.Stop()
    










