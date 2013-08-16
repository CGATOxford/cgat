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
plot_data.py - create scatter plot from tabular data
====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts creates a scatter plot from tabular data.

Reads data in tabular format with observations in rows in
different samples in columns. The first column provides the
bins.

Reads data from stdin unless infile is given.

Usage
-----

Example::

   python plot_data.py --help

Type::

   python plot_data.py --help

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
import optparse

import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import CGAT.MatlabTools as MatlabTools
import numpy
import numpy.ma as MA

def readTable( file,
               separator = "\t",
               numeric_type = numpy.float,
               take = "all",
               has_headers = True,
               truncate = None,
               missing_value = -99999,
               color = None,
               ):
    """read a matrix.

    There probably is a routine for this in Numpy, which I haven't found yet.
    Uses the MA module to deal with missing values.
    """

    lines = filter( lambda x: x[0] != "#", file.readlines())

    if len(lines) == 0:
        raise ValueError( "no data" )

    if take == "all":
        num_cols = len(string.split(lines[0][:-1], "\t"))
        take = range( 0, num_cols)
    else:
        num_cols = len(take)
        
    colors = None
    if color != None:
        take = [ x for x in take if x != color ]
        num_cols = len(take)
        colors = []

    if has_headers:
        headers = lines[0][:-1].split("\t")
        headers = map( lambda x: headers[x], take )        
        del lines[0]
    else:
        headers = None

    num_rows = len(lines)

    matrix = numpy.zeros( (num_rows, num_cols), numeric_type )

    nrow = 0
    if truncate:
        min_val, max_val = truncate
        
    for l in lines:
        data = l[:-1].split("\t")

        if color != None: 
            colors.append( float(data[color]) )

        try:
            data = map( lambda x: data[x], take)
        except IndexError:
            continue

        ## try conversion. Unparseable fields set to missing_value
        for x in range(len(data)):
            try:
                data[x] = float(data[x])
            except ValueError:
                data[x] = missing_value

        if truncate:
            if data[0] < min_val or data[0] > max_val: continue
            
        matrix[nrow] = data
        nrow += 1

    if truncate != None:
        # truncate matrix
        matrix = matrix[ 0:nrow+1, 0:num_cols]

    return MA.masked_values( matrix, missing_value ), headers, colors

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: plot_data.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-f", "--file", dest="filenames", type="string",
                      help="files[s] to take data from,stdin = -."  )
    parser.add_option("-l", "--legend", dest="legend", type="string",
                      help="legend for plot. If not given, the legend is build from the column titles [default=%Default]."  )
    parser.add_option("--no-legend", dest="no_legend", action="store_true",
                      help="do not output a legend [default=%Default]"  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="title for plot."  )
    parser.add_option("-p", "--hardcopy", dest="hardcopy", type="string",
                      help="filename for hardcopy.", metavar = "FILE"  )
    parser.add_option("", "--xrange", dest="xrange", type="string",
                      help="xrange."  )
    parser.add_option("", "--yrange", dest="yrange", type="string",
                      help="yrange."  )
    parser.add_option("-o", "--logscale", dest="logscale", type="string",
                      help="log-transform one or both axes [default=%Default]."  )
    parser.add_option("-x", "--xtitle", dest="xtitle", type="string",
                      help="title for x axis"  )
    parser.add_option("-y", "--ytitle", dest="ytitle", type="string",
                      help="title for y axis"  )
    parser.add_option("-n", "--normalize", dest="normalize", action = "store_true",
                      help="normalize histograms"  )
    parser.add_option("", "--legend-location", dest="legend_location", type="string",
                      help="location of legend"  )
    parser.add_option("", "--backend", dest="backend", type="string",
                      help="backend to use [Agg|SVG|PS]"  )
    parser.add_option("", "--symbols", dest="symbols", type="string",
                      help="symbols to use for each histogram [steps|...]."  )
    parser.add_option("", "--dump", dest="dump", action="store_true",
                      help="dump data for debug purposes."  )
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for plotting [default=%Default].")
    parser.add_option( "--color", "--colour", dest="color", type="int",
                      help="column with field to use for coloring [default=%Default].")
    parser.add_option("--transpose", dest="transpose", action="store_true",
                      help="transpose values (swop x and y axis).")
    parser.add_option("--truncate", dest="truncate", action="store_true",
                      help="truncate histogram before plotting to xrange.")
    parser.add_option("--layout", dest="layout", type="string",
                      help="layout.")
    parser.add_option("--scatter", dest="scatter", action="store_true",
                      help="do a scatter plot.")
    parser.add_option("--point-size", dest="point_size", type="int",
                      help="point-size.")
    parser.add_option("--add-function", dest="function", type="string",
                      help="add a function to the plot." )
    parser.add_option("--add-error-bars", dest="error_bars", type="choice",
                      choices=("interleaved", "block"),
                      help="add error bars." )
    parser.add_option( "--no-titles", dest ="has_headers", action="store_false",
                       help="no column titles given." )
    parser.add_option( "--palette", dest="palette", type="choice",
                       choices=("rainbow", "gray", "blue-white-red",
                                "autumn", "bone", "cool", "copper", "flag", "gray", "hot", "hsv", "jet", "pink", "prism",
                                "spring", "summer", "winter", "spectral",
                                "RdBu", "RdGy", "BrBG", "BuGn", "Blues", "Greens", "Reds", "Oranges", "Greys" ),
                       help="colour palette [default=%Default]")
    parser.add_option( "--reverse-palette", dest="reverse_palette", action="store_true",
                       help="reverse colour palette [default=%Default].")

    parser.set_defaults(
        filenames = "-",
        legend = None,
        title = None,
        hardcopy = None,
        logscale = None,
        xtitle = None,
        ytitle = None,
        xrange = None,
        yrange = None,        
        normalize = None,
        columns = "all",
        legend_location = "upper right",
        backend = "Agg",
        symbols = "gD,bh,r+,c+,m+,y+,ko,g^,b<,r>,cD,mh",
        dump = False,
        truncate = False,
        layout = "1,1",
        scatter = False,
        point_size=0,
        transpose = False,
        function = None,
        error_bars = None,
        has_headers = True,
        color = None,
        palette = None,
        reverse_palette = False,
        no_legend = False,
        )

    (options, args) = E.Start( parser )

    options.symbols=options.symbols.split(",")

    # import matplotlib/pylab. Has to be done here
    # for batch scripts without GUI.
    import matplotlib
    if options.hardcopy: matplotlib.use( "cairo" )
    import pylab

    plot_kwargs = {}
    
    if options.xrange: options.xrange = map(float, options.xrange.split(","))
    if options.yrange: options.yrange = map(float, options.yrange.split(","))
    if options.legend: options.legend = options.legend.split(",")  
    if options.layout: options.layout = tuple(map(int, options.layout.split(",")))
    if options.color: 
        options.color -= 1
        plot_kwargs["edgecolor"] = 'none'
    if options.columns != "all": options.columns = map(lambda x: int(x)-1, options.columns.split(","))

    if options.palette:
        if options.palette == "gray":
            color_scheme = colors_gray
        else:
            if options.reverse_palette:
                color_scheme = eval( "pylab.cm.%s_r" % options.palette)                    
            else:
                color_scheme = eval( "pylab.cm.%s" % options.palette)
        plot_kwargs["cmap"] = color_scheme
    else:
        color_scheme = None
    
    options.filenames = options.filenames.split(",")
    
    if len(args) > 0:
        options.filenames = args

    plotno = 0
    for filename in options.filenames:
        if filename == "-":
            infile = sys.stdin
        else:
            infile = open(filename)

        try:
            data, legend, colors = readTable( infile,
                                              numeric_type=numpy.float32,
                                              take=options.columns,
                                              has_headers = options.has_headers,
                                              truncate= options.xrange,
                                              color = options.color )
        except ValueError, msg:
            E.warn( "parsing error in %s: %s" % (filename, msg) )
            continue

        if filename != "-":
            infile.close()

        if not options.has_headers:
            legend = ["col%i" % i for i in range( 0, len(data[0]) ) ]

        if options.legend: legend = options.legend
        
        plotno += 1
        pylab.subplot( "".join(map(str, options.layout + (plotno,))))

        if options.title:
            pylab.title( options.title )

        if options.xtitle:
            xlab = options.xtitle      
        else:
            xlab = legend[0] 
    
        if options.ytitle:
            ylab = options.ytitle
        else:
            ylab = None

        if options.transpose:
            xlab, ylab = ylab, xlab

        if xlab: pylab.xlabel( xlab )
        if ylab: pylab.ylabel( ylab )
        
        if options.normalize:
            m = [0] * len(legend)
            for d in data:
                for x in range(1, len(legend)):
                    m[x] = max(m[x], d[x])
            for d in data:
                for x in range(1, len(legend)):
                    if m[x] > 0:
                        d[x] = d[x] / m[x]

        if options.dump:
            for d in data:
                print d

        lines = []
        
        for x in range(1, len(legend)):
            s = options.symbols[x % len(options.symbols)]

            if not options.transpose:
                first, second = 0, x
            else:
                first, second = x, 0

            yvals = data[:,second]
            xvals = MA.masked_array(data[:,first], MA.getmask(yvals))
            
            xvals = xvals.compressed()
            yvals = yvals.compressed()            

            if options.scatter:
                if options.point_size:
                    lines.append (pylab.scatter( xvals, yvals, options.point_size, c = colors, **plot_kwargs ))
                else:
                    # lines.append ( pylab.scatter( xvals, yvals, c = colors ) ) # , vmin=min(colors),vmax=max(colors), facetted=False))
                    lines.append( pylab.scatter( xvals, yvals, c = colors, vmin=min(colors),vmax=max(colors), **plot_kwargs ) )
            else:
                if options.error_bars:
                    lines.append(pylab.errorbar( xvals, yvals,
                                                 yerr = None,
                                                 fmt=s ) )
                else:
                    assert len(xvals) == len(yvals), "unequal number of values: %i != %i" % (len(xvals),len(yvals))
                    lines.append(pylab.plot( xvals, yvals, s ) )

        if not options.no_legend:
            pylab.legend( lines, legend[1:], options.legend_location)

        if options.logscale:
            if "x" in options.logscale:
                pylab.gca().set_xscale('log')
            if "y" in options.logscale:
                pylab.gca().set_yscale('log')

        if options.xrange:
            pylab.xlim( options.xrange )
            
        if options.yrange:
            pylab.ylim( options.yrange )

        if options.color != None:
            pylab.colorbar()

    if options.function:
        xstart, xend = pylab.gca().get_xlim()
        increment = (xend - xstart) / 100.0
        exec("f = lambda x: %s" % options.function )
        xvals, yvals = [], []
        for x in range(0,100):
            xvals.append( xstart )
            yvals.append( f(xstart) )
            xstart += increment
        xvals.append( xstart )
        yvals.append( f(xstart) )
        
        pylab.plot( xvals, yvals )
            
    if options.hardcopy:
        pylab.savefig( os.path.expanduser(options.hardcopy) )
    else:
        pylab.show()

    E.Stop()
