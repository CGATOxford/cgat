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
plot_histogram.py - create line plot of tabular data
====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Reads data in tabular format with observations in rows in
different samples in columns. The first column provides the
bins.

Missing values (na) are ignored.

Reads data from stdin unless infile is given.

Usage
-----

Example::

   python plot_histogram.py --help

Type::

   python plot_histogram.py --help

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
import array
import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import CGAT.MatlabTools as MatlabTools
import numpy
import CGAT.IOTools as IOTools
import numpy.ma

def main():

    parser = E.OptionParser( version = "%prog version: $Id: plot_histogram.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-l", "--legend", dest="legend", type="string",
                      help="legend for plot [default=%default]."  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="title for plot [default=%default]."  )
    parser.add_option("-p", "--hardcopy", dest="hardcopy", type="string",
                      help="filename for hardcopy of plot. The extension defines the format. Known extensions are: 'emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg, svgz' [default=%default].", metavar = "FILE"  )
    parser.add_option("", "--xrange", dest="xrange", type="string",
                      help="x viewing range of plot [default=%default]."  )
    parser.add_option("", "--yrange", dest="yrange", type="string",
                      help="y viewing range of plot[default=%default]."  )
    parser.add_option("-o", "--logscale", dest="logscale", type="string",
                      help="use logscale on x, y or xy [default=%default]"  )
    parser.add_option("-x", "--xtitle", dest="xtitle", type="string",
                      help="title for x axis [default=%default]"  )
    parser.add_option("-y", "--ytitle", dest="ytitle", type="string",
                      help="title for y axis [default=%default]"  )
    parser.add_option("-d", "--dpi", dest="dpi", type="int",
                      help="dpi of images [default=%default]"  )
    parser.add_option("-n", "--normalize", dest="normalize", action = "store_true",
                      help="normalize histograms [default=%default]"  )
    parser.add_option( "--cumulate", dest="cumulate", action = "store_true",
                      help="calculate cumulative histogram [default=%default]."  )
    parser.add_option( "--reverse-cumulate", dest="reverse_cumulate", action = "store_true",
                      help="calculate cumulative histogram in reverse order [default=%default]."  )
    parser.add_option( "--legend-location", dest="legend_location", type="choice",
                       choices=("upper left", "upper right", "lower left", "lower right", "center", "center right", "center left", "none"),
                       help="location of legend [default=%default]"  )
    parser.add_option( "--backend", dest="backend", type="string",
                      help="backend to use [Agg|SVG|PS] [default=%default]"  )
    parser.add_option( "--symbols", dest="symbols", type="string",
                      help="symbols to use for each histogram [steps|...] [default=%default]."  )
    parser.add_option( "--dump", dest="dump", action="store_true",
                      help="dump data for debug purposes [default=%default]."  )
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for plotting [default=%default].")
    parser.add_option("--truncate", dest="truncate", action="store_true",
                      help="truncate date within x range. If not set, xrange is simply a viewing range [default=%default].")
    parser.add_option("--as-lines", dest="as_lines", action="store_true",
                      help="plot only lines, no symbols [default=%default].")
    parser.add_option("--noheaders", dest="headers", action="store_false",
                      help="do not take first input line as header [default=%default].")
    parser.add_option("--stacked", dest="stacked", action="store_true",
                      help="do a stacked plot [default=%default].")
    parser.add_option("--add-function", dest="function", type="string",
                      help="add a function to the plot [default=%default]." )
    parser.add_option("--add-error-bars", dest="error_bars", type="choice",
                      choices=("interleaved", "blocked"),
                      help="add error bars. The input format is 'interleaved' or 'blocked'. In the interleaved format the error follows each column. I the blocked format first the data, then the errors in the same order [default=%default]." )

    parser.set_defaults(
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
        headers = True,
        legend_location = "upper right",
        backend = "cairo",
        symbols = "g-D,b-h,r-+,c-+,m-+,y-+,k-o,g-^,b-<,r->,c-D,m-h",
        dump = False,
        truncate = False,
        cumulate = False,
        reverse_cumulate = False,
        function = None,
        add_error_bars = None,
        as_lines = False,
        stacked = False,
        dpi = 80,
        )

    (options, args) = E.Start( parser )

    # import matplotlib/pylab. Has to be done here
    # for batch scripts without GUI.
    import matplotlib
    if options.hardcopy: matplotlib.use( "cairo" )
    import pylab

    # put this method here (because it requires pylab)
    def doStackedPlot( data, legend ):

        colors = ["red",
                  "blue",
                  "green",
                  "cyan",
                  "magenta",
                  "yellow",
                  "brown",
                  "silver",
                  "purple",
                  "lightyellow",
                  "black",
                  "ivory",
                  "pink",
                  "orange",
                  "gray",
                  "teal"]

        ax = data[:,0]
        xvals = numpy.concatenate((ax,ax[::-1]))
        y_top = numpy.zeros(len(ax))

        min_y = min(data[:,1:].flat)
        max_y = min_y
        new_legend, dummy_lines = [], []

        for i in range(1,len(legend)):
            new_y_top = y_top + data[:,i]
            yvals = numpy.concatenate((new_y_top,y_top[::-1]))
            p = pylab.fill(xvals,
                           yvals,
                           colors[i % len(colors)])

            y_top = new_y_top
            max_y = max(y_top)

            dummy_lines.append(pylab.plot( xvals,
                                           yvals,
                                           colors[i % len(colors)] ) )

            new_legend.append( legend[i] )

        if not options.xrange:
            options.xrange=min(data[:,0]),max(data[:,0])

        if not options.yrange:
            options.yrange=0,max_y

        return dummy_lines, new_legend


    if options.as_lines:
        options.symbols = []
        for y in ("-",":","--"):
            for x in "gbrcmyk":
                options.symbols.append( y+x )
    else:
        options.symbols = options.symbols.split(",")
    
    if options.xrange: options.xrange = map(float, options.xrange.split(","))
    if options.yrange: options.yrange = map(float, options.yrange.split(","))    

    # Added support for (inclusive) range format: "1,3,5,7-100"  (Gerton 13/12/06)
    if options.columns != "all":
        cols = []
        for d in options.columns.split(','):
            colopts = d.split('-')
            if len(colopts)==2: cols += range( int(colopts[0]),int(colopts[1])+1 )
            else: cols += [int(d)-1]
        options.columns = cols

    if args:
        if args[0] == "-":
            infile = sys.stdin
        else:
            infile = open(args[0],"r")
    else:
        infile = sys.stdin
        
    if options.truncate:
        xr = options.xrange
    else:
        xr = None

    data, legend = IOTools.readTable( infile,
                                      numeric_type=numpy.float,
                                      take=options.columns,
                                      headers = options.headers,
                                      truncate= xr )
    
    if infile != sys.stdin: infile.close()
    if len(data) == 0: # or data == None:
        E.info("empty table: no plot" )
        E.Stop()
        return

    nrows, ncols = data.shape

    ## note: because of MA, iteration makes copy of slices
    ## Solution: inplace edits.
    if options.cumulate:
        if options.add_error_bars:
            raise "can not add error bars to cumulative histogram."
        if data.mask.any():
            # cumsum does not work with masked arrays, so do it manually
            for y in range(1, ncols):
                c = 0
                for x in range(0, nrows):
                    if not data.mask[x,y]:
                        data[x,y] += c
                        c = data[x,y]
        else:
            for x in range(1,ncols):
                data[:,x] = data[:,x].cumsum()

    elif options.reverse_cumulate:
        if options.add_error_bars:
            raise "can not add error bars to cumulative histogram."
        if data.mask.any():
            l = [0] * ncols
            for x in range(nrows-1,-1,-1):
                for y in range(1, ncols):
                    if not data.mask[x,y]:
                        data[x,y] += l[y]
                        l[y] = data[x,y]
        else:
            l = [0] * ncols
            for x in range(nrows-1,-1,-1):
                for y in range(1, ncols):
                    data[x,y] += l[y]
                    l[y] = data[x,y]

    if options.normalize:
        if options.add_error_bars:
            raise "can not add error bars to normalized histogram."
        if data.mask.any():
            m = [0] * ncols
            for x in range(nrows):
                for y in range(1, ncols):
                    if not data.mask[x,y]:
                        m[y] = max( m[y], float(data[x,y]) )

            for y in range(1, ncols):
                if m[y] == 0: m[y] = 1.0

            for x in range(nrows):
                for y in range(1, ncols):
                    data[x,y] = data[x,y] / m[y]
        else:
            for x in range(1,ncols):
                m = float(data[:,x].max())
                data[:,x] /= m

    if options.legend:
        legend = options.legend.split(",")
        
    if options.dump:
        for d in data:
            print d

    if options.title:
        pylab.title( options.title )

    if options.xtitle:
        pylab.xlabel( options.xtitle )
    else:
        pylab.xlabel( legend[0] )

    if options.ytitle:
        pylab.ylabel( options.ytitle )
        
    lines = []
    # use dummy_lines to workaround a bug in errorbars that
    # causes the line styles to be set incorrectly.
    dummy_lines = []
    new_legend = []

    if options.error_bars:
        if options.error_bars == "interleaved":
            step_size = 2
            max_size = len(legend)
        elif options.error_bars == "blocked":
            step_size = 1
            max_size = (len(legend) -1 ) / 2
    else:
        step_size = 1
        max_size= len(legend)

    if options.stacked:
        dummy_lines, new_legend = doStackedPlot( data, legend )
    else:
        nplotted = 0
        nskipped = 0
        for x in range(1, max_size, step_size):

            s = options.symbols[nplotted % len(options.symbols)]

            yvals = data[:,x]

            xvals = numpy.ma.masked_array(data[:,0], numpy.ma.getmask(yvals))

            xvals = xvals.compressed()
            yvals = yvals.compressed()            

            if len(xvals) == 0:
                E.warn( "skipped empty column %i: %s" % (x,legend[x] ) )

            if options.error_bars == "interleaved":
                yerr = data[:,x+1]
                yerr = yerr.compressed()
            else:
                yerr = None

            lines.append(pylab.errorbar( xvals,
                                         yvals,
                                         yerr = yerr,
                                         fmt = s ) )

            dummy_lines.append(pylab.plot( xvals,
                                           yvals,
                                           s ) )

            new_legend.append( legend[x] )

            nplotted += 1

        E.info("nplotted=%i, nskipped=%i" % (nplotted, nskipped))

    if len(lines) == 0:
        E.Stop()
        return 

    if options.legend_location != "none":
        pylab.figlegend( dummy_lines,
                         new_legend,
                         options.legend_location )

    if options.logscale:
        if "x" in options.logscale:
            pylab.gca().set_xscale('log')
        if "y" in options.logscale:
            pylab.gca().set_yscale('log')

    if options.xrange:
        pylab.xlim( options.xrange )
        
    if options.yrange:
        pylab.ylim( options.yrange )

    if options.function:
        xstart, xend = pylab.gca().get_xlim()
        increment = (xend - xstart) / 100.0
        exec("f = lambda x: %s" % options.function ) in locals()
        xvals, yvals = [], []
        for x in range(0,100):
            xvals.append( xstart )
            yvals.append( f(xstart) )
            xstart += increment
        xvals.append( xstart )
        yvals.append( f(xstart) )
        
        pylab.plot( xvals, yvals )

    if options.hardcopy:
        pylab.savefig( os.path.expanduser(options.hardcopy), dpi = options.dpi )
    else:
        pylab.show()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())
