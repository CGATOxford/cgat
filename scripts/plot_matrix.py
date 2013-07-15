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
plot_matrix.py - plot a matrix (imshow)
=======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script plots a matrix given in tabular format.

Usage
-----

Example::

   python plot_matrix.py --help

Type::

   python plot_matrix.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import os
import string
import re
import time
import optparse


import scipy
import numpy
import CGAT.MatlabTools as MatlabTools
import CGAT.Experiment as E
from array import array
import math
import numpy

param_grid_size_x = 6

param_hardcopy = "test.ps"
param_chunk_size = 2000
param_background = 255

def GetRange( matrix, r ):
    """get ranges from a range string.
    can be two numbers separated by comma, or min/max + a number
    """
    
    rr = r.split(",")
    
    if len(rr) != 2:
        raise "please supply two values for range separated by a comma."

    vv = []
    for x in rr:
        if x == "min": v = min(matrix.flat)
        elif x == "max": v = max(matrix.flat)
        else: v = float(x)
            
        vv.append(v)
        
    return vv

def plotMatrix( matrix, color_scheme, row_headers, col_headers, vmin, vmax, options ):

    pylab.imshow(matrix,
                 cmap=color_scheme,
                 origin='lower',
                 vmax = vmax,
                 vmin = vmin,
                 interpolation='nearest')

    # offset=0: x=center,y=center
    # offset=0.5: y=top/x=right
    offset = 0.0

    if options.xticks:
        pylab.xticks( [ offset + x for x in range(len(options.xticks)) ],
                      options.xticks,
                      rotation="vertical",
                      fontsize="8" )
    else:
        if col_headers and len(col_headers) < 100:
            pylab.xticks( [ offset + x for x in range(len(col_headers)) ],
                          col_headers,
                          rotation="vertical",
                          fontsize="8" )

    if options.yticks:
        pylab.yticks( [ offset + y for y in range(len(options.yticks)) ],
                      options.yticks,
                      fontsize="8" )
    else:
        if row_headers and len(row_headers) < 100:
            pylab.yticks( [ offset + y for y in range(len(row_headers)) ],
                          row_headers,
                          fontsize="8" )


if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: plot_matrix.py 2782 2009-09-10 11:40:29Z andreas $")
    
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take from table." )

    parser.add_option("-a", "--hardcopy", dest="hardcopy", type="string",
                      help="write hardcopy to file.", metavar = "FILE" )

    parser.add_option("-f", "--file", dest="input_filename", type="string",
                      help="filename with table data.",
                      metavar = "FILE")

    parser.add_option("-p", "--plot", dest="plot", type="string",
                      help="plots to plot.",
                      action = "append")

    parser.add_option("-t", "--threshold", dest="threshold", type="float",
                      help="min threshold to use for counting method.")

    parser.add_option("-o", "--colours", dest="colours", type="int",
                      help="column with colour information.")

    parser.add_option("-l", "--labels", dest="labels", type="string",
                      help="column labels for x and y in matched plots.")

    parser.add_option("-e", "--headers", dest="headers", action="store_true",
                      help="headers are supplied in matrix.")

    parser.add_option( "--no-headers", dest="headers", action="store_false",
                      help="headers are not supplied in matrix.")

    parser.add_option( "--normalize", dest="normalize", action="store_true",
                      help="normalize matrix.")

    parser.add_option( "--palette", dest="palette", type="choice",
                       choices=("rainbow", "gray", "blue-white-red",
                                "autumn", "bone", "cool", "copper", "flag", "gray", "hot", "hsv", "jet", "pink", "prism",
                                "spring", "summer", "winter", "spectral",
                                "RdBu", "RdGy", "BrBG", "BuGn", "Blues", "Greens", "Reds", "Oranges", "Greys" ),
                       help="colour palette [default=%Default]")

    parser.add_option( "--reverse-palette", dest="reverse_palette", action="store_true",
                      help="reverse the palette [default=%default].")

    parser.add_option("", "--xrange", dest="xrange", type="string",
                      help="xrange."  )

    parser.add_option("", "--yrange", dest="yrange", type="string",
                      help="yrange."  )

    parser.add_option("", "--zrange", dest="zrange", type="string",
                      help="zrange."  )

    parser.add_option("", "--xticks", dest="xticks", type="string",
                      help="xticks."  )
    
    parser.add_option("", "--yticks", dest="yticks", type="string",
                      help="yticks."  )

    parser.add_option("--bar-format", dest="bar_format", type="string",
                      help="format for ticks on colourbar."  )

    parser.add_option("--title", dest="title", type="string",
                      help="title to use."  )
    
    parser.add_option("--missing", dest="missing", type="float",
                      help="value to use for missing data."  )

    parser.add_option("--subplots", dest="subplots", type="string",
                      help="split matrix into several subplots. Supply number of rows and columns separated by a comma."  )
    
    parser.set_defaults( \
        hardcopy = None,
        input_filename = "-",
        columns = "all",
        statistics = [],
        plot=[],
        threshold=0.0,
        labels = "x,y",
        colours= None,
        xrange = None,
        yrange = None,
        zrange = None,
        palette = None,
        reverse_palette = False,
        xticks = None,
        yticks = None,
        normalize = False,
        bar_format = "%1.1f",
        headers = True,
        missing = None,
        title = None,
        subplots = None)

    (options, args) = E.Start( parser )

    # import matplotlib/pylab. Has to be done here
    # for batch scripts without GUI.
    import matplotlib
    if options.hardcopy: matplotlib.use( "cairo" )
    import pylab

    if len(args) > 0:
        options.input_filename = ",".join(args)

    if options.xticks: options.xticks = options.xticks.split(",")
    if options.yticks: options.yticks = options.yticks.split(",")

    if options.xrange: options.xrange = map(float, options.xrange.split(","))
    if options.yrange: options.yrange = map(float, options.yrange.split(","))

    if options.columns != "all":
        options.columns = map( lambda x: int(x)-1, options.columns.split(","))

    filenames = options.input_filename.split(",")

    if len(filenames) > 1:
        nsubrows = (len(filenames) / 3) + 1
        nsubcols = 3
    elif options.subplots:
        nsubrows, nsubcols = [ int(x) for x in options.subplots.split( ",") ]
    else:
        nsubrows, nsubcols = 1, 1

    nsubplots = nsubrows * nsubcols

    ## Setting up color maps
    if options.palette:
        if options.palette == "gray":
            _gray_data =  {'red':   ((0., 1, 1), (1., 0, 0)),
                           'green': ((0., 1, 1), (1., 0, 0)),
                           'blue':  ((0., 1, 1), (1., 0, 0))}

            LUTSIZE = pylab.rcParams['image.lut']
            colors_gray   = matplotlib.colors.LinearSegmentedColormap('gray',   _gray_data, LUTSIZE)

    plot_id = 0
    for filename in filenames:

        plot_id += 1
        pylab.subplot( nsubrows, nsubcols, plot_id)
        
        if filename == "-":
            infile = sys.stdin
        else:
            infile = open(filename, "r")
            
        matrix,row_headers,col_headers = MatlabTools.readMatrix( infile,
                                                                 numeric_type=numpy.float32,
                                                                 take=options.columns,
                                                                 headers = options.headers,
                                                                 missing = options.missing )

        if min(matrix.flat) == max(matrix.flat):
            options.stderr.write( "matrix is uniform - no plotting done.\n")
            sys.exit(0)

        if options.normalize:
            v = max(matrix.flat)
            matrix = matrix / v
        
        if options.zrange: options.zrange = GetRange( matrix, options.zrange )

        nrows, ncols = matrix.shape
            
        if options.palette:
            if options.palette == "gray":
                color_scheme = colors_gray
            else:
                if options.reverse_palette:
                    color_scheme = eval( "pylab.cm.%s_r" % options.palette)                    
                else:
                    color_scheme = eval( "pylab.cm.%s" % options.palette)
        else:
            color_scheme = None

        if options.zrange:
            vmin, vmax = options.zrange
            for x in range(nrows):
                for y in range(ncols):
                    if matrix[x,y] < options.zrange[0]:
                        matrix[x,y] = options.zrange[0]
                    if matrix[x,y] > options.zrange[1]:
                        matrix[x,y] = options.zrange[1]
        else:
            vmin, vmax = None, None
        
        if options.subplots:

            if nsubcols > 1:
                increment_x = int(float(nrows+1) / nsubcols)
                increment_y = nrows
                
                x = 0
                y = 0
                for n in range(nsubplots):
                    pylab.subplot( nsubrows, nsubcols, plot_id)
                    plot_id +=1 

                    print n, "rows=",nsubrows, "cols=", nsubcols, y, y+ increment_y, x, x + increment_x
                    print matrix[y:y+increment_y,x:x+increment_x].shape
                    print matrix.shape
                    plotMatrix( matrix[y:y+increment_y,x:x+increment_x], 
                                color_scheme, 
                                row_headers[y:y+increment_y], 
                                col_headers[x:x+increment_x], 
                                0, 100, options )
                
                x += increment_x

            elif nsubrows > 1:
                increment_x = int(float(ncols+1) / nsubrows)
                
                x = 0
                for n in range(nsubplots):
                    pylab.subplot( nsubrows, nsubcols, plot_id)
                    plot_id +=1 
                    plotMatrix( matrix[0:nrows,x:x+increment_x], 
                                color_scheme, 
                                row_headers, 
                                col_headers[x:x+increment_x], 
                                vmin, vmax, options )
                
                    x += increment_x
        else:
            plotMatrix( matrix, color_scheme, row_headers, col_headers, vmin, vmax, options )

        if options.xrange:
            pylab.xlim( options.xrange )
            
        if options.yrange:
            pylab.ylim( options.yrange )

        if options.labels:
            xlabel, ylabel = options.labels.split(",")
            pylab.xlabel(xlabel)
            pylab.ylabel(ylabel)

        if not options.subplots:
            pylab.colorbar( format = options.bar_format)

        if options.title == None or options.title != "":
            pylab.title( filename )
        
    if options.hardcopy:
        pylab.savefig( os.path.expanduser(options.hardcopy) )
    else:
        pylab.show()

    E.Stop()
