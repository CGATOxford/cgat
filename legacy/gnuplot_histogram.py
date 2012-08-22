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
gnuplot_histogram.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gnuplot_histogram.py --help

Type::

   python gnuplot_histogram.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, re, string, os, getopt, time, optparse

USAGE = """python %s < stdin > stdout

plot a histogram (or a set of histograms).

'#' at start of line is a comment

""" % sys.argv[0]

import Experiment
import Histogram
import Gnuplot

parser = optparse.OptionParser( version = "%prog version: $Id: gnuplot_histogram.py 2782 2009-09-10 11:40:29Z andreas $")


if __name__ == "__main__":

    parser.add_option("-l", "--legend", dest="legend", type="string",
                      help="legend for plot."  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="title for plot."  )
    parser.add_option("-c", "--hardcopy", dest="hardcopy", type="string",
                      help="filename for hardcopy.", metavar = "FILE"  )
    parser.add_option("-w", "--with", dest="with", type="string",
                      help="gnuplot: with statemenat."  )
    parser.add_option("-r", "--xrange", dest="xrange", type="string",
                      help="gnuplot: xrange."  )
    parser.add_option("-o", "--logscale", dest="logscale", type="string",
                      help="gnuplot: logscale"  )
    parser.add_option("-x", "--xtitle", dest="xtitle", type="string",
                      help="gnuplot: title for x axis"  )
    parser.add_option("-y", "--ytitle", dest="ytitle", type="string",
                      help="gnuplot: title for y axis"  )
    parser.add_option("-n", "--normalize", dest="normalize", action = "store_true",
                      help="normalize histograms"  )

    parser.set_defaults(
        legend = None,
        title = None,
        hardcopy = None,
        terminal = "postscript",
        xwith = "linespoints",
        logscale = None,
        xtitle = None,
        ytitle = None,
        xrange = None,
        normalize = None)

    (options, args) = Experiment.Start( parser )
    
    if options.hardcopy:
        if re.search( "\.png$", options.hardcopy):
            options.terminal = "png"

    if options.xrange:
        options.xrange = map(float, options.xrange.split(","))

    if len(args) > 0:
        print USAGE, "no arguments needed."
        sys.exit(1)

    if not options.hardcopy:
        Gnuplot.GnuplotOpts.gnuplot_command = 'gnuplot -persist'

    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    if options.legend:
        legend = options.legend.split(",")
    else:
        legend = string.split(lines[0][:-1], "\t")
        del lines[0]
    
    data = []
    for l in lines:
        try:
            v = map(string.atof, string.split(l[:-1], "\t"))
        except ValueError:
            continue
        data.append(v)

    if options.normalize:
        m = [0] * len(legend)
        for d in data:
            for x in range(1, len(legend)):
                m[x] = max(m[x], d[x])
        for d in data:
            for x in range(1, len(legend)):
                d[x] = d[x] / m[x]
        
    g = Gnuplot.Gnuplot(debug=1)
    g.clear()

    if options.title:
        g.title( options.title )

    g.xlabel( legend[0] )
    g.ylabel( legend[1] )

    if options.xrange:
        g("set xrange [%f:%f]" % (options.xrange[0], options.xrange[1]) )
    
    if options.logscale:
        g("set logscale %s" % options.logscale)
        
    for x in range(1, len(legend)):
        g.replot( Gnuplot.Data( data, cols=(0, x),
                                xwith=options.xwith,
                                title=legend[x]) )

    g.replot()

    if options.hardcopy:
        g.hardcopy( options.hardcopy,
                    terminal = options.terminal )
        
    Experiment.Stop()
