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
append_histogram.py - append columns to histogram
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts appends additional columns like cumulative sum, normalized density,
etc, to an histogram. A histogram is a tab-separated table with the
bins in the first column and the counts in the second.

Usage
-----

Example::

   python append_histogram.py --help

Type::

   python append_histogram.py --help

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


import CGAT.Experiment as Experiment
import CGAT.Histogram as Histogram

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: append_histogram.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--is-int", dest="is_ints", action="store_true",
                      help="categories are integers."  )

    parser.set_defaults(
        is_ints = False
        )

    (options, args) = Experiment.Start( parser )

    vals = []
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    # check if first line contains a header
    d = string.split(lines[0][:-1], "\t")[0]
    try:
        if options.is_ints:
            value = int(d)
        else:
            value = float( d )
    except ValueError:
        print string.join( (d, "counts", "frequency",
                            "cumulative counts", "increasing cumulative frequency",
                            "cumulative counts", "decreasing cumulative frequency"), "\t")
        del lines[0]
        
    data = map( lambda x: map(float, string.split(x[:-1], "\t")), lines)

    if len(data) == 0:
        raise "No data found."
        
    total = float(reduce( lambda x,y: x+y, map( lambda x: x[1], data)))

    cumul_down = int(total)
    cumul_up = 0

    if options.is_ints:
        form = "%i\t%i\t%6.4f\t%i\t%6.4f\t%i\t%6.4f"         
    else:
        form = "%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f" 
        
    for bin,val in data:
        percent     = float(val) / total
        cumul_up   += val
        percent_cumul_up = float(cumul_up) / total
        percent_cumul_down = float(cumul_down) / total        
        
        print form % \
              (bin, val, percent, cumul_up, percent_cumul_up, cumul_down, percent_cumul_down)

        cumul_down -= val

    Experiment.Stop()









