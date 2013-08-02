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
calculate_histogram_2D.py - calculate 2D histogram
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes a 2D histogram. The data is a tab-separated
with ``x``, ``y`` in the first two columns. The columns can be
changed by the ``-1`` and ``-2`` options.

Usage
-----

Example::

   python calculate_histogram_2D.py --help

Type::

   python calculate_histogram_2D.py --help

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
import CGAT.Histogram2D as Histogram2D

USAGE = """python calculate_histogram.py < stdin > stdout

-1, --column1    column to take [default = 0]
-2, --column2    column to take [default = 1]
# at start of line is a comment
"""

param_column1 = 0
param_column2 = 1
param_bin_size1 = 1
param_bin_size2 = 1
param_titles = True

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: calculate_histogram_2D.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-t", "--titles", dest="titles", action="store_true",
                      help="input data has title in first row [default=%default]." )

    parser.add_option( "--no-titles", dest="titles", action="store_false",
                      help="input data has no title in first row [default=%default]." )

    parser.add_option("-1", "--column1", dest="column1", type="int",
                      help="first column to use [default=%default]." )
 
    parser.add_option("-2", "--column2", dest="column2", type="int",
                      help="second column to use [default=%default]." )
 
    parser.add_option("--bin-size1", dest="bin_size1", type="float",
                      help="bin size for first column [default=%default]." )

    parser.add_option("--bin-size2", dest="bin_size2", type="float",
                      help="bin size for second column [default=%default]." )

    parser.set_defaults( 
        column1 = 1,
        column2 = 2,
        bin_size1 = 1.0,
        bin_size2 = 1.0,
        titles = True )


    (options, args) = E.Start( parser )                                  
    options.column1 -= 1
    options.column2 -= 1

    histograms = []

    vals = []
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    if options.titles:
        data = string.split(lines[0][:-1], "\t")        
        print "\t".join((data[options.column1], data[options.column2], "counts"))
        del lines[0]

    ninput, noutput, nskipped = 0, 0, 0

    for l in lines:
        ninput += 1
        data = string.split(l[:-1], "\t")

        try:
            val = map( string.atof, (data[options.column1],data[options.column2]))
        except IndexError:
            nskipped += 1
            continue
        except ValueError:
            nskipped += 1
            continue

        vals.append( val )
        noutput += 1

    lines = None
    
    h = Histogram2D.Calculate( vals, bin_function=lambda x:(int(x[0] / options.bin_size1), int(x[1] / options.bin_size2)) )

    Histogram2D.Print(h, bin_function=lambda x: (x[0] * options.bin_size1, x[1] * options.bin_size2, x[2]))

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput,noutput,nskipped) )

    E.Stop()

