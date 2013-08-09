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
histograms2kl.py - compute Kullback-Leibler distance
====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Output the Kullback-Leibler distance between two or
more distributions given as histograms.

Usage
-----

Example::

   python histograms2kl.py --help

Type::

   python histograms2kl.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import time
import random

import math
import numpy
import CGAT.Experiment as E
import pgdb
import CGAT.IOTools as IOTools

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: compare_histograms.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       help="method to use [kl=kullback-leibler]",
                       choices=("kl",) )
    parser.add_option( "-n", "--no-normalize", dest="normalize", action="store_false",
                       help="do not normalize data" )
    parser.add_option( "-p", "--pseudocounts", dest="pseudocounts", type="int",
                       help="pseudocounts to add." )
    parser.add_option( "-f", "--number-format", dest="number_format", type="string",
                       help="number format." )

    parser.set_defaults(
        method = "kl",
        columns = "all",
        headers = True,
        xrange = None,
        pseudocounts = 1,
        normalize = True,
        number_format = "%6.4f"
        )
    
    (options, args) = E.Start( parser,
                                        add_pipe_options = True )

    if options.xrange: options.xrange = map(float, options.xrange.split(","))

    data, legend = IOTools.readTable( sys.stdin,
                                      numeric_type=numpy.float32,
                                      take=options.columns,
                                      headers = options.headers,
                                      truncate= options.xrange )

    nrows, ncols = data.shape

    ## first: normalize rows
    for y in range(1, ncols):
        for x in range(nrows): data[x,y] = data[x,y] + float(options.pseudocounts)
        if options.normalize:
            t = numpy.sum( data[:,y] )            
            for x in range(nrows): data[x, y] = data[x,y] / t

    for x in range(1, len(legend) - 1):
        for y in range(x+1, len(legend)):

            if options.method == "kl":
                d1 = 0.0
                d2 = 0.0
                for bin in range(nrows):
                    p = data[bin,x]
                    q = data[bin,y]
                    d1 += p * math.log( p / q)
                    d2 += q * math.log( q / p)
                    
                options.stdout.write( "%s\t%s\t%s\n" %
                                      (legend[x], legend[y],
                                       options.number_format % d1) )
                options.stdout.write( "%s\t%s\t%s\n" %
                                      (legend[y], legend[x],
                                       options.number_format % d2) )
                            
    E.Stop()
