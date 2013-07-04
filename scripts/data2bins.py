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
data2bins.py - split table according to values in a column
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Separate data in table according to a column.

Missing values (na) are ignored.

Reads data from stdin unless infile is given.

Usage
-----

Example::

   python data2bins.py --help

Type::

   python data2bins.py --help

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
import math
import bisect
import CGAT.Experiment as E
import CGAT.CSV as CSV

class Outputter:
    def __init__(self, filename, headers = None):
        self.mFilename = filename
        self.mOutfile = open(filename,"w")
        self.mCounts = 0
        if headers:
            self.mOutfile.write( "\t".join(headers) + "\n" )

    def write(self, data ):
        self.mOutfile.write( "\t".join(map(str,data)) + "\n" )
        self.mCounts += 1

    def __del__(self):
        self.mOutfile.close()
        if self.mCounts == 0:
            os.remove( self.mFilename )

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: data2bins.py 2782 2009-09-10 11:40:29Z andreas $", 
                                    usage = globals()["__doc__"])

    parser.add_option("--column", dest="column", type="int",
                      help="column to split on."  )

    parser.add_option("--num-bins", dest="num_bins", type="int",
                      help="number of bins to create."  )

    parser.add_option("--method", dest="method", type="choice",
                      choices=("equal-sized-bins" ,),
                      help="method to use to bin data."  )

    parser.add_option("--no-headers", dest="has_headers", action="store_false",
                      help="matrix has no row/column headers."  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename with histogram information on aggregate coverages [%default].")


    parser.set_defaults(
        has_headers = True,
        method = "equal-sized-bins",
        column = 1,
        num_bins = 4,
        output_filename_pattern = "bin%i",
        )

    (options, args) = E.Start( parser )
    options.column -= 1

    if args:
        if args[0] == "-":
            infile = sys.stdin
        else:
            infile = open(args[0],"r")
    else:
        infile = sys.stdin

    fields, data = CSV.ReadTable( infile )

    c = options.column 
    values = [ float(x[c]) for x in data ]

    bins = []

    if options.method == "equal-sized-bins":
        increment = int(math.floor( float(len(values)) / options.num_bins ) )
        indices = range( 0, len(values) )
        indices.sort( key = lambda x: values[x] )
        for x in xrange( len(values)): values[indices[x]] = x
        bins = range( 0, len(values)-increment, increment )

    elif options.method == "pass":
        pass
    
    E.debug( "bins=%s" % str(bins))

    outputters = []
    for x in xrange(0,len(bins)):
        outputters.append( Outputter(options.output_filename_pattern % x, fields ) )
        
    ## output tables
    for x in xrange(0, len(data)):
        bin = bisect.bisect( bins, values[x]) - 1
        outputters[bin].write(data[x])
        
    ## stats
    if options.loglevel >= 1:
        options.stdlog.write("# bin\tstart\tcounts\tfilename\n" )
        for x in xrange(0,len(bins)):
            options.stdlog.write("# %i\t%f\t%i\t%s\n" % (x, bins[x], outputters[x].mCounts, outputters[x].mFilename) )
            
    E.info( "ninput=%i, noutput=%i" % (len(data), sum( (x.mCounts for x in outputters )) ))
    
    E.Stop()
