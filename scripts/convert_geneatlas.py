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
convert_geneatlas.py - 
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

   python convert_geneatlas.py --help

Type::

   python convert_geneatlas.py --help

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
import tempfile
import subprocess
import optparse
import math

USAGE = """Convert gene atlas table into a convenient format.

Input:
   * the gene atlas table.
   * a map between probe sets and identifiers.

The map need not be unique.

"""

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.MatlabTools as MatlabTools
import scipy
import numpy

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: convert_geneatlas.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-m", "--filename-map", dest="filename_map", type="string",
                      help="filename map." )

    parser.add_option("-i", "--filename-info", dest="filename_info", type="string",
                      help="output filename with mapping information." )

    parser.add_option("-t", "--filename-tissues", dest="filename_tissues", type="string",
                      help="output filename with tissue list - only keep these columns." )

    parser.add_option("-a", "--aggregate", dest="filename_map", type="choice",
                      choices=("mean", "max", "min", "median" ),
                      help="aggregation function." )

    parser.set_defaults( filename_map = None,
                         filename_info = None,
                         filename_tissues = None,
                         headers = True,
                         aggregate = "mean",
                         value_format = "%5.2f",
                         method="counts")
    
    (options, args) = E.Start( parser )

    if not options.filename_map:
        raise "please supply filename mapping probesets to identifiers."
    
    map_probe2locus = IOTools.ReadMap( open(options.filename_map, "r") )

    matrix, row_headers, col_headers = MatlabTools.readMatrix( sys.stdin,
                                                               format="full", 
                                                               headers = options.headers )

    if options.filename_tissues:
        tissues, nerrors = IOTools.ReadList( open(options.filename_tissues, "r") )
        tissues = set(tissues)
        columns = []
        for x in range(len(col_headers)):
            if col_headers[x] in tissues:
                columns.append( x )
    else:
        columns = range(len(col_headers))
        
    nrows, ncols = len(row_headers), len(col_headers)
    
    ninput, noutput, nkept = 0, 0, 0
    
    no_map = []
    degenerate = []

    map_old2new = {}
    new_row_headers = []
    
    for x in range(len(row_headers)):

        ninput += 1
        try:
            new_header = map_probe2locus[row_headers[x]]
        except KeyError:
            no_map.append( row_headers[x] )
            continue

        if "," in new_header:
            degenerate.append( new_header )
            continue

        if new_header not in map_old2new:
            new_row_headers.append( new_header )
            map_old2new[new_header] = []

        map_old2new[new_header].append( x )
            
        nkept += 1
    
    ## output - aggregate values
    options.stdout.write( "locuslink\t" + "\t".join( [ col_headers[x] for x in columns ] ) + "\n" )

    if options.aggregate == "mean":
        f = numpy.mean
    elif options.aggregate == "median":
        f = numpy.median
    elif options.aggregate == "min":
        f = max
    elif options.aggregate == "max":
        f = min

    for x in range(len(new_row_headers)):

        new_values = []
        row_header = new_row_headers[x]
        old_rows = map_old2new[row_header]
        
        for y in columns:
            new_values.append( f( [ matrix[x][y] for x in old_rows ] ) )
            
        options.stdout.write("%s\t%s\n" % (row_header, "\t".join( map( lambda x: options.value_format % x, new_values) ) ) )
        noutput += 1

    if options.filename_info:
        outfile = open(options.filename_info, "w")
        outfile.write( "locuslink\tnprobesets\tprobesets\n" )

        for x in range(len(new_row_headers)):

            row_header = new_row_headers[x]
            old_rows = map_old2new[row_header]
            outfile.write( "%s\t%i\t%s\n" % (row_header, len(old_rows), "\t".join( [ row_headers[x] for x in old_rows] ) ) )

        for row_header in no_map:
            outfile.write( "%s\t0\tnomap\n" % row_header )
            
        for row_header in degenerate:
            outfile.write( "%s\t0\tdegenerate\n" % row_header )

        outfile.close()
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, nkept=%i, noutput=%i, nno_map=%i, ndegenerate=%i\n" % (len(row_headers), nkept, noutput, len(no_map), len(degenerate)))

    E.Stop()
