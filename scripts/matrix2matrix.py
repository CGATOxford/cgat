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
matrix2matrix.py - operate on matrices
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

   * full: full matrix with row and column headers (unless --no-headers is given.)
   * sparse: sparse matrix
   * phylip: phylip formatted matrix, but using tabs as separators and long names.

Methods:

sort-rows
   sort rows by order in --filename-rows

sort-columns
   sort columns by order in --filename-columns

mask-rows
   set rows matching ids in --filename-rows to --value

mask-columns
   set columns matching ids in --filename-columns to --value

mask-rows-and-columns
   set rows and columns matching ids in --filename-columns to --value (and)

Usage
-----

Example::

   python matrix2matrix.py --help

Type::

   python matrix2matrix.py --help

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
import math
import StringIO

import numpy
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.CorrespondenceAnalysis as CorrespondenceAnalysis
import CGAT.MatlabTools as MatlabTools
import scipy

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: matrix2matrix.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("normalize-by-min-diagonal", "normalize-by-column",
                               "log", "ln", "negzero2value", 
                               "set-diagonal", 
                               "subtract-matrix", "mix-matrix", "normalize-by-matrix", 
                               "normalize-by-column-max", "normalize-by-row-max", 
                               "normalize-by-column-min", "normalize-by-row-min", 
                               "normalize-by-column-median", "normalize-by-row-median", 
                               "normalize-by-column-mean", "normalize-by-row-mean", 
                               "normalize-by-column-total", "normalize-by-row-total", 
                               "correspondence-analysis",
                               "normalize-by-value", 
                               "add-value", 
                               "sort-rows", "sort-columns",
                               "transpose", 
                               "upper-bound", "lower-bound", 
                               "subtract-first-col", "multiply-by-value", "divide-by-value",
                               "mask-rows", "mask-columns", "mask-rows-and-columns",
                               "symmetrize-mean", "symmetrize-max", "symmetrize-min",
                               ),
                      help="""method to use [default=%default]"""  )
    
    parser.add_option("-s", "--scale", dest="scale", type="float",
                      help="factor to scale matrix by [default=%default]."  )
    
    parser.add_option("-f", "--format", dest="format", type="string",
                      help="output number format [default=%default]."  )

    parser.add_option( "--filename-rows", dest="filename_rows", type="string",
                      help="filename with rows to mask [default=%default]."  )

    parser.add_option( "--filename-columns", dest="filename_columns", type="string",
                      help="filename with columns to mask [default=%default]."  )

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="Parameters for various functions."  )

    parser.add_option("-t", "--headers", dest="headers", action="store_true",
                      help="matrix has row/column headers."  )

    parser.add_option("--no-headers", dest="headers", action="store_false",
                      help="matrix has no row/column headers."  )

    parser.add_option("-a", "--value", dest="value", type="float",
                      help="value to use for various algorithms."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("full", "sparse", "phylip" ),
                      help="""input format for matrix."""  )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("full", "sparse", "phylip" ),
                      help="""output format for matrix."""  )
    
    parser.add_option( "--missing", dest="missing", type="float",
                      help="value to use for missing values. If not set, missing values will cause the script to fail [default=%default]."  )
    
    parser.set_defaults(
        methods = [],
        scale = 1.0,
        headers = True,
        format = "%6.4f",
        output_format = "full",
        input_format = "full",        
        value = 0.0,
        parameters = "",
        write_separators = True,
        filename_rows = None,
        filename_columns = None,
        missing = None,
        )

    (options, args) = E.Start( parser )

    options.parameters = options.parameters.split(",")

    lines = filter( lambda x: x[0] != "#", sys.stdin.readlines())
    
    if len(lines) == 0: raise IOError("no input")

    chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)) )

    if not chunks:
        options.write_separators = False
        chunks = [-1]
        
    chunks.append( len(lines) )

    if options.filename_rows:
        row_names, n = IOTools.ReadList( open( options.filename_rows, "r") )
    if options.filename_columns:
        column_names, n = IOTools.ReadList( open( options.filename_columns, "r") )

    for chunk in range(len(chunks) -1 ):

        try:
            raw_matrix, row_headers, col_headers = MatlabTools.readMatrix( StringIO.StringIO("".join(lines[chunks[chunk]+1:chunks[chunk+1]])),
                                                                           format=options.input_format, 
                                                                           headers = options.headers,
                                                                           missing = options.missing )
        except ValueError, msg:
            E.warn( "matrix could not be read: %s" % msg)
            continue

        nrows, ncols = raw_matrix.shape

        E.debug("read matrix: %i x %i, %i row titles, %i colum titles" %\
                    (nrows, ncols, len(row_headers), len(col_headers)))

        parameter = 0

        for method in options.methods:

            matrix = numpy.reshape( numpy.array(raw_matrix), raw_matrix.shape )

            if method in ("normalize-by-matrix", "subtract-matrix", "mix-matrix", "add-matrix"):

                other_matrix, other_row_headers, other_col_headers = MatlabTools.ReadMatrix( open(options.parameters[parameter], "r"),
                                                                                             headers = options.headers )

                other_nrows, other_ncols = other_matrix.shape

                if options.loglevel >= 2:
                    options.stdlog.write( "# read second matrix from %s: %i x %i, %i row titles, %i colum titles.\n" %\
                                          (options.parameters[parameter],
                                          other_nrows, other_ncols, len(other_row_headers), len(other_col_headers)))

                parameter += 1


            elif method == "normalize-by-min-diagonal":
                for x in range(nrows) :
                    for y in range(ncols):
                        m = min(raw_matrix[x,x], raw_matrix[y,y])
                        if m > 0:
                            matrix[x,y] = raw_matrix[x,y] / m

            elif method == "normalize-by-column":
                if nrows != ncols:
                    raise "only supported for symmeric matrices."

                for x in range(nrows) :
                    for y in range(ncols):
                        if raw_matrix[y,y] > 0:
                            matrix[x,y] = raw_matrix[x,y] / raw_matrix[y,y]

            elif method == "normalize-by-value":
                matrix = raw_matrix / float(options.parameters[parameter])
                parameter += 1

            elif method == "normalize-by-row":
                if nrows != ncols:
                    raise "only supported for symmeric matrices."

                for x in range(nrows) :
                    for y in range(ncols):
                        if raw_matrix[y,y] > 0:
                            matrix[x,y] = raw_matrix[x,y] / raw_matrix[x,x]

            elif method == "subtract-first-col":
                for x in range(nrows) :
                    for y in range(ncols):
                        matrix[x,y] -= raw_matrix[x,0]
                        
            elif method.startswith( "normalize-by-column"):
                if method.endswith( "max" ): f = max
                elif method.endswith( "min" ): f = min
                elif method.endswith( "median" ): f = scipy.median
                elif method.endswith( "mean" ): f = scipy.mean
                elif method.endswith( "total" ): f = sum

                for y in range(ncols):
                    m = f(matrix[:,y])
                    if m != 0:
                        for x in range(nrows) :
                            matrix[x,y] = matrix[x,y] / m

            elif method.startswith( "normalize-by-row"):
                if method.endswith( "max" ): f = max
                elif method.endswith( "min" ): f = min
                elif method.endswith( "median" ): f = scipy.median
                elif method.endswith( "mean" ): f = scipy.mean
                elif method.endswith( "total" ): f = sum

                for x in range(nrows) :
                    m = f(matrix[x,:])
                    if m != 0:
                        for y in range(ncols):
                            matrix[x,y] = raw_matrix[x,y] / m

            elif method == "negzero2value":
                ## set zero/negative values to a value
                for x in range(nrows) :
                    for y in range(ncols):
                        if matrix[x,y] <= 0:
                            matrix[x,y] = options.value

            elif method == "minmax":
                ## set zero/negative values to a value
                for x in range(nrows) :
                    for y in range(ncols):
                        matrix[x,y], matrix[y,x] = \
                                     min(matrix[x,y], matrix[y,x]), \
                                     max(matrix[x,y], matrix[y,x])

            elif method == "log":
                ## apply log to all values.
                for x in range(nrows) :
                    for y in range(ncols):
                        if matrix[x,y] > 0:
                            matrix[x,y] = math.log10(matrix[x,y])

            elif method == "ln":
                for x in range(nrows) :
                    for y in range(ncols):
                        if matrix[x,y] > 0:
                            matrix[x,y] = math.log(matrix[x,y])

            elif method == "transpose":
                matrix = numpy.transpose(matrix)
                row_headers, col_headers = col_headers, row_headers
                nrows, ncols = ncols, nrows

            elif method == "mul":
                matrix = numpy.dot( matrix, numpy.transpose(matrix))
                col_headers = row_headers

            elif method == "multiply-by-value":
                matrix *= options.value

            elif method == "divide-by-value":
                matrix /= options.value

            elif method == "add-value":
                matrix += options.value

            elif method == "angle":
                ## write angles between col vectors
                v1 = numpy.sqrt( numpy.sum( numpy.power( matrix, 2), 0) )
                matrix = numpy.dot( numpy.transpose(matrix), matrix)
                row_headers = col_headers
                nrows = ncols
                for x in range(nrows):
                    for y in range(ncols):
                        matrix[x,y] /= v1[x] * v1[y]

            elif method == "euclid":
                ## convert to euclidean distance matrix
                matrix = numpy.zeros( (ncols, ncols), numpy.float )
                for c1 in range(0, ncols-1):
                    for c2 in range(c1+1, ncols):
                        for r in range(0, nrows):
                            d = raw_matrix[r][c1] - raw_matrix[r][c2]
                            matrix[c1,c2] += (d * d)
                        matrix[c2,c1] = matrix[c1,c2]
                matrix = numpy.sqrt(matrix)
                row_headers = col_headers
                nrows = ncols

            elif method.startswith("symmetrize"):
                f = method.split( "-" )[1]
                if f == "max": f = max
                elif f == "min": f = min
                elif f == "mean": f = lambda x,y: float(x+y)/2
                
                if nrows != ncols: 
                    raise ValueError("symmetrize only available for symmetric matrices")
                if row_headers != col_headers: 
                    raise ValueError("symmetrize not available for permuted matrices")
                for x in range(nrows) :
                    for y in range(ncols):
                        matrix[x,y] = matrix[y,x] = f( matrix[x,y], matrix[y,x] )
            elif method == "sub":
                matrix = options.value - matrix

            elif method in ("lower-bound", "upper-bound"):

                boundary = float(options.parameters[parameter])
                new_value = float(options.parameters[parameter+1])
                parameter += 2
                if method == "upper-bound":
                    for x in range(nrows):
                        for y in range(ncols):
                            if matrix[x,y] > boundary: matrix[x,y] = new_value
                else:
                    for x in range(nrows):
                        for y in range(ncols):
                            if matrix[x,y] < boundary: matrix[x,y] = new_value

            elif method == "subtract-matrix":
                matrix = matrix - other_matrix

            elif method == "add-matrix":
                matrix = matrix + other_matrix

            elif method == "normalize-by-matrix":

                ## set 0s to 1 in the other matrix
                for x in range(nrows):
                    for y in range(ncols):
                        if other_matrix[x,y] == 0: other_matrix[x,y] = 1.0

                matrix = matrix / other_matrix 

            elif method == "mix-matrix":
                for x in range(len(other_row_headers)-1):
                    for y in range(x+1,len(other_col_headers)):
                        matrix[x,y] = other_matrix[x,y] 

            elif method == "set-diagonal":
                value = float(options.parameters[parameter])
                for x in range(min(nrows, ncols)):
                    matrix[x,x] = value
                parameter += 1

            elif method == "transpose":
                matrix = numpy.transpose(raw_matrix)
                row_headers, col_headers = col_headers, row_headers

            elif method == "correspondence-analysis":
                row_indices, col_indices =  CorrespondenceAnalysis.GetIndices( raw_matrix )
                map_row_new2old = numpy.argsort(row_indices)
                map_col_new2old = numpy.argsort(col_indices)

                matrix, row_headers, col_headers = CorrespondenceAnalysis.GetPermutatedMatrix( raw_matrix,
                                                                                               map_row_new2old,
                                                                                               map_col_new2old,
                                                                                               row_headers = row_headers,
                                                                                               col_headers = col_headers)

            elif method == "mask-rows":
                r = set(row_names)
                for x in range(len(row_headers)):
                    if row_headers[x] in r:
                        matrix[x,:] = options.value

            elif method == "mask-columns":
                r = set(column_names)
                for x in range(len(col_headers)):
                    if col_headers[x] in r:
                        matrix[:,x] = options.value

            elif method == "mask-rows-and-columns":

                r = set(row_names)
                c = set(column_names)
                for x in range(len(row_headers)):
                    for y in range(len(col_headers)):
                        if row_headers[x] in r and col_headers[y] in c:
                            matrix[x,y] = options.value

            raw_matrix = numpy.reshape( numpy.array(matrix), matrix.shape )

        else:
            # for simple re-formatting jobs
            matrix = raw_matrix

        if options.write_separators:
            options.stdout.write( lines[chunks[chunk]] )

        MatlabTools.writeMatrix( sys.stdout, matrix,
                                 value_format = options.format,
                                 format = options.output_format,
                                 row_headers = row_headers,
                                 col_headers = col_headers )

    E.Stop()
