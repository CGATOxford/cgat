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
MatlabTools.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys, string, re
import numpy


def WriteMatrixOld( matrix, separator = "\t"):
    # this is clumsy

    lines = []
    for x in range(0, matrix.shape[0]):
        lines.append( string.join( map(str, matrix[x,]), separator))
        
    return string.join(lines, "\n" )

def WriteMatrix( matrix, outfile = sys.stdout, separator = "\t", format="%f",
                 row_headers = None, col_headers = None):
    """write matrix to stream.
    """
    if col_headers:
        outfile.write( separator + separator.join( col_headers ) + "\n" )
        
    for x in range(0, matrix.shape[0]):
        if row_headers: outfile.write(row_headers[x] + separator )
        outfile.write( string.join( map(lambda x: format % x, matrix[x,]), separator) + "\n")

def ReadMatrix( file,
                separator = "\t",
                numeric_type = numpy.float,
                take = "all",
                headers = False
                ):
    """read a matrix. There probably is a routine for this in Numpy, which
    I haven't found yet.
    """
    
    lines = filter( lambda x: x[0] != "#", file.readlines())

    row_headers, col_headers = [], []

    if headers:
        col_headers = lines[0][:-1].split("\t")[1:]
        del lines[0]

    num_rows = len(lines)
    
    if take != "all":
        num_cols = len(take)
    else:
        l = len(string.split(lines[0][:-1], "\t"))
        if headers:
            take = range( 1, l)
        else:
            take = range( 0, l)            

    num_cols = len(take)
            
    matrix = numpy.zeros( (num_rows, num_cols), numeric_type )
    
    nrow = 0
    for l in lines:
        data = l[:-1].split("\t")
        if headers: row_headers.append( data[0] )
        
        try:
            data = map( lambda x: float(data[x]), take)
        except ValueError:
            print "error parsing data", data
            raise

        matrix[nrow] = data
        nrow += 1

    return matrix, row_headers, col_headers
        

def ReadSparseMatrix( filename,
                      separator = "\t",
                      numeric_type = numpy.float,
                      is_symmetric = None):
    """read sparse matrix."""

    lines = map(lambda x: string.split(x[:-1], separator)[:3],
                filter(lambda x: x[0] != "#", open(filename,"r").readlines()))

    data = map(lambda x: (string.atoi(x[0]), string.atoi(x[1]), string.atof(x[2])), lines)

    num_rows = max( map(lambda x: x[0], data))
    num_cols = max( map(lambda x: x[1], data))   

    if (is_symmetric):
        num_rows= max( num_rows, num_cols)
        num_cols= max( num_rows, num_cols)

    matrix = numpy.zeros( (num_rows, num_cols), numeric_type )

    for row, col, weight in data:
        matrix[row-1,col-1] = weight

    if is_symmetric:
        for col, row, weight in data:
            matrix[row-1,col-1] = weight
            
    return matrix


def ReadBinarySparseMatrix( filename,
                            separator = "\t",
                            numeric_type = numpy.float,
                            is_symmetric = None):
    """read sparse matrix."""

    lines = map(lambda x: string.split(x[:-1], separator)[:2],
                filter(lambda x: x[0] != "#", open(filename,"r").readlines()))

    data = map(lambda x: (string.atoi(x[0]), string.atoi(x[1])), lines)

    num_rows = max( map(lambda x: x[0], data))
    num_cols = max( map(lambda x: x[1], data))   

    if (is_symmetric):
        num_rows= max( num_rows, num_cols)
        num_cols= max( num_rows, num_cols)
    matrix = numpy.zeros( (num_rows, num_cols), numeric_type )


    for row, col in data:
        matrix[row-1,col-1] = 1

    if is_symmetric:
        for col, row in data:
            matrix[row-1,col-1] = 1
            
    return matrix


def readMatrix( infile,
                format = "full",
                separator = "\t",
                numeric_type = numpy.float,
                take = "all",
                headers = True,
                missing = None,
                ):
    """read a matrix from file ane return a numpy matrix.

    formats accepted are:
    * full
    * sparse
    * phylip
    """
    
    row_headers, col_headers = [], []

    lines = filter( lambda x: x[0] != "#", infile.readlines())

    if len(lines) == 0: raise IOError("no input")

    if format == "full":

        if headers:
            col_headers = lines[0][:-1].split("\t")[1:]
            del lines[0]

        num_rows = len(lines)

        if take != "all":
            num_cols = len(take)
        else:
            l = len(string.split(lines[0][:-1], "\t"))
            if headers:
                take = range( 1, l)
            else:
                take = range( 0, l)            

        num_cols = len(take)

        matrix = numpy.zeros( (num_rows, num_cols), numeric_type )

        nrow = 0
        for l in lines:
            data = l[:-1].split("\t")
            if headers: row_headers.append( data[0] )

            if missing == None:
                try:
                    data = map( lambda x: float(data[x]), take)
                except ValueError, msg:
                    raise ValueError( "error %s: data=%s" % (msg, str( data )))
                except IndexError, msg:
                    raise IndexError( "error %s: data=%s" % (msg, str( data )))

            else:
                d = []
                for x in take:
                    try:
                        d.append( float(data[x]) )
                    except ValueError:
                        d.append( missing )
                    except IndexError, msg:
                        raise IndexError( "error %s: data=%s" % (msg, str( data )))

                data = d

            matrix[nrow] = data

            nrow += 1

    elif format == "phylip":
        ## read in symmetric phylip matrices
        ## note: they can wrap around
        if take != "all":
            raise "phylip matrix does not support take - only full matrices are processed."

        if not headers:
            raise "phylip matrix always has headers."
            
        num_rows = int(lines[0].strip())
        
        num_cols = num_rows

        matrix = numpy.zeros( (num_rows, num_cols), numeric_type )
        take = range( 1, num_rows)

        nrow = 0
        ncol = 0        
        for l in lines[1:]:

            data = re.split("\s+", l[:-1])
            
            if ncol == 0:
                row_headers.append( data[0] )
            
            try:
                data = map( float, data[1:len(data)] )
            except ValueError:
                raise "parsing error in conversion to float in line %s" % l

            for x in range(len(data)):
                matrix[nrow][ncol] = data[x]
                ncol += 1

            ## deal with wrapping around
            if ncol == num_cols:
                ncol = 0
                nrow += 1
        
        col_headers = row_headers
        
    return matrix, row_headers, col_headers
        
def writeMatrix( outfile, matrix,
                 format= "full",
                 separator = "\t",
                 value_format="%f",
                 row_headers = None,
                 col_headers = None):
    """write matrix to stream.
    """

    if format == "full":
        if col_headers:
            outfile.write( separator + separator.join( col_headers ) + "\n" )

        for x in range(0, matrix.shape[0]):
            if row_headers: outfile.write(row_headers[x] + separator )
            outfile.write( string.join( map(lambda x: value_format % x, matrix[x,]), separator) + "\n")
            
    elif format == "phylip":
        if not row_headers:
            raise "phylip output requires row headers."

        nrows = len(row_headers)
        outfile.write( "%i\n" % nrows )

        for x in range(0, nrows):
            outfile.write(row_headers[x] + separator )
            outfile.write( separator.join( [value_format % y for y in matrix[x,]]) + "\n" )
            
    else:
        raise "unknown output format %s" % output_format





