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
CorrespondenceAnalysis.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import numpy, numpy.linalg, numpy.linalg.linalg

##---------------------------------------------------------------------
def GetIndices( matrix ):
    """return order (1st eigenvector) of row and column indicies.

    This procedure fails if there are row or columns with a sum of 0.
    """

    nrows, ncols = matrix.shape

    # calculate row and column sums
    row_sums = numpy.sum( matrix, 1 )
    col_sums = numpy.sum( matrix, 0 )    

    # check for empty rows/columns
    # return the original permutation
    if 0 in row_sums or 0 in col_sums:
        return range(nrows), range(ncols)
    
    a = numpy.zeros( (nrows, nrows), numpy.float)
    for x in range( 0, nrows):
        a[x,x] = 1.0 / float(row_sums[x])

    b = numpy.zeros( (ncols, ncols), numpy.float)    
    for x in range( 0, ncols):
        b[x,x] = 1.0 / float(col_sums[x])
        
    M = numpy.dot( \
        a, numpy.dot( \
        matrix, numpy.dot( \
        b, numpy.transpose( matrix ))))

    try:
        row_eigenvector = numpy.linalg.eig(M)[1][:,1]
    except numpy.linalg.linalg.LinAlgError, msg:
        raise ValueError( msg )
 
    M = numpy.dot( \
        b, numpy.dot( \
            numpy.transpose(matrix), numpy.dot( \
                a, matrix )))

    try:
        col_eigenvector = numpy.linalg.eig(M)[1][:,1]
    except numpy.linalg.linalg.LinAlgError, msg:
        raise ValueError( msg )
         
    ## insert columns ignored at the computation and give them the lowest
    ## eigenvalue
    row_eigenvector = row_eigenvector.astype(numpy.float) 
    col_eigenvector = col_eigenvector.astype(numpy.float)

    return row_eigenvector, col_eigenvector

##---------------------------------------------------------------------
def GetPermutatedMatrix( matrix, 
                         map_row_new2old, map_col_new2old,
                         row_headers = None, col_headers = None):
    """return a permuted matrix. Note, that currently this is very
    inefficient, as I do not know how to do this in numpy.
    """

    nrows, ncols = matrix.shape

    result = numpy.zeros( (nrows, ncols), matrix.dtype)
    for r in range(0, nrows):
        for c in range(0,ncols):
            result[r,c] = matrix[map_row_new2old[r], map_col_new2old[c]]

    if not row_headers or not col_headers:
        return result

    rows = []
    for x in map_row_new2old:
        rows.append( row_headers[x] )
        
    cols = []
    for x in map_col_new2old:
        cols.append( col_headers[x] )

    return result, rows, cols
        
##---------------------------------------------------------------------
def PermuteRows( matrix ):
    pass
    

if __name__ == "__main__":

    num_rows = 6
    num_cols = 5
    matrix = numpy.zeros( (num_rows,num_cols), numpy.int)

    matrix[0,2] = 1
    matrix[0,3] = 1
    matrix[1,0] = 1
    matrix[1,1] = 1
    matrix[1,4] = 1
    matrix[2,1:5] = 1
    matrix[3,2:4] = 1
    matrix[4,0] = 1
    matrix[4,4] = 1
    matrix[5,0] = 1
    matrix[5,2:5] = 1


    print "matrix=", matrix


    row_indices, col_indices =  GetIndices( matrix )

    map_row_new2old = numpy.argsort(row_indices)
    map_col_new2old = numpy.argsort(col_indices)

    print map_row_new2old
    print map_col_new2old

    print GetPermutatedMatrix( matrix, map_row_new2old, map_col_new2old)

    

