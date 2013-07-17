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
sparse2full.py - convert a sparse matrix to adjacency matrix
============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a sparse matrix to adjacency matrix and vice versa.

Usage
-----

Example::

   python sparse2full.py --help

Type::

   python sparse2full.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import optparse
import string
import re
import random

import CGAT.Experiment as E

def CountElements( matrix, default_value ):
    """count elements that are of default value and those that aren't."""

    ndefault = 0
    for x in matrix:
        for y in x:
            if default_value == y:
                ndefault += 1
                
    nfound = reduce( lambda x,y : x + y, [ len(x) for x in matrix] ) - ndefault
    return ndefault, nfound

def Sparse2Matrix( outfile, matrix_id, lines, options, in_map_token2row = {}, in_map_token2col = {}):

    # remove comments
    lines = filter(lambda x:x[0] != "#" and len(x[:-1]) > 0, lines)

    if len(lines) == 0: raise IOError("no input")

    # forget about titles
    v = lines[0][:-1].split("\t")[2]
    if v not in ("na", "NaN"):
        try:
            v = float(v)
        except ValueError:
            del lines[0]

    if in_map_token2row:
        map_token2row = in_map_token2row
    else:
        map_token2row = {}

    if in_map_token2col:
        map_token2col = in_map_token2col
    else:
        map_token2col = {}

    row_converter, col_converter = str, str
    if options.format == "string":

        has_row_names = len(map_token2row) > 0
        has_col_names = len(map_token2col) > 0

        ## if either row/column names are not given:
        if not map_token2row or not map_token2col:
            row_tokens = map(lambda x: string.split(x[:-1], "\t")[0], lines )
            col_tokens = map(lambda x: string.split(x[:-1], "\t")[1], lines )

            if options.input_format == "row-col-weight-weight":
                # merge row and col tokens
                row_tokens.extend(col_tokens)
                col_tokens = row_tokens
            
            if options.is_numeric:
                try:
                    row_tokens = map(int, row_tokens)
                    row_converter = int
                except ValueError:
                    row_tokens = map(float, row_tokens)
                    row_converter = float

                try:
                    col_tokens = map(int, col_tokens)                    
                    col_converter = int
                except ValueError:
                    col_tokens = map(float, col_tokens)                    
                    col_converter = float

            row_tokens.sort()
            col_tokens.sort()

            if not has_row_names:
                for row_token in row_tokens:
                    if row_token not in map_token2row:
                        map_token2row[row_token] = len(map_token2row)
            if not has_col_names:
                for col_token in col_tokens:
                    if col_token not in map_token2col:
                        map_token2col[col_token] = len(map_token2col)

        if not options.asymmetric:
            for col_token in map_token2col.keys():
                if col_token not in map_token2row:
                    map_token2row[col_token] = len(map_token2row)            
            map_token2col = map_token2row

            

        matrix = [ [ options.default for j in range(len(map_token2col))] for i in range(len(map_token2row)) ]

        if len(map_token2col) == len(map_token2row):
            for j in range(len(map_token2col)):
                matrix[j][j] = options.default_diagonal

        ## replicate counts
        replicates = [ [ 0 for j in range(len(map_token2col))] for i in range(len(map_token2row)) ]
        
        if options.input_format == "row-col-weight":
            for line in lines:
                row_token, col_token, weight = string.split(line[:-1], "\t")[:3]
                row_token, col_token = row_converter( row_token ), col_converter(col_token)
                matrix[map_token2row[row_token]][map_token2col[col_token]] = weight
                replicates[map_token2row[row_token]][map_token2col[col_token]] += 1                
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][map_token2row[row_token]] = weight
                    replicates[map_token2col[col_token]][map_token2row[row_token]] += 1
                    
        elif options.input_format == "row-col-weight-replicates":
            for line in lines:
                row_token, col_token, weight, n = string.split(line[:-1], "\t")[:4]
                matrix[map_token2row[row_token]][map_token2col[col_token]] = weight
                replicates[map_token2row[row_token]][map_token2col[col_token]] = int(n)
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][map_token2row[row_token]] = weight
                    replicates[map_token2col[col_token]][map_token2row[row_token]] = int(n)                     

        elif options.input_format == "row-col-weight-weight":
            for line in lines:
                row_token, col_token, weight1, weight2 = string.split(line[:-1], "\t")[:4]
                matrix[map_token2row[row_token]][map_token2col[col_token]] = weight1
                matrix[map_token2col[col_token]][map_token2row[row_token]] = weight2
                replicates[map_token2row[row_token]][map_token2col[col_token]] += 1                
                replicates[map_token2row[col_token]][map_token2col[row_token]] += 1                
                    
        col_tokens = map_token2col.items()
        col_tokens.sort( lambda x,y: cmp(x[1],y[1]))
        row_tokens = map_token2row.items()
        row_tokens.sort( lambda x,y: cmp(x[1],y[1]))

        ndefault, nfound = CountElements( matrix, options.default )

        ## Apply filtering to decide whether matrix should be output

        ## 1. Filter by number of rows/columns
        if (options.filter_numrows and len(row_tokens) != options.filter_numrows) or\
               (options.filter_numcols and len(col_tokens) != options.filter_numcols):
            if options.loglevel >= 1:
                options.stdlog.write( "# id=%s, nrows=%i, ncols=%i, nfound=%i, ndefault=%i\n" % ( matrix_id, len(row_tokens), len(col_tokens), nfound, ndefault) )
            return False

        if options.output_format == "square":

            for col_token, index in col_tokens:
                outfile.write("\t%s" % col_token)
            outfile.write("\n")

            for row in range(len(matrix)):
                outfile.write("%s\t%s\n" % (row_tokens[row][0], string.join(matrix[row], "\t")))

        elif options.output_format == "phylip":

            if len(row_tokens) != len(col_tokens):
                raise ValueError, "phylip needs symmetric matrices."
            
            outfile.write ("%i\n" % len(row_tokens) )
            
            for row in range(len(matrix)):
                outfile.write("%-10s\t%s\n" % (row_tokens[row][0], " ".join(map(lambda x: "  %10s" % x, matrix[row] ))))

        elif options.output_format == "phylip-replicates":

            if len(row_tokens) != len(col_tokens):
                raise ValueError, "phylip needs symmetric matrices."
            
            outfile.write ("%i\n" % len(row_tokens) )

            for row in range(len(matrix)):
                outfile.write("%-10s" % row_tokens[row][0])
                rr = matrix[row]
                re = replicates[row]
                for c in range(len(matrix[row])):
                    outfile.write(" %10s %i" % (rr[c], re[c]))
                outfile.write("\n")
    else:
        if not options.row_names or not options.col_names:
            raise "Please specify row and column range."

        row_range = eval(options.row_names)
        col_range = eval(options.col_names)

        map_row = {}
        map_col = {}
        for x in range(len(row_range)): map_row[row_range[x]] = x
        for x in range(len(col_range)): map_col[col_range[x]] = x

        matrix = [ [ options.default for j in col_range] for i in row_range ]
        for line in lines:
            row_token, col_token, weight = string.split(line[:-1], "\t")[:3]

            row_pos = map_row[int(float(row_token))]
            col_pos = map_col[int(float(col_token))]
            # col_pos = int(float(col_token)) - col_range[0]

            if row_pos < 0 or row_pos > row_range[-1]: continue
            if col_pos < 0 or col_pos > col_range[-1]: continue                

            matrix[row_pos][col_pos] = weight

            if not options.asymmetric:
                matrix[col_pos][row_pos] = weight

        for col in col_range:
            outfile.write("\t%i" % col)
        outfile.write("\n")

        ndefault, nfound = CountElements( matrix, options.default )

        for row in range(len(matrix)):
            outfile.write("%i\t%s\n" % (row_range[row], "\t".join(map(str, matrix[row]))))
            
    if options.loglevel >= 1:
        options.stdlog.write( "# id=%s, nrows=%i, ncols=%i, nfound=%i, ndefault=%i\n" % ( matrix_id, len(row_tokens), len(col_tokens), nfound, ndefault) )

    return True

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: sparse2full.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="""method to use [normalize-by-min-diagonal|normalize-by-column|
                      log|ln|negzero2value|normalize-by-matrix|subtract-matrix|mix-matrix|
                      normalize-by-column-max|normalize-by-row-max|correspondence-analysis|
                      transpose|upper-bound|lower-bound
                      """  )
    
    parser.add_option("-f", "--format", dest="format", type="string",
                      help="format."  )
    
    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("square", "phylip", "phylip-replicates"),
                      help="output format."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("row-col-weight", "row-col-weight-replicates", "row-col-weight-weight" ),
                      help="input format."  )
    
    parser.add_option("-a", "--asymmetric", dest="asymmetric", action="store_true",
                      help="matrix is asymmetric."  )

    parser.add_option("-d", "--default", dest="default", type="string",
                      help="default value for missing values."  )

    parser.add_option("-D", "--default-diagonal", dest="default_diagonal", type="string",
                      help="default value for missing values on diagonal."  )

    parser.add_option("-r", "--row-names", dest="row_names", type="string",
                      help="list of row names separated by ','."  )
    
    parser.add_option("-c", "--col-names", dest="col_names", type="string",
                      help="list of col names separated by ','."  )

    parser.add_option( "--file-row-names", dest="file_row_names", type="string",
                      help="filename with row names."  )
    
    parser.add_option( "--file-col-names", dest="file_col_names", type="string",
                      help="filename with col names."  )

    parser.add_option("--full2sparse", dest="full2sparse", action="store_true",
                      help="convert full to sparse matrix."  )

    parser.add_option("--skip-separators", dest="write_separators", action="store_false",
                      help="do not echo separators (starting with >)")

    parser.add_option("-n", "--numeric", dest="is_numeric", action="store_true",
                      help="row and column titles are numeric (for sorting)."  )

    parser.add_option( "--filter-numcols", dest="filter_numcols", type="int",
                      help="only output matrices with # columns."  )

    parser.add_option( "--filter-numrows", dest="filter_numrows", type="int",
                      help="only output matrices with # rows."  )

    parser.add_option( "--filename-map", dest="filename_map", type="string",
                      help="filename with mapping between input and output chunks."  )

    parser.set_defaults(
        default = "0",
        default_diagonal = "0",
        asymmetric = False,
        row_names = None,
        col_names = None,
        file_row_names = None,
        file_col_names = None,
        full2sparse = False,
        format = "string",
        is_numeric = False,
        output_format = "square",
        input_format = "row-col-weight",
        filter_numcols = 0,
        filter_numrows = 0,
        filename_map = None,
        write_separators = True,
        )

    (options, args) = E.Start( parser )

    if options.full2sparse:
        # convert a full matrix to a sparse matrix

        options.stdout.write("row\tcol\tvalue\n" )

        if options.row_names:
            row_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.row_names, "r").readlines())
        else:
            row_tokens = None
            
        col_tokens = None
        
        row = 0
        for line in sys.stdin:
            if line[0] == "#": continue

            data = line[:-1].split("\t")

            row += 1
            if row == 1:
                if not col_tokens:
                    if options.col_names:
                        col_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.col_names, "r").readlines())
                    else:
                        if not row_tokens: del data[0]
                        col_tokens = map(str, data) 
                    continue
                
            if row_tokens:
                row_token = row_tokens[row]
            else:
                row_token = data[0]
                del data[0]
            
            for x in range(len(data)):
                if data[x] != options.default:
                    options.stdout.write("%s\t%s\t%s\n" % (row_token, col_tokens[x], data[x] ) )

    else:
        nskipped = 0
        
        if options.filename_map:
            outfile_map = open(options.filename_map, "w")
            outfile_map.write("## Map between matrices in input file and output matrices.\nnew\told\n" )
            
        # convert a sparse matrix to a full matrix        
        lines = filter( lambda x: x[0] != "#", sys.stdin.readlines())

        chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)) )

        if not chunks:
            chunks = [-1]
            options.write_separators = False
            
        chunks.append( len(lines) )

        if options.loglevel >= 2:
            print "# processing chunks:", chunks

        map_token2row = {}
        map_token2col = {}

        if options.file_row_names:
            row_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.file_row_names, "r").readlines())
            for row_token in row_tokens:
                map_token2row[row_token] = len(map_token2row)

        if options.row_names:
            for x in options.row_names.split(","):
                map_token2row[x] = len(map_token2row)

        if options.file_col_names:
            col_tokens = map( lambda x: string.split(x[:-1], "\t")[0], open(options.file_col_names, "r").readlines())
            for col_token in col_tokens:
                map_token2col[col_token] = len(map_token2col)

        if options.col_names:
            for x in options.col_names.split(","):
                map_token2col[x] = len(map_token2col)

        if not options.asymmetric:
            if map_token2row and not map_token2col:
                map_token2col = map_token2row
            elif map_token2col and not map_token2row:
                map_token2row = map_token2col

        noutput = 0
        
        for x in range(len(chunks) -1 ):

            if options.write_separators:
                options.stdout.write( lines[chunks[x]] )
                    
            output = Sparse2Matrix( options.stdout, str(x+1), lines[chunks[x]+1:chunks[x+1]],
                                    options,
                                    map_token2row, map_token2col )

            if not output:
                nskipped += 1
            else:
                noutput += 1
                if options.filename_map:
                    outfile_map.write( "%i\t%i\n" % (noutput, x+1))

        if options.filename_map:
            outfile_map.close()
                
    E.Stop()
