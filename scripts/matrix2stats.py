################################################################################
#   Gene prediction pipeline 
#
#   $Id: matrix2stats.py 2795 2009-09-16 15:29:23Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
matrix2stats.py - compute statistics on matrices
================================================

:Author: Andreas Heger
:Release: $Id: matrix2stats.py 2795 2009-09-16 15:29:23Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

compute statistics on matrices.

Chi-Squared test
++++++++++++++++

Compute the Chi-squared test on a contingency table.
If --pairwise is given, the pairwise

The command line option *iteration* denotes what values will
be computed:

   pairwise
      compute all pairs assuming symmetry: `((1,2), (1,3), (2,3))`

   all-vs-all
      compute all-vs-all pairwise explicetely: `((1,2), (1,3), (2,1), (2,3), (3,1), (3,2))`

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

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
import CGAT.MatlabTools as MatlabTools
import CGAT.Stats as Stats
import scipy

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: matrix2stats.py 2795 2009-09-16 15:29:23Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("chi-squared", "pearson-chi-squared" ),
                      help="statistical methods to apply." )
    
    parser.add_option("-t", "--headers", dest="headers", action="store_true",
                      help="matrix has row/column headers."  )

    parser.add_option("--no-headers", dest="headers", action="store_false",
                      help="matrix has no row/column headers."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("full", "sparse", "phylip" ),
                      help="""input format for matrix."""  )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("full", "sparse", "phylip" ),
                      help="""output format for matrix."""  )
    
    parser.add_option("-p", "--parameters", dest="parameters", action="append", type="string",
                      help="parameters for various functions." )


    parser.add_option("-a", "--iteration", dest="iteration", type="choice",
                      choices=("pairwise", "all-vs-all" ),
                      help="""how to compute stats [%default]."""  )

    parser.set_defaults(
        method = "chi-squared",
        headers = True,
        value_format = "%6.4f",
        pvalue_format = "%6.4e",
        input_format = "full",        
        write_separators = True,
        parameters = [],
        iteration = None,
        )

    (options, args) = E.Start( parser )
    
    lines = filter( lambda x: x[0] != "#", sys.stdin.readlines())

    chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)) )

    if not chunks:
        options.write_separators = False
        chunks = [-1]
        
    chunks.append( len(lines) )

    ninput, noutput, nskipped = 0, 0, 0

    if options.write_separators:
        options.stdout.write( "test\t" )

    header_prefix = ""

    if options.method == "chi-squared":
        header_prefix = "observed\texpected"
        options.stdout.write( "\t".join( ( header_prefix, "n", "min", "max", "chi", "df", "P", "passed", "phi" ) ) + "\n" )

    elif options.method in ("pearson-chi-squared",):
        options.stdout.write("column\t")
        options.stdout.write( "\t".join( ( header_prefix, "n", "prob", "obs", "exp", "chi", "df", "P", "passed", "phi" ) ) + "\n" )

        if len(options.parameters) == 0:
            raise "out of parameters - please supply probability or filename with probabilities."

        param = options.parameters[0]
        del options.parameters[0]

        if options.write_separators:
            probabilities = IOTools.ReadMap( open( param, "r" ), map_functions = (str,float) )
        else:
            probability = float( param )


    for x in range(len(chunks) -1 ):
        ninput += 1
        matrix, row_headers, col_headers = MatlabTools.readMatrix( StringIO.StringIO("".join(lines[chunks[x]+1:chunks[x+1]])),
                                                                       format=options.input_format, 
                                                                       headers = options.headers )
        nrows, ncols = matrix.shape

        if options.loglevel >= 2:
            options.stdlog.write( "# read matrix: %i x %i, %i row titles, %i colum titles.\n" %\
                                      (nrows, ncols, len(row_headers), len(col_headers)))

        if options.write_separators:
            options.stdout.write( lines[chunks[x]][1:-1] + "\t" )

        pairs = []
        if options.iteration == "pairwise":
            pairs = []
            for row1 in range( 0, len(row_headers) ):
                for row2 in range( row1+1, len(row_headers) ):
                    pairs.append( (row1, row2) )
        elif options.iteration == "all-vs-all":
            pairs = []
            for row1 in range( 0, len(row_headers) ):
                for row2 in range( 0, len(row_headers) ):
                    if row1 == row2: continue
                    pairs.append( (row1, row2) )
    
        if options.method == "chi-squared":
            
            for row1, row2 in pairs:
                row_header1 = row_headers[row1]
                row_header2 = row_headers[row2]
                try:
                    result = Stats.doChiSquaredTest( numpy.vstack( (matrix[row1], matrix[row2] ) ) )
                except ValueError:
                    nskipped += 1
                    continue

                noutput += 1
                options.stdout.write( "\t".join( ( "%s" % row_header1,
                                                   "%s" % row_header2,
                                                   "%i" % result.mSampleSize,
                                                   "%i" % min(matrix.flat),
                                                   "%i" % max(matrix.flat),
                                                   options.value_format % result.mChiSquaredValue,
                                                   "%i" % result.mDegreesFreedom,
                                                   options.pvalue_format % result.mProbability,
                                                   "%s" % result.mSignificance,
                                                   options.value_format % result.mPhi ) ) + "\n" )
            

        elif options.method == "pearson-chi-squared":

            if nrows != 2:
                raise "only implemented for 2xn table"
            
            if options.write_separators:
                id = re.match( "(\S+)", lines[chunks[x]][1:-1] ).groups()[0]
                probability = probabilities[id]
                
            for col in range(ncols):
                options.stdout.write( "%s\t" % col_headers[col] )
                result = Stats.doPearsonChiSquaredTest( probability , sum(matrix[:,col]), matrix[0,col] )
                options.stdout.write( "\t".join( ( "%i" % result.mSampleSize,
                                                   "%f" % probability,
                                                   "%i" % result.mObserved,
                                                   "%f" % result.mExpected,
                                                   options.value_format % result.mChiSquaredValue,
                                                   "%i" % result.mDegreesFreedom,
                                                   options.pvalue_format % result.mProbability,
                                                   "%s" % result.mSignificance,
                                                   options.value_format % result.mPhi ) )  )
                if col < ncols - 1:
                    options.stdout.write("\n")
                    if options.write_separators:
                        options.stdout.write( lines[chunks[x]][1:-1] + "\t" )

            options.stdout.write( "\n" )

    E.info( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    E.Stop()
