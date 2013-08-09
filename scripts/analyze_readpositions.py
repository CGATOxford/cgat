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
analyze_readpositions.py - compute read coverage on transcripts
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


input = stream of counts on a transcript

Usage
-----

Example::

   python analyze_readpositions.py --help

Type::

   python analyze_readpositions.py --help

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
import getopt
import time
import optparse
import math
import tempfile

USAGE=""" program $Id: analyze_readpositions.py 2781 2009-09-10 11:33:14Z andreas $


"""
import CGAT.Experiment as E
import CGAT.CSV as CSV
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats
import numpy

def main():

    parser = E.OptionParser( version = "%prog version: $Id: analyze_readpositions.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="pattern for additional output files [%default]."  )

    parser.set_defaults(
        length=1000,
        minimum_coverage=0.90,
        maximum_reads = [1,10,20,50,100],
        output_filename_pattern = "%s",
        normalize = True,
        )

    (options, args) = E.Start( parser, add_csv_options = True )

    fields, table = CSV.ReadTable( sys.stdin, dictreader=CSV.DictReaderLarge )
    
    map_fields2column = {}
    for x in fields: map_fields2column[x] = len(map_fields2column)
    
    coverage_5prime = numpy.zeros( options.length, numpy.float )
    coverage_3prime = numpy.zeros( options.length, numpy.float )

    coverage_maxreads5prime = numpy.zeros( options.length, numpy.float )
    coverage_maxreads3prime = numpy.zeros( options.length, numpy.float )

    coverage_full5prime = numpy.zeros( options.length, numpy.float )
    coverage_full3prime = numpy.zeros( options.length, numpy.float )

    coverage_min5prime = numpy.zeros( options.length, numpy.float )
    coverage_min3prime = numpy.zeros( options.length, numpy.float )
    
    histograms = []
    for x in range(len( options.maximum_reads) ):
        histograms.append( 
            [ numpy.zeros( options.length, numpy.float ),
              numpy.zeros( options.length, numpy.float ),
              0 ] )

    ninput, noutput, nfull, nmincov, nskipped, nlength, nmaxreads = 0, 0, 0, 0, 0, 0, 0
    for row in table:
        length, covered, meancov, data, nreads = (int(row["cov_nval"]), 
                                                          float(row["cov_covered"]), 
                                                          float(row["cov_mean"]), 
                                                          row["cov_values"],
                                                          int(row["nover2"]) )
        ninput += 1
        if length < options.length: 
            nlength += 1
            continue
        
        if data == "na":
            nskipped += 1
            continue

        noutput += 1
        mincov = covered / length
        values = map( float, data.split(";") )
        m = max(values)
        values = [ x / m for x in values ]
        coverage_5prime += values[0:1000]
        coverage_3prime += values[-1000:]
        
        if mincov >= 1.0:
            coverage_full5prime += values[0:1000]
            coverage_full3prime += values[-1000:]
            nfull += 1

        if meancov >= options.minimum_coverage:
            coverage_min5prime += values[0:1000]
            coverage_min3prime += values[-1000:]
            nmincov += 1

        for maxreads in range( len(options.maximum_reads) ):
            if nreads <= options.maximum_reads[maxreads]:
                histograms[maxreads][0] += values[0:1000]
                histograms[maxreads][1] += values[-1000:]
                histograms[maxreads][2] += 1

    if options.normalize:
        for x5,x3 in ((coverage_5prime,
                       coverage_3prime),
                      (coverage_min5prime,
                       coverage_min3prime),
                      (coverage_full5prime,
                       coverage_full3prime) ):
            m = max( (max(x5), max(x3)) )
            x3 /= m
            x5 /= m

        for x5,x3,c in histograms:
            m = max( (max(x5), max(x3) ) )
            x5 /= m
            x3 /= m
            
    outfile = options.stdout
    outfile.write( "\t".join( ("distance", "minlen-5'", "minlen-3'", "mincov-5'", "mincov-3'", "full-5'", "full-3'" ) ) + "\n" )

    for x in range(0, options.length):
        outfile.write( "\t".join( [ "%6.4f" % x for x in \
                                        (x, 
                                         coverage_5prime[x],
                                         coverage_3prime[x],
                                         coverage_min5prime[x],
                                         coverage_min3prime[x],
                                         coverage_full5prime[x],
                                         coverage_full3prime[x] ) ] ) + "\n" )


    outfile5 = open(options.output_filename_pattern % "reads5", "w")
    outfile3 = open(options.output_filename_pattern % "reads3", "w")

    outfile5.write( "\t".join( ["distance",] + ["reads%i" % options.maximum_reads[y] for y in range(len(options.maximum_reads) ) ] ) + "\n" )
    outfile3.write( "\t".join( ["distance",] + ["reads%i" % options.maximum_reads[y] for y in range(len(options.maximum_reads) ) ] ) + "\n" )
    for x in range(0, options.length):
        outfile5.write( "%i\t%s\n" % (x, "\t".join( [ "%6.4f" % histograms[y][0][x] for y in range(len(options.maximum_reads) ) ] )))
        outfile3.write( "%i\t%s\n" % (x, "\t".join( [ "%6.4f" % histograms[y][1][x] for y in range(len(options.maximum_reads) ) ] )))

    E.info( "ninput=%i, noutput=%i, nmaxreads=%i, nfull=%i, nmincov=%i, nskipped=%i, nlength=%i" %\
                (ninput, noutput, nmaxreads, nfull, nmincov, nskipped, nlength) )
    
    E.Stop()

if __name__ == "__main__":
    sys.exit( main() )
