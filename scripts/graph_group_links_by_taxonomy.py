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
graph_group_links_by_taxonomy.py - 
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

   python graph_group_links_by_taxonomy.py --help

Type::

   python graph_group_links_by_taxonomy.py --help

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
import optparse

import CGAT.Experiment as E
import scipy
import scipy.stats

USAGE="""python %s [OPTIONS]

Aggregate Evalues according to taxonomic ranges.

Taxonomic ranges are given by: name:from-to

The default input is a tab-separated table, sorted by query with the following fields:
   1. Query
   2. Taxonomic node_id
   3. Weight (blast E-value, etc.)
   4. Identifier of sequence
"""



def printResult( outfile, id, ranges, values_within, values_without, f, options):
    """print results in one line."""

    outfile.write( id )

    if options.write_id:
        
        ## reporting ids: sort by function, instead of selecting value
        
        for header, a, b in ranges:
            
            v = values_within[header]
            if len( v ) == 0:
                f1 = "na\tna"
            else:
                v.sort()
                if options.function == "max": v.reverse()
                f1 = "%s\t%s" % (v[0][1], options.value_format % v[0][0])

            v = values_without[header]
            if len( v ) == 0:
                f2 = "na\tna"
            else:
                v.sort()
                if options.function == "max": v.reverse()
                f2 = "%s\t%s" % (v[0][1], options.value_format % v[0][0])

            outfile.write( "\t%s\t%s" % (f1, f2) )
        
    else:
    
        for header, a, b in ranges:

            v = values_within[header]            
            if len( v ) == 0:
                f1 = "na"
            else:
                f1 =options.value_format % f(v)

            v = values_without[header]
            if len( v ) == 0:
                f2 = "na"
            else:
                f2 =options.value_format % f(v)

            outfile.write( "\t%s\t%s" % (f1, f2) )
        
    outfile.write( "\n" )

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: graph_group_links_by_taxonomy.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms." )

    parser.add_option("-r", "--range", dest="aranges", type="string", action="append",
                      help="taxonomic range." )

    parser.add_option("-f", "--function", dest="function", type="choice",
                      choices=("count", "sum", "min", "max" ),
                      help="choice of aggregate function." )

    parser.add_option("-i", "--write-id", dest="write_id", action="store_true",
                      help="report id and value." )

    parser.set_defaults(
        aranges = [],
        function="count",
        value_format = "%6.4f",
        write_id = False,
        )
    
    (options, args) = E.Start( parser )

    if options.write_id:
        if options.function not in ("min", "max"):
            raise "only min and max allowed with --write-id option"
        f = None
    else:
        if options.function == "count":
            f = len
        elif options.function == "min":
            f = min
        elif options.function == "max":        
            f = max
        elif options.function == "mean":
            f = scipy.mean

    options.aranges += args

    ranges = []
    for r in options.aranges:

        a, b = r.split(":")
        first, last = map(int, b.split("-"))

        ranges.append( (a, first, last ) )

    ## print header
    options.stdout.write("id" )
    for header, first, last in ranges:
        if options.write_id:
            options.stdout.write("\t%s_within_id\t%s_within\t%s_without_id\t%s_without" % (header, header, header, header) )
        else:
            options.stdout.write("\t%s_within\t%s_without" % (header, header ) )
    options.stdout.write("\n")
    
    values_within = {}
    values_without = {}
    for header, first, last in ranges:
        values_within[header] = []
        values_without[header] = []

    last_id = None

    for line in sys.stdin:
        
        if line[0] == "#": continue
        
        data = line[:-1].split( "\t" )

        if data[0] != last_id:
            if last_id:
                printResult( options.stdout, last_id, ranges, values_within, values_without, f,
                             options = options )
                
            values_within = {}
            values_without = {}
            
            for header, first, last in ranges:
                values_within[header] = []
                values_without[header] = []
                
            last_id = data[0]
            
        val = float(data[2])

        for header, first, last in ranges:
            within = first <= int(data[1]) <= last

            if options.write_id:
                xval = ( val, data[3] )
            else:
                xval = val
                
            if within:
                values_within[header].append( xval )
            else:
                values_without[header].append( xval )
                
    printResult( options.stdout, last_id, ranges, values_within, values_without, f,
                 options = options )

    E.Stop()
