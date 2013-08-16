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
graph_check.py - 
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

   python graph_check.py --help

Type::

   python graph_check.py --help

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

## patch for old python installations
if sys.version_info < (2,4):
    from sets import *
    set = Set

""" program $Id: graph_check.py 2782 2009-09-10 11:40:29Z andreas $

python graph_check.py < graph.in

Check graph for completeness.

"""
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

def writeSet( outfile, outset ):
    """write set to file."""
    outlist = list(outset)
    outlist.sort()
    for x in outlist: outfile.write( "%s\n" % x )

def writeInfo( outfile, vertices, nlinks, nlines, nerrors, ncomments, is_sorted ):
    
    all = set(vertices.keys())
    queries = set(filter( lambda x: vertices[x] & 1, all) )
    sbjcts  = set(filter( lambda x: vertices[x] & 2, all) )
    ## count only those as missed self, that do have queries.
    missed_self = set(filter( lambda x: vertices[x] & 4 == 0, all) ).intersection( queries )

    missed_queries = all.difference( queries)
    missed_sbjcts =  all.difference( sbjcts )

    if is_sorted:
        sorted= "yes"
    else:
        sorted = "no"

    outfile.write( "\t".join( map(str, (len(queries), len(sbjcts),
                                               len(queries.union( sbjcts ) ),
                                               nlinks,
                                               nlines, nerrors, ncomments,
                                               sorted,
                                               len(all),
                                               len(missed_queries),
                                               len(missed_sbjcts),
                                               len(missed_self) ) ) ) + "\n") 
    outfile.flush()

    return missed_queries, missed_sbjcts, missed_self
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph_check.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("--filename-missing", dest="filename_missing", type="string",
                      help="missing entries.")
    parser.add_option("--filename-found", dest="filename_found", type="string",
                      help="found entries.")
    parser.add_option("--report-step1", dest="report_step1", type="int",
                      help="report interval for input.")
    parser.add_option("--report-step2", dest="report_step2", type="int",
                      help="report interval for processing.")
    parser.add_option("-n", "--filename-vertices", dest="filename_vertices", type="string",
                      help="filename with vertices.")
    parser.add_option("-u", "--num-fields", dest="num_fields", type="int",
                      help="number of fields to expect.")
    parser.add_option("-o", "--filename-output-pattern", dest="filename_output_pattern", type="string",
                      help="filenames for output (should contain one %s for one section).")
    parser.add_option("-s", "--sort-order", dest="sort_order", type="choice",
                      choices=("numeric", "alphanumeric" ),
                      help="sort order - if numeric, vertices are cast to int.")

    parser.set_defaults(
        filename_vertices = None,
        report_step1 = 100000,
        report_step2 = 10000,
        filename_output_pattern = "%s",
        subsets = False,
        num_fields = 11,
        sort_order = "alphanumeric",
        )

    (options, args) = E.Start( parser )

    if options.loglevel >= 1:
        options.stdlog.write("# output goes to:\n" )
        options.stdlog.write("# errors: %s\n" % options.filename_output_pattern % "errors" )
        options.stdlog.write("# missed query: %s\n" % options.filename_output_pattern % "missed_queries" )
        options.stdlog.write("# missed sbjct: %s\n" % options.filename_output_pattern % "missed_sbjcts" )
        options.stdlog.write("# missed self: %s\n" % options.filename_output_pattern % "missed_self" )                
        
    outfile_errors = open( options.filename_output_pattern % "errors", "w" )

    if options.sort_order == "numeric":
        f = int
    else:
        f= str

    if options.filename_vertices:
        vv, errors = IOTools.ReadList( open( options.filename_vertices, "r" ), map_function = f )
        vertices = {}
        ## use flags for vertices
        ## 1st bit: is query: 1
        ## 2nd bit: is sbjct: 2
        ## 3rd bit: has self: 4
        for v in vv:
            vertices[v] = 0
    else:
        raise "for the time being, specify a vertex file."

    options.stdout.write( "nqueries\tnsbjcts\tnvertices\tnlinks\tnlines\tnerrors\tncomments\tis_sorted\tnexpected\tnmissed_queries\tnmissed_sbjcts\tnmissed_self\n" )

    ncomments, nlinks, nerrors, nlines = 0, 0, 0, 0

    is_sorted = True

    last = None

    
    for line in sys.stdin:

        nlines += 1
        if line[0] == "#":
            ncomments += 1
            continue

        nlinks += 1

        data = line[:-1].split("\t")

        if len(data) != options.num_fields:
            nerrors += 1
            outfile_errors.write( line )
            outfile_errors.flush()
            continue

        q, s = f(data[0]), f(data[1])

        if q == s: vertices[q] |= 4
        vertices[q] |= 1
        vertices[s] |= 2

        if last and last > q:
            is_sorted = False
            outfile_errors.write("# sort inconsistency between %s and %s at line %i\n" % ( last, q, nlines) )
            outfile_errors.flush()
            if options.loglevel >= 1:
                options.stdlog.write("# sort inconsistency between %s and %s at line %i\n" % ( last, q, nlines) )
                options.stdlog.flush()
                
        if options.report_step1 and nlines % options.report_step1 == 0:
            writeInfo( options.stdlog, vertices, nlinks, nlines, nerrors, ncomments, is_sorted )
            
        last = q

    missed_queries, missed_sbjcts, missed_self = writeInfo( options.stdout, vertices, nlinks, nlines, nerrors, ncomments, is_sorted )

    if nerrors == 0:
        os.remove( options.filename_output_pattern % "errors" )

    if missed_queries:
        writeSet( open( options.filename_output_pattern % "missed_queries", "w" ), missed_queries )

    if missed_sbjcts:
        writeSet( open( options.filename_output_pattern % "missed_sbjcts", "w" ), missed_sbjcts )

    if missed_self:        
        writeSet( open( options.filename_output_pattern % "missed_self", "w" ), missed_self )

    E.Stop()
