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
graph_check_transitivity.py - 
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

   python graph_check_transitivity.py --help

Type::

   python graph_check_transitivity.py --help

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

""" program $Id: graph_check_transitivity.py 2782 2009-09-10 11:40:29Z andreas $

python graph_check_transitivity < graph.in

check whether all edges in a graph are transitive, i.e.,
for every two edges A->B and B->C check whether A->C exists.

Edges are taken to be undirected.
"""
import CGAT.Experiment as E
import CGAT.Histogram as Histogram

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph_check_transitivity.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("--filename-missing", dest="filename_missing", type="string",
                      help="missing entries.")
    parser.add_option("--filename-found", dest="filename_found", type="string",
                      help="found entries.")
    parser.add_option("--report-step1", dest="report_step1", type="int",
                      help="report interval for input.")
    parser.add_option("--report-step2", dest="report_step2", type="int",
                      help="report interval for processing.")
    parser.add_option("--use-subsets", dest="subsets", action="store_true",
                      help="do subset calculation. Third field contains a redundancy code.")
    
    parser.set_defaults(
        filename_missing = None,
        filename_found = None,
        report_step1 = 100000,
        report_step2 = 10000,
        subsets = False,
        )

    (options, args) = E.Start( parser )

    # retrieve data
    vals = {}
    niterations = 0
    ninput = 0

    for line in sys.stdin:
        if line[0] == "#": continue
        
        niterations += 1
        
        if options.loglevel >= 1 and (niterations % options.report_step1 == 0):
            options.stdlog.write( "# input: %i\n" % (niterations))
            options.stdlog.flush()

        v1, v2, w = line[:-1].split("\t")[:3]
        
        if v1 == v2: continue
        
        if v1 not in vals: vals[v1] = []
        if v2 not in vals: vals[v2] = []
        
        if not options.subsets:
            w = ninput
            
        vals[v1].append( (v2, w) )
        vals[v2].append( (v1, w) )
            
        ninput += 1

    ## make everything unique
    for key, v1 in vals.items():
        vals[key] = tuple(set(v1))

    keys = vals.keys()
    keys.sort()
    niterations = 0
    nkeys = len(keys)
    missing = []
    ntotal = 0
    nfound = 0
    counted = {}
    nremoved = 0

    if options.filename_found:
        outfile_found = open(options.filename_found, "w")
        
    for v1 in keys:

        niterations += 1
        
        if options.loglevel >= 1 and (niterations % options.report_step2 == 0):
            options.stdlog.write( "# loop: %i\n" % (niterations))
            options.stdlog.flush()

        for v2, c2 in vals[v1]:

            ## only to half-symmetric test
            
            for v3, c3 in vals[v2]:
                if (c2, c3) in counted:
                    nremoved += 1
                    # print "v1=", v1, "v2=", v2, "v3=", v3, "c2=", c2, "c3=", c3, "removed"                    
                    continue

                ## do not do self-comparisons
                if v1 == v3: continue
                if c2 == c3: continue
                
                counted[(c2,c3)] = True
                ntotal += 1
                if v3 in map(lambda x: x[0], vals[v1]) or v1 in map(lambda x: x[0], vals[v3]):
                    nfound += 1
                    if options.filename_found:
                        outfile_found.write( "\t".join( (v1, v2, v3) ) + "\n" )
                    # print "v1=", v1, "v2=", v2, "v3=", v3, "c2=", c2, "c3=", c3, "found"
                else:
                    missing.append( (v1, v2, v3) )
                    # print "v1=", v1, "v2=", v2, "v3=", v3, "c2=", c2, "c3=", c3, "missing"
        
    nmissing = len(missing)

    options.stdout.write( "number of egdes\t%i\n" % ninput)
    options.stdout.write( "number of vertices\t%i\n" % nkeys)
    options.stdout.write( "number of removed triplets\t%i\n" % nremoved)
    options.stdout.write( "number of tested triplets\t%i\t%6.4f\n" % (ntotal, float(ntotal) / float(ntotal)))
    options.stdout.write( "number of realized triplets\t%i\t%6.4f\n" % (nfound, float(nfound) / float(ntotal)))
    options.stdout.write( "number of incomplete triplets\t%i\t%6.4f\n" % (nmissing, float(nmissing) / float(ntotal)))

    if options.filename_missing:
        outfile = open(options.filename_missing, "w")
        for v1, v2, v3 in missing:
            outfile.write( "\t".join( (v1, v2, v3) ) + "\n")
        outfile.close()

    if options.filename_found:
        outfile_found.close()
        
    E.Stop()
