'''
modify_table.py - 
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

   python modify_table.py --help

Type::

   python modify_table.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import tempfile
import subprocess
import optparse
import time
import math

"""read in data and append columns to a density histogram

-> relative frequencies
-> cumulative counts and frequencies in both directions

'#' at start of line is a comment
"""

import CGAT.Experiment as E
import numpy
import scipy

parser = E.OptionParser( version = "%prog version: $Id: modify_table.py 2782 2009-09-10 11:40:29Z andreas $")

##---------------------------------------------------------------------------------------------------------        

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take from table." )

    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="methods to apply to columns.",
                      action="append")

    parser.add_option("-e", "--echo", dest="echo",
                      action="store_true",
                      help="echo columns not taken." )

    parser.add_option("-r", "--replace", dest="replace",
                      action="store_true",
                      help="replace orginial values." )

    parser.set_defaults( \
        columns="1",
        echo = False,
        replace = False,
        format="%5.2f",
        methods=[])

    (options, args) = E.Start( parser )
    
    options.columns = map( lambda x: int(x)-1, options.columns.split(","))

    print E.GetHeader()
    print E.GetParams()

    vals = []
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    headers = lines[0][:-1].split("\t")
    del lines[0]

    notcolumns = filter( lambda x: x not in options.columns, range( len(headers)))
    
    data = [ [] for x in range(len(headers))]


    for l in lines:
        d = l[:-1].split("\t")
        for c in options.columns:
            data[c].append( float(d[c]) )
        for c in notcolumns:
            data[c].append( d[c] )
            
    if len(data) == 0:
        raise "No data found."
        
    totals = [0] * len(headers)
    
    for c in options.columns:
        totals[c] = reduce( lambda x,y: x+y, data[c] )

    new_columns = []
    new_headers = []

    if options.echo:
        for c in notcolumns:
            new_headers.append( headers[c] )
            new_columns.append( data[c] )
        
    for c in options.columns:
        if not options.replace:
            new_columns.append( data[c] )
            new_headers.append( headers[c] )
            
        for method in options.methods:
            if method == "normalize":
                new_columns.append( map( lambda d: d / totals[c], data[c] ) )
                new_headers.append( "normalized" )

    print string.join( new_headers, "\t")

    for d in zip( *new_columns ):
        print string.join( map(str, d), "\t")
    
    E.Stop()








if __name__ == "__main__":
    sys.exit( main( sys.argv) )

