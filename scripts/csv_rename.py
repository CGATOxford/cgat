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
csv_rename.py - rename columns in a table
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

rename columns in a csv file

Usage
-----

Example::

   csv_rename.py gene=id < stdin

Type::

   python csv_rename.py --help

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

import CGAT.Experiment as E
import csv

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: csv_rename.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-r", "--remove", dest="remove", action="store_true",
                       help="remove specified columns, keep all others." )

    parser.add_option( "-u", "--unique", dest="unique", action="store_true",
                       help="output rows are uniq." )

    parser.add_option( "-f", "--filename-fields", dest="filename_fields", type="string",
                       help="filename with field information.")

    parser.set_defaults(
        filename_fields = None,
        )

    (options, args) = E.Start( parser,
                                        add_csv_options  = True )
    mapper = {}
    for x in args:
        a,b = x.split("=")
        mapper[a.strip()] = b.strip()
    
    while 1:
        line = sys.stdin.readline()
        
        if not line:
            E.Stop()
            sys.exit(0)
        
        if line[0] == "#": 
            sys.stdout.write( line )
            continue
        
        break

    header = []
    nreplaced = 0
    for x in line[:-1].split():
        if x in mapper:
            nreplaced += 1
            header.append( mapper[x] )
        else:
            header.append( x )
    
    options.stdout.write( "\t".join( header ) + "\n" )
    nlines = 0
    for line in sys.stdin:
        nlines += 1
        sys.stdout.write( line )
            
    if options.loglevel >= 1:
        ninput = len(header)
        noutput = ninput
        options.stdout.write( "# ninput=%i, noutput=%i, nreplaced=%i, nlines=%i\n" % (ninput,noutput,nreplaced,nlines) )

    E.Stop()
