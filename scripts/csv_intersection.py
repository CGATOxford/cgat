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
csv_intersection.py - intersect two tables
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

   python csv_intersection.py --help

Type::

   python csv_intersection.py --help

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
import CGAT.CSV as CSV
import csv
import hashlib



class UniqueBuffer:
    mKeys = {}
    def __init__(self, outfile):
        self.mOutfile = outfile
    def write( self, out ):
        key = hashlib.md5(out).digest()
        if key not in self.mKeys:
            self.mKeys[key] = True
            self.mOutfile.write(out)
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: csv_intersection.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-u", "--unique", dest="unique", action="store_true",
                       help="output rows are uniq." )

    parser.set_defaults(
        remove = False,
        unique = False,
        )

    (options, args) = E.Start( parser, add_csv_options  = True)

    if len(args) != 2:
        raise "please specify two files to join."
        
    options.filename1, options.filename2 = args

    table1 = CSV.ReadTable( open(options.filename1, "r") )
    table2 = CSV.ReadTable( open(options.filename2, "r") )    

    if options.unique:
        outfile = UniqueBuffer(sys.stdout)
    else:
        outfile = options.stdout

    ## build new field list
    new_fields = []

    for x in options.join_fields1:
        new_fields.append( x )

    for x in fields1:
        if x not in options.join_fields1:
            new_fields.append( x )
        if x not in options.join_fields2:
            new_fields.append( x )
            
        writer = csv.DictWriter( outfile,
                                 fields,
                                 dialect=options.csv_dialect,
                                 lineterminator = options.csv_lineterminator,
                                 extrasaction = 'ignore' )
        

    
    if len(lines) > 0:
        
        old_fields = lines[0][:-1].split("\t")

        if options.remove:
            fields = []
            for x in old_fields:
                if x not in input_fields:
                    fields.append( x )
        else:
            fields = input_fields

        reader = csv.DictReader( lines,
                                 dialect=options.csv_dialect )
        
        print "\t".join(fields)        

        first_row = True
        for row in reader:
            row = CSV.ConvertDictionary( row )
            writer.writerow(row)

    E.Stop()
