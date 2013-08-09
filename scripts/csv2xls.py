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
csv2xls.py - convert table to excel format
==========================================

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

   python csv2xls.py --help

Type::

   python csv2xls.py --help

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
import CGAT.CSV as CSV

import openpyxl

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: csv2xls.py 2782 2009-09-10 11:40:29Z andreas $")
    
    parser.add_option( "-o", "--outfile=", dest="output_filename", type="string",
                       help="write to output filename." )

    parser.set_defaults(
        output_filename=None,
        )

    (options, args) = E.Start( parser, add_csv_options  = True)

    if not options.output_filename:
        raise ValueError("please specify an output filename.")

    w = openpyxl.Workbook( optimized_write = True )

    ## create styles
    header_style = GetHeaderStyle()
    data_style = GetDataStyle()

    for filename in args:

        lines = filter( lambda x: x[0] != "#", open(filename, "r").readlines())

        if len(lines) == 0: continue

        if options.loglevel >= 2:
            print "# read %i rows" % len(lines)
            sys.stdout.flush()

        headers = lines[0][:-1].split("\t")

        ws = w.add_sheet(os.path.basename(filename))
        
        cur_row = 0
        
        ws.append( headers )
            
        cur_row += 1
        
        reader = csv.DictReader( lines, dialect=options.csv_dialect )

        for row in reader:
            row = CSV.ConvertDictionary( row )

            data = [ row.get( headers[x], "") for x in range(len(headers))]
            ws.append( data )
                
            cur_row += 1

    
    w.save(options.output_filename)
        
    E.Stop()


