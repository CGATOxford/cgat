#! /bin/env python
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
csv2db.py - upload table to database
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

create a table from a csv separated file and load data into it.

This module supports backends for postgres and sqlite3. Column types are
auto-detected.

Read a table from stdin and create an sqlite3 database. By default, the database
will reside in a file called csvdb and in a table csv.

  
.. todo::

   Use file import where appropriate to speed up loading. Currently, this is
   not always the case.

Usage
-----

Example::

   python csv2db.py -b sqlite < stdin 

Type::

   python csv2db.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import CGAT.Experiment as E
import CGAT.CSV2DB as CSV2DB
import CGAT.CSV as CSV

import csv
import sqlite3

csv.field_size_limit(sys.maxint)

def main( argv = sys.argv ):

    parser = CSV2DB.buildParser()

    (options, args) = E.Start( parser, argv = argv, add_psql_options = True )

    if options.from_zipped:
        import gzip
        infile = gzip.GzipFile( fileobj= options.stdin, mode='r')

    else:
        infile = options.stdin

    CSV2DB.run( infile, options )

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())
