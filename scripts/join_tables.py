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
join_tables.py - join tables
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

   python join_tables.py --help

Type::

   python join_tables.py --help

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

""" program $Id: join_tables.py 2782 2009-09-10 11:40:29Z andreas $

"""
import CGAT.Experiment as E

USAGE = """python %s < stdin > stdout

OPTIONS:

-t, --titles=           column titles
-m, --missing=          missing value
-h, --headers=          add headers for files
-s, --sort=             sort by column titles (given by sort order)
'#' at start of line is a comment
""" % sys.argv[0]

def GetData( f ):

    while 1:
        line = f.readline()
        if not line: return None
        if line[0] == "#": continue
        break
    
    return line[:-1].split("\t")
        
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: join_tables.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to do join on. Files have to sorted alphabeticlly and incrementally on these columns." )
    parser.add_option("-m", "--missing", dest="missing_value", type="string",
                      help="value for missing entries." )

    parser.set_defaults(
        headers = False,
        pattern_filename = None,
        title = "",
        footer = "",
        columns = "1",
        missing_value = "",
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) < 2:
        raise "there have to be at least two tables."

    if options.columns: options.columns = map(lambda x: int(x) - 1, options.columns.split(","))
    
    ## open all files
    files = []
    for filename in args:

        if os.path.exists( filename ):
            files.append( open(filename, "r") )

    if len(files) <= 1:
        raise "less than two files opened."

    nfiles = len(files)
    
    ## go through all files and iteratively find smallest entry
    data = []
    entries = []
    lengths = []
    takes = []
    
    for f in files:
        d = GetData( f )
        l = len(d)
        t = []
        for x in range( l ):
            if x not in options.columns:
                t.append(x)
                
        takes.append(t)
        data.append( map( lambda x: d[x], t) )
        entries.append( map( lambda x: d[x], options.columns) )
        lengths.append( len(t) )

    activa = list(range( nfiles ))
    ## take only first entry for duplicate entries
    ## columns need to be sorted incrementally on all columns
    while len(activa) > 0:

        line = []
        
        min_field = min([entries[x] for x in activa])
        
        for f in range(nfiles):

            if files[f] and entries[f] == min_field:
                
                line += data[f]

                while entries[f] == min_field:
                    
                    d = GetData(files[f])

                    if not d:
                        files[f].close()
                        files[f] = None
                        activa.remove( f )
                        break
                    
                    data[f] = map( lambda x: d[x], takes[f])
                    entries[f] = map( lambda x: d[x], options.columns)
                
            else:
                line += [options.missing_value for x in range(lengths[f])]
                
        options.stdout.write( "\t".join(min_field) + "\t" + "\t".join(line) + "\n" )

    ## close all files
    for f in files:
        if f:
            f.close()

    E.Stop()
