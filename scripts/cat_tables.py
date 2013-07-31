################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2781 2009-09-10 11:33:14Z andreas $
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
"""
cat_tables.py - concatenate tables
==================================

:Author: Andreas Heger
:Release: $Id: cgat_script_template.py 2781 2009-09-10 11:33:14Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

concatenate tables. Headers of subsequent files are ignored.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse
import fileinput

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.set_defaults(
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) == 0 or (len(args) == 1 and args[0] == "-"):
        infile = options.stdin
    else:
        infile = fileinput.FileInput( args )
        
    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    header = False

    for line in infile:
        ninput += 1
        if line.startswith( "#" ):
            pass
        elif not header:
            header = line
        elif line == header:
            nskipped += 1
            continue
            
        options.stdout.write( line )
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

