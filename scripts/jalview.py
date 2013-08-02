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
jalview.py - build annotation for viewing in jalview
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will build annotations for a multiple alignment
such that both can be visualized with :file:`jalview`.

Usage
-----

Example::

   python jalview.py --help

Type::

   python jalview.py --help

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
import tempfile
import subprocess
import optparse
import time
import math

import CGAT.Experiment as E
import CGAT.Mali as Mali

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: jalview.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice", 
                      choices=("list2annotation", ),
                      help="methods.")

    parser.add_option("--filename-mali", dest="filename_mali", type="string",
                      help="filename with multiple alignment used for calculating sites - used for filtering" )
    
    parser.add_option("--jalview-title", dest="jalview_title", type="string",
                      help="title for jalview annotation." )
    

    parser.set_defaults(
        method = None,
        jalview_symbol = "*",
        jalview_title = "anno",
        filename_mali = None,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if not options.filename_mali:
        raise "please specify a multiple alignment."
    
    mali = Mali.Mali()
    mali.readFromFile( open(options.filename_mali, "r") )

    if options.method == "list2annotation":

        options.stdout.write("JALVIEW_ANNOTATION\n" )
        options.stdout.write("# Created: %s\n\n" % (time.asctime(time.localtime(time.time()))))

        codes = [""] * mali.getWidth()

        first = True
        for line in sys.stdin:
            if line[0] == "#": continue
            if first:
                first= False
                continue

            position = int(line[:-1].split("\t")[0])
            codes[position-1] = options.jalview_symbol
            
        options.stdout.write("NO_GRAPH\t%s\t%s\n" % (options.jalview_title, "|".join( codes ) ))
            
    E.Stop()
