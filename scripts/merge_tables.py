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
merge_tables.py - 
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

   python merge_tables.py --help

Type::

   python merge_tables.py --help

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

"""merge two tables with the same number of lines
"""

import CGAT.Experiment as E

import CGAT.WrapperAdaptiveCAI as WrapperAdaptiveCAI
import numpy

parser = E.OptionParser( version = "%prog version: $Id: merge_tables.py 2782 2009-09-10 11:40:29Z andreas $")

if __name__ == "__main__":

    parser.add_option("-t", "--table", dest="tables", type="string",
                      help="tables to merge.",
                      action="append")


    parser.set_defaults( \
        tables = [])

    (options, args) = E.Start( parser )

    if len(options.tables) < 1:
        raise "please specify at least one table."

    files = []
    for t in options.tables:
        files.append( open( t, "r"))

    while 1:
        frags = []
        
        stop = False
        
        for f in files:
            l = f.readline()
            if not l:
                stop = True
                break
            frags.append( l[:-1] )

        if stop: break

        print string.join( frags, "\t" )
    
    E.Stop()
