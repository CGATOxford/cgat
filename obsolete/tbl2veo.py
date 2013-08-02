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
tbl2veo.py - 
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

   python tbl2veo.py --help

Type::

   python tbl2veo.py --help

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

"""convert tab-separated table to VidaExpert formatted data file.
"""

import CGAT.Experiment as E

parser = E.OptionParser( version = "%prog version: $Id: tbl2veo.py 2782 2009-09-10 11:40:29Z andreas $")

if __name__ == "__main__":


    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take from table." )


    (options, args) = E.Start( parser )

    lines = map( string.strip, filter( lambda x: x[0] != "#", sys.stdin.readlines()))

    num_vals = len(lines)
    header = lines[0].split("\t")
    del lines[0]
    num_cols = len(header)

    print "%i %i" % (num_cols, num_vals)

    for h in header:
        print h, "FLOAT"


    for l in lines:
        print l.replace( "\t", " ")

    E.Stop()
