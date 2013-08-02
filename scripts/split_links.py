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
split_links.py - 
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

   python split_links.py --help

Type::

   python split_links.py --help

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
import CGAT.WrapperAdaptiveCAI as WrapperAdaptiveCAI
import numpy

open_files = {}

def WriteLine( a, b, line, prefix="%s-%s" ):

    key1 = prefix % (a,b)
    key2 = prefix % (b,a)    
    if key1 in open_files or os.path.exists( key1 ):
        key = key1
    else:
        key = key2
    
    if key not in open_files:
        open_files[key] = open( key, "a" )

    f = open_files[key]
    f.write( line )
    f.flush()
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: split_links.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method for splitting." )

    parser.add_option("-r", "--regex", dest="regex", type="string",
                      help="regex to find prefix." )

    parser.add_option("-o", "--output", dest="output", type="string",
                      help="output filename." )

    parser.add_option("-t", "--targets", dest="targets", type="string",
                      help="output filename." )

    parser.set_defaults()

    (options, args) = E.Start( parser )

    if options.targets: options.targets = options.targets.split(",")

    nsame = 0
    ndiff = 0
    
    if options.method=="prefix":

        for line in sys.stdin:
            if line[0] == "#": continue

            data = line[:-1].split("\t")

            g1 = re.search( options.regex, data[0] ).groups()[0]
            g2 = re.search( options.regex, data[1] ).groups()[0]            

            if g1 == g2:
                for t in options.targets:
                    if g1 == t: continue
                    WriteLine( g1, t, line, options.output )
                    nsame += 1
            else:
                WriteLine( g1, g2, line, options.output )
                ndiff += 1

    print "nsame=%i, ndiff=%i" % (nsame, ndiff )

    E.Stop()
            
        
    
    
