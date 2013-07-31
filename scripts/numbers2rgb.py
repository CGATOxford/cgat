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
numbers2rgb.py - map numbers to RGB values
==========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Reads numbers from stdin and maps to random RGB values.

Usage
-----

Example::

   python numbers2rgb.py --help

Type::

   python numbers2rgb.py --help

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
import optparse
import time
import random

import CGAT.Experiment as E

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: numbers2rgb.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-m", "--method", dest="method", type="string" ,
                       help="method to use.")

    parser.set_defaults(
        method = "random_rgb",
        )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    numbers = map( lambda x: int(x[:-1].split("\t")[0]),
                   filter( lambda x: x[0] != "#", sys.stdin.readlines()) )

    if options.method == "random_rgb":
        f = lambda x: "%i,%i,%i" % (random.randint(0,256),
                                    random.randint(0,256),
                                    random.randint(0,256) )

    map_number2output = {}
    for x in numbers:
        if x not in map_number2output:
            map_number2output[x] = f(x)

    for x, y in map_number2output.items():
        options.stdout.write( "%i\t%s\n" % (x, y) )

    E.Stop()
    
