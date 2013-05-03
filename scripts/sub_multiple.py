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
sub_multiple.py - substitute multiple patterns per row
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

   python sub_multiple.py --help

Type::

   python sub_multiple.py --help

for command line help.

Documentation
-------------

Code
----

'''
import string
import re
import sys

if __name__ == "__main__":
    param_filename_patterns = sys.argv[1]

    infile = open(param_filename_patterns, "r")

    patterns = []
    for line in infile:
        (old, new) = line[:-1].split("\t")
        patterns.append( (re.compile( "%s" % old), new) )

    infile.close()

    for line in sys.stdin:

        for pattern, new in patterns:
            line = pattern.sub( new, line )

        print line[:-1]
