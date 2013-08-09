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
mali_phylip2fasta.py - 
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

   python mali_phylip2fasta.py --help

Type::

   python mali_phylip2fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    line = sys.stdin.readline()
    num_lines, width = map(int, re.search("(\d+)\s+(\d+)", line[:-1]).groups())

    alignment = []

    for x in range(num_lines):
        line = sys.stdin.readline()
        id, ali = re.search( "^(\S+)\s+(.+)", line[:-1]).groups()
        alignment.append( (id, [ali]) )

    while 1:
        line = sys.stdin.readline()
        if not line: break
        for x in range(num_lines):
            line = sys.stdin.readline()            
            alignment[x][1].append( line[:-1] )

    for x in range(num_lines):
        print ">%s\n%s" % (alignment[x][0], re.sub("\s", "", string.join(alignment[x][1], "")))

    ## write footer and output benchmark information.
    E.Stop()

