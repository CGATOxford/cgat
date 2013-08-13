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
cgat_log2wiki.py - convert a log file into a wiki page
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

   python cgat_log2wiki.py --help

Type::

   python cgat_log2wiki.py --help

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
import CGAT.Experiment as E

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: cgat_log2wiki.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-s", "--start", dest="start", type="string",
                      help="start of section." )

    parser.add_option("-e", "--end", dest="end", type="string",
                      help="end of section." )

    parser.add_option("-l", "--level", dest="level", type="int",
                      help="depth of sections." )

    parser.set_defaults(
        level = 2,
        start = None,
        end = None)

    (options, args) = E.Start( parser )

    if not options.start:
        keep = True
    else:
        keep = False

    section_id1 = 0
    section_id2 = 0
    last_l = None
    
    for line in sys.stdin:

        if options.start and not keep and re.search( options.start, line ):
            keep = True

        if not keep: continue
            
        if options.end and re.search( options.end, line ):
            break

        if line[0] != "#":
            l  = "| " + re.sub( "\t", " | ", line[:-1] ) + " | "
        else:
            
            if last_l:

                if re.match( "# ", line):
                    header = line[2:-1]
                    print last_l % header
                    continue
                
            last_l = None
                
            if re.match( "#----------------------------", line):
                section_id2 += 1
                last_l = "\n---%s subsection %i: %%s" % ("+" * (options.level+1), section_id2)
                l = None
            elif re.match( "#===========================", line):
                section_id1 += 1
                section_id2 = 0
                last_l = "\n---%s section %i: %%s" % ("+" * (options.level), section_id1)
                l = None
            else:
                l = line[:-1]

        if l:
            print l


    E.Stop()
