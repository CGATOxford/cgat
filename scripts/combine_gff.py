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
combine_gff.py - merge overlapping intervals
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

combine gff features as overlapping regions.

Usage
-----

Example::

   python combine_gff.py --help

Type::

   python combine_gff.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse
import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: combine_gff.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--format", dest="format",
                      help="output format.", type="choice", choices=("flat", "full", "first") )

    parser.set_defaults(
        format="flat"
        )

    (options, args) = E.Start( parser )

    last_e = None
    for line in sys.stdin:
        if line[0] == "#" : continue
        if options.format in ("full", "first" ):
            last_e = GTF.Entry()
        else:
            last_e = GTF.Entry()
        last_e.Read(line)
        break
    
    for line in sys.stdin:
        
        if line[0] == "#" : continue

        if options.format in ("full", "first" ):
            e = GTF.Entry()
        else:
            e = GTF.Entry()
        e.Read(line)

        if not GTF.Overlap( last_e, e):
            print str(last_e)
            last_e = e
        else:
            last_e.start = min( last_e.start, e.start )
            last_e.end = max( last_e.end, e.end )
            if options.format == "full":
                last_e.mInfo += " ; " + e.mInfo
            continue
        
    print str(last_e)

    E.Stop()

        
