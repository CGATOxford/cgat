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
set_diff.py - compare contents of several files
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

pairwise comparison of multiple sets. 

The script outputs a tab-separated table with the following fields:

set1, set2: labels of sets that are compared
n1, n2: number of elements in the two sets
union: size of union of both sets
inter: size of the intersection of two sets
unique1: number of unique elements in set1
unique2: number of unique elements in set2

If the options --add-percent is given, the following columns are added:
pinter: intersection / union
punique1: unique1 / n1
punique1: unique2 / n2
pcov1: coverage of set1 by set2 = intersection / n1
pcov2: coverage of set2 by set1 intersection / n2
pcovmax: maximum coverage of either set = max(pcov1, pcov2)

Usage
-----

Example::

   python set_diff.py --help

Type::

   python set_diff.py --help

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
import CGAT.IOTools as IOTools

if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: set_diff.py 2782 2009-09-10 11:40:29Z andreas $" )

    parser.add_option("-p", "--add-percent", dest="add_percent", action="store_true",
                      help="add percentage information to each line." )

    parser.add_option("-t", "--headers", dest="headers", type="string",
                      help="comma separated list of headers. If empty or set to '-', filenames are used." )

    parser.add_option("--skip-header", dest="add_header", action="store_false",
                      help="do not add header to flat format."  )

    parser.add_option("--write-header", dest="write_header", action="store_true",
                      help="write header and exit."  )

    parser.add_option("--with-title", dest="with_title", action="store_true",
                      help="use column titles in input data [%default]."  )

    parser.add_option("--no-title", dest="with_title", action="store_false",
                      help="there are no titles in input data [%default]."  )

    parser.set_defaults(
        add_percent = False,
        percent_format = "%5.2f",
        headers = None,
        add_header = True,
        write_header = False,
        with_title = True,
        )

    (options, args) = E.Start( parser )

    if options.add_header:
        options.stdout.write("set1\tset2\tn1\tn2\tunion\tinter\tunique1\tunique2" )
        if options.add_percent:
            options.stdout.write( "\tpinter\tpunique1\tpunique2\tpcov1\tpcov2\tpcovmax" )
        options.stdout.write("\n")

        if options.write_header:
            sys.exit(0)
            
    if len(args) < 2:
        raise ValueError( "please supply at least two filenames.")

    headers, titles, sets = [], [], []

    if options.headers:
        if options.headers == "-":
            headers=args
        else:
            headers=options.headers.split(",")
            if len(headers) != len(args):
                raise ValueError ("please supply the same number of headers as there are filenames." )

    for f in args:
        if options.with_title:
            title, data = IOTools.readList( open(f,"r"), with_title = options.with_title )
            titles.append( title )
        else:
            data = IOTools.readList( open(f,"r") )
        sets.append( set( data ))
        
    if not headers and titles:
        headers = titles
    else:
        headers = args

    for x in range(len(sets)-1):
        set1=sets[x]

        for y in range(x+1, len(sets)):
            set2=sets[y]
            l1,l2 = len(set1), len(set2)
            options.stdout.write("%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i" % (headers[x],headers[y],
                                                                     l1,l2,
                                                                     len(set1.union(set2)),                                                                   
                                                                     len(set1.intersection(set2)),
                                                                     len(set1.difference(set2)),
                                                                     len(set2.difference(set1))))



            
            if options.add_percent:
                if len(set1) == 0:
                    ri, r1, r2 = 0, 1, 0
                    c1, c2, cm = 1, 0, 0
                elif len(set2) == 0:
                    ri, r1, r2 = 0, 0, 1
                    c1, c2, cm = 0, 1, 0
                else:
                    i = len(set1.intersection(set2))
                    ri, r1, r2 = (
                        i / float(len(set1.union(set2))),
                        len(set1.difference(set2)) / float(l1),
                        len(set2.difference(set1)) / float(l2))
                    c1, c2 = (i / float(l1), i / float(l2))
                    cm = max( c1, c2)

                options.stdout.write( "\t" + ("\t".join([options.percent_format for z in range(6)] )) % (ri, r1, r2, c1, c2, cm))
            
            options.stdout.write("\n")
    
    E.Stop()
