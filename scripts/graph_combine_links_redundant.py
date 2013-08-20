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
graph_combine_links_redundant.py - remove redundant links
=========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

remove redundant links from a sorted graph.

Usage
-----

Example::

   python graph_combine_links_redundant.py --help

Type::

   python graph_combine_links_redundant.py --help

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
import getopt

import CGAT.Experiment as E

param_long_options=["verbose=", "help", "version" ]
param_short_options="v:h"

param_loglevel = 1
##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ( "-h", "--help" ):
            print globals()["__doc__"]
            sys.exit(0)
            
    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(1)

    print E.GetHeader()
    print E.GetParams()

    print "token1\ttoken2\tmin\tmax\tcount\ttotal\tavg"

    last = None
    n = 0
    
    for line in sys.stdin:
        
        if line[0] == "#": continue
        
        data = line[:-1].split("\t")[:3]
        # ignore header lines
        try:
            data[2] = float(data[2])
        except ValueError:
            continue

        if not last:
            last = data
            mmin, mmax, mtotal, n  = data[2], data[2], data[2], 1
            continue
        
        if last[0] != data[0] or last[1] != data[1]:
            print string.join( map( str, (last[0], last[1], mmin, mmax, n, mtotal, mtotal/n)), "\t" )
            last = data
            mmin, mmax, mtotal, n  = data[2], data[2], 0, 0
            
        mtotal += data[2]
        mmin = min(mmin, data[2] )
        mmax = max(mmax, data[2] )
        n += 1

    if last:
        print string.join( map(str, (last[0], last[1], mmin, mmax, n, mtotal, mtotal/n)), "\t" )        
            
            
    print E.GetFooter()        
