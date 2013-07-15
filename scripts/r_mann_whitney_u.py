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
r_mann_whitney_u.py - Mann-Whitney U test
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Input: two sets of values. These can be either given
as values directly or as categories, which will be mapped to
values.

Usage
-----

Example::

   python r_mann_whitney_u.py --help

Type::

   python r_mann_whitney_u.py --help

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
import rpy2
from rpy2.robjects import r as R
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: r_mann_whitney_u.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-m", "--method", dest="method", type="string",
                       help="method to use [ks=Kolmogorov-Smirnov,mwu=Mann-WhitneyU]")
    parser.add_option( "-a", "--hardcopy", dest="hardcopy", type="string",
                      help="write hardcopy to file.", metavar = "FILE" )
    parser.add_option( "-1", "--infile1", dest="filename_input1", type="string" ,
                       help="input filename for distribution 1.")
    parser.add_option( "-2", "--infile2", dest="filename_input2", type="string" ,
                       help="input filename for distribution 2.")
    parser.add_option( "-p", "--infile-map", dest="filename_input_map", type="string" ,
                       help="input filename for mapping categories to values.")

    parser.set_defaults(
        method = "ks",
        filename_input1 = None,
        filename_input2 = None,
        filename_input_map = None,
        )
    
    (options, args) = E.Start( parser,
                               add_pipe_options = True,
                               add_psql_options = True,)


    map_category2value = {}
    if options.filename_input_map:
        map_category2value = IOTools.ReadMap( open(options.filename_input_map, "r"),
                                              map_functions=(str,float))
    
    values1, errors1 = IOTools.ReadList( open(options.filename_input1, "r"),
                                         map_category=map_category2value )
    values2, errors2 = IOTools.ReadList( open(options.filename_input2, "r"),
                                         map_category=map_category2value )    
    
    E.info( "ninput1=%i, nerrors1=%i, ninput2=%i, nerrors2=%i" % (len(values1), len(errors1),
                                                                  len(values2), len(errors2)) )

    if options.hardcopy:
        R.png(options.hardcopy, width=1024, height=768)

    if options.method == "ks":
        result = R.ks_test( values1, values2 )
    elif options.method == "mwu":
        result = R.wilcox_test( values1, values2, paired=False)

    R.assign("v1", values1)
    R.assign("v2", values2)    

    R.layout(R.matrix((1,2,3,4), 2, 2, byrow = True))
        
    R.boxplot( values1, values2, col=('white','red'), main="Boxplot" )

    R("""qqplot( v1, v2, main ='Quantile-quantile plot' ); lines( c(0,1), c(0,1) );""")

    R("""hist( v1, freq=FALSE, width=0.5, density=10, main='Relative frequency histogram')""")
    R("""hist( v2, freq=FALSE, add=TRUE,   width=0.5, col='red', offset=0.5, density=20, angle=135)""")
    R("""hist( v1, freq=TRUE,  width=0.5, density=10, main='Absolute frequency histogram')""")
    R("""hist( v2, freq=TRUE,  add=TRUE,   width=0.5, col='red', offset=0.5, density=20, angle=135)""")

    print "## Results for %s" % result['method']
    for x in ['p.value', 'statistic', 'alternative', 'method']:
        print x, result[x]


    E.Stop()
