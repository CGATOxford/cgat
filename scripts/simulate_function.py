'''
simulate_function.py - create datapoints for a function
=======================================================

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

   python simulate_function.py --help

Type::

   python simulate_function.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import string
import os
import getopt
import time
import optparse

import CGAT.Experiment as E
from scipy.stats import *

USAGE="""simulate_function.py [options] [infile] < stdin

Create data points for a function. All the scipy stats function
are imported.

Examples:

--function="norm,cdf( x, loc=0.0150, scale=0.1225 )" --xrange=-3,3,0.01
"""

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: simulate_function.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-x", "--xrange", dest="xrange", type="string",
                      help="xrange."  )
    parser.add_option("-y", "--yrange", dest="yrange", type="string",
                      help="yrange."  )
    parser.add_option("-f", "--function", dest="function", type="string",
                      help="function."  )


    parser.set_defaults(
        xtitle = "x",
        ytitle = "y",
        xrange = "0,1,0.1",
        yrange = None,
        function = "x",
        value_format = "%6.4f" )

    (options, args) = E.Start( parser )

    if options.xrange: options.xrange = map(float, options.xrange.split(","))
    if options.yrange: options.yrange = map(float, options.yrange.split(","))

    exec("f = lambda x: %s" % options.function )
        
    options.stdout.write( "x\ty\n" )
    x = options.xrange[0]
    while x <= options.xrange[1]:
        options.stdout.write( "\t".join( map( lambda x: options.value_format % x, 
                                              (x, f(x)) )) + "\n" )
        x += options.xrange[2]

    E.Stop()
    
        
