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
r_compare_distributions.py - statistical test for distributions
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares two distributions from two samples
of data. If only one sample is given, the sample can be
compared the normal distribution.

The tests are implemented in ``R``.

Usage
-----

Example::

   python r_compare_distributions.py --help

Type::

   python r_compare_distributions.py --help

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

"""Tests to compare distribution of values

Input: two sets of values. These can be either given
as values directly or as categories, which will be mapped to
values.

Question: is the distribution of values different?

Implementation: uses R 

See http://www.graphpad.com/library/BiostatsSpecial/article_197.htm
for an interesting article about normality tests.


"""

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: r_compare_distributions.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       help="method to use: ks=Kolmogorov-Smirnov, mwu=Mann-WhitneyU, shapiro=Shapiro-Wilk, paired-mwu=paired Mann-WhitneyU, paired-t=paired t-test [default=%default]",
                       choices=("ks", "mwu", "shapiro", "paired-mwu", "paired-t" ) )
    parser.add_option( "-a", "--hardcopy", dest="hardcopy", type="string",
                      help="write hardcopy to file.", metavar = "FILE" )
    parser.add_option( "-1", "--infile1", dest="filename_input1", type="string" ,
                       help="input filename for distribution 1.")
    parser.add_option( "-2", "--infile2", dest="filename_input2", type="string" ,
                       help="input filename for distribution 2.")
    parser.add_option( "--legend", dest="legend", type="string" ,
                       help="legend for histograms.""")
    parser.add_option( "-f", "--infile-map", dest="filename_input_map", type="string" ,
                       help="input filename for mapping categories to values.")
    parser.add_option( "-n", "--norm-test", dest="norm_test", action="store_true",
                       help="""test if a set of values is normally distributed. Mean and variance
                       are calculated from the data.""")
    parser.add_option( "-b", "--num-bins", dest="num_bins", type="int",
                       help="""number of bins (for plotting purposes only).""")
    parser.add_option( "--bin-size", dest="bin_size", type="float",
                       help="""bin size for plot.""")
    parser.add_option( "--min-value", dest="min_value", type="float",
                       help="""minimum_value for plot.""")
    parser.add_option( "--max-value", dest="max_value", type="float",
                       help="""maximum_value for plot.""")
    parser.add_option( "--skip-plot", dest="plot", action="store_false",
                       help="""skipping plotting.""")
    parser.add_option( "--header", dest="header", type="string",
                       help="""header of value column [default=%default].""")
    parser.add_option( "--title", dest="title", type="string",
                       help="""plot title [default=%default].""")


    parser.set_defaults(
        method = "ks",
        filename_input1 = None,
        filename_input2 = None,
        filename_input_map = None,
        legend = None,
        norm_test = False,
        num_bins = 0,
        legend_range = "2,2",
        bin_size = None,
        min_value = None,
        plot = True,
        header = "value",
        title = None,
        )
    
    (options, args) = E.Start( parser,
                               add_pipe_options = True,
                               add_psql_options = True,)

    kwargs = {}
    xargs = []
    for arg in args:
        if "=" in arg:
            key,value = arg.split("=")
            kwargs[key] = value
        else:
            xargs.append( arg )

    if options.legend: options.legend = options.legend.split(",")
    
    map_category2value = {}
    if options.filename_input_map:
        map_category2value = IOTools.ReadMap( open(options.filename_input_map, "r"),
                                              map_functions=(str,float))
        f = str
    else:
        f = float

    if options.filename_input1:
        infile1 = open(options.filename_input1, "r")
    else:
        infile1 = sys.stdin

    values1, errors1 = IOTools.ReadList( infile1,
                                         map_function=f,
                                         map_category=map_category2value )

    if options.filename_input1:
        infile1.close()

    if errors1 and options.loglevel >= 3:
        options.stdlog.write("# errors in input1: %s\n" % ";".join(map(str, errors1)))

    if options.norm_test:
        mean = R.mean(values1)
        stddev = R.sd(values1)
        options.stdlog.write( "# creating %i samples from normal distribution with mean %f and stddev %f\n" % (len(values1), mean, stddev))
        
        values2 = R.rnorm( len(values1), mean, stddev )
        errors2 = ()
    else:
        values2, errors2 = IOTools.ReadList( open(options.filename_input2, "r"),
                                             map_function=f,                                             
                                             map_category=map_category2value )    

    if errors2 and options.loglevel >= 3:
        options.stdlog.write("# errors in input2: %s\n" % ";".join(map(str, errors2)))

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput1=%i, nerrors1=%i, ninput2=%i, nerrors2=%i\n" % (len(values1), len(errors1),
                                                                                        len(values2), len(errors2)) )

    if options.method in ("paired-mwu", "paired-t"):
        if len(values1) != len(values2):
            raise ValueError("number of values must be equal for paired tests.")

    if options.hardcopy:
        R.png(options.hardcopy, width=1024, height=768)

    if options.method == "ks":
        result = R.ks_test( values1, values2, *xargs, **kwargs )
    elif options.method == "mwu":
        result = R.wilcox_test( values1, values2, paired=False, correct = True, *xargs, **kwargs)
    elif options.method == "paired-mwu":
        result = R.wilcox_test( values1, values2, paired=True, correct = True, *xargs, **kwargs)
    elif options.method == "paired-t":
        result = R.t_test( values1, values2, paired=True, *xargs, **kwargs )
    elif options.method == "shapiro":
        if len(values1) > 5000:
            E.warn("shapiro-wilk test only accepts < 5000 values, a random sample has been created." )
            values1 = random.sample( values1, 5000 )
        result = R.shapiro_test( values1, *xargs, **kwargs )

    if options.plot:
        R.assign("v1", values1)
        R.assign("v2", values2)    

        if options.title:
            # set the size of the outer margins - the title needs to be added at the end
            # after plots have been created
            R.par(oma=R.c(0,0,4,0) )                     

        R.layout(R.matrix((1,2,3,4), 2, 2, byrow = True))

        R.boxplot( values1, values2, col=('white','red'), main="Boxplot" )
        R("""qqplot( v1, v2, main ='Quantile-quantile plot' ); lines( c(0,1), c(0,1) );""")

        ## compute breaks:

        min_value = min( min(values1), min(values2) )
        if options.min_value != None:
            min_value = min( min_value, options.min_value)

        max_value = max( max(values1), max(values2) )
        if options.max_value != None:
            max_value = max( max_value, options.max_value)

        extra_options = ""
        if options.num_bins and not (options.min_value or options.max_value):
            extra_options += ", breaks=%i" % options.num_bins

        elif options.num_bins and (options.min_value or options.max_value):
            bin_size = float((max_value - min_value )) / (options.num_bins + 1)
            breaks = [ min_value + x * bin_size for x in range(options.num_bins) ]
            extra_options += ", breaks=c(%s)" % ",".join(map( str, breaks) )

        elif options.bin_size != None:
            num_bins = int(((max_value - min_value) / options.bin_size)) + 1
            breaks = [ min_value + x * options.bin_size for x in range(num_bins+1) ]
            extra_options += ", breaks=c(%s)" % ",".join(map( str, breaks) )

        R("""h1 <- hist( v1, freq=FALSE,           density=20, main='Relative frequency histogram' %s)""" % extra_options )
        R("""h2 <- hist( v2, freq=FALSE, add=TRUE, density=20, col='red', offset=0.5, angle=135 %s)""" % extra_options)
        if options.legend:
            R("""legend( ( max(c(h1$breaks[-1], h2$breaks[-1])) - min(c(h1$breaks[1], h2$breaks[1]) ) ) / 2,
            max( max(h1$density), max(h2$density)) / 2, c('%s'), fill=c('white','red'))""" % (
                "','".join(options.legend)) )

        R("""h1 <- hist( v1, freq=TRUE,            density=20, main='Absolute frequency histogram' %s)""" % extra_options )
        R("""h2 <- hist( v2, freq=TRUE,  add=TRUE, density=20, col='red', offset=0.5, angle=135 %s )""" % extra_options)
        if options.legend:    
            R("""legend( ( max(c(h1$breaks[-1], h2$breaks[-1])) - min(c(h1$breaks[1], h2$breaks[1]) ) ) / 2,
            max( max(h1$counts), max(h2$counts)) / 2, c('%s'), fill=c('white','red'))""" % (
                "','".join(options.legend)) )

        if options.title:
            R.mtext( options.title,3,outer=True,line=1,cex=1.5)  

    if options.loglevel >= 1:
        options.stdout.write("## Results for %s\n" % result['method'])

    options.stdout.write( "%s\t%s\n" % ("key", options.header) )

    for key in result.keys():
        if key == "data.name": continue
        options.stdout.write("\t".join( (key, str(result[key]))) + "\n")

    stat = Stats.Summary(values1)
    for key, value in stat.items():
        options.stdout.write( "%s1\t%s\n" % (str(key), str(value)))

    stat = Stats.Summary(values2)
    for key, value in stat.items():
        options.stdout.write( "%s2\t%s\n" % (str(key), str(value)))

    if options.plot:
        if options.hardcopy:
            R.dev_off()

    E.Stop()
