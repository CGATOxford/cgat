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
rates2rates.py - operate on rates
=================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

input = tab separated table ids, rates and G+C

This script corrects a substitution rate for mutational
biases. The mutational bias is estimated from a given
distribution of rates and G+C content in neutrally 
evolving regions, for example ancestral repeats.

Usage
-----

Example::

   python rates2rates.py --help

Type::

   python rates2rates.py --help

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
import time
import optparse
import math
import tempfile
import CGAT.Experiment as E
import CGAT.CSV as CSV
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats
import numpy

def readRates( infile ):
    """read rates and G+C from a tab-separated file."""
    from rpy import r as R
    import rpy

    handle, name = tempfile.mkstemp()
    os.close(handle)
    outfile = open(name, "w")

    first = True
    headers = []
    for line in infile:
        if line[0] == "#": continue
        data = line[:-1].split("\t")
        if first:
            headers = data
            first = False
            continue

        outfile.write( line )
    outfile.close()

    assert len(headers) == 3, "malformatted file of rates, please supply id, g+c, rate"

    rpy.set_default_mode(rpy.NO_CONVERSION)
    matrix = R.read_table( name, na_string = ("NA", "na"), col_names=headers ) 
    rpy.set_default_mode(rpy.BASIC_CONVERSION)
    os.remove( name )
    return matrix, headers

def main():

    parser = E.OptionParser( version = "%prog version: $Id: rates2rates.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="pattern for additional output files [%default]."  )

    parser.add_option( "--input-filename-neutral", dest="input_filename_neutral", type="string",
                      help="a tab-separated file with rates and G+C content in neutrally evolving regions [%default]."  )

    parser.set_defaults(
        input_filename_neutral = None,
        output_filename_pattern = "%s",
        normalize = True,
        hardcopy = None,
        )

    (options, args) = E.Start( parser, add_csv_options = True )

    if not options.input_filename_neutral:
        raise ValueError( "please supply a file with neutral rates." )

    lines = options.stdin.readlines()
    if len(lines) == 0:
        raise IOError ( "no input" )

    from rpy import r as R
    import rpy

    R.png( options.output_filename_pattern % "fit" + ".png", width=1024, height=768, type="cairo")
    matrix, headers = readRates( open( options.input_filename_neutral, "r" ) )
    R.assign("matrix", matrix)
    R.assign("headers", headers)
    nref = R( """length( matrix[,1] )""" )

    dat = R("""dat <- data.frame(x = matrix[,2], y = matrix[,3])""")
    mod = R("""mod <- lm( y ~ x, dat)""")

    R("""plot( matrix[,2], matrix[,3], cex=%s, col="blue", pch="o", xlab="%s", ylab="%s" %s)""" % (0.3, headers[1], headers[2], "") )
    R("""new <- data.frame(x = seq( min(matrix[,2]), max(matrix[,2]), (max(matrix[,2]) - min(matrix[,2])) / 100))""")
    R("""predict(mod, new, se.fit = TRUE)""")
    R("""pred.w.plim <- predict(mod, new, interval="prediction")""")
    R("""pred.w.clim <- predict(mod, new, interval="confidence")""")
    R("""matpoints(new$x,cbind(pred.w.clim, pred.w.plim[,-1]), lty=c(1,2,2,3,3), type="l")""")
    R.mtext(
        "y = %f * x + %f, r=%6.4f, n=%i" % (mod["coefficients"]["x"], 
                                            mod["coefficients"]["(Intercept)"], 
                                            R("""cor( dat )[2]"""), 
                                            nref ),
        3,
        cex = 1.0)

    R("""mean_rate <- mean( matrix[,3] )""")

    data_matrix, data_headers = readRates( lines )
    R.assign("data_matrix", data_matrix)
    R.assign("data_headers", data_headers)
    ndata = R( """length( data_matrix[,1] )""" )
    
    R("""points( data_matrix[,2], data_matrix[,3], cex=%s, col="red", pch="o" %s)""" % (0.3, "") )
    R("""topred <- data.frame( x = data_matrix[,2] )""")
    R("""corrected_rates <- predict( mod, topred, se.fit = TRUE )""")
    uncorrected = R("""uncorrected <- data_matrix[,3] / mean_rate """) 
    corrected = R("""corrected <- as.vector(data_matrix[,3] / corrected_rates$fit)""")
    R.dev_off()
    
    R.png( options.output_filename_pattern % "correction" + ".png", width=1024, height=768, type="cairo")
    R("""plot( uncorrected, corrected, cex=%s, col="blue", pch="o", xlab="uncorrected rate", ylab="corrected rate" %s)""" % (0.3, "") )
    R.dev_off()

    E.Stop()

if __name__ == "__main__":
    sys.exit( main() )
