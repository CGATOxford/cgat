################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
coverage2stats.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The output from bedtools genomecov is converted to a set of statistics
for each contig that is output. Statistics output are the mean and sd
of coverage for each contig.

Usage
-----

Example::

   python coverage2stats.py --help

Type::

   python coverage2stats.py --help

for command line help.

Documentation
-------------

The script assumes the input is from stdin and outputs the results to stdout.
The input is expected to be in tab delimited text format (no header line)

          <contig><base-position><coverage>

Code
----

'''

import os
import sys
import re
import optparse
import numpy as np
import collections
import itertools

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    inf = options.stdin
    
    coverage_result = collections.defaultdict(list)
    E.info("reading in coverage data")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contig, coverage = data[0], data[2]
        coverage_result[contig].append(coverage)
    E.info("read %i contigs" % len(coverage_result.keys()))
        
    options.stdout.write("contig\tcov_mean\tcov_sd\n")
    for contig, coverage in coverage_result.iteritems():
        coverage = map(int, coverage)
        options.stdout.write("%s\t%s\t%s\n" % (contig, str(np.mean(coverage)), str(np.std(coverage))))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
