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
nucmer2bed.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes as input the output from nucmer 'show-coords' and outputs
a bed file with alignment intervals.

Usage
-----

Example::

   python nucmer2bed.py --help

Type::

   python nucmer2bed.py --help

for command line help.

Documentation
-------------

Read from stdin and output to stdout. Requires that the -T option in 
show-coords was used to produce tabular output.

Code
----

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-t", "--type", dest="type", type="choice",
                      choices = ("bed", "bed4"), help="output bed type: bed or bed4"  )
    
    parser.set_defaults(type="bed")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )



    # read in header lines
    inf = sys.stdin
    input_files = inf.readline()
    program = inf.readline()
    spacer = inf.readline()
    header = inf.readline()
    
    # contig is the reference tag


    # get the coordinates for the alignment with respect to
    # the reference
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contig, start, end = data[-2], data[0], data[1]
        if options.type == "bed4":
            name = data[-1]
            options.stdout.write("\t".join(map(str,[contig, start, end, name])) + "\n")
        else:
            options.stdout.write("\t".join(map(str,[contig, start, end])) + "\n")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


#################################################################
# convert nucmer output coords file to tab delimited text file
#################################################################
import sys
import os
import re

inf = sys.stdin

input = inf.readline()
program = inf.readline()
spacer = inf.readline()

header = inf.readline()
headers = "\t".join(["s1", "e1", "s2", "e2", "len1", "len2", "%id", "len_r", "len_q", "ref_tag", "query_tag"])
sys.stdout.write("%s\n" % headers)
equals_spacer = inf.readline()
for line in inf.readlines():
    data = [x for x in line[:-1].split(" ") if x and x != "|"]
    sys.stdout.write("\t".join(map(str,data)) + "\n")
