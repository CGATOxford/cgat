################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
solexa2stats.py - compute stats from solexa export files
========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script returns the genomic insert size distances between mate-paired reads mapping to the same chromosome.  Hence, these values could in theory be negative if the fragment was so short that sequenced ends overlap.  It assumes it is piped an Illumina .export file::

   zcat whatever.export.gz | $RMAAPATH/map/return_insert_sizes > file.for.further.manipulation

The export file is assumed to be organized with mate pairs next to each other as follows:

#. Read 1234/1
#. Read 1234/2
#. Read 1635/1
#. Read 1635/2
#. etc.

Inner-distances are returned on individual lines as follows:

#. 101
#. 120
#. 115
#. etc.

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import collections

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-i", "--test-option", dest="test_option", type="string",
                      help="test option [default=%default]."  )

    parser.set_defaults(
        test_option = "test"
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    infile = options.stdin
    outfile = options.stdout

    line=infile.readline()
    line2=infile.readline()
    seq_len=len(line.split("\t")[8])
   
    counts = collections.defaultdict( int )

    while 1:

       data = line.split("\t")
       data2 = line2.split("\t")

       try:
          if data[10]==data[10] and 'chr' in data[10]:
             size = abs(int(data[12])-int(data2[12]))-seq_len
             counts[size] += 1
       except (ValueError, IndexError):
          pass

       line=infile.readline()
       line2=infile.readline()

       if not line or not line2: break

    outfile.write( "size\tcounts\n" )
    for key in sorted(counts.keys()):
       outfile.write( "%i\t%i\n" % (key, counts[key]) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


