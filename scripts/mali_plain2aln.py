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
mali_plain2aln.py - 
======================================================

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

   python mali_plain2aln.py --help

Type::

   python mali_plain2aln.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import os
import string
import optparse

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # Length of each alignment segment
    LINE_WIDTH = 60

    # Get file data from stdin
    aliLines = sys.stdin.read().split("\n")
    while "" in aliLines: aliLines.remove("")
    while len(aliLines) > 2: del aliLines[2]

    # Sort data and get alignment length
    qArray = aliLines[0].split("\t")
    queryStart, queryAli, queryEnd, queryId = int(qArray[0]), qArray[1], int(qArray[2]), qArray[3]
    queryAli = queryAli.replace("-",".")
    sArray = aliLines[1].split("\t")
    sbjctStart, sbjctAli, sbjctEnd, sbjctId = int(sArray[0]), sArray[1], int(sArray[2]), sArray[3]
    sbjctAli = sbjctAli.replace("-",".")
    length = len(queryAli)

    CURR_POS = 0
    qFrom, qTo, sFrom, sTo = queryStart, 0, sbjctStart, 0
    while (CURR_POS < length):
	    if length - CURR_POS < 60: LINE_WIDTH = length - CURR_POS
	    qGaps, sGaps = 0, 0
	    # Get alignment slices
	    qAli = queryAli[CURR_POS:CURR_POS+LINE_WIDTH]
	    sAli = sbjctAli[CURR_POS:CURR_POS+LINE_WIDTH]
	    # Count gaps
	    for x in range(0,len(qAli),1):
		    if qAli[x] == ".": qGaps += 1
		    if sAli[x] == ".": sGaps += 1
	    qTo = qFrom + (LINE_WIDTH-1) - qGaps
	    sTo = sFrom + (LINE_WIDTH-1) - sGaps
	    print "\n"
	    print "%-20s %-10d%s\t%d" %(queryId, qFrom, qAli, qTo)
	    print "%-20s %-10d%s\t%d" %(sbjctId, sFrom, sAli, sTo)

	    qFrom = qTo+1
	    sFrom = sTo+1
	    CURR_POS += LINE_WIDTH

    print "\n"

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

