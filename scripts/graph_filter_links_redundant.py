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
graph_filter_links_redundant.py -
=============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

remove from a sorted list of links those which are redundant.

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

import CGAT.Experiment as E

class Map:

    def __init__(self ):
        (self.mQueryToken, self.mSbjctToken, self.score,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli ) =\
         ("", "", 0, 0, 0, "", 0, 0, "" )
        self.mIsExpanded = False
        self.mMapQuery2Sbjct = None

    def Read( self, line ):
        
        (self.mQueryToken, self.mSbjctToken, self.score,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.score = float(self.score)

    def __str__( self ):

        return string.join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.score,
            self.mQueryFrom, self.mQueryTo, self.mQueryAli,
            self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli)), "\t")

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help"  )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    last_m = None

    for line in options.stdin:
        
        if line[0] == "#": continue

        m = Map()
        m.Read( line )

        if not last_m:
            last_m = m
            continue
        
        if m.mQueryToken != last_m.mQueryToken or \
           m.mSbjctToken != last_m.mSbjctToken or \
           m.mQueryFrom > last_m.mQueryTo:
            print str(last_m)
            last_m = m
        else:
            if m.mQueryTo - m.mQueryFrom > last_m.mQueryTo - last_m.mQueryFrom:
                last_m = m

    print str(last_m)
        
    ## write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )

