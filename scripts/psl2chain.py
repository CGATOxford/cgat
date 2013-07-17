################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: psl2chain.py 2901 2010-04-13 14:38:07Z andreas $
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
"""
psl2chain.py - convert psl formatted file to a chain formatted file
===================================================================

:Author: Andreas Heger
:Release: $Id: psl2chain.py 2901 2010-04-13 14:38:07Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

convert a `psl <http://genome.ucsc.edu/FAQ/FAQformat.html#format2>`_
formatted file to a `chain <http://www.breyer.com/ucsc/htdocs/goldenPath/help/chain.html>`_ 
formatted file for use in liftover.

The nomenclature the UCSC uses is :file:`TargetToQuery.chain` for mapping ``target``
to ``query`` (according to the UCSC documentation, ``target`` is the first entry
in ``chain`` files). I have been using the nomenclature ``QueryToTarget.psl``. Hence,
the correct way to converting a psl file is::

   python psl2chain.py < queryToTarget.psl > targetToQuery.chain

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.Blat as Blat

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: psl2chain.py 2901 2010-04-13 14:38:07Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    for psl in Blat.iterator(options.stdin):
        ninput += 1
        if psl.strand == "-":
            qstart, qend = psl.mQueryLength - psl.mQueryTo, psl.mQueryLength - psl.mQueryFrom
        else:
            qstart, qend = psl.mQueryFrom, psl.mQueryTo

        options.stdout.write( "chain %i %s %i %s %i %i %s %i %s %i %i %i\n" % \
                                  (psl.mNMatches, 
                                   psl.mSbjctId,
                                   psl.mSbjctLength,
                                   "+",
                                   psl.mSbjctFrom,
                                   psl.mSbjctTo,
                                   psl.mQueryId,
                                   psl.mQueryLength,
                                   psl.strand,
                                   qstart,
                                   qend,
                                   ninput ) )
            
        size, tend, qend = 0, None, None
        for qstart, tstart, size in psl.getBlocks():
            if tend != None:
                options.stdout.write( "\t%i\t%i\n" % (tstart - tend, qstart - qend) )
            qend, tend = qstart + size, tstart + size
            options.stdout.write( "%i" % (size,))
        options.stdout.write( "\n" )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
