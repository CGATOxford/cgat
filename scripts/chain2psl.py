################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: chain2psl.py 2899 2010-04-13 14:37:37Z andreas $
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
chain2psl.py - convert a chain file to a psl file
=================================================

:Author: Andreas Heger
:Release: $Id: chain2psl.py 2899 2010-04-13 14:37:37Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

convert a `chain <http://www.breyer.com/ucsc/htdocs/goldenPath/help/chain.html>`_ 
formatted file to a `psl <http://genome.ucsc.edu/FAQ/FAQformat.html#format2>`_
formatted file.

This tool is equivalent to chainToPsl except that it will not compute the number
of matching, mismatching, etc. bases and thus does not require the sequences.

The nomenclature the UCSC uses is for its chain files is :file:`targetToQuery.chain` for 
mapping ``query`` to ``target`` (reference). According to the UCSC documentation, 
``target`` is the first entry in ``chain`` files. 

I have been using the nomenclature ``QueryToTarget.psl``. In following this convention,
the correct way to converting a psl file is::

   python chain2psl.py < targetToQuery.chain > QueryToTarget.psl

If you would like to keep the TargetToQuery convention, you will need to add a pslSwap::

   python chain2psl.py < targetToQuery.chain | pslSwap stdin stdout > targetToQuery.psl

Usage
-----

Type::

   python chain2psl.py --help

for command line help.

Command line options
--------------------

""" 

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.Blat as Blat
import alignlib

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: chain2psl.py 2899 2010-04-13 14:37:37Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    psl = None

    def chain_iterator( infile ):
        lines = []
        for line in options.stdin:
            
            if line.startswith("#"): continue
            if line.strip() == "": continue
            if line.startswith("chain"):
                if lines: yield lines
                lines = []
            lines.append( line )
            
        yield lines

    for lines in chain_iterator(options.stdin):
        
        ninput += 1
        psl = Blat.Match()

        ( _, 
          _, 
          psl.mSbjctId,
          target_length,
          target_strand,
          target_start,
          target_end,
          psl.mQueryId,
          query_length,
          query_strand,
          query_start, 
          query_end,
          alignment_id ) = lines[0][:-1].split()
        
        ( psl.mQueryStart, psl.mQueryEnd, psl.mQueryLength,
          psl.mSbjctStart, psl.mSbjctEnd, psl.mSbjctLength ) = \
        [ int(x) for x in 
          (query_start, 
           query_end,
           query_length,
           target_start, 
           target_end,
           target_length) ]

        map_query2target = alignlib.makeAlignmentBlocks()
        
        qstart, tstart = psl.mQueryStart, psl.mSbjctStart
        
        for line in lines[1:-1]:
            size, dt, dq = [int(x) for x in line[:-1].split() ]
            map_query2target.addDiagonal( qstart,
                                          qstart + size,
                                          tstart - qstart )
            qstart += size + dq
            tstart += size + dt

        size = int(lines[-1][:-1])

        map_query2target.addDiagonal( qstart,
                                      qstart + size,
                                      tstart - qstart )

        psl.fromMap( map_query2target )

        # sort out strand
        # target_strand is always positive
        assert( target_strand == "+" )

        # if query strand is negative
        if query_strand == "-": 
            # invert both query and target
            psl.switchTargetStrand()
            # manually invert the query coordinates
            psl.mQueryFrom, psl.mQueryTo = psl.mQueryLength - psl.mQueryTo, psl.mQueryLength - psl.mQueryFrom

        options.stdout.write("%s\n" % psl )
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
