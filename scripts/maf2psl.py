################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: maf2psl.py 2879 2010-04-06 14:44:34Z andreas $
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
map2psql.py - convert maf formatted file to a psl formatted file
================================================================

:Author: Andreas Heger
:Release: $Id: maf2psl.py 2879 2010-04-06 14:44:34Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

convert a maf file to a psl file.

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
from bx.align import maf
from bx.align.tools import get_components_for_species

import CGAT.Blat as Blat

def threaditer( reader, species ):
    '''iterate over reader and return components for species.'''
    for m in reader:
        components = get_components_for_species( m, species ) 
        if components != None:
            yield components
            
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: maf2psl.py 2879 2010-04-06 14:44:34Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-q", "--query", dest="query", type="string",
                      help="sequence to use for query [default=%default]."  )

    parser.add_option("-t", "--target", dest="target", type="string",
                      help="sequence to use for target [default=%default]."  )

    parser.set_defaults(
        query = None,
        target = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if options.query == None or options.target == None:
        if len(args) != 2:
            raise ValueError( "please supply two sequence identifiers for query and target" )
        options.query, options.target = args
        
    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    reader = maf.Reader(options.stdin)
    
    psl = Blat.Match()
    for cc in threaditer( reader, (options.query, options.target) ):
        
        ninput += 1
        query, target = cc

        # treat identfiers like Hsap.GL000223.1
        try:
            data = query.src.split(".")
            qs, qcontig = data[0], ".".join(data[1:])
        except ValueError, msg:
            raise ValueError( "error: could not parse query %s: msg=%s" % (query.src, msg) )
                              
        try:        
            data = target.src.split(".")
            ts, tcontig = data[0], ".".join(data[1:])
        except ValueError, msg:
            raise ValueError( "error: could not parse target %s: msg=%s" % (target.src, msg) )

        assert qs == options.query
        assert ts == options.target
        psl.mQueryId = qcontig
        psl.mSbjctId = tcontig

        psl.fromPair( query.start, query.src_size, query.strand, query.text.upper(),
                       target.start, target.src_size, target.strand, target.text.upper() )
        
        E.debug( "%s\t%s\t%i\t%i\t%s\t%s" % \
                     (qs, qcontig, query.start, query.src_size, query.strand, query.text ) )
        E.debug( "%s\t%s\t%i\t%i\t%s\t%s" % \
                     (ts, tcontig, target.start, target.src_size, target.strand, target.text ) )
        options.stdout.write( "%s\n" % str(psl) )
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

