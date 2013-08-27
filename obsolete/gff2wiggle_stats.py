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
gff2wiggle_stats.py - compute coverage of intervals
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script builds a :term:`wiggle` file from a :term:`gff` formatted
file by stacking the intervals in the :term:`gff` formatted file on
top of each other. The value in the :term:`wiggle` file is the height
of the stack.

Usage
-----

Example::

   python gff2wiggle_stats.py --help

Type::

   python gff2wiggle_stats.py --help

for command line help.

Documentation
-------------

Code
----

'''

import sys
import re
import string
import optparse
import time
import os
import glob

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.GTF as GTF

import CGAT.Wiggle as Wiggle
import alignlib

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: gff2wiggle_stats.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("--wiggle-files", dest="wiggle_files", type="string",
                      help="glob expression for wiggle files [%default]."  )

    parser.add_option("--prefix", dest="prefix", type="string",
                      help="prefix to add to contig names before lookup [%default]."  )

    parser.add_option( "-a", "--as-gtf", dest="as_gtf", action="store_true",
                       help="input is in gtf format.")

    parser.add_option( "--test", dest="test", type="int",
                       help="test - stop after # rows of parsing [%default]." )

    parser.add_option( "--with-values", dest="with_values", action="store_true",
                       help="output values in last column [%default]." )

    parser.set_defaults( wiggle_files = "*.data.bz2",
                         prefix = "",
                         with_values = False,
                         test = None )

    (options, args) = E.Start( parser, add_pipe_options = True )

    # open indexed access to wiggles
    wiggle_files = glob.glob( options.wiggle_files )
    if not wiggle_files:
        raise IOError( "could not find wiggle files with '%s'" % options.wiggle_files )

    index = Wiggle.WiggleMultiIndexedAccess( wiggle_files,
                                             keep_open = True,
                                             use_cache = False )

    if options.as_gtf:
        iterator = GTF.flat_gene_iterator( GTF.iterator( sys.stdin) )
        id = "gene_id"
    else:
        iterator = GTF.chunk_iterator( GTF.iterator( sys.stdin ) )
        id = "query"

    ninput, noutput, nskipped = 0, 0, 0

    options.stdout.write( "%s\tnali\t%s" % (id, "\t".join( Stats.DistributionalParameters().getHeaders() )))
    if options.with_values:
        options.stdout.write( "\tvalues" )
    options.stdout.write("\n" )

    for gffs in iterator:

        if options.test and ninput >= options.test:
            break

        ninput += 1
        
        if options.loglevel >= 2:
            for gff in gffs:
                options.stdlog.write( str(gff) + "\n" )

        map_genome2query = alignlib.makeAlignmentBlocks()
        for gff in gffs:
            map_genome2query.addDiagonal( gff.start, gff.end, 0 )

        data = index.get( options.prefix + gff.contig,
                          map_genome2query.getRowFrom(),
                          map_genome2query.getRowTo() )

        values = []
        for x, vv in data:
            for v in vv:
                if map_genome2query.mapRowToCol( x ) >= 0:
                    values.append( v )
                x += 1
        if len(values) == 0:
            nskipped += 1
            continue
            
        noutput += 1

        if options.loglevel >= 2:
            options.stdlog.write( "# %s\n" % ",".join( [ "%5.3f" % v for v in values] ))

        s = Stats.DistributionalParameters( values )
        if options.as_gtf:
            id = gff.gene_id
        else:
            id = gff.mAttributes
        options.stdout.write( "%s\t%i\t%s" % ( id,
                                               map_genome2query.getNumAligned(),
                                               str( s ) ) )
        
        if options.with_values:
            options.stdout.write( "\t%s" % ",".join( [ "%5.3f" % v for v in values] ))
            
        options.stdout.write("\n" )    

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    E.Stop()
