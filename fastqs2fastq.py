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
fastqs2fastq.py - merge two paired fastq files
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes two fastq files and outputs a single fastq file.

Usage
-----

Example::

   python script_template.py --help

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

This script assumes that the fastq files are sorted by read
id.

Code
----

'''

import os, sys, re, optparse, math, random

import IOTools
import Experiment as E
import Fastq

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ('reconcile' ),
                      help="method to apply [default=%default]."  )

    parser.add_option("-o", "--outfile-pattern", dest="outfile_pattern", type="string",
                      help="pattern for output files [default=%default]."  )

    parser.set_defaults(
        method = "reconcile",
        outfile_pattern = "%i.fastq.gz",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    c = E.Counter()
    
    if options.method == "reconcile":

        def getIds( infile ):
            '''return ids in infile.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l: break
                yield l[0].split()[0]

        def write( outfile, infile, take ):
            '''filter fastq files with ids in take.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l: break
                if (l[0].split()[0]) not in take: continue
                outfile.write("\n".join(l) + "\n" )

        inf1 = IOTools.openFile( fn1 )
        ids1 = set( getIds( inf1 ) )
        inf2 = IOTools.openFile( fn2 )
        ids2 = set( getIds( inf2 ) )
        take = ids1.intersect( ids2 )
        
        E.info( "first pair: %i reads, second pair: %i reads, shared: %i reads" % \
                    ( len(ids1),
                      len(ids2),
                      len(take)) )
 
        with IOTools.openFile( outfile_pattern % 1, "w" ) as outf:
            inf = IOTools.openFile( fn1 )
            E.info( "writing first in pair" )
            write( outf, inf, take )

        with IOTools.openFile( outfile_pattern % 2, "w" ) as outf:
            inf = IOTools.openFile( fn2 )
            E.info( "writing second in pair" )
            write( outf, inf, take )
        
    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
