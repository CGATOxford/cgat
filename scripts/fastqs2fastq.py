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

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

This script assumes that the fastq files are sorted by read
id.

Code
----

'''

import os
import sys
import re
import optparse
import math
import random

import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ('reconcile', ),
                      help="method to apply [default=%default]."  )
                      
    parser.add_option("-c", "--chop", dest="chop", action="store_true",
                      help="whether or not to trim last character of sequence name (sometimes used to indicte end) [default=%default]."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="pattern for output files [default=%default]."  )

    parser.set_defaults(
        method = "reconcile",
        chop = False,
        output_pattern = "%i.fastq.gz",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError("please supply at least two fastq files on the commandline")

    fn1, fn2 = args
    c = E.Counter()
    
    if options.method == "reconcile":

        def getIds( infile ):
            '''return ids in infile.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]: break
                r = l[0].split()[0]
                # decide if to chop read number off
                if options.chop: yield r[:-1]
                else: yield r

        def write( outfile, infile, take ):
            '''filter fastq files with ids in take.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]: break
                r = l[0].split()[0]
                if options.chop: r = r[:-1]
                if r not in take: continue
                outfile.write("\n".join(l) + "\n" )

        E.info( "reading first in pair" )
        inf1 = IOTools.openFile( fn1 )
        ids1 = set( getIds( inf1 ) )

        E.info( "reading second in pair" )
        inf2 = IOTools.openFile( fn2 )
        ids2 = set( getIds( inf2 ) )

        take = ids1.intersection( ids2 )
        
        E.info( "first pair: %i reads, second pair: %i reads, shared: %i reads" % \
                    ( len(ids1),
                      len(ids2),
                      len(take)) )
 
        with IOTools.openFile( options.output_pattern % 1, "w" ) as outf:
            inf = IOTools.openFile( fn1 )
            E.info( "writing first in pair" )
            write( outf, inf, take )

        with IOTools.openFile( options.output_pattern % 2, "w" ) as outf:
            inf = IOTools.openFile( fn2 )
            E.info( "writing second in pair" )
            write( outf, inf, take )
        
    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
