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
fastq2fastq.py - manipulate fastq files
=============================================

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

   python script_template.py --help

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os, sys, re, optparse, math, random

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

    parser.add_option("-f", "--change-format", dest="change_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'integer' ),
                      help="guess quality score format and set quality scores to format [default=%default]."  )

    parser.add_option( "--guess-format", dest="guess_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'integer' ),
                      help="quality score format to assume if ambiguous [default=%default]."  )

    parser.add_option( "--sample", dest="sample", type="float",
                       help="sample a proportion of reads [default=%default]."  )

    parser.add_option( "--trim3", dest="trim3", type="int",
                       help="trim # bases from 3' end [default=%default]."  )

    parser.set_defaults(
        change_format = None,
        guess_format = None,
        sample = None,
        trim3 = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    c = E.Counter()

    if options.change_format:
        for record in Fastq.iterate_convert( options.stdin, 
                                             format = options.change_format,
                                             guess = options.guess_format ):
            c.input += 1
            options.stdout.write( "%s\n" % record )
            c.output += 1

    elif options.sample:
        sample_threshold = min( 1.0, options.sample)
        
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if random.random() <= sample_threshold:
                c.output += 1
                options.stdout.write( "%s\n" % record )

    elif options.trim3:
        trim3 = options.trim3
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            record.trim( trim3 )
            options.stdout.write( "%s\n" % record )
            c.output += 1

    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
