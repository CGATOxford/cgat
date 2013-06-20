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

import os
import sys
import re
import optparse
import math
import random
import itertools

import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq

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

    parser.add_option( "--pair", dest="pair", type="string",
                       help="if data is paired, filename with second pair. "
                       "Implemented for sampling [default=%default]."  )

    parser.add_option( "--outfile-pair", dest="outfile_pair", type="string",
                       help="if data is paired, filename for second pair. "
                       "Implemented for sampling [default=%default]."  )

    parser.add_option( "--uniq", dest="uniq", action="store_true",
                       help="remove duplicate reads (by name) [default=%default]."  )

    parser.add_option( "--apply", dest="apply", type="string",
                       help="apply a filter to fastq file (taking only reads in filename) [default=%default]."  )

    parser.add_option( "--trim3", dest="trim3", type="int",
                       help="trim # bases from 3' end [default=%default]."  )
    parser.add_option( "--sort", dest="sort", action="store_true",
                       help="sort fastq by sequence id [default=%default]."  )


    parser.set_defaults(
        change_format = None,
        guess_format = None,
        sample = None,
        trim3 = None,
        pair = None,
        apply = None,
        uniq = False,
        outfile_pair = None,
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
        
        if options.pair:
            if not options.outfile_pair:
                raise ValueError( "please specify output filename for second pair (--outfile-pair)")

            outfile1 = options.stdout
            outfile2 = IOTools.openFile( options.outfile_pair, "w" )
            
            for record1, record2 in itertools.izip( Fastq.iterate( options.stdin ), Fastq.iterate( IOTools.openFile( options.pair) ) ):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    outfile1.write( "%s\n" % record1 )
                    outfile2.write( "%s\n" % record2 )

        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if random.random() <= sample_threshold:
                c.output += 1
                options.stdout.write( "%s\n" % record )

    elif options.apply:
        ids = set(IOTools.readList( IOTools.openFile( options.apply ) ))
        
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if re.sub(" .*", "", record.identifier).strip() in ids:
                c.output += 1
                options.stdout.write( "%s\n" % record )

    elif options.trim3:
        trim3 = options.trim3
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            record.trim( trim3 )
            options.stdout.write( "%s\n" % record )
            c.output += 1
            
    elif options.uniq:
        keys = set()
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if record.identifier in keys: continue
            else: keys.add( record.identifier )
            options.stdout.write( "%s\n" % record )
            c.output += 1
        
    elif options.sort:
        statement = "paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n'"
        os.system(statement)

    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
