'''
fastqs2fastqs.py - reconcile reads from a pair of fastq files
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS FASTQ FASTQ Manipulation

Purpose
-------

This script takes two fastq files and outputs two fastq files. The aim
of the script is to reconcile reads from a pair of fastq files. For example,
if reads from each file are independently filtered based on quality scores,
the read representation will vary across files. This script will output reads
that are common to both files.

The two files must be sorted by read identifier.

Usage
-----

Example::

   python fastqs2fastqs.py myReads.1.fastq.gz myReads.2.fastq.gz --reconcile --output-pattern myReads_reconciled.%i.fastq

In this example we take a pair of fastq files, reconcile by read identifier and output 2 new
fastq files named myReads_reconciled.1.fastq.gz myReads_reconciled.2.fastq.gz

Type::

   python fastqs2fastqs.py --help

for command line help.

Command line options
--------------------

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
                      help="whether or not to trim last character of sequence name. For example sometimes ids in the first \
                            file in the pair will end with \1 and teh second with \2. If --chop is not specified \
                            then the results will be wrong [default=%default]."  )

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
