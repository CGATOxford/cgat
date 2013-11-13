'''
fastq2table.py - compute stats on reads in fastq files
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Sequences FASTQ Annotation

Purpose
-------

This script iterates over a fastq file and outputs
summary statistics for each read.

The output is a tab-delimited text file with the following columns:

+----------+----------------------------------------+
|*Column*  |*Content*                               |
+----------+----------------------------------------+
|read      |read identifier present in input fastq  |
|          |file                                    |
+----------+----------------------------------------+
|nfailed   |number of reads that fall below Q10     |
+----------+----------------------------------------+
|nN        |number of ambiguous base calls (N)      |
+----------+----------------------------------------+
|nval      |number of bases in the read             |
+----------+----------------------------------------+
|min       |minimum base quality score for the read |
+----------+----------------------------------------+
|max       |maximum base quality for the read       |
+----------+----------------------------------------+
|mean      |mean base quality for the read          |
+----------+----------------------------------------+
|median    |median base quality for the read        |
+----------+----------------------------------------+
|stddev    |standard devitation of quality scores   |
|          |for the read                            |
+----------+----------------------------------------+
|sum       |sum of quality scores for the read      |
+----------+----------------------------------------+
|q1        |25th percentile of quality scores for   |
|          |the read                                |
+----------+----------------------------------------+
|q3        |25th percentile of quality scores for   |
|          |the read                                |
+----------+----------------------------------------+

Usage
-----

Example::

   python fastq2table.py --guess-format=sanger < in.fastq > out.tsv

In this example we know that our data have quality scores formatted as sanger. Given that
illumina-1.8 quality scores are highly overlapping with sanger, this option defaults to
sanger qualities. In default mode the script may not be able to distinguish
highly overlapping sets of quality scores.

Type::

   python fastq2table.py --help

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

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.Fastq as Fastq

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "--guess-format", dest="guess_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'illumina-1.8', 'integer' ),
                      help="The default behaviour of the script is to guess the quality format of the input fastq file. The user can specify \
                            the quality format of the input file using the --format option. The script will use this format if the \
                            sequence qualities are ambiguous.[default=%default]."  )

    parser.add_option("-f", "--change-format", dest="change_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'illumina-1.8', 'integer' ),
                      help="The script will guess the quality format of the input file and convert \
                            quality scores to the destination format unless --format is specified [default=%default]."  )

    parser.set_defaults(
        change_format = None,
        guess_format = None,
        min_quality = 10,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    c = E.Counter()
    
    if options.change_format:
        iterator = Fastq.iterate_convert( options.stdin, 
                                          format = options.change_format,
                                          guess = options.guess_format )
    else:
        iterator = Fastq.iterate_guess( options.stdin,
                                        guess = options.guess_format )

    options.stdout.write( "read\tnfailed\tnN\t%s\n" % ("\t".join(Stats.Summary().getHeaders()) ) )

    min_quality = options.min_quality

    for record in iterator:
        c.input += 1
        quals = record.toPhred()
        nfailed = len( [x for x in quals if x < min_quality ] )
        nns = record.seq.count("N") + record.seq.count(".")
        options.stdout.write( "%s\t%i\t%i\t%s\n" % (record.identifier,
                                                    nfailed,
                                                    nns,
                                                    str(Stats.Summary(quals))
                                                    ) )
        c.output += 1


    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
