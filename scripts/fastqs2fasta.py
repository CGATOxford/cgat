'''
fastqs2fasta.py - interleave two fastq files
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS FASTQ FASTA Conversion

Purpose
-------

This script is used to interleave two :term:`fastq`-formatted files
(paired data) into a single :term:`fasta`-formatted file. Read1 is
followed by read2 in the resultant file.

:term:`fastq` files MUST be sorted by read identifier.

Usage
-----

Example::

   python fastqs2fasta.py --fastq1 in.fastq.1.gz --fastq2 in.fastq.2.gz > out.fasta

Type::

   python fastqs2fasta.py --help

for command line help.


Command line options
--------------------

'''

import os
import sys
import re
import optparse
import itertools
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq

import CGAT.Experiment as E

class PairedReadError(Exception):
    '''
    exception raised when reads aren't paired -
    could be not sorted or files of different lengths
    '''
    pass

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                             usage = globals()["__doc__"] )

    parser.add_option("-a", "--fastq1", dest="fastq1", type="string",
                      help="supply read1 fastq file"  )
    parser.add_option("-b", "--fastq2", dest="fastq2", type="string",
                      help="supply read2 fastq file"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    fastq1 = IOTools.openFile(options.fastq1)
    fastq2 = IOTools.openFile(options.fastq2)

    E.info("iterating over fastq files")
    f1_count = 0
    for f1, f2 in itertools.izip_longest(Fastq.iterate(fastq1), Fastq.iterate(fastq2)):
        if not (f1 and f2) or (not f2 and f1):
            try:
                raise PairedReadError("unpaired reads detected. Are files sorted? are files of equal length?")
            except PairedReadError, e:
                raise PairedReadError(e), None, sys.exc_info()[2]
        else:
            assert f1.identifier.endswith("/1") and f2.identifier.endswith("/2"), "Reads in file 1 must end with /1 and reads in file 2 with /2"
            options.stdout.write(">%s\n%s\n>%s\n%s\n" % (f1.identifier, f1.seq, f2.identifier, f2.seq))
            f1_count += 1

    E.info("output: %i pairs" % f1_count)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

