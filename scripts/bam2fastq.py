################################################################################
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
bam2fastq.py - 
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Convert a BAM file to a FASTQ files.

Usage
-----

Example::

   python bam2fastq.py in.bam out.1.fastq out.2.fastq

This command converts the BAM file in.bam into fastq files containing forward reads (out.1.fastq) and reverse reads 
(out.2.fastq). 

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

This tool converts a BAM file containing paired-end reads into forward and reverse fastq files.

For example::

   python bam2fastq.py in.bam out.1.fastq out.2.fastq 

Code
----

'''

import os
import sys
import re
import optparse
import gzip

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

import pysam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-i", "--test-option", dest="test_option", type="string",
                      help="test option [default=%default]."  )

    parser.set_defaults(
        test_option = "test"
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    assert len(args) == 3, "expected three command line arguments" 

    samfile = pysam.Samfile( args[0], "rb" )
    outstream1 = IOTools.openFile( args[1], "w" )
    outstream2 = IOTools.openFile( args[2], "w" )

    c = E.Counter()
    for read in samfile.fetch():
        c.input += 1
        if read.mate_is_unmapped:
            c.is_unmapped += 1
            if read.is_read1:
                seq1, seq2 = read.seq, "N" * read.qlen 
                qual1, qual2 = read.qual, "B" * read.qlen
            else:
                seq1, seq2 = "N" * read.qlen, read.seq
                qual1, qual2 = "B" * read.qlen, read.qual
        else:
            if read.is_read2: continue
            try:
                mate = samfile.mate( read )
                c.found += 1
            except ValueError, msg:
                mate = None
                c.failed += 1
                
            if mate:
                seq1, seq2 = read.seq, mate.seq
                qual1, qual2 = read.qual, mate.qual
            else:
                seq1, seq2 = read.seq, "N" * read.qlen 
                qual1, qual2 = read.qual, "B" * read.qlen

        c.output += 1
        outstream1.write( "@%s\n%s\n+\n%s\n" % (read.qname, seq1, qual1 ) )
        outstream2.write( "@%s\n%s\n+\n%s\n" % (read.qname, seq2, qual2 ) )
        
    E.info( "%s" % str(c) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
