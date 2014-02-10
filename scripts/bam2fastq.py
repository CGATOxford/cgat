'''
bam2fastq.py - output fastq files from a bam-file
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Sequences BAM FASTQ Conversion

Purpose
-------

Convert a BAM file to a FASTQ files. This script ouputs
fastq records from a aligned reads in a bam file.

Usage
-----

Example::

   python bam2fastq.py in.bam out.1.fastq out.2.fastq

This command converts the BAM file in.bam into fastq files containing
forward reads (out.1.fastq) and reverse reads (out.2.fastq).

Type::

   python bam2fastq.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import gzip
import tempfile

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

    parser.set_defaults(
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    assert len(args) == 3, "expected three command line arguments" 

    fastqfile1, fastqfile2 = args[1], args[2]
    
    # only output compressed data
    if not fastqfile1.endswith(".gz"): fastqfile1 += ".gz"
    if not fastqfile2.endswith(".gz"): fastqfile2 += ".gz"

    samfile = pysam.Samfile( args[0], "rb" )

    tmpdir = tempfile.mkdtemp( )
    
    outtemp1 = os.path.join( tmpdir, "pair1.gz" )
    outtemp2 = os.path.join( tmpdir, "pair2.gz" )

    outstream1 = IOTools.openFile( outtemp1, "w" )
    outstream2 = IOTools.openFile( outtemp2, "w" )

    found1, found2 = set(), set()
    read1_qlen, read2_qlen = 0, 0

    c = E.Counter()
    for read in samfile.fetch():
        c.input += 1
        if read.is_read1:
            if read.qname not in found1:
                outstream1.write( "\t".join( (read.qname, read.seq, read.qual) ) + "\n" )
                found1.add(read.qname)
                if not read1_qlen: read1_qlen = read.qlen
                c.output1 += 1
        elif read.is_read2:
            if read.qname not in found2:
                outstream2.write( "\t".join( (read.qname, read.seq, read.qual) ) + "\n" )
                found2.add(read.qname)
                if not read2_qlen: read2_qlen = read.qlen
                c.output2 += 1
            
    for qname in found2.difference( found1):
        outstream1.write( "\t".join( (qname, "N" * read1_qlen, "B" * read1_qlen) ) + "\n" )
        c.extra1 += 1

    for qname in found1.difference( found2):
        outstream2.write( "\t".join( (qname, "N" * read2_qlen, "B" * read2_qlen) ) + "\n" )
        c.extra2 += 1

    E.info( "%s" % str(c) )

    outstream1.close()
    outstream2.close()

    E.info( "sorting fastq files" )
    statement = '''zcat %s 
                   | sort -k1,1 
                   | awk '{printf("@%%s\\n%%s\\n+\\n%%s\\n", $1,$2,$3)}' 
                   | gzip > %s'''
        
    E.run( statement % (outtemp1, fastqfile1) )
    E.run( statement % (outtemp2, fastqfile2) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
