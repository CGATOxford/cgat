"""
bam2bed.py - convert bam formatted file to bed formatted file
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Intervals BAM BED Conversion

Purpose
-------

Convert BAM files into BED files.

Usage
-----

Example::

   python bam2bed.py in.bam > out.bed

This command converts the BAM file in.bam into a BED file named out.bed.

Type::

   python bam2bed.py --help

for command line help.

Documentation
-------------

This tool converts BAM files into BED files supplying the intervals for each read in the BAM file.  BAM files must
have a corresponding index file ie. example.bam and example.bam.bai

For example::

   samtools view example.bam

   READ1     163     1       13040   15      76M     =       13183   219     ...
   READ1     83      1       13183   7       76M     =       13040   -219    ...
   READ2     147     1       13207   0       76M     =       13120   -163    ...

   python bam2bed.py example.bam 

   1       13039   13115   READ1     15      +
   1       13119   13195   READ2     0       +
   1       13182   13258   READ1     7       -
   1       13206   13282   READ2     0       -

Options
^^^^^^^
+-------------------+--------------------------------------------------------+
|--region, -r       |output read intervals that overlap a specified region   |
+-------------------+--------------------------------------------------------+
|--merge-pairs, -m  |merge paired-end reads and output interval for entire   |
|                   |fragment                                                |
+-------------------+--------------------------------------------------------+
|--max-insert-size  |only merge if insert size is less than specified no. of |
|                   |bases                                                   |
+-------------------+--------------------------------------------------------+
|--min-insert-size  |only merge if insert size is greater than specified     |
|                   |no. of bases                                            |
+-------------------+--------------------------------------------------------+

For example,

To output read intervals that overlap chromosome 1, coordinates 13000-13100::

   samtools view example.bam

   READ1     163     1       13040   15      76M     =       13183   219 ...
   READ2     99      1       13120   0       76M     =       13207   163 ...
   READ1     83      1       13183   7       76M     =       13040   -219...
   READ2     147     1       13207   0       76M     =       13120   -163...

   python bam2bed.py example.bam --region '1:13000:13100'

   1       13039   13115   READ1     15      +

To merge paired-end reads and output fragment interval ie. leftmost mapped base to rightmost mapped base::

   python bam2bed.py example.bam --merge-pairs

   1       13119   13282   READ2     0       +
   1       13039   13258   READ1     7       +
   
Command line options
--------------------

""" 

import os
import sys
import re
import optparse
import pysam

import CGAT.Experiment as E

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2bed
except ImportError:
    import CGAT._bam2bed as _bam2bed

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )


    parser.add_option("-r", "--region", dest="region", type="string",
                      help = "output read intervals that overlap samtools region string [default=%default]. " )

    parser.add_option("-m", "--merge-pairs", dest="merge_pairs", action="store_true",
                      help = "merge paired-ended reads and output interval for entire fragment [default=%default]. " )

    parser.add_option( "--max-insert-size", dest="max_insert_size", type = "int",
                      help = "only merge paired-end reads if they are less than # bases apart. " 
                             " 0 turns off this filter [default=%default]. " )

    parser.add_option( "--min-insert-size", dest="min_insert_size", type = "int",
                       help = "only merge paired-end reads if they are at least # bases apart. "
                              " 0 turns off this filter [default=%default]. " )

    parser.add_option( "--bed-format", dest="bed_format", type = "choice", 
                       choices = ('3','4','5','6'),
                       help = "bed format to output. "
                              " [default=%default]" )

    parser.set_defaults(
        region = None,
        call_peaks = None,
        merge_pairs = None,
        min_insert_size = 0,
        max_insert_size = 0,
        bed_format = '6',
        )

    (options, args) = E.Start( parser, argv = argv )

    if len(args) == 0:
        args.append( "-" )
    
    samfile = pysam.Samfile( args[0], "rb" )

    options.bed_format = int(options.bed_format)

    if options.merge_pairs != None:
        counter = _bam2bed.merge_pairs( samfile, 
                                        options.stdout,
                                        min_insert_size = options.min_insert_size,
                                        max_insert_size = options.max_insert_size,
                                        bed_format = options.bed_format )

        options.stdlog.write( "category\tcounts\n%s\n" % counter.asTable() )

    else:
        if options.region != None:
            if args[0] == "-":
                raise ValueError("can't use region with a file from stdin" )
            it = samfile.fetch( region = options.region )
        else:
            # use until_eof. Files from stdin have no index
            it = samfile.fetch( until_eof = True )

        # more comfortable cigar parsing will
        # come with the next pysam release
        BAM_CMATCH = 0
        BAM_CDEL = 2
        BAM_CREF_SKIP = 3
        take = (BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP)
        outfile = options.stdout

        for read in it:
            if read.is_unmapped: continue

            t = 0
            for op, l in read.cigar:
                if op in take: t += l

            start = read.pos
            if read.is_reverse: strand = "-"
            else: strand = "+"
            # IMS: converted rname to reference name
            outfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %\
                          ( samfile.getrname(read.rname),
                            read.pos,
                            read.pos+t,
                            read.qname,
                            read.mapq,
                            strand) )            

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

