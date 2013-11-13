"""
bam2wiggle.py - convert bam to wig/bigwig file
==============================================

:Author: Andreas Heger
:Release: $Id: bam2wiggle.py 2832 2009-11-24 16:11:06Z andreas $
:Date: |today|
:Tags: Genomics NGS Intervals Conversion BAM WIGGLE

Purpose
-------

convert a bam file to a bigwig or bedgraph file.   

The script requires the executables :file:`wigToBigWig` and 
:file:`bedToBigBed` to be in the user's PATH.

If no --shift or --extend option are given, the coverage is computed directly on reads. 
Counting can be performed at a certain resolution.
The counting currently is not aware of spliced reads, i.e., an inserted intron will be 
included in the coverage.

If --shift or --extend are given, the coverage is computed by shifting read
alignment positions upstream for positive strand reads or downstream for negative
strand reads and extend them by a fixed amount. 

For RNASEQ data it might be best to run genomeCoverageBed directly on the 
bam file.

Usage
-----

Type::

   python bam2wiggle.py --output-format=bigwig in.bam out.bigwig

to convert the :term:`bam` file file:`in.bam` to :term:`bigwig` format 
and save the result in :file:`out.bigwig`.

Type::

   python bam2wiggle.py --help

for command line help.

Command line options
--------------------

""" 

import os
import sys
import re
import optparse
import itertools
import tempfile
import shutil
import subprocess

import CGAT.Experiment as E
import pysam
import CGAT.IOTools as IOTools
# for merging pairs
import bam2bed

class SpanWriter(object):
    '''output values within spans.
    
    values are collected according to a span and an average is
    output.
    '''
    def __init__(self,span):
        self.span = span

        self.laststart = 0
        self.lastend = None
        self.val = 0

        self.lastout = None

    def __call__( self, outfile, contig, start, end, val ):
        
        # deal with previous window
        if self.lastend:
            last_window_start = self.lastend - self.lastend % self.span
            last_window_end = last_window_start + self.span
        else:
            last_window_start = start - start % self.span
            last_window_end = start + self.span

        # print start, end, val, last_window_start, last_window_end

        if self.lastend and start > last_window_end:
            # no overlap, output previous span
            assert self.lastout != last_window_start, \
                "start=%i, end=%i, laststart=%i, lastend=%i, last_window_start=%i" % \
                (start, end, self.laststart, self.lastend, last_window_start)
            self.lastout = last_window_start
            v = self.val / float(self.span)
            outfile.write( "%i\t%f\n" % (last_window_start, v ) )
            self.val = 0

            last_window_start = start - start % self.span
            last_window_end = start + self.span

        if end < last_window_end:
            # window too small to output, simply add values
            self.val += val * (end - start)
            self.lastend = max( end, self.lastend )
        else:
            # output first window
            v = self.val + val * (self.span - start % self.span) / float(self.span)
            
            s = last_window_start
            assert self.lastout != s, \
                    "start=%i, end=%i, laststart=%i, lastend=%i, s=%i" % \
                    (start, end, self.laststart, self.lastend, s)
            outfile.write( "%i\t%f\n" % (s, v ) )
            self.lastout = s
            self.val = 0

            # Output middle windows
            for x in range( start + self.span - start % self.span, end - self.span, self.span ):
                assert self.lastout != x
                outfile.write( "%i\t%f\n" % (x, val ) )
                self.lastout = x

            if end % self.span:
                # save rest
                self.lastend = end 
                self.val = val * end % self.span
            elif end - self.span != last_window_start:
                # special case, end ends on window
                assert self.lastout != end - self.span, \
                    "start=%i, end=%i, laststart=%i, lastend=%i" % \
                    (start, end, self.laststart, self.lastend)
                self.lastout = end - self.span
                outfile.write( "%i\t%f\n" % (end - self.span, val ) )
                self.lastend = None
            else:
                # special case, end ends on window and only single window - already
                # output as start
                self.lastend = None


    def flush( self, outfile ):
        if self.lastend:
            outfile.write( "%i\t%f\n" % (self.lastend - self.lastend % self.span, 
                                         self.val / (self.lastend % self.span) ) )

def main( argv = None ):
    """script main.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: bam2wiggle.py 2832 2009-11-24 16:11:06Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("bedgraph", "wiggle", "bigbed", "bigwig", "bed"),
                      help="output format [default=%default]" )

    parser.add_option("-b", "--output-filename", dest="output_filename", type="string",
                      help="filename for output [default=%default]" )

    parser.add_option( "-s", "--shift", dest="shift", type = "int",
                       help = "shift reads by a certain amount (ChIP-Seq) [%default]" )

    parser.add_option( "-e", "--extend", dest="extend", type = "int",
                       help = "extend reads by a certain amount (ChIP-Seq) [%default]" )

    parser.add_option( "-p", "--span", dest="span", type = "int",
                       help = "span of a window in wiggle tracks [%default]" )

    parser.add_option("-m", "--merge-pairs", dest="merge_pairs", action="store_true",
                      help="merge paired-ended reads into a single bed interval [default=%default]. " )

    parser.add_option( "--max-insert-size", dest="max_insert_size", type = "int",
                      help = "only merge if insert size less that # bases. 0 turns of this filter [default=%default]."  )

    parser.add_option( "--min-insert-size", dest="min_insert_size", type = "int",
                       help = "only merge paired-end reads if they are at least # bases apart. "
                              " 0 turns of this filter. [default=%default]" )
    parser.set_defaults(
        samfile = None,
        output_format = "wiggle",
        output_filename = None,
        shift = 0,
        extend = 0,
        span = 1,
        merge_pairs = None,
        min_insert_size = 0,
        max_insert_size = 0,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) >= 1:
        options.samfile = args[0]

    if not options.samfile:
        raise ValueError("please provide a bam file")

    if len(args) == 2:
        options.output_filename = args[1]

    samfile = pysam.Samfile( options.samfile, "rb" )

    contig_sizes = dict( zip( samfile.references, samfile.lengths) )

    if options.shift or options.extend:
        if options.output_format != "bigwig":
            raise ValueError( "shift and extend only available for bigwig output" )

    if options.output_format in ("bigwig", "bigbed"):
        
        if not options.output_filename:
            raise ValueError("please output file for bigwig/bigbed computation.")

        if options.output_format == "bigwig":
            executable_name = "wigToBigWig"
        elif options.output_format == "bigbed":
            executable_name = "bedToBigBed"
        else:
            raise ValueError("unknown output format `%s`" % options.output_format)

        executable = IOTools.which( executable_name )

        if not executable:
            raise OSError( "could not find %s in path." % executable_name )

        tmpdir = tempfile.mkdtemp()
        E.debug( "temporary files are in %s" % tmpdir)

        tmpfile_wig = os.path.join( tmpdir, "wig" )
        tmpfile_sizes = os.path.join( tmpdir, "sizes" )

        # write contig sizes
        outfile_size = open( tmpfile_sizes, "w")
        for contig, size in contig_sizes.items():
            outfile_size.write("%s\t%s\n" % (contig, size) )
        outfile_size.close()    
        
        outfile = open( tmpfile_wig, "w" )
        E.info( "starting output to %s" % tmpfile_wig )

    else:
        outfile = options.stdout
        E.info( "starting output to stdout" )        

    if options.output_format in ("wiggle", "bigwig"):
        # wiggle is one-based, so add 1, also step-size is 1, so need to output all bases
        if options.span == 1:
            outf = lambda outfile, contig, start, end, val: \
                outfile.write( "".join( [ "%i\t%i\n" % (x, val) for x in xrange(start+1,end+1) ] ) )
        else:
            outf = SpanWriter( options.span )

    elif options.output_format in ("bed", "bigbed"):
        # bed is 0-based, open-closed
        outf = lambda outfile, contig, start, end, val: \
            outfile.write("%s\t%i\t%i\t%i\n" % (contig, start, end,val))

    ninput, nskipped, ncontigs = 0, 0, 0

    output_filename = options.output_filename 
    if output_filename:
        output_filename = os.path.abspath( output_filename )

    if options.shift > 0 or options.extend > 0 or options.merge_pairs:

        if options.merge_pairs:
            E.info( "merging pairs to temporary file" )
            counter = bam2bed._bam2bed.merge_pairs( samfile, 
                                                    outfile,
                                                    min_insert_size = options.min_insert_size,
                                                    max_insert_size = options.max_insert_size,
                                                    bed_format = 3)
        else:
            # create bed file with shifted tags
            shift, extend = options.shift, options.extend
            shift_extend = shift + extend

            for contig in samfile.references:
                E.debug("output for %s" % contig )
                lcontig = contig_sizes[contig]

                for read in samfile.fetch( contig ):
                    pos = read.pos
                    if read.is_reverse:
                        start = max(0, read.pos + read.alen - shift_extend )
                    else: 
                        start = max(0, read.pos + shift)

                    # intervals extending beyond contig are removed
                    if start >= lcontig: continue

                    end = min( lcontig, start + extend )
                    outfile.write( "%s\t%i\t%i\n" % (contig, start, end))
                
        outfile.close()

        tmpfile_bed = os.path.join( tmpdir, "bed" )
        E.info("computing coverage")
        # calculate coverage - format is bedgraph
        statement = "genomeCoverageBed -bg -i %(tmpfile_wig)s -g %(tmpfile_sizes)s > %(tmpfile_bed)s" % locals()
        E.run( statement )
        
        E.info("converting to bigwig" )

        tmpfile_sorted = os.path.join( tmpdir, "sorted" )
        statement = "sort -k 1,1 -k2,2n %(tmpfile_bed)s > %(tmpfile_sorted)s; bedGraphToBigWig %(tmpfile_sorted)s %(tmpfile_sizes)s %(output_filename)s" % locals()
        E.run( statement)

        shutil.rmtree( tmpdir )

    else:

        def column_iter( iterator ):
            
            start = None
            end = 0
            n = None
            for t in iterator:
                if t.pos - end > 1 or n != t.n: 
                    if start != None: yield start, end, n
                    start = t.pos
                    end = t.pos
                    n = t.n
                end = t.pos
            yield start, end, n

        for contig in samfile.references:
            # if contig != "chrX": continue
            E.debug("output for %s" % contig )
            lcontig = contig_sizes[contig]

            if options.output_format in ("wiggle", "bigwig"):
                outfile.write( "variableStep chrom=%s span=%i\n" % (contig, options.span) )
                
            for start, end, val in column_iter( samfile.pileup(contig)):

                # patch: there was a problem with bam files and reads overextending at the end.
                # These are usually Ns, but need to check as otherwise wigToBigWig fails.
                if lcontig <= end: 
                    E.warn( "read extending beyond contig: %s: %i > %i" % (contig, end, lcontig))
                    end = lcontig 
                    if start >= end: continue

                if val > 0: 
                    outf( outfile, contig, start, end, val)
            ncontigs += 1

        if type(outf) == type(SpanWriter):
            outf.flush(outfile)

        E.info( "finished output" )

        E.info( "ninput=%i, ncontigs=%i, nskipped=%i" % (ninput, ncontigs, nskipped) )

        if options.output_format in ("bigwig", "bigbed"):
            outfile.close()

            E.info( "starting %s conversion" % executable )
            try:
                retcode = subprocess.call( " ".join( (executable,
                                                      tmpfile_wig,
                                                      tmpfile_sizes,
                                                      output_filename )),
                                           shell=True)
                if retcode != 0:
                    E.warn( "%s terminated with signal: %i" % (executable, -retcode))
                    return -retcode
            except OSError, msg:
                E.warn( "Error while executing bigwig: %s" % e)
                return 1

            shutil.rmtree( tmpdir )

            E.info( "finished bigwig conversion" )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )

