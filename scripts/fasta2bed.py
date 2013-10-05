"""
fasta2bed.py - segment sequences
================================

:Author: Andreas Heger
:Release: $Id: fasta2bed.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Genomics Sequences Intervals FASTA BED Conversion

Purpose
-------

This script takes a genomic sequence in :term:`fasta` format
and applies various segmentation algorithms.

The methods implemented (``--methods``) are:

GCProfile
   use `GCProfile <http://tubic.tju.edu.cn/GC-Profile/>`_ to segment by
   G+C content. GCProfile needs to be installed and in the PATH.

cpg
   output all locations of cpg in the genome

fixed-width-windows-gc
   output fixed width windows of a certain size adding their
   G+C content as score

gaps
   ouput all locations of assembly gaps (blocks of `N`)
   in the genomic sequences

ungapped
   output ungapped locations in the genomic sequences

Usage
-----

Type::

   python fasta2bed.py --method=gap < in.fasta > out.bed


Type::

   python fasta2bed.py --help

for command line help.

Command line options
--------------------

""" 

import os
import sys
import re
import optparse
import tempfile
import subprocess
import glob
import shutil

import CGAT.Experiment as E

EXECUTABLE="GCProfile"

import CGAT.FastaIterator as FastaIterator

def segmentWithCpG( infile, options ):
    '''segment a fasta file, output locations of CpG.'''

    ninput, nskipped, noutput = 0, 0, 0
    
    iterator = FastaIterator.FastaIterator( infile )

    segments = []
    
    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        if cur_record is None: break
        ninput += 1
        contig = re.sub("\s.*", "", cur_record.title )
        last = None
        for pos, this in enumerate( cur_record.sequence.upper()):
            if last == "C" and this == "G":
                segments.append( (contig, pos - 1, pos + 1, 1.0))
            last = this

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    return segments

def segmentWithGCProfile( infile, options ):
    '''segment a fasta file with GCProfile.'''

    ninput, nskipped, noutput = 0, 0, 0
    
    iterator = FastaIterator.FastaIterator( infile )

    tmpdir = tempfile.mkdtemp()
    E.info( "working in %s" % tmpdir )
    
    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        if cur_record is None: break
        ninput += 1
        
        contig = re.sub("\s.*", "", cur_record.title )
        filename = os.path.join( tmpdir, contig) + ".fasta" 
        with open( filename, "w") as f:
            f.write( ">%s\n%s\n" % (contig, cur_record.sequence ) )
        
        E.info( "running %s on %s" % (EXECUTABLE, contig ) )

        statement = "%s %s -t %i -g %f -i %i " % \
                    (EXECUTABLE, filename,
                     options.gcprofile_halting_parameter,
                     options.gcprofile_gap_size,
                     options.gcprofile_min_length )
        
        E.debug( "statement = %s" % statement )

        process = subprocess.Popen(  statement,
                                     cwd = tmpdir,
                                     shell = True,
                                     stdin = None,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE )


        stdout, stderr = process.communicate()

        E.debug( "GProfile: stdout=%s" % stdout )
        E.debug( "GProfile: stderr=%s" % stderr )
                
        # process.stdin.close()
        if process.returncode != 0:
            raise OSError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))

        noutput += 1

    segments = []

    for f in glob.glob( os.path.join( tmpdir, "*.SegGC" ) ):
        contig = os.path.basename(f[:-len(".SegGC")])
        with open(f,"r") as infile:
            start = None
            for line in infile:
                pos, gc = line[:-1].split("\t")
                try: gc = float(gc)
                except ValueError: gc = None
                end = int(pos)
                # ignore the double lines
                if gc != None and start != None and end - start > 2:
                    segments.append( (contig, start, end, gc ) )
                start = int(pos)-1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    # shutil.rmtree( tmpdir )
    
    return segments

def segmentFixedWidthWindows( infile, window_size ):
    """return a list of fixed contig sizes."""

    ninput, nskipped, noutput = 0, 0, 0
    
    iterator = FastaIterator.FastaIterator( infile )
    window_increment = window_size
    # at most 50% can be gap
    gap_cutoff = int( window_size // 2 )
    segments = []
            
    while 1:
        ninput += 1
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        if cur_record is None: break
        contig = re.sub("\s.*", "", cur_record.title )
        seq = cur_record.sequence
        size = len(cur_record.sequence)
        
        for x in range(0, size, window_increment):
            s = seq[x:x+window_size].upper()
            gc, at = 0,0
            for c in s:
                if c in "GC": gc += 1
                elif c in "AT": at += 1

            # skip segments containing mostly gaps
            if window_size - (gc + at) > gap_cutoff:
                nskipped += 1
                continue
            
            segments.append( (contig, x, x+window_size, float(gc) / (gc+at) ) )
        noutput += 1
        
    E.info( "ninput=%i, noutput=%i, nskipped_windows=%i" % (ninput, noutput,nskipped) )

    return segments

def gapped_regions( seq, gap_chars ):
    '''iterator yielding gapped regions in seq.'''
    is_gap = seq[0] in gap_chars
    last = 0
    size = len(seq)
    for x,c in enumerate(seq):
        if c in gap_chars:
            if not is_gap:
                last = x
                is_gap = True
        else:
            if is_gap:
                yield( last, x )
                last = x
                is_gap = False
    if is_gap:
        yield last, size

def segmentGaps( infile, gap_char ):

    iterator = FastaIterator.FastaIterator( infile )

    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        if cur_record is None: break
        contig = re.sub("\s.*", "", cur_record.title )

        for start, end in gapped_regions(cur_record.sequence, gap_char):
            yield( contig, start, end, 0)

def segmentUngapped( infile, gap_char, min_gap_size = 0 ):

    iterator = FastaIterator.FastaIterator( infile )

    while 1:
        try:
            cur_record = iterator.next()
        except StopIteration:
            break

        if cur_record is None: break
        contig = re.sub("\s.*", "", cur_record.title )
        size = len(cur_record.sequence)

        last_end = 0
        for start, end in gapped_regions(cur_record.sequence, gap_char):
            if end - start < min_gap_size: continue
            
            if last_end != 0:
                yield( contig,last_end,start, 0)
            last_end = end

        if last_end < size:
            yield( contig,last_end,size, 0)

            
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: fasta2bed.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("--gcprofile-gap-size", dest="gcprofile_gap_size", type="float",
                      help="GCProfile: Choose n as the gap size to be filtered "
                      "(If n > 1, n bp is set as the gap size to be filtered. "
                      "If 0 < n < 1, for example, n = 0.01, means gaps less than 1% "
                      "of the input sequence length will be filtered "
                      "[default=%default]."  )

    parser.add_option("--gcprofile-min-length", dest="gcprofile_min_length", type="int",
                      help="GCProfile: minimum length [default=%default]" )

    parser.add_option( "--gcprofile-halting-parameter", dest="gcprofile_halting_parameter", type="int",
                       help="GCProfile: halting parameter [default=%default]" )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ("GCProfile", 
                                 "fixed-width-windows-gc",
                                 "cpg",
                                 "gaps",
                                 "ungapped",),
                      help="Method to use for segmentation [default=%default]" )

    parser.add_option( "-w", "--window-size=", dest="window_size", type="int",
                       help="window size for fixed-width windows [default=%default]." )
    
    parser.add_option("--min-length", dest="min_length", type="int",
                      help="minimum length for ungapped regions [default=%default]" )

    parser.set_defaults(
        gcprofile_halting_parameter = 1000,
        gcprofile_gap_size = 0.01,
        gcprofile_min_length = 3000,
        window_size = 10000,
        method = "GCProfile",
        gap_char = "NnXx",
        min_length = 0,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth

    if options.method == "GCProfile":
        segments = segmentWithGCProfile( options.stdin, options )
    elif options.method == "cpg":
        segments = segmentWithCpG( options.stdin, options )
    elif options.method == "Isoplotter":
        segments = segmentWithIsoplotter( options.stdin, options )
    elif options.method == "fixed-width-windows-gc":
        segments = segmentFixedWidthWindows( options.stdin, options.window_size )
    elif options.method == "gaps":
        segments = segmentGaps( options.stdin, options.gap_char )
    elif options.method == "ungapped":
        segments = segmentUngapped( options.stdin, options.gap_char, options.min_length )
    else:
        raise ValueError("unknown method %s" % (method))
    x = 0
    for contig, start, end, gc in segments:
        x += 1
        options.stdout.write( "%s\n" % "\t".join( \
            (contig, str(start), str(end), str(x), "%6.4f" % (100.0 * gc) ) ) )
        
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
