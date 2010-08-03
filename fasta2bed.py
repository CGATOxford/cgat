################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: fasta2bed.py 2861 2010-02-23 17:36:32Z andreas $
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
"""
fasta2bed.py - compute G+C profile for a DNA sequence
=====================================================

:Author: Andreas Heger
:Release: $Id: fasta2bed.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Uses GCProfile

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse, tempfile, subprocess, glob, shutil

import Experiment as E

EXECUTABLE="GCProfile"

import FastaIterator

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

        statement = "%s %s -t %i -g %f -i %i" % \
                    (EXECUTABLE, filename,
                     options.segmentation_threshold,
                     options.gap_size,
                     options.min_length )
        
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
            raise PipelineError( "Child was terminated by signal %i: \nThe stderr was: \n%s\n%s\n" % (-process.returncode, stderr, statement ))

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
                
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: fasta2bed.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-t", "--segmentation-threshold", dest="segmentation_threshold", type="int",
                      help="Choose n as the threshold for segmentation. [default=%default]."  )

    parser.add_option("-g", "--gap-size", dest="gap_size", type="float",
                      help="Choose n as the gap size to be filtered (If n > 1, n bp is set as the gap size to be filtered. If 0 < n < 1, for example, n = 0.01, means gaps less than 1% of the input sequence length will be filtered [default=%default]."  )

    parser.add_option("-i", "--min-length", dest="min_length", type="int",
                      help="Choose n as minimum length [default=%default]" )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ("GCProfile", "fixed-width-windows" ),
                      help="Method to use for segmentation [default=%default]" )

    parser.add_option( "-w", "--window-size=", dest="window_size", type="int",
                       help="window size for fixed-width windows [default=%default]." )

    parser.set_defaults(
        segmentation_threshold = 1000,
        gap_size = 0.01,
        min_length = 3000,
        window_size = 10000,
        method = "GCProfile",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth

    if options.method == "GCProfile":
        segments = segmentWithGCProfile( options.stdin, options )
    elif options.method == "fixed-width-windows":
        segments = segmentFixedWidthWindows( options.stdin, options.window_size )
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
