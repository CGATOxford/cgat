'''
gff2chunks.py - split a gff file into chunks
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Genesets GFF Manipulation

Purpose
-------

This scripts splits a gff file into chunks. 
The input file needs to be sorted appropriately.

The gff file is split into chunks ensuring that 
overlapping intervals are part of the same chunk

Usage
-----

Example::

   python gff2chunk.py < in.gff 

Type::

   python gff2chunk.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import os
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Experiment as E

class OutputChunk:
    def __init__(self, options):
        self.nchunk = 0
        self.options = options

    def createOpen( self, mode = "w" , header = None):
        """open file. Check first, if directory exists.
        """

        self.nchunk += 1
        filename = self.options.output_pattern % self.nchunk

        if self.options.dry_run:
            E.info( "opening file %s" % filename )
            return open("/dev/null", mode)

        if mode in ("w", "a"):
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists( dirname ):
                os.makedirs( dirname )

        if os.path.exists( filename ):
            existed = True
        else:
            existed = False
        
        f = IOTools.openFile( filename, mode )

        if header and not existed:
            f.write( header + "\n" )

        return f

    def __call__( self, chunk ):
        """output a chunk into a new file."""
        outfile = self.createOpen()
        for c in chunk: outfile.write( str(c) + "\n" )
        outfile.close()
        return len(chunk)

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2chunks.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option( "-p", "--output-pattern", dest="output_pattern", type="string",
                       help="output pattern for filenames. Should contain a '%i' [default=%default]." )

    parser.add_option( "-i", "--min-chunk-size", dest="min_chunk_size", type="int",
                       help="minimum chunk size [default=%default]." )

    parser.add_option( "-n", "--dry-run", dest="dry_run", action="store_true",
                       help="do not create any files [default=%default]." )

    parser.set_defaults(
        method = "overlap",
        dry_run = False,
        output_pattern = "%06i.chunk",
        min_chunk_size = 1,
        )

    (options, args) = E.Start( parser )

    gffs = GTF.iterator( sys.stdin )
    
    ninput, noutput, nchunks = 0, 0, 0

    outputChunk = OutputChunk( options )

    if options.method == "overlap":
        
        last_contig, last_to = None, None
        chunk = []
        for gff in gffs:
            ninput += 1
            if len(chunk) > options.min_chunk_size and \
                    (gff.contig != last_contig or \
                         gff.start > last_to):
                noutput += outputChunk( chunk )
                nchunks += 1
                chunk = []
                last_contig, last_to = gff.contig, gff.end
            
            chunk.append( gff )
            last_to = max( gff.end, last_to )
            
        noutput += outputChunk( chunk )
        nchunks += 1

    E.info( "ninput=%i, noutput=%i, nchunks=%i" % (ninput, noutput, nchunks ) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

