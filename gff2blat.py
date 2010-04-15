####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@helsinki.fi>
##
## $Id: gff2blat.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####

USAGE="""python blat2gff.py [OPTIONS] < input > output

convert gff output to blat output.

Todo:
- count gaps
- count mismatches
"""

import sys, re, string, optparse, time, os, glob

import Experiment
import IndexedFasta
import Blat
import Genomics
import GFF, GTF
import alignlib
import Intervals

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gff2blat.py 2781 2009-09-10 11:33:14Z andreas $", usage=USAGE )

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf."  )

    parser.add_option( "--no-header", dest="with_header", action="store_false",
                      help="do not output BLAT header [default=%default]."  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("--input-filename-queries", dest="input_filename_queries", type="string",
                      help="fasta filename with queries [default=%default]."  )

    parser.add_option("--allow-duplicates", dest="allow_duplicates", action="store_true",
                      help="""permit duplicate entries. Adjacent exons of a transcript will still be merged [default=%default]."""  )

    parser.set_defaults( is_gtf = False,
                         genome_file = None,
                         with_header = True,
                         allow_duplicates = False,
                         test = None )
    
    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        genome_fasta = None

    if options.input_filename_queries:
        queries_fasta =  IndexedFasta.IndexedFasta( options.input_filename_queries )
    else:
        queries_fasta = None

    ninput, noutput, nskipped = 0, 0, 0

    if options.is_gtf:
        iterator = GTF.transcript_iterator( GTF.iterator_filtered( GTF.iterator( sys.stdin ), feature="exon" ), 
                                            strict = not options.allow_duplicates )
    else:
        iterator = GFF.joined_iterator( GFF.iterator(sys.stdin) )

    if options.with_header:
        options.stdout.write( Blat.Match().getHeader() + "\n" )

    for gffs in iterator:
        
        if options.test and ninput >= options.test:
            break

        ninput += 1

        result = alignlib.makeAlignmentBlocks()

        xstart = 0

        intervals = Intervals.combine( [ (gff.start, gff.end) for gff in gffs ] )

        for start, end in intervals:
            xend = xstart + end - start
            
            result.addDiagonal( xstart, xend,
                                start - xstart )
            xstart = xend

        entry = Blat.Match()
        entry.mQueryId = gff.transcript_id
        entry.mSbjctId = gff.contig
        entry.strand  = gff.strand

        if genome_fasta:
            if entry.mSbjctId in genome_fasta:
                entry.mSbjctLength = genome_fasta.getLength( entry.mSbjctId )
            else:
                entry.mSbjctLength = result.getColTo()

        if queries_fasta:
            if entry.mQueryId in queries_fasta:
                entry.mQueryLength = queries_fasta.getLength( entry.mQueryId )
        else:
            entry.mQueryLength = result.getRowTo()

        entry.fromMap( result )
        
        options.stdout.write (str(entry) + "\n" )
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    Experiment.Stop()
