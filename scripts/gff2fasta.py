'''
gff2fasta.py - output sequences from genomic features
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Sequences GFF Fasta Transformation

Purpose
-------

This script outputs the genomic sequences for intervals within
a :term:`gff` formatted file.

The ouput can be optionall masked and filtered.

Usage
-----

Example::

   python gff2fasta.py --genome-file=hg19 < features.gff > features.fasta

Type::

   python gff2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import re
import optparse
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Genomics as Genomics
import CGAT.AGP as AGP
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Intervals as Intervals

import bx.intervals.io
import bx.intervals.intersection

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2fasta.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf instead of gff."  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-m", "--merge", dest="merge", action="store_true",
                      help="merge adjacent intervals with the same attributes. "
                      "[default=%default]")

    parser.add_option("-e", "--feature", dest="feature", type = "string",
                      help="filter by a feature, for example 'exon', 'CDS'. If "
                      "set to the empty string, all entries are output [%default]."  )

    parser.add_option("-f", "--filename-masks", dest="filename_masks", type = "string",
                      metavar = "gff",
                      help="mask sequences with regions given in gff file [%default]."  )

    parser.add_option( "--remove-masked-regions", dest="remove_masked_regions", action="store_true",
                      help="remove regions instead of masking [%default]."  )

    parser.add_option( "--min-length", dest="min_length", type="int",
                       help="set minimum length for sequences output [%default]" )

    parser.add_option( "--max-length", dest="max_length", type="int",
                       help="set maximum length for sequences output [%default]" )

    parser.add_option( "--extend-at", dest="extend_at", type="choice",
                       choices=("none", "3", "5", "both", "3only", "5only" ),
                       help="extend at no end, 3', 5' or both ends. If "
                       "3only or 5only are set, only the added sequence "
                       "is returned [default=%default]" )

    parser.add_option( "--extend-by", dest="extend_by", type="int",
                       help="extend by # bases [default=%default]" )

    parser.set_defaults(
        is_gtf = False,
        genome_file = None,
        merge = False,
        feature = None,
        filename_masks = None,
        remove_masked_regions = False,
        min_length = 0,
        max_length = 0,
        extend_at = None,
        extend_by = 100,
        )

    (options, args) = E.Start( parser )

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = fasta.getContigSizes()

    if options.is_gtf:
        iterator = GTF.transcript_iterator( GTF.iterator( sys.stdin ) )
    else:
        gffs = GTF.iterator( sys.stdin )
        if options.merge:
            iterator = GTF.joined_iterator( gffs )
        else:
            iterator =  GTF.chunk_iterator( gffs )

    masks = None
    if options.filename_masks:
        masks = {}
        with open( options.filename_masks, "r") as infile:
            e = GTF.readAsIntervals( GFF.iterator( infile ) )

        # convert intervals to intersectors
        for contig in e.keys():
            intersector = bx.intervals.intersection.Intersecter()
            for start, end in e[contig]:
                intersector.add_interval( bx.intervals.Interval(start,end) )
            masks[contig] = intersector

    ninput, noutput, nmasked, nskipped_masked = 0, 0, 0, 0
    nskipped_length = 0
    nskipped_noexons = 0

    feature = options.feature

#    for item in iterator:
#	print len(item) # 3, 2
#	for i in item:
#	   print len(i) # 9, 9, 9, 9, 9
#	   print i.contig
#	   print i.strand
#	   print i.transcript_id

    # iterator is a list containing groups (lists) of features.
    # Each group of features have in common the same transcript ID, in case of GTF files.
    for ichunk in iterator:

        ninput += 1

        if feature:
            chunk = filter( lambda x: x.feature == feature, ichunk )
        else:
            chunk = ichunk

        if len(chunk) == 0:
            nskipped_noexons += 1
            E.info("no features in entry from %s:%i..%i - %s" % (ichunk[0].contig,
                                                                 ichunk[0].start,
                                                                 ichunk[0].end,
                                                                 str(ichunk[0])))
            continue

        contig, strand = chunk[0].contig, chunk[0].strand 
        if options.is_gtf:
            name = chunk[0].transcript_id
        else:
            name = str(chunk[0].attributes)

        lcontig = contigs[contig]
        positive = Genomics.IsPositiveStrand( strand )
        intervals = [ (x.start, x.end) for x in chunk ]
        intervals.sort()

        if masks:
            if contig in masks:
                masked_regions = []
                for start, end in intervals:
                    masked_regions += [ (x.start, x.end) for x in masks[contig].find( start, end ) ]
                
                masked_regions = Intervals.combine( masked_regions )
                if len(masked_regions): nmasked += 1

                if options.remove_masked_regions:
                    intervals = Intervals.truncate( intervals, masked_regions )
                else:
                    raise "unimplemented"

                if len(intervals) == 0:
                    nskipped_masked += 1
                    if options.loglevel >= 1:
                        options.stdlog.write( "# skipped because fully masked: %s: regions=%s masks=%s\n" %\
                                                  (name, str([ (x.start, x.end) for x in chunk ]), masked_regions) )
                    continue

        out = intervals

        if options.extend_at:
            if options.extend_at == "5only":
                intervals = [ (max( 0, intervals[0][0] - options.extend_by ), intervals[0][0]) ]            
            elif options.extend_at == "3only":
                intervals = [ (intervals[-1][1], min( lcontig, intervals[-1][1] + options.extend_by)) ]
            else:
                if options.extend_at in ("5", "both"): intervals[0] = (max( 0, intervals[0][0] - options.extend_by ), intervals[0][1])
                if options.extend_at in ("3", "both"): intervals[-1] = (intervals[-1][0], min( lcontig, intervals[-1][1] + options.extend_by))
            
        if not positive:
            intervals = [ (lcontig - x[1], lcontig - x[0]) for x in intervals[::-1] ]
            out.reverse()

        s = [ fasta.getSequence( contig, strand, start, end ) for start,end in intervals ]
        l = sum( [len(x) for x in s ] )
        if l < options.min_length or (options.max_length and l > options.max_length):
            nskipped_length += 1
            if options.loglevel >= 1:
                options.stdlog.write( "# skipped because length out of bounds %s: regions=%s len=%i\n" %\
                                          (name, str(intervals), l) )
            continue

        options.stdout.write( ">%s %s:%s:%s\n%s\n" % (name, 
                                                      contig,
                                                      strand,
                                                      ";".join( ["%i-%i" % x for x in out] ),
                                                      "\n".join( s )) )
        
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nmasked=%i, nskipped_noexons=%i, nskipped_masked=%i, nskipped_length=%i" %\
                (ninput, noutput, nmasked, nskipped_noexons, nskipped_masked, nskipped_length ) )

    E.Stop()
    
if __name__ == "__main__":
    sys.exit( main( sys.argv) )
