"""
beds2beds.py - decompose bed files
==================================

:Author: Andreas Heger
:Release: $Id: diff_bed.py 2866 2010-03-03 10:18:49Z andreas $
:Date: |today|
:Tags: Genomics Intervals BED Manipulation

Purpose
-------

decompose a collection of bed-files into a collection unions/intersections.

merged-combinations
   merge intervals across :term:`bed` files and only report those
   that appear in every file.

unmerged-combinations
   for each :term:`bed` file, report intervals that overlap with intervals
   in every other :term:`bed` file.

If the ``--exclusive`` option is set, report exclusive overlap. Only intervals 
will be reported that overlap in a pairwise comparison but do not overlap with 
intervals in any of the other sets.

This script requires bed files indexed by tabix_.

Usage
-----

For example, you have ChIP-Seq data for PolII and two transcription factors tf1 and tf2.
The following statement will output four :term:`bed` files::

   python beds2beds.py polii.bed.gz tf1.bed.gz tf2.bed.gz

The four files are containing intervals, that

1. have PolII and tf1 present,
2. have PolII and tf2 present,
3. have tf1 and tf2 present, or
4. have PolII and tf1 and tf2 present.

If the --exclusive option is set, three sets will be output with intervals that

1. have PolII and tf1 present but no tf2,
2. have PolII and tf2 present but no tf1,
3. have tf1 and tf2 present bu no PolII.

Type::

   python beds2beds.py --help

for command line help.

Command line options
--------------------

""" 

import os
import sys
import re
import optparse
import itertools
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import numpy
import pysam
import CGAT.Intervals as Intervals

def isContainedInAll( contig, start, end, bedfiles ):

    for bedfile in bedfiles:
        try:
            if len(list(bedfile.fetch( contig, start, end ))) == 0:
                return False
        except ValueError:
            return False

    return True

def isContainedInOne( contig, start, end, bedfiles ):

    for bedfile in bedfiles:
        try:
            if len(list(bedfile.fetch( contig, start, end ))) > 0:
                return True
        except ValueError:
            pass

    return False

def combineMergedIntervals( bedfiles ):
    '''combine intervals in a collection of bed files.
    
    Overlapping intervals between tracks are merged.

    Algorithm:

    1. collect all intervals in all tracks into a single track
    2. merge overlapping intervals 
    3. report all intervals that overlap with an interval in each track.

    '''
    
    # get all intervals
    data_per_contig = collections.defaultdict( list )
    for bedfile in bedfiles:
        for contig in bedfile.contigs:
            i = []
            for bed in bedfile.fetch( contig, parser = pysam.asBed() ):
                i.append( (bed.start, bed.end) )
            data_per_contig[contig].extend( i )
    
    # merge intervals
    for contig in data_per_contig.keys():
        data_per_contig[contig] = Intervals.combine( data_per_contig[contig] )

    # filter intervals - take only those present in all bedfiles
    for contig, data in data_per_contig.iteritems():
        for start, end in data:
            if isContainedInAll( contig, start, end, bedfiles ):
                yield contig, start, end

def combineUnmergedIntervals( foreground, background ):
    '''combine intervals in a collection of bed files.
    
    Only intervals in the first track are reported.

    Algorithm:

    1. report all intervals in the first track that overlap with an interval in every other track.
    '''

    intervals = []
    c = 0
    for bed in foreground.fetch( parser = pysam.asBed() ):
        c += 1
        if isContainedInAll( bed.contig, bed.start, bed.end, background ):
            yield bed

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: diff_bed.py 2866 2010-03-03 10:18:49Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-e", "--exclusive", dest="exclusive", action="store_true",
                      help="Intervals reported will be merged across the positive set"
                           " and do not overlap any interval in any of the other sets"
                           " [default=%default]." )

    parser.add_option("-p", "--pattern-id", dest="pattern_id", type="string",
                      help="pattern to convert a filename to an id [default=%default]." )

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       choices = ("merged-combinations", 
                                  "unmerged-combinations" ),
                       help = "method to perform [default=%default]" )
    
    parser.set_defaults(
        pattern_id = "(.*).bed.gz",
        exclusive = False,
        method = "merged-combinations",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if len(args) < 2:
        raise ValueError( "at least two arguments required" )

    tags, bedfiles = [], []
    for infile in args:
        bedfiles.append( pysam.Tabixfile( infile, "r" ) )
        tags.append( re.search( options.pattern_id, infile ).groups()[0] )

    indices = range(len(bedfiles))
    is_exclusive = options.exclusive

    if options.method == "merged-combinations":

        if is_exclusive: start = 1    
        else: start = 2

        options.stdout.write( "combination\twithout\tcounts\n" )

        for ncombinants in range( start, len(bedfiles) + 1):
            for combination in itertools.combinations( indices, ncombinants ):
                other = [ x for x in indices if x not in combination ]
                tag = ":".join( [ tags[x] for x in combination ] )
                E.debug( "combination %s started" % tag)
                E.debug( "other: %s" % ":".join( [tags[x] for x in other]) )
                other_bed = [ bedfiles[x] for x in other ]
                outf = IOTools.openFile( E.getOutputFile( tag ), "w", create_dir = True)
                c = E.Counter()
                for contig, start, end in combineMergedIntervals( [bedfiles[x] for x in combination ] ):
                    c.found += 1
                    if is_exclusive and isContainedInOne( contig, start, end, other_bed ):
                        c.removed += 1
                        continue
                    c.output += 1
                    outf.write( "%s\t%i\t%i\n" % (contig, start, end ) )

                outf.close()
                E.info( "combination %s finished: %s" % (tag, c))
                
                options.stdout.write( "%s\t%s\t%i\n" % (
                        ":".join( [tags[x] for x in combination] ), 
                        ":".join( [tags[x] for x in other] ), 
                        c.output ) )

    elif options.method == "unmerged-combinations":
        options.stdout.write( "track\tcombination\twithout\tcounts\n" )

        for foreground in indices:

            start = 0

            background = [x for x in indices if x != foreground ]            
            for ncombinants in range( 0, len(background) + 1 ):
                for combination in itertools.combinations( background, ncombinants ):
                    other = [ x for x in background if x not in combination ]
                    combination_bed = [ bedfiles[x] for x in combination ] 
                    other_bed = [ bedfiles[x] for x in other ]
                    tag = ":".join( [tags[foreground]] + [ tags[x] for x in combination ] )

                    E.debug( "fg=%i, combination=%s, other=%s" % (foreground, combination, other) )
                    E.debug( "combination %s started" % tag)
                    E.debug( "other: %s" % ":".join( [tags[x] for x in other]) )

                    outf = IOTools.openFile( E.getOutputFile( tag ), "w", create_dir = True)
                    c = E.Counter()
                    for bed in combineUnmergedIntervals(
                        bedfiles[foreground],
                        combination_bed ):
                        c.found += 1
                        if is_exclusive and isContainedInOne( bed.contig, bed.start, bed.end, other_bed ):
                            c.removed += 1
                            continue
                        c.output += 1
                        outf.write( "%s\n" % str(bed))

                    outf.close()
                    E.info( "combination %s finished: %s" % (tag, c))

                    options.stdout.write( "%s\t%s\t%s\t%i\n" % (
                            tags[foreground],
                            ":".join( [tags[x] for x in combination] ), 
                            ":".join( [tags[x] for x in other] ), 
                            c.output ) )


    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )
