####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: psl2map.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####

'''
psl2map.py - build a mappping from blat alignments
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts reads :term:`psl` formatted alignments and builds a map of queries 
to targets. The mapping can be restricted by different measures of uniqueness.

If polyA processing is turned on, the non-overlapping terminus of each read will 
be checked if they are mostly A. If they are, the query coverage will be adjusted
appropriately and the read flagged in the section polyA.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

'''

import sys
import re
import string
import optparse
import time
import os
import tempfile
import shutil

import CGAT.Intervals as Intervals
import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import CGAT.Blat as Blat
import CGAT.IndexedFasta as IndexedFasta

def printHistogram( values, section, options, min_value = 0, increment = 1.0 ):

    if len(values) == 0:
        if options.loglevel >= 1:
            options.stdlog.write("# no histogram data for section %s\n" % (section))
        return

    outfile = open(options.output_filename_pattern % section, "w" )
    h = Histogram.Calculate( values, no_empty_bins = True, min_value = 0, increment = 1.0 )

    outfile.write("bin\t%s\n" % section)
    for bin, val in h:
        outfile.write( "%5.2f\t%i\n" % (bin, val) )
    outfile.close()

def printMatched( query_ids, section, options ):

    outfile = open(options.output_filename_pattern % section, "w" )

    for query_id in query_ids:
        outfile.write( "%s\n" % (query_id) )
    outfile.close()


def detectPolyA( query_id, matches, options, queries_fasta = None ):
    """detect PolyA tails and adjust coverage.

    polyA detection. If there are ambiguous matches, the location
    of the polyA tail is not straight-forward as suboptimal matches
    might be declared to be a tail.
    
    1. collect all matches with polyA tail - the remainder is unchanged
    2. declare polyA tail by the shortest! unaligned segment and compute
       coverage for each match appropriately. 
    
    The method checks whether the tails are consistent (always at the same end).
    If not, an AssertionError is thrown
    """
    max_total = matches[0].mQueryLength
    best_start, best_end, best_pA, best_pT, best_tail = 0, 0, 0, 0, ""

    tail_matches = []
    new_matches = []
    for match in matches:
        missing_start = match.mQueryFrom
        missing_end = match.mQueryLength - match.mQueryTo
        if missing_start < missing_end:
            smaller = missing_start
            larger = missing_end
            start,end = match.mQueryTo,match.mQueryLength
        else:
            smaller = missing_end
            larger = missing_start
            start,end = 0, match.mQueryFrom

        # check if tail is at least polyA_min_aligned and at most polyA_max_unaligned
        # are missing from the other end.
        if not(smaller < options.polyA_max_unaligned and larger > options.polyA_min_unaligned):
            new_matches.append( match )
            continue

        tail = queries_fasta.getSequence( query_id )[start:end]

        counts = {"A": 0, "T": 0, "N": 0}
        for c in tail.upper(): counts[c] = counts.get(c, 0) + 1
        total = end-start
        pA = 100.0 * (counts["A"] + counts["N"]) / total
        pT = 100.0 * (counts["T"] + counts["N"]) / total

        if options.loglevel >= 5:
            options.stdlog.write( "# polyA detection: %s:%i-%i pA=%5.2f pT=%5.2f tail=%s\n" % (query_id, start, end, pA, pT, tail))

        if max(pA,pT) < options.polyA_min_percent:
            new_matches.append( match )
            continue

        if total < max_total:
            max_total = total
            best_start, best_end, best_pA, best_pT, best_tail = start, end, pA, pT, tail

        if not( best_start == start or best_end == end ):
            if options.loglevel >= 1:
                options.stdlog.write("# inconsistent polyA tails for %s: %i-%i, %i-%i - analysis skipped\n" %\
                                         (query_id, best_start, best_end, start,end) )
            return matches

        tail_matches.append(match)

    if tail_matches:
        for match in tail_matches:
            match.mQueryCoverage += 100.0 * float(len(best_tail)) / match.mQueryLength
            assert match.mQueryCoverage <= 100.0, "%s: coverage %f > 100.0: incr=%f\n%s" % (query_id, match.mQueryCoverage, float(len(best_tail)) / match.mQueryLength, str(match))
            new_matches.append( match )

        options.outfile_polyA.write( "%s\t%i\t%i\t%i\t%5.2f\t%5.2f\t%s\n" % \
                                         (query_id, 
                                          len(tail_matches),
                                          best_start, best_end,
                                          best_pA, best_pT, best_tail ) )

    assert len(new_matches) == len(matches)

    return new_matches

def selectMatches( query_id, matches, options, queries_fasta = None ):
    """find the best match."""

    if options.loglevel >= 2:
        options.stdlog.write("# attempting to select best match for %s\n" % query_id )
        
        if options.loglevel >= 3:
            for match in matches:
                options.stdlog.write("# match=%s\n" % str(match) )

    new_matches = []

    if options.polyA:
        matches = detectPolyA( query_id, matches, options, queries_fasta )

    if options.matching_mode == "all":
        return matches, None

    elif options.matching_mode in ("best-coverage", "best-query-coverage", "best-sbjct-coverage",
                                   "best-pid", 
                                   "best-covpid", "best-query-covpid", "best-sbjct-covpid",
                                   "best-min-covpid", "best-query-min-covpid", "best-sbjct-min-covpid",):
        if options.matching_mode == "best-coverage":
            f = lambda match: min(match.mQueryCoverage, match.mSbjctCoverage)
        elif options.matching_mode == "best-query-coverage":
            f = lambda match: match.mQueryCoverage
        elif options.matching_mode == "best-sbjct-coverage":
            f = lambda match: match.mSbjctCoverage
        elif options.matching_mode == "best-pid":
            f = lambda match: match.mPid
        elif options.matching_mode == "best-covpid":
            f = lambda match: min(match.mQueryCoverage, match.mSbjctCoverage) * match.mPid
        elif options.matching_mode == "best-query-covpid":
            f = lambda match: match.mQueryCoverage * match.mPid
        elif options.matching_mode == "best-sbjct-covpid":
            f = lambda match: match.mSbjctCoverage * match.mPid
        elif options.matching_mode == "best-min-covpid":
            f = lambda match: min( (match.mQueryCoverage, match.mSbjctCoverage, match.mPid) )
        elif options.matching_mode == "best-query-min-covpid":
            f = lambda match: min(match.mQueryCoverage, match.mPid)
        elif options.matching_mode == "best-sbjct-min-covpid":
            f = lambda match: min(match.mSbjctCoverage, match.mPid)

        for match in matches:
            match.mMatchScore = f(match)

        # collect "significant" matches
        # this filter removes matches out of contention, i.e., matches
        # with very low score are not considered when assessing the uniqueness
        # of the highest scoring match
        matches.sort( lambda x,y: cmp(x.mMatchScore, y.mMatchScore) )
        matches.reverse()
        best_score = min( matches[0].mMatchScore * options.collection_threshold, 
                          matches[0].mMatchScore - options.collection_distance ) 

        for match in matches:
            ## stop when matchscore drops below best score
            if match.mMatchScore < best_score: break
            new_matches.append( match )

        if not options.keep_all_best:

            if len( new_matches ) > 1:
            
                if len(new_matches) == 2:
                    # accept matches against chrX and chrX_random (or vice versa)
                    if new_matches[0].mSbjctId == "%s_random" % new_matches[1].mSbjctId:
                        return new_matches[1:], None
                    elif new_matches[1].mSbjctId == "%s_random" % new_matches[0].mSbjctId:
                        return new_matches[:1], None
                    # or against chrUn or chrU.
                    else:
                        new_matches = [ x for x in new_matches if not( x.mSbjctId.endswith( "Un" ) or x.mSbjctId.endswith("chrU" )) ]
                        if len(new_matches) == 1:
                            return new_matches, None

                if options.ignore_all_random:
                    new_matches = [ x for x in new_matches if not( x.mSbjctId.endswith( "_random" ) or x.mSbjctId.endswith( "Un" ) or x.mSbjctId.endswith("chrU" ) )]
                    if len(new_matches) == 1:
                        return new_matches, None

                return [], "not unique: %s" % (" ".join( map(lambda x: str(x.mMatchScore), matches) ) )
        
    elif options.matching_mode == "unique":
        ## only return matches if they are "unique", i.e. no other match
        if len(matches) == 1:
            new_matches.append( matches[0] )
        else:
            return [], "not unique: %s" % (" ".join( map(lambda x: str(x.mMatchScore), matches) ) )

    matches = new_matches

    if options.best_per_sbjct:
        new_matches = []
        sbjcts = set()
        for match in matches:
            if match.mSbjctId in sbjcts: continue
            new_matches.append( match )
            sbjcts.add(match.mSbjctId)

        matches = new_matches

    return matches, None

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: psl2map.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("--input-filename-queries", dest="input_filename_queries", type="string",
                      help="fasta filename with queries - required for polyA analysis [%default]."  )

    parser.add_option( "--polyA", dest="polyA", action="store_true",
                       help="detect polyA tails [%default].")

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename with histogram information on aggregate coverages [%default].")

    parser.add_option( "--output-filename-empty", dest="output_filename_empty", type="string" ,
                       help="OUTPUT filename with queries for which all matches have been discarded [%default].")

    parser.add_option( "-o", "--output-format", dest="output_format", type="choice" ,
                       choices=("map", "psl" ),
                       help="output format to choose [%default].")

    parser.add_option( "-z", "--from-zipped", dest="from_zipped", action="store_true",
                       help="input is zipped.")

    parser.add_option( "--threshold-min-pid", dest="threshold_min_pid", type="float",
                       help="minimum thresholds for pid [%default]." )

    parser.add_option( "--threshold-min-matches", dest="threshold_min_matches", type="int",
                       help="minimum threshold for number of matching residues [%default]." )

    parser.add_option( "--threshold-max-error-rate", dest="threshold_max_error_rate", type="float",
                       help="maximum threshold for error of aligned part [%default]." )

    parser.add_option( "--threshold-good-query-coverage", dest="threshold_good_query_coverage", type="float",
                       help="minimum query coverage for segments to be counted as good [%default]." )

    parser.add_option( "--threshold-min-query-coverage", dest="threshold_min_query_coverage", type="float",
                       help="minimum query coverage for segments to be accepted [%default]." )

    parser.add_option( "--threshold-max-query-gapchars", dest="threshold_max_query_gapchars", type="int",
                       help="maximum number of gap characters  in query[%default]." )

    parser.add_option( "--threshold-max-query-gaps", dest="threshold_max_query_gaps", type="int",
                       help="maximum number of gaps  in query[%default]." )

    parser.add_option( "--threshold-max-sbjct-gapchars", dest="threshold_max_sbjct_gapchars", type="int",
                       help="maximum number of gap characters  in sbjct[%default]." )

    parser.add_option( "--keep-unique-matches", dest="keep_unique_matches", action="store_true",
                       help="ignore filters for unique matches [%default]." )

    parser.add_option( "--keep-all-best", dest="keep_all_best", action="store_true",
                       help="when sorting matches, keep all matches within the collection threshold [%default]." )

    parser.add_option( "--best-per-sbjct", dest="best_per_sbjct", action="store_true",
                       help="keep only the best entry per sbjct (for transcript mapping) [%default]." )

    parser.add_option( "--threshold-max-sbjct-gaps", dest="threshold_max_sbjct_gaps", type="int",
                       help="maximum number of gaps  in sbjct[%default]." )

    parser.add_option( "--test", dest="test", type="int",
                       help="test - stop after # rows of parsing[%default]." )

    parser.add_option( "-m", "--matching-mode", dest="matching_mode", type="choice",
                       choices=("best-coverage", "best-query-coverage", "best-sbjct-coverage", 
                                "best-pid", "best-covpid", "best-query-covpid", "best-sbjct-covpid", 
                                "best-min-covpid", "best-query-min-covpid", "best-sbjct-min-covpid", 
                                "unique", "all" ),
                       help="determines how to selecte the best match [%default]." )

    parser.add_option( "--filename-filter-sbjct", dest="filename_filter_sbjct", type="string",
                       help="gff file for filtering sbjct matches. Matches overlapping these regions are discarded, but see --keep-forbidden [%default]." )

    parser.add_option( "--keep-forbidden", dest="keep_forbidden", action="store_true",
                       help="if set, keep only matches that overlap the regions supplied with --filename-filter-sbjct [%default]." )

    parser.add_option( "--query-forward-coordinates", dest="query_forward_coordinates", action="store_true",
                       help="use forward coordinates for query, strand will refer to sbjct [%default]." )

    parser.add_option( "--ignore-all-random", dest="ignore_all_random", action="store_true",
                       help="if there are multiple best matches, ignore all those to chrUn and _random [%default]." )

    parser.add_option( "--collection-threshold", dest="collection_threshold", type="float",
                       help="threshold for collecting matches, percent of best score [%default]." )

    parser.add_option( "--collection-distance", dest="collection_distance", type="float",
                       help="threshold for collecting matches, difference to best score [%default]." )

    parser.set_defaults( input_filename_domains = None,
                         input_filename_queries = None,
                         threshold_good_query_coverage = 90.0,
                         threshold_min_pid = 30.0,
                         threshold_min_matches = 0,
                         threshold_max_error_rate = None,
                         output_filename_pattern = "%s",
                         keep_unique_matches = False,
                         output_format = "map",
                         print_matched = ["full", "partial", "good" ],
                         from_zipped = False,
                         combine_overlaps = True,
                         min_length_domain = 30,
                         threshold_min_query_coverage = 50,
                         min_length_singletons = 30,
                         new_family_id = 10000000,
                         add_singletons = False,
                         matching_mode = "best-coverage",
                         best_per_sbjct = False,
                         threshold_max_query_gapchars = None,
                         threshold_max_query_gaps = None,
                         threshold_max_sbjct_gapchars = None,
                         threshold_max_sbjct_gaps = None,
                         filename_filter_sbjct = None,
                         keep_forbidden = False,
                         keep_all_best = False,
                         test = None,
                         query_forward_coordinates=False,
                         output_filename_empty = None,
                         collection_threshold = 1.0, 
                         collection_distance = 0,
                         polyA = False,
                         # max residues missing from non polyA end
                         polyA_max_unaligned = 3,
                         # min residues in tail
                         polyA_min_unaligned = 10,
                         # min percent residues that are A/T in tail
                         polyA_min_percent = 70.0,
                         ## ignore duplicate matches if they are on Un or _random
                         ignore_all_random = False,
                         )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) == 1:
        if options.from_zipped or args[0][-3:] == ".gz":
            import gzip
            infile = gzip.open( args[0], "r" )
        else:
            infile = open(args[0], "r" )
    else:
        infile = sys.stdin

    if options.input_filename_queries:
        queries_fasta =  IndexedFasta.IndexedFasta( options.input_filename_queries )
    else:
        queries_fasta = None
            
    if options.filename_filter_sbjct:

        try:
            import bx.intervals.io
            import bx.intervals.intersection
        except ImportError:
            raise "filtering for intervals requires the bx tools."

        intervals = GTF.readGFFFromFileAsIntervals( open( options.filename_filter_sbjct, "r" ) )

        intersectors = {}

        for contig, values in intervals.items():
            intersector = bx.intervals.intersection.Intersecter()
            for start, end in values:
                intersector.add_interval( bx.intervals.Interval(start,end) )
            intersectors[contig] = intersector

        if options.loglevel >= 1:
            options.stdlog.write("# read %i intervals for %i contigs.\n" %\
                                 (sum( [ len(x) for x in intervals.values() ] ),
                                  len( intersectors ) ))
    else:
        intersectors = None

    ################################################
    ################################################
    ################################################
    ## processing of a chunk (matches of same query)
    ################################################        
    ninput, noutput, nskipped = 0, 0, 0

    ## number of sequences with full/partial/good matches
    nfull_matches, npartial_matches, ngood_matches = 0, 0, 0
    ## number of sequences which are fully/good/partially matched
    ## i.e., after combining all aligned regions
    nfully_matched, npartially_matched, nwell_matched = 0, 0, 0
        
    nremoved_pid, nremoved_query_coverage, nempty = 0, 0, 0
    nremoved_gaps, nremoved_nmatches = 0, 0
    nremoved_regions = 0
    nqueries_removed_region = 0
    
    aggregate_coverages = []
    mapped_coverages = []
    fully_matched = []
    well_matched = []
    partially_matched = []
    new_family_id = options.new_family_id

    if options.output_filename_empty:
        outfile_empty = open( options.output_filename_empty, "w" )
        outfile_empty.write( "read_id\tcomment\n" )
    else:
        outfile_empty = None

    if options.polyA:
        options.outfile_polyA = open( options.output_filename_pattern % "polyA", "w")
        options.outfile_polyA.write( "query_id\tstart\tend\tpA+N\tpT+N\ttail\n" )

    def processChunk( query_id, matches ):
        """process a set of matches from query_id"""

        global ninput, noutput, nskipped
        global nfull_matches, npartial_matches, ngood_matches
        global nremoved_pid, nremoved_query_coverage, nempty, nremoved_gaps, nremoved_nmatches
        global nremoved_regions, nqueries_removed_region
        global outfile_empty
        ninput += 1

        full_matches = []
        good_matches = []
        partial_matches = []

        x_nremoved_pid, x_nquery_coverage, x_nremoved_gaps, x_nremoved_nmatches = 0, 0, 0, 0
        nmatches = len(matches)

        new_matches = []

        # absolute filters applicable to non-fragmentory matches

        for match in matches:

            if match.mPid < options.threshold_min_pid:
                nremoved_pid += 1
                continue
                
            if match.mNMatches < options.threshold_min_matches:
                nremoved_nmatches += 1
                continue

            if options.threshold_max_error_rate:
                r = 100.0 * math.power( options.threshold_max_error_rate, match.mNMatches + match.mNMismatches)
                if match.mPid < r:
                    nremoved_pid += 1
                    x_nremoved_pid += 1
                    continue
            
            new_matches.append(match)

        matches = new_matches

        # filter matches        
        if len(matches) == 0:
            if outfile_empty:
                outfile_empty.write( "%s\tall matches removed after applying thresholds: before=%i, npid=%i, nqcoverage=%i, ngaps=%i, nmatches=%i\n" %\
                                     (query_id, nmatches, x_nremoved_pid, x_nquery_coverage, x_nremoved_gaps, x_nremoved_nmatches ) )
            nskipped += 1
            return
        
        if options.keep_unique_matches and len(matches) == 1:
            pass
        else:
            new_matches = []

            for match in matches:

                if match.mQueryCoverage < options.threshold_min_query_coverage:
                    nremoved_query_coverage += 1
                    x_nquery_coverage += 1
                    continue

                if options.threshold_max_query_gaps and options.threshold_max_query_gaps > match.mQueryNGapsCounts:
                    nremoved_gaps += 1
                    x_nremoved_gaps += 1
                    continue

                if options.threshold_max_query_gapchars and options.threshold_max_query_gapchars > match.mQueryNGapsBases:
                    nremoved_gaps += 1
                    x_nremoved_gaps += 1
                    continue

                if options.threshold_max_sbjct_gaps and options.threshold_max_sbjct_gaps > match.mSbjctNGapsCounts:
                    nremoved_gaps += 1
                    x_nremoved_gaps += 1
                    continue

                if options.threshold_max_sbjct_gapchars and options.threshold_max_sbjct_gapchars > match.mSbjctNGapsBases:
                    nremoved_gaps += 1
                    x_nremoved_gaps += 1
                    continue
                
                new_matches.append( match )
            matches = new_matches

        if len(matches) == 0:
            if outfile_empty:
                outfile_empty.write( "%s\tall matches removed after applying thresholds: before=%i, npid=%i, nqcoverage=%i, ngaps=%i, nmatches=%i\n" %\
                                     (query_id, nmatches, x_nremoved_pid, x_nquery_coverage, x_nremoved_gaps, x_nremoved_nmatches ) )
            nskipped += 1
            return

        ## Remove queries matching to a forbidden region. This section
        ## will remove the full query if any of its matches matches in a
        ## forbidden region.
        keep = True
        for match in matches:
            if intersectors and match.mSbjctId in intersectors:
                found = intersectors[match.mSbjctId].find( match.mSbjctFrom, match.mSbjctTo )
                if found and not options.keep_forbidden or (found and not options.keep_forbidden):
                    nremoved_regions += 1
                    keep = False
                    continue

        if not keep:
            nqueries_removed_region += 1
            if outfile_empty:
                outfile_empty.write( "%s\toverlap with forbidden region\n" % query_id )
            return 

        ## check for full length matches
        for match in matches:
            if match.mQueryCoverage >= 99.9:
                full_matches.append(match)
            if match.mQueryCoverage > options.threshold_good_query_coverage:
                good_matches.append(match)
            else:
                partial_matches.append(match)
            
        if full_matches:
            nfull_matches += 1
        elif good_matches:
            ngood_matches += 1
        elif partial_matches:
            npartial_matches += 1

        ## compute coverage of sequence with matches
        intervals = []
        for match in full_matches + good_matches + partial_matches:
            intervals.append( (match.mQueryFrom, match.mQueryTo) )
        
        rest = Intervals.complement( intervals, 0, match.mQueryLength )
        
        query_coverage = 100.0 * (match.mQueryLength - sum( map( lambda x: x[1] - x[0], rest) ) ) / match.mQueryLength

        if query_coverage >= 99.9:
            fully_matched.append( query_id )
        elif  query_coverage > options.threshold_good_query_coverage:
            well_matched.append( query_id )
        else:
            partially_matched.append( query_id )

        aggregate_coverages.append( query_coverage )

        ## select matches to output
        matches, msg = selectMatches( query_id, matches, options, queries_fasta )

        if len(matches) > 0:
            for match in matches:
                if options.query_forward_coordinates:
                    match.convertCoordinates()

                if options.output_format == "map":
                    options.stdout.write( "%s\n" %\
                                              "\t".join( map(str, (
                                match.mQueryId, match.mSbjctId, 
                                match.strand,
                                "%5.2f" % match.mQueryCoverage,
                                "%5.2f" % match.mSbjctCoverage,
                                "%5.2f" % match.mPid,
                                match.mQueryLength,
                                match.mSbjctLength,
                                match.mQueryFrom, match.mQueryTo,
                                match.mSbjctFrom, match.mSbjctTo,
                                ",".join( map(str,match.mBlockSizes) ),
                                ",".join( map(str,match.mQueryBlockStarts)),
                                ",".join( map(str,match.mSbjctBlockStarts)), 
                                ))))
                elif options.output_format == "psl":
                    options.stdout.write( str(match) + "\n" )

            noutput += 1
        else:
            if outfile_empty:
                outfile_empty.write( "%s\tno matches selected: %s\n" % (query_id, msg) )
            nempty += 1
            
    if options.output_format == "map":
        options.stdout.write( "\t".join( ("query_id", "sbjct_id", "sstrand", "qcoverage","scoverage", "pid", "qlen", "slen", "qfrom", "qto", "sfrom", "sto", "blocks", "qstarts", "sstarts" ) ) + "\n" )
    elif options.output_format == "psl":
        options.stdout.write( Blat.Match().getHeader() + "\n" )

    ################################################
    ################################################
    ################################################
    ## main loop
    ################################################        
    nfully_covered = None
    matches = []
    last_query_id = None
    is_complete = True
    ninput_lines = 0

    skip = 0

    iterator = Blat.BlatIterator( infile )

    while 1:

        try:
            match = iterator.next()
        except Blat.ParsingError:
            iterator = Blat.BlatIterator( infile )
            continue
        
        if match == None: break
        
        ninput_lines += 1

        if options.test and ninput_lines > options.test:
            break
        
        if match.mQueryId != last_query_id:
            if last_query_id:
                processChunk( last_query_id, matches )
            matches = []
            last_query_id = match.mQueryId

        matches.append(match)

    processChunk( last_query_id, matches )

    printHistogram( aggregate_coverages, "aggregate", options )
    
    printHistogram( mapped_coverages, "mapped", options )

    if "full" in options.print_matched:
        printMatched( fully_matched, "full", options )

    if "good" in options.print_matched:
        printMatched( well_matched, "good", options )

    if "partial" in options.print_matched:
        printMatched( partially_matched, "partial", options )

    if options.loglevel >= 1:
        options.stdlog.write("# alignments: ninput=%i, is_complete=%s\n" % (ninput_lines, str(is_complete)))
        options.stdlog.write("# queries: ninput=%i, noutput=%i\n" % (ninput, noutput ))
        options.stdlog.write("# individual coverage: full=%i, good=%i, partial=%i\n" % (nfull_matches, ngood_matches, npartial_matches ) )
        options.stdlog.write("# aggregate  coverage: full=%i, good=%i, partial=%i\n" % (len(fully_matched), len(well_matched), len(partially_matched) ) )
        options.stdlog.write("# omitted queries: total=%i, thresholds=%i, regions=%i, selection=%i\n" %\
                             (nskipped+nqueries_removed_region+nempty,
                              nskipped, nqueries_removed_region, nempty))
        options.stdlog.write("# omitted matches: pid=%i, query_coverage=%i, gaps=%i, regions=%i, nmatches=%i\n" % (nremoved_pid, nremoved_query_coverage, nremoved_gaps, nremoved_regions, nremoved_nmatches ))

    E.Stop()
