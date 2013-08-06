####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@helsinki.fi>
##
## $Id: psl2psl.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####
'''
psl2psl.py - manipulate psl files
===================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads alignments in :term:`psl` format,
applies one ore more operations on them and then
outputs them again in :term:`psl`.

The following operators are available:

merge
   merges psl alignments between overlapping entries of query/target
        in a file. Merging is stopped on query/target and/or strand change.
        This script will check for overlaps and will take the longest path
        through a set of overlapping paths.

map 
   map alignments. Similar to mapPsl, but will apply a map to both query
   and target. This method requires the two options filter-query and filter-target.
   Intervals are given in gff or gtf format. In the latter case, the alignment is 
   threaded through the exons. 

filter-keep
   filter alignments. Only parts of the alignments are kept that are part of the intervals.

filter-remove
   filter alignments. Only parts of the alignments are kept that are not part of the intervals.

complement
   complement a psl sequence. For example this will convert exons to introns.
   If an entry contains only a single exon, it is omitted.

select-longest
   select only the longest transcript (nmatches)

rename-query
   rename matches (query) consecutively.

test
   go through file and remove all non-parseable lines

sanitize
   change query/target names according to fasta sequence database

remove-overlapping-query
   remove all alignments that are not unique with respect to query segment
   (requires alignments to be sorted by query)

remove-overlapping-target
   remove all alignments that are not unique with respect to a target segment
   (requires alignments to be sorted by target)

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
import collections
import warnings

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.Genomics as Genomics
import CGAT.Blat as Blat
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome

## ignore warnings from networkx/matplotlib that a display
## can not be found
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import networkx

import alignlib
import bx.intervals.io
import bx.intervals.intersection

def readIntervals( infile, options ):

    ninput = 0
    t = time.time()

    if options.format == "gtf":

        index = IndexedGenome.IndexedGenome()

        for gffs in GTF.transcript_iterator( GTF.iterator(infile) ):
        
            ali = alignlib.makeAlignmentBlocks()
            for gff in gffs:
                if gff.feature != "exon": continue
                ali.addDiagonal( gff.start, gff.end, 0 )

            index.add( min( [x.start for x in gffs]), 
                       max( [x.end for x in gffs] ),
                       ali )
            ninput += 1
        
            if ninput % options.report_step == 0:
                E.info("reading intervals - progress: ninput=%i, time=%i, avg=%f" % (ninput, time.time() -t, float(time.time() -t) / ninput))
                
    elif options.format == "gff":

        index = IndexedGenome.Simple()

        for g in GTF.iterator(infile):
        
            index.add( g.contig, g.start, g.end )
            ninput += 1

            if ninput % options.report_step == 0:
                E.info("reading intervals - progress: ninput=%i, time=%i, avg=%f" % (ninput, time.time() -t, float(time.time() -t) / ninput))

    E.info("read intervals: %i contigs, %i intervals" % (len(index), ninput))
    return index

def iterator_psl_intervals( options ):
    """iterate over psl file yield an entry together with overlapping entries.

    returns tuples of (match, list(query_intervals), list(target_intervals))
    """

    if options.filename_filter_query:
        intervals_query = readIntervals( open(options.filename_filter_query, "r" ), options )
    else:
        intervals_query = None

    if options.filename_filter_target:
        intervals_target = readIntervals( open(options.filename_filter_target, "r" ), options )
    else:
        intervals_target = None
        
    iterator = Blat.BlatIterator( options.stdin )

    ninput = 0

    while 1:
        
        match = iterator.next()
        if not match: break

        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if options.loglevel >= 1 and ninput % options.report_step == 0:
            options.stdlog.write( "# progress: ninput=%i\n" % (ninput))
            options.stdlog.flush()

        qx, tx = None, None
        if intervals_query:
            try:
                qx = list(intervals_query.get( match.mQueryId, match.mQueryFrom, match.mQueryTo ))
            except KeyError:
                qx = []

        if intervals_target:
            try:
                tx = list(intervals_target.get( match.mSbjctId, match.mSbjctFrom, match.mSbjctTo ))
            except KeyError:
                tx = []

        if options.loglevel >= 2:
            options.stdlog.write( "###################################################\n" )
            options.stdlog.write( "# testing match %s\n" % (str(match)))
            options.stdlog.write( "###################################################\n" )

        yield match, qx, tx

def pslFilter( options, keep = True ):
    """filter psl entries. 

    Only positions contained in intervals are kept, unless remove
    is set, in which case only positions not in intervals are kept.

    """

    ## read as gff
    options.format = "gff"

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0

    if keep:
        raise "not implemented"
    else:
        for match, qx, tx in iterator_psl_intervals( options ):

            map_query2target = match.getMapQuery2Target()

            if not qx and not tx:
                options.stdout.write( str(match) + "\n" )
                noutput += 1
                continue
            
            if qx:
                for qstart, qend, v in qx:
                    if match.strand == "-":
                        qstart, qend = match.mQueryLength - qend, match.mQueryLength - qstart
                    map_query2target.removeRowRegion( qstart, qend )
        
            if tx:
                for tstart, tend, v in tx:
                    map_query2target.removeColRegion( tstart, tend )
            
            if map_query2target.getNumAligned() > options.min_aligned:
                match.fromMap( map_query2target, use_strand = True )
                options.stdout.write( str(match) + "\n" )
                noutput += 1
            else:
                ndiscarded += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def pslMap( options ):
    """thread psl alignments using intervals.

    """

    if options.format == "gtf":
        use_copy = False
    else:
        use_copy = True

    ninput, noutput, ndiscarded, nskipped, nskipped_small_queries = 0, 0, 0, 0, 0

    min_length = options.min_aligned

    for match, qx, tx in iterator_psl_intervals( options ):

        map_query2target = match.getMapQuery2Target()

        ninput += 1

        ## if no filter on qx or tx, use full segment
        if qx == None:
            qx = [ (match.mQueryFrom,match.mQueryTo,0) ]
        elif tx == None:
            tx = [ (match.mSbjctFrom,match.mSbjctTo,0) ]

        ## if no overlap: return
        if not qx or not tx: 
            nskipped += 1
            continue

        for query in qx:

            qstart, qend, qval = query

            # skip elements that are too small
            if qend - qstart < min_length: 
                E.debug( "query too small - skipped at %s:%i-%i" % (match.mQueryId, qstart, qend) )
                nskipped_small_queries += 1
                continue

            E.debug( "working on query %s:%i-%i" % (match.mQueryId, qstart, qend) )

            mqstart, mqend = ( map_query2target.mapRowToCol(qstart, 
                                                            alignlib.RIGHT), 
                               map_query2target.mapRowToCol(qend, 
                                                            alignlib.LEFT) )
                        
                
            if match.strand == "-":
                qstart, qend = match.mQueryLength - qend, match.mQueryLength - qstart

            for target in tx:

                tstart, tend, tval = target
                if tstart >= mqend or tend <= mqstart: continue
                if tend - tstart < min_length: continue

                new = alignlib.makeAlignmentBlocks()
                    
                if use_copy:
                    # do copy with range filter
                    if options.loglevel >= 3:

                        mtstart, mtend = map_query2target.mapColToRow(tstart), map_query2target.mapColToRow(tend) 
                        E.debug( "query: %i-%i (len=%i)-> %i-%i(len=%i); target: %i-%i (len=%i)-> %i-%i (len=%i)" % \
                                     (qstart, qend,
                                      qend - qstart,
                                      mqstart, mqend,
                                      mqend - mqstart,
                                      tstart, tend,
                                      tend - tstart,
                                      mtstart, mtend,
                                      mtend - mtstart ) )
                                     
                    alignlib.copyAlignment( 
                        new, 
                        map_query2target,
                        qstart, qend,
                        tstart, tend )
                else:
                    # do copy with alignment filter
                    map_query = qval
                    if map_query:
                        tmp = alignlib.makeAlignmentBlocks()                        
                        alignlib.copyAlignment( tmp, map_query2target, map_query, alignlib.RR )
                        if options.loglevel >= 5:
                            options.stdlog.write( "######## mapping query ###########\n" )
                            options.stdlog.write( "# %s\n" % str(alignlib.AlignmentFormatEmissions( map_query2target ) ))
                            options.stdlog.write( "# %s\n" % str(alignlib.AlignmentFormatEmissions( map_query ) ))
                            options.stdlog.write( "# %s\n" % str(alignlib.AlignmentFormatEmissions( tmp ) ))
                    else:
                        tmp = map_query2target
                        
                    map_target = tval
                    if map_target:
                        new = alignlib.makeAlignmentBlocks()
                        alignlib.copyAlignment( new, tmp, map_target, alignlib.CR )                        
                        if options.loglevel >= 5:
                            options.stdlog.write( "######## mapping target ###########\n" )
                            options.stdlog.write( "# before: %s\n" % str(alignlib.AlignmentFormatEmissions( tmp ) ))
                            options.stdlog.write( "# map   : %s\n" % str(alignlib.AlignmentFormatEmissions( map_target ) ))
                            options.stdlog.write( "# after : %s\n" % str(alignlib.AlignmentFormatEmissions( new ) ))
                    else:
                        new = tmp

                if options.loglevel >= 4:
                    E.debug("putative match with intervals: %s and %s: %i-%i" % \
                                (str(query), str(target), qstart, qend ))
                    if options.loglevel >= 5:
                        E.debug( "input : %s" % str(alignlib.AlignmentFormatEmissions( map_query2target ) ))
                        E.debug( "final : %s" % str(alignlib.AlignmentFormatEmissions( new ) ) )

                    if new.getLength() > 0:
                        n = match.copy()
                        n.fromMap( new, use_strand = True )
                        E.info( "match : %s" % (str(n)))

                if new.getNumAligned() > options.min_aligned:
                    n = match.copy()
                    n.fromMap( new, use_strand = True )
                    options.stdout.write( str(n) + "\n" )
                    noutput += 1
                else:
                    ndiscarded += 1

    E.info( "map: ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i, nsmall_queries=%i" % \
                (ninput, noutput, nskipped, ndiscarded, nskipped_small_queries) )

def pslMerge( options ):
    """merge psl alignments.
    """

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0

    last_query = None
    last_target = None
    last_strand = None

    def process( matches ):

        new = matches[0].copy()

        map_query2target = alignlib.makeAlignmentBlocks()

        graph = networkx.DiGraph()
        graph.add_nodes_from(xrange(len(matches)+2))

        matches.sort( key = lambda x: x.mQueryFrom )

        if Genomics.IsPositiveStrand( matches[0].strand ):
            f = lambda x,y: x.mSbjctTo < y.mSbjctFrom
        else:
            f = lambda x,y: x.mSbjctFrom > y.mSbjctTo

        for x in range(0,len(matches)):

            xx = matches[x]
            if options.loglevel >= 6:
                options.stdlog.write( "# graph: %2i %s\n" % (x, str(xx)))

            for y in range(x+1,len(matches)):
                yy = matches[y]
                d = min(xx.mQueryTo,yy.mQueryTo) - max(xx.mQueryFrom,yy.mQueryFrom)
                if d > 0 or not f(xx,yy):
                    continue
                else:
                    graph.add_edge( x, y, { 'weight' : -d } )

        source = len(matches)
        target = len(matches) + 1
        for x in range(len(matches)):
            xx = matches[x]
            graph.add_edge( source, x, { 'weight': xx.mQueryFrom } ) 
            graph.add_edge( x, target, { 'weight' : xx.mQueryLength - xx.mQueryTo } )

        if options.loglevel >= 6:
            networkx.write_edgelist( graph, options.stdlog )

        path = networkx.dijkstra_path( graph, source, target ) 

        if options.loglevel >= 6:
            options.stdlog.write( "# path: %s\n" % (str(path)))

        new_matches = [ matches[x] for x in path[1:-1] ]

        if len(matches) != len(new_matches):
            E.warn( "query=%s, target=%s, strand=%s: removed overlapping/out-of-order segments: before=%i, after=%i" %\
                                          (matches[0].mQueryId,
                                           matches[0].mSbjctId,
                                           matches[0].strand,
                                           len(matches),
                                           len(new_matches) ) )

        matches = new_matches

        for match in matches:
            m = match.getMapQuery2Target()
            alignlib.addAlignment2Alignment( map_query2target, m )

        new.fromMap( map_query2target, use_strand = True )
        
        options.stdout.write( str( new ) + "\n" )
        options.stdout.flush()
        return 1

    while 1:
        
        match = iterator.next()
        if not match: break

        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if options.loglevel >= 10:
            options.stdlog.write( "# input: %s\n" % (str(match)))

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        if match.mQueryId != last_query or match.strand != last_strand or match.mSbjctId != last_target:
            if last_query:
                noutput += process( matches )
            matches = []
            last_query, last_target, last_strand = match.mQueryId, match.mSbjctId, match.strand
        
        matches.append( match )
            
    if last_query:
        noutput += process( matches )

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def pslAddSequence( query_fasta, sbjct_fasta, options ):

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0


    while 1:
        
        match = iterator.next()
        if not match: break

        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        new = Blat.MatchPSLX()
        new.fromPSL( match,
                     query_fasta.getSequence( match.mQueryId, "+", match.mQueryFrom, match.mQueryTo ),
                     sbjct_fasta.getSequence( match.mSbjctId, "+", match.mSbjctFrom, match.mSbjctTo ) )
        
        options.stdout.write( str( new ) + "\n" )
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def pslComplement( query_fasta, target_fasta, options) :
    """complenment psl entries.
    """

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0

    border = options.complement_border
    min_length = options.complement_min_length

    while 1:
        
        match = iterator.next()
        if not match: break

        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        if match.mNBlocks <= 1: 
            nskipped += 1
            continue

        pairs = []
        for qstart, tstart, size in match.getBlocks():

            qend = qstart + size - border
            qstart += border
            
            if qend - qstart < options.complement_min_length:
                continue

            tend = tstart + size - border
            tstart += border

            if tend - tstart < options.complement_min_length:
                continue

            query_sequence = query_fasta.getSequence( match.mQueryId, match.strand, qstart, qend )
            sbjct_sequence = sbjct_fasta.getSequence( match.mSbjctId, "+", tstart, tend )
            

        ndiscarded += 1


        options.stdout.write( str( new ) + "\n" )
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def pslComplementQuery( options) :
    """complement psl entries. 

    Fill the regions from a second psl file. 
    """

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0

    border = options.complement_border
    min_length = options.complement_min_length

    while 1:
        
        match = iterator.next()
        if not match: break

        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        if match.mNBlocks <= 1: 
            nskipped += 1
            continue

        pairs = []
        for qstart, tstart, size in match.getBlocks():

            qend = qstart + size - border
            qstart += border
            
            if qend - qstart < options.complement_min_length:
                continue

            tend = tstart + size - border
            tstart += border

            if tend - tstart < options.complement_min_length:
                continue

        ndiscarded += 1


        options.stdout.write( str( new ) + "\n" )
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def pslSelectQuery( options ):

    ninput, noutput, ndiscarded, nskipped = 0, 0, 0, 0

    value, field = options.select.split( "-" )

    if field == "nmatches":
        f = lambda x: x.mNMatches
    elif field == "nmismatches":
        f = lambda x: x.mNMisMatches

    for data in Blat.iterator_per_query( Blat.iterator( options.stdin ) ):
        
        ninput += 1        
        if options.test and ninput >= options.test:
            break

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        data.sort( key = f )

        if value == "most":
            options.stdout.write( "%s\n" % str(data[-1]) )
        elif value == "least":
            options.stdout.write( "%s\n" % str(data[0]) )
            
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, ndiscarded=%i" % (ninput, noutput, nskipped, ndiscarded) )

def iterator_filter_overlapping_query( psls, options ):
    '''remove alignments that overlap on query.

    If multiple alignments overlap, the one with the highest number
    of matching nucleotides is chosen.
    '''

    # note: only takes the full ranges, but does not check for individual overlap of blocks
    # use connected components and hasAlignmentOverlap
    ninput, noutput, ndiscarded = 0, 0, 0    

    last_contig = None

    for block in Blat.iterator_query_overlap( psls, options.threshold_merge_distance ):

        # commented code is for base-level filtering, which is very slow
        # disabled for now
        # if block[0].mQueryId != last_contig:
        #     last_contig = block[0].mQueryId
        #     E.info( "processing %s" % last_contig )

        l = len(block)
        ninput += l
        if l > 1: 
            ndiscarded += l
            # components = Blat.getComponents( block, by_query = True )
            # for component in components:
            #     m = [ block[x] for x in component ]
            #     m.sort( key = lambda x: -x.mNMatches )
            #     ndiscarded += len(m) - 1
            #     yield m[0]
            #     noutput += 1
        else:
            yield block[0]
            noutput += 1

    E.info( "iterator_filter_overlapping_query: ninput=%i, noutput=%i, ndiscarded=%i" % (ninput, noutput,ndiscarded) )

def iterator_filter_overlapping_target( psls, options ):

    ninput, noutput, ndiscarded = 0, 0, 0    
    for block in Blat.iterator_target_overlap( psls, options.threshold_merge_distance ):
        l = len(block)
        ninput += l
        if l > 1: 
            ndiscarded += l
        else:
            yield block[0]
            noutput += 1

    E.info( "iterator_filter_overlapping_target: ninput=%i, noutput=%i, ndiscarded=%i" % (ninput, noutput,ndiscarded) )

def iterator_rename_query( infile, options ):

    ninput, noutput, nerrors = 0, 0, 0

    map_old2new = {}
    x = 1
    while 1:

        match = infile.next()

        if not match: break
        ninput += 1

        if match.mQueryId not in map_old2new or options.unique:
            new = options.id_format % x
            map_old2new[ match.mQueryId ] = new
            x += 1
        else:
            new = map_old2new[match.mQueryId]
            match.mQueryId = new
        yield match

    if options.output_filename_map:
        outfile = open(options.output_filename_map, "w")
        outfile.write( "%s\t%s\n" % ("old", "new") )
        for old,new in map_old2new.iteritems():
            outfile.write( "%s\t%s\n" % (old, new) )
        outfile.close()

    E.info( "ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput,nerrors) )

def iterator_sanitize( infile, query_fasta, sbjct_fasta, options ):

    fq = query_fasta.getToken
    ft = sbjct_fasta.getToken

    ninput, noutput, nerrors = 0, 0, 0

    while 1:
        try:
            x = infile.next()
        except Blat.ParsingError, msg:
            nerrors += 1
            ninput += 1
            E.warn( str(msg) )
            continue
        except StopIteration:
            break

        if not x: break

        ninput += 1        

        if query_fasta: x.mQueryId = fq( x.mQueryId )
        if sbjct_fasta: x.mSbjctId = ft( x.mSbjctId )

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        yield x

    E.info( "ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput,nerrors) )

def iterator_filter_fasta( infile, query_fasta, sbjct_fasta, options ):
    
    ninput, noutput, nerrors = 0, 0, 0
    
    qmissing, smissing = collections.defaultdict(int), collections.defaultdict(int),
    
    while 1:
        try:
            x = infile.next()
        except Blat.ParsingError, msg:
            nerrors += 1
            ninput += 1
            E.warn( str(msg) )
            continue
        except StopIteration:
                    break

        if not x: break

        ninput += 1        
        skip = False
        if query_fasta: 
            if x.mQueryId not in query_fasta:
                qmissing[x.mQueryId] += 1
                skip = True
        if sbjct_fasta: 
            if x.mSbjctId not in sbjct_fasta:
                smissing[x.mSbjctId] += 1
                skip = True
        if skip:
            nerrors += 1
            continue

        if ninput % options.report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        yield x
        noutput += 1

    if qmissing:
        E.warn( "query: %i ids missing: %s" % (len(qmissing), str(qmissing)))
    if smissing:
        E.warn( "target: %i ids missing: %s" % (len(smissing), str(smissing)))

    E.info( "ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput,nerrors) )

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: psl2psl.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("--filter-query", dest="filename_filter_query", type="string",
                      help="filename with intervals in the query to filter (in gff format) [default=%default]."  )

    parser.add_option("--filter-target", dest="filename_filter_target", type="string",
                      help="filename with intervals in the target to filter (in gff format) [default=%default]."  )


    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("map", "merge", 
                               "add-sequence", "complement", 
                               "select-query", "test", 
                               "filter-keep", "filter-remove", 
                               "rename-query",
                               "sanitize",
                               "filter-fasta",
                               "remove-overlapping-query",
                               "remove-overlapping-target"),
                      help="""action to perform [default=%default].""" )

    parser.add_option("--select", dest="select", type="choice",
                      choices=("most-nmatches", "least-nmatches", "most-nmismatches", "least-nmismatches" ),
                      help="entry to select [default=%default]."  )

    parser.add_option("--header", dest="header", type="choice",
                      choices=( "none", "table", "full" ),
                      help="output psl header [default=%default]."  )

    parser.add_option("--format", dest="format", type="choice",
                      choices=("gff", "gtf" ),
                      help="format of intervals [default=%default]."  )

    parser.add_option("--filename-queries", dest="filename_queries", type="string",
                      help="fasta filename with queries."  )

    parser.add_option("--filename-target", dest="filename_sbjcts", type="string",
                      help="fasta filename with sbjct [default=%default]."  )

    parser.add_option("--id-format", dest="id_format", type="string",
                      help="format of new identifiers for the rename function [default=%default]."  )

    parser.add_option("--unique", dest="unique", action="store_true",
                      help="in the rename function, make each match unique [default=%default]."  )

    parser.add_option("--output-filename-map", dest="output_filename_map", type="string",
                      help="filename with map of old to new labels for rename function [default=%default]."  )

    parser.add_option("--complement-min-length", dest="complement_min_length", type="int",
                      help="minimum length for complemented blocks [default=%default]."  )

    parser.add_option("--complement-border", dest="complement_border", type="int",
                      help="number of residues to exclude before alignment at either end [default=%default]."  )

    parser.add_option("--complement-aligner", dest="complement_aligner", type="choice",
                      choices=("clustal", "dba", "dialign", "dialign-lgs" ),
                      help="aligner for complemented segments [default=%default]."  )

    parser.add_option( "--threshold-merge-distance", dest="threshold_merge_distance", type="int",
                       help="distance in nucleotides at which two adjacent reads shall be merged even if they are not overlapping [%default]." )

    parser.set_defaults( filename_filter_target = None,
                         filename_filter_query = None,
                         filename_queries = None,
                         filename_sbjcts = None,
                         threshold_merge_distance = 0,
                         report_step = 100000,
                         min_aligned = 100,
                         methods = [],
                         format = "gff",
                         select = "most-nmatches",
                         id_format = "%06i",
                         unique = False,
                         output_filename_map = None,
                         header = None,
                         test = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    if options.filename_queries:
        query_fasta =  IndexedFasta.IndexedFasta( options.filename_queries )
    else:
        query_fasta = None

    if options.filename_sbjcts:
        sbjct_fasta =  IndexedFasta.IndexedFasta( options.filename_sbjcts )
    else:
        sbjct_fasta = None

    if "add-sequence" in options.methods and (sbjct_fasta == None or query_fasta == None):
        raise ValueError( "please supply both indexed query and target/genome sequence data." )

    iterator = Blat.iterator( options.stdin )

    if options.header != None or options.header != "none":
        if options.header == "table":
            options.stdout.write( "\t".join( Blat.FIELDS ) + "\n" )
        elif options.header == "full":
            options.stdout.write( Blat.HEADER + "\n" )

    for method in options.methods:
    
        if "map" == method:
            pslMap( options )
            break
        elif "filter-keep" == method:
            pslFilter( options, keep = True )
            break
        elif "filter-remove" == method:
            pslFilter( options, keep = False )
            break
        elif "merge" == method:
            pslMerge( options )
            break
        elif "add-sequence" == method:
            pslAddSequence( query_fasta, sbjct_fasta, options )
            break
        elif "complement" == method:
            pslComplement( query_fasta, sbjct_fasta, options )
            break
        elif "select-query" == method:
            pslSelectQuery( options )
            break
        elif "test" == method:
            iterator = Blat.iterator_test( iterator, options.report_step )
        elif "rename-query" == method:
            iterator = iterator_rename_query( iterator, options )
        elif "sanitize" == method:
            iterator = iterator_sanitize( iterator, query_fasta, sbjct_fasta, options )
        elif "filter-fasta" == method:
            iterator = iterator_filter_fasta( iterator, query_fasta, sbjct_fasta, options )
        elif "remove-overlapping-query" == method:
            iterator = iterator_filter_overlapping_query( iterator, options )
        elif "remove-overlapping-target" == method:
            iterator = iterator_filter_overlapping_target( iterator, options )

    for psl in iterator:
        options.stdout.write( "%s\n" % str( psl ) )

    E.Stop()
