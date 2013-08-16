####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $
##
##
####
####
'''
bed2bed.py - manipulate bed files
=================================

Purpose
-------

manipulate bed-formatted files.

The script currently implements the following methods:

1. merge: merge adjacent intervals. The option ``--merge-distance``
   permits the merging of segments that are up to a certain distance apart.

2. block: build blocked bed file (bed12 format) from individual blocks.

3. sanitize-genome: remove all empty intervals and intervals on unknown contigs. 
   Intervals extending beyond a contig a truncated. 
   
4. filter-genome: remove all intervals on unknown contigs or extending beyond 
   contigs.
5. extend: extend the start and end of each interval by a specified distance.
   Note: currently this method does not check that new interval is outside the confines of the contig

Usage
-----

Type::

   python bed2bed.py --help

for command line usage.

Methods
-------

This script provides several methods:

merge

filter-genome

bins
   merge adjacent bins by score. Several merging regimes are possible.
   equal-base, equal-intervals.
   This options requires the fifth field of the bed input file to be
   present.


Code
----
'''

import sys
import re
import string
import optparse
import time
import os
import itertools
import tempfile
import subprocess
import shutil

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import pysam

def merge( iterator, max_distance = 0, by_name = False, min_intervals = 1 ):
    """iterator for merging adjacent bed entries.

    *max_distance* > 0 permits merging of intervals that are
    not directly adjacent.

    If *by_name = True*, only entries with the same name are merged.
    
    The score gives the number of intervals that have been merged.
    """

    def iterate_chunks( iterator ):
        last = iterator.next()
        max_end = last.end
        to_join = [last]

        for bed in iterator:
            d = bed.start - max_end
            if bed.contig == last.contig:
                assert bed.start >= last.start, "input file should be sorted by contig and position: d=%i:\n%s\n%s\n" % (d, last, bed)
            
            if bed.contig != last.contig or \
                    d > max_distance or \
                    (by_name and last.name != bed.name) :
                yield to_join
                to_join = []
                max_end = 0

            last = bed
            max_end = max( last.end, max_end )
            to_join.append( bed )

        if to_join: yield to_join
        raise StopIteration

    ninput, noutput, nskipped_min_intervals = 0, 0, 0

    for to_join in iterate_chunks(iterator):

        ninput += 1
        if len(to_join) < min_intervals: 
            nskipped_min_intervals += 1
            continue

        a = to_join[0]
        a.end = to_join[-1].end
        a.score = len(to_join)
        yield a
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped_min_intervals=%i" % (ninput, noutput, nskipped_min_intervals) )

def filterGenome( iterator, contigs ):
    """remove bed intervals that are outside of contigs.
    
    contigs is a dictionary of contig sizes."""
    
    ninput, noutput = 0, 0
    nskipped_contig, nskipped_range, nskipped_endzero = 0, 0, 0

    for bed in iterator:
        ninput += 1
        if bed.contig not in contigs: 
            nskipped_contig += 1
            continue
        if bed.end >= contigs[bed.contig]:
            nskipped_range += 1
            continue
        if bed.end == 0:
            nskipped_endzero += 1
            continue    
        noutput += 1
        yield bed

    E.info( "ninput=%i, noutput=%i, nskipped_contig=%i, nskipped_range=%i, nskipped_endzero=%i" % \
                (ninput, noutput, nskipped_contig, nskipped_range, nskipped_endzero) )

def sanitizeGenome( iterator, contigs ):
    """truncate bed intervals that extend beyond contigs.

    removes empty intervals (start == end).

    throws an error if start > end.
    """
    
    ninput, noutput = 0, 0
    ntruncated_contig, nskipped_contig, nskipped_empty = 0, 0, 0

    for bed in iterator:
        ninput += 1
        if bed.contig not in contigs: 
            nskipped_contig += 1
            continue
        if bed.end >= contigs[bed.contig]:
            bed.end = contigs[bed.contig]
            ntruncated_contig += 1
        if bed.start < 0:
            bed.start = 0
            ntruncated_contig += 1
        if bed.start == bed.end:
            nskipped_empty += 1
            continue
        elif bed.start > bed.end:
            raise ValueError( "invalid interval: start > end for %s" % str(bed))

        noutput += 1
        yield bed

    E.info( "ninput=%i, noutput=%i, nskipped_contig=%i, ntruncated=%i, nskipped_empty=%i" % \
                (ninput, noutput, nskipped_contig, ntruncated_contig, nskipped_empty) )

def shiftIntervals( iterator, contigs, offset ):
    """shift intervals by a certain offset and ensure size is maintaned even id contig end reached.
    
    contigs is a dictionary of contig sizes."""
    
    ninput, noutput = 0, 0
    nskipped_contig, nskipped_range = 0, 0

    for bed in iterator:
        ninput += 1
        if bed.contig not in contigs: 
            nskipped_contig += 1
            continue
        if bed.end >= contigs[bed.contig]:
            nskipped_range += 1
            continue
        noutput += 1

        # add offset to each start and end, and adjust for contig length
        l = bed.end - bed.start
        newstart = bed.start + offset
        newend = bed.end + offset
        if newstart <0:
            newstart = 0
            newend = l
        if newend > contigs[bed.contig]:
            newstart = contigs[bed.contig] - l
            newend = contigs[bed.contig]

        bed.start = newstart
        bed.end = newend

        yield bed

    E.info( "ninput=%i, noutput=%i, nskipped_contig=%i, nskipped_range=%i" % \
                (ninput, noutput, nskipped_contig, nskipped_range) )

#IMS: new method: extend intervals by set amount
def extendInterval(iterator, distance):

    ninput, noutput, nskipped = 0, 0, 0
    for bed in iterator:
        ninput += 1
        newstart = bed.start - distance
        newend = bed.end + distance

        if newstart < 0:
            nskipped += 1
            continue

        bed.start = newstart
        bed.end = newend

        noutput += 1
        yield bed

    E.info("ninput = %i, noutput=%i, nskipped=%i" % (ninput,noutput,nskipped))


def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    #IMS: new method: extend intervals by set amount
    parser.add_option( "-m", "--method", dest="methods", type="choice", action="append",
                       choices=("merge", "filter-genome", "bins", "block", "sanitize-genome", "shift", "extend" ),
                       help="method to apply [default=%default]"  )

    parser.add_option( "--num-bins", dest="num_bins", type="int",
                      help="number of bins into which to merge (used for method `bins) [default=%default]"  )

    parser.add_option( "--bin-edges", dest="bin_edges", type="string",
                      help="bin_edges for binning method [default=%default]"  )

    parser.add_option( "--binning-method", dest="binning_method", type="choice",
                       choices = ("equal-bases", "equal-intervals", "equal-range" ),
                       help="method used for binning (used for method `bins` if no bin_edges is given) [default=%default]"  )

    parser.add_option( "--merge-distance", dest="merge_distance", type="int",
                      help="distance in bases over which to merge that are not directly adjacent [default=%default]"  )

    parser.add_option( "--merge-min-intervals", dest="merge_min_intervals", type="int",
                      help="only output merged intervals that are build from at least x intervals [default=%default]"  )

    parser.add_option( "--merge-by-name", dest="merge_by_name", action="store_true",
                      help="only merge intervals with the same name [default=%default]"  )

    parser.add_option( "--offset", dest="offset",  type="int",
                      help="offset for shifting intervals [default=%default]"  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-b", "--bam-file", dest="bam_file", type="string",
                      help="bam-formatted filename with genome."  )
    #IMS: new method: extend intervals by set amount
    parser.add_option("--extend-distance",dest="extend_distance", type="int",
                      help = "distance by which to extend intervals on either side for method=extend [default=%default]")

    parser.set_defaults( methods = [],
                         merge_distance = 0,
                         binning_method = "equal-bases",
                         genome_file = None,
                         bam_file = None,
                         num_bins = 5,
                         merge_min_intervals = 1,
                         bin_edges = None,
                         offset = 10000,
                         test = None,
                         extend_distance=1000)
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    contigs = None

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = genome_fasta.getContigSizes()

    if options.bam_file:
        samfile = pysam.Samfile( options.bam_file )
        contigs = dict( zip(samfile.references, samfile.lengths) )

    processor = Bed.iterator( options.stdin )

    for method in options.methods:
        if method ==  "filter-genome":
            if not contigs: raise ValueError("please supply contig sizes" )
            processor = filterGenome( processor, contigs )
        elif method == "sanitize-genome":
            if not contigs: raise ValueError("please supply contig sizes" )
            processor = sanitizeGenome( processor, contigs )
        elif method == "merge":
            processor = merge( processor, options.merge_distance,
                               by_name = options.merge_by_name,
                               min_intervals = options.merge_min_intervals )
        elif method == "bins":
            if options.bin_edges: bin_edges = map(float, options.bin_edges.split(","))
            else: bin_edges = None
            processor, bin_edges = Bed.binIntervals( processor,
                                                     num_bins = options.num_bins,
                                                     method = options.binning_method,
                                                     bin_edges = bin_edges)
            E.info("# split bed: bin_edges=%s" % (str(bin_edges)))
            
        elif method == "block":
            processor = Bed.blocked_iterator( processor )
        elif method == "shift":
            processor = shiftIntervals( processor, contigs, offset=options.offset )
        #IMS: new method: extend intervals by set amount
        elif method == "extend":
            processor = extendInterval( processor, options.extend_distance )
            
    noutput = 0
    for bed in processor:
        options.stdout.write( str(bed) + "\n" )
        noutput += 1

    E.info( "noutput=%i" % (noutput) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
