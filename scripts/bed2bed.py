'''
bed2bed.py - manipulate bed files
=================================

Purpose 
-------

This script provides various methods for merging (by position, by name
or by score), filtering and moving bed formatted intervals and
outputting the results as a bed file

Usage 
-----

   python bed2bed.py --method=[METHOD] [options]

This script provides several methods:

merge 
+++++

Merge together overlapping or adjacent intervals. The basic
functionality is similar to bedtools merge, but with some additions: 

* Merging by name: specifying the --merge-by-name option will mean
  that only overlaping (or adjacent intervals) with the same value in
  the 4th column of the bed will be merged 

* Removing overlapping intervals with inconsistent names: set the
   ``--remove-inconsistent`` option.

.. caution:: 
   Intervals of the same name will only be merged if they 
   are consecutive in the bed file.

* Only output merged intervals: By specifiying the --merge-min-intervals=n
  options, only those intervals that were created by merging at least n
  intervals together will be output

bins 
++++

Merges together overlapping or adjecent intervals only if they have
"similar" scores. Score similarity is assessed by creating a number of
score bins and assigning each interval to a bin. If too adjacent
intervals are in the same bin, the intervals are merged. Note that in
contrast to merge-by-name above, two intervals do not need to be
overlapping or within a certain distance to be merged.

There are several methods to create the bins: 

* equal-bases: Bins are created to that they contain the same number of bases.  Specified
  by passing "equal-bases" to --binning-method. This is the default.  

* equal-intervals: Score bins are create so that each bin contains the
  same number of intervals. Specified by passing "equal-intervals" to
  --binning-method.  

* equal-range: Score bins are created so that
  each bin covers the same fraction of the total range of
  scores. Specified by passing "equal-range" to --binning-method.  

* bin-edges: Score binds can be specified by manually passing a comma
  seperated list of bin edges to --bin-edges.

The number of bins is specified by the --num-bins options, and the
default is 5.

block 
+++++

Creates blocked bed12 outputs from a bed6, where intervals with the
same name are merged together to create a single bed12 entry.

.. Caution:: Input must be sorted so that entries of the same
name are together.

filter-genome 
+++++++++++++

Removes intervals that are on unknown contigs or extend off the 3' or
5' end of the contig.  Requires a tab seperated input file to -g which
lists the contigs in the genome, plus their lengths.

sanitize-genome 
+++++++++++++++

As above, but instead of removing intervals overlapping the ends of
contigs, truncates them.  Also removes empty intervals.

shift 
+++++

Moves intervals by the specified amount, but will not allow them to be
shifted off the end of contigs. Thus if a shift will shift the start
of end of the contig, the interval is only moved as much as is
possible without doing this.

Command line options
--------------------
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

def merge( iterator, 
           max_distance = 0, 
           by_name = False, 
           min_intervals = 1,
           remove_inconsistent = False ):
    """iterator for merging adjacent bed entries.

    *max_distance* > 0 permits merging of intervals that are
    not directly adjacent.

    If *by_name = True*, only entries with the same name are merged.
    
    If *remove_inconsistent*, overlapping intervals where the names
    are inconsistent will be removed.
    
    The score gives the number of intervals that have been merged.
    """

    if remove_inconsistent and by_name:
        assert ValueError( "using both remove_inconsistent and by_name makes no sense" )

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

    c = E.Counter()

    for to_join in iterate_chunks(iterator):

        c.input += 1
        if len(to_join) < min_intervals: 
            c.skipped_min_intervals += 1
            continue

        if remove_inconsistent:
            names = set([ x.name for x in to_join ] )
            if len(names) > 1: 
                c.skipped_inconcistent_intervals += 1
                continue

        a = to_join[0]
        a.end = to_join[-1].end
        a.score = len(to_join)
        yield a
        c.output += 1

    E.info( str(c) )

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
        #IMS: add filtering for filtering <0 co-ordinates
        if bed.start < 0 or bed.end < 0:
            nskipped_range +=1
            continue
        # should this not be just >, as co-ordinates are half-closed, so 
        # if end = contigs[bed.contig], then interval ends on last base?
        if bed.end > contigs[bed.contig]:
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
        #IMS: changing >= to > in if statement: next line sets bed.end = contigs[bed.contig]
        # this shouldn't count as a truncation.
        if bed.end > contigs[bed.contig]:
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
        #IMS: if we skip intervals off the end of the contig we should skipp ones
        # off the start as well
        if bed.start < 0 or bed.end < 0:
            nskipped_range += 1
            continue
        #IMS: changing >= to > as bed is half-open
        if bed.end > contigs[bed.contig]:
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

        if bed.contig not in contigs:
            nskipped_contig += 1
            continue
        if bed.start <0 or bed.end < 0:
            nskipped_range += 1
            continue
        if bed.end > contigs[bed.contig]:
            nskipped_range += 1
            continue
        
        newstart = bed.start - distance
        newend = bed.end + distance

        if newstart < 0:
            newstart = 0

        if newend > contigs[bed.contig]:
            newend = contigs[bed.contig]

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

    parser.add_option( "--remove-inconsistent", dest="remove_inconsistent", action="store_true",
                      help="when merging, do not output intervals where the names of overlapping intervals "
                      "do not match [default=%default]"  )

    parser.add_option( "--offset", dest="offset",  type="int",
                      help="offset for shifting intervals [default=%default]"  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-b", "--bam-file", dest="bam_file", type="string",
                      help="bam-formatted filename with genome."  )
   
    parser.set_defaults( methods = [],
                         merge_distance = 0,
                         binning_method = "equal-bases",
                         merge_by_name = False,
                         genome_file = None,
                         bam_file = None,
                         num_bins = 5,
                         merge_min_intervals = 1,
                         bin_edges = None,
                         offset = 10000,
                         test = None,
                         extend_distance=1000,
                         remove_inconsistent = False)
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    contigs = None

    ## Why provide full indexed genome, when a tsv of contig sizes would do? 
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
            processor = merge( processor, 
                               options.merge_distance,
                               by_name = options.merge_by_name,
                               min_intervals = options.merge_min_intervals,
                               remove_inconsistent = options.remove_inconsistent )
        elif method == "bins":
            if options.bin_edges:
                bin_edges = map(float, options.bin_edges.split(","))
                #IMS: check bin edges are valid
                if not(len(bin_edges) == options.num_bins +1): raise ValueError(
                    "Number of bin edge must be one more than number of bins")
            else: bin_edges = None
            processor, bin_edges = Bed.binIntervals( processor,
                                                     num_bins = options.num_bins,
                                                     method = options.binning_method,
                                                     bin_edges = bin_edges)
            E.info("# split bed: bin_edges=%s" % (str(bin_edges)))
            
        elif method == "block":
            processor = Bed.blocked_iterator( processor )
        elif method == "shift":
            #IMS: test that contig sizes are availible
            if not contigs: raise ValueError("please supply genome file")
            processor = shiftIntervals( processor, contigs, offset=options.offset )
        #IMS: new method: extend intervals by set amount
        elif method == "extend":
            if not contigs: raise ValueError("please supply genome file")
            processor = extendInterval( processor, contigs, options.offset )
            
    noutput = 0
    for bed in processor:
        options.stdout.write( str(bed) + "\n" )
        noutput += 1

    E.info( "noutput=%i" % (noutput) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
