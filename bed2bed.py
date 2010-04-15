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

import sys, re, string, optparse, time, os, itertools, tempfile, subprocess, shutil

import Experiment as E
import alignlib
import Stats
import GFF, GTF
import IndexedFasta, IOTools
import Bed

def merge( iterator, max_distance = 0 ):
    """iterator for merging adjacent bed entries.

    *max_distance* > 0 permits merging of intervals that are
    not directly adjacent.
    """

    def iterate_chunks( iterator ):
        last = iterator.next()
        to_join = [last]

        for bed in iterator:
            d = bed.start - last.end
            if bed.contig == last.contig:
                assert bed.start >= last.start, "input file should be sorted by contig and position: d=%i:\n%s\n%s\n" % (d, last, bed)
            
            if bed.contig != last.contig or \
                    d > max_distance:
                yield to_join
                to_join = []

            last = bed
            to_join.append( bed )

        if to_join: yield to_join
        raise StopIteration

    ninput, noutput = 0, 0

    for to_join in iterate_chunks(iterator):

        ninput += 1
        a = to_join[0]
        a.end = to_join[-1].end
        yield a
        noutput += 1

    E.info( "ninput=%i, noutput=%i\n" % (ninput, noutput) )


def filterGenome( iterator, contigs ):
    """remove bed intervals that are outside of contigs.
    
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
        yield bed

    E.info( "ninput=%i, noutput=%i, nskipped_contig=%i, nskipped_range=%i" % \
                (ninput, noutput, nskipped_contig, nskipped_range) )

    
def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--method", dest="methods", type="choice", action="append",
                       choices=("merge", "filter-genome", "bins" ),
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

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.set_defaults( methods = [],
                         merge_distance = 0,
                         binning_method = "equal-bases",
                         num_bins = 5,
                         bin_edges = None,
                         test = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = genome_fasta.getContigSizes()


    processor = Bed.iterator( options.stdin )

    for method in options.methods:
        if method ==  "filter-genome":
            processor = filterGenome( processor, contigs )
        elif method == "merge":
            processor = merge( processor, options.merge_distance )
        elif method == "bins":
            if options.bin_edges: bin_edges = map(float, options.bin_edges.split(","))
            else: bin_edges = None
            processor, bin_edges = Bed.binIntervals( processor,
                                                     num_bins = options.num_bins,
                                                     method = options.binning_method,
                                                     bin_edges = bin_edges)
            E.info("# split bed: bin_edges=%s" % (str(bin_edges)))
                
    noutput = 0
    for bed in processor:
        options.stdout.write( str(bed) + "\n" )
        noutput += 1

    E.info( "noutput=%i\n" % (noutput) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
