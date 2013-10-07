'''
gff2coverage.py - compute genomic coverage of gff intervals
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Summary GFF

Purpose
-------

This script computes the genomic coverage of intervals 
in a :term:`gff` formatted file. The coverage is computed
per feature.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Command line options
--------------------

'''

import sys
import string
import re
import optparse
import time
import os
import shutil
import tempfile
import math
import collections

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF

def printValues( contig, max_size, window_size, values, options ):
    """output values."""

    outfile = E.openOutputFile( contig, "w" )

    outfile.write( "abs_pos\trel_pos" )

    for feature in options.features:
        outfile.write("\tabs_%s\trel_%s" % (feature,feature) )
    outfile.write("\n")

    max_vv = []

    for f in range(len(options.features)):
        max_vv.append( float(max(map(lambda x: x[f], values ) )))

    bin = 0
    for vv in values:
        outfile.write( "%i\t" % bin )
        outfile.write( options.value_format % (float(bin) / max_size) )
        
        for x in range(len(options.features)):
            outfile.write( "\t%i\t%s" % (vv[x] ,
                                         options.value_format % (vv[x] / max_vv[x]) ) )
        outfile.write("\n" )            
        bin += window_size

    outfile.close()
    
def processChunk( contig, chunk, options, fasta = None ):
    """
    This function requires segments to be non-overlapping.
    """
    
    if len(chunk) == 0: return

    # check whether there are overlapping features or not
    checked = []
    for feature in chunk:
        checked.append(feature)
	others = [x for x in chunk if x not in checked]
	for otherFeature in others:
	   if GTF.Overlap(feature, otherFeature):
		raise ValueError( " Histogram could not be created since the file contains overlapping features! \n%s\n%s  " % (feature, otherFeature) )
    # clear auxiliary list
    del checked[:]

    # compute min_coordinate, max_coordinate for the histogram
    min_coordinate = 1
    max_coordinate = max( map( lambda x: x.end, chunk) )
    ## compute window size
    if options.window_size:
        window_size = options.window_size
        num_bins = int(math.ceil((float(max_coordinate) / window_size)))
    elif options.num_bins and fasta:
        contig_length = fasta.getLength( contig )
        assert max_coordinate <= contig_length, "maximum coordinate (%i) larger than contig size (%i) for contig %s" % (max_coordinate, contig_length, contig )
        max_coordinate = contig_length
        window_size = int(math.floor( float(contig_length) / options.num_bins))
        num_bins = options.num_bins
    else:
        raise ValueError("please specify a window size of provide genomic sequence with number of bins.")

    values = [ [] for x in range(num_bins) ]

    ## do several parses for each feature, slow, but easier to code
    ## alternatively: sort by feature and location.
    for feature in options.features:
        total = 0
        bin = 0
        end = window_size
        for entry in chunk:
            if entry.feature != feature: continue

            while end < entry.start:
                values[bin].append( total )
                bin += 1
                end += window_size

            while entry.end > end:
                seg_start = max(entry.start, end - window_size)
                seg_end = min(entry.end, end)
                total += seg_end - seg_start
                values[bin].append( total )                
                end += window_size
                bin += 1
            else:
                seg_start = max(entry.start, end - window_size)
                seg_end = min(entry.end, end)
                total += seg_end - seg_start

        while bin < num_bins:
            values[bin].append(total)
            bin += 1

    printValues( contig, max_coordinate, window_size, values, options )
    
def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2coverage.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome [default=%default]"  )

    parser.add_option( "-f", "--features", dest="features", type="string", action="append",
                       help="features to collect " 
                       "[default=%default]" )

    parser.add_option( "-w", "--window-size", dest="window_size", type="int", 
                       help="window size in bp for histogram computation. "
                       "Determines the bin size.  "
                       "[default=%default]" )

    parser.add_option( "-b", "--num-bins", dest="num_bins", type="int", 
                       help="number of bins for histogram computation "
                       "if window size is not given. "
                       "[default=%default]" )

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       choices = ("genomic", "histogram", ),
                       help="methods to apply. "
                       "[default=%default]"  )

    parser.set_defaults(
        genome_file = None,
        window_size = None,
        num_bins = 1000,
        value_format = "%6.4f",
        features = [],
        method = "genomic",
        )

    (options, args) = E.Start( parser, add_output_options = True )

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None
        
    if options.method == "histogram":

        gff = GTF.readFromFile( sys.stdin )

        gff.sort( lambda x,y: cmp( (x.contig, x.start), (y.contig, y.start) ) )

        chunk = []
        last_contig = None

        for entry in gff:

            if last_contig != entry.contig:
                processChunk( last_contig, chunk, options, fasta )
                last_contig = entry.contig
                chunk = []

            chunk.append( entry )

        processChunk( last_contig, chunk, options, fasta )

    elif options.method == "genomic":
        intervals = collections.defaultdict( int )
        bases = collections.defaultdict( int )
        total = 0
        for entry in GTF.iterator( sys.stdin ):
            intervals[ (entry.contig, entry.source, entry.feature) ] += 1
            bases[ (entry.contig, entry.source, entry.feature) ] += entry.end - entry.start
            total += entry.end - entry.start

        options.stdout.write( "contig\tsource\tfeature\tintervals\tbases" )
        if fasta:
            options.stdout.write( "\tpercent_coverage\ttotal_percent_coverage\n" )
        else:
            options.stdout.write( "\n" )

        total_genome_size = sum( fasta.getContigSizes( with_synonyms = False ).values() )

        for key in sorted (intervals.keys() ):
            nbases = bases[key]
            nintervals = intervals[key]
            contig, source, feature = key
            options.stdout.write( "\t".join( ( "\t".join(key),
                                     str(nintervals),
                                     str(nbases) ) ) )
            if fasta: 
                options.stdout.write( "\t%f" % ( 100.0 * float(nbases) / fasta.getLength( contig )))
                options.stdout.write( "\t%f\n" % ( 100.0 * float(nbases) / total_genome_size ))
            else: options.stdout.write( "\n" )
                                     
            
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


