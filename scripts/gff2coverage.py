################################################################################
#
#   $Id: gff2coverage.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
gff2coverage.py - compute genomic coverage of gff intervals
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Intervals Sequences Statistics Summary


Purpose
-------

This script computes the genomic coverage of intervals 
in a :term:`gff` formatted file. The coverage is computed
per feature.


Usage
-----

Genomic method
++++++++++++++

Compute how many times a feature appears in a gff file and the total number of bases:
$ python gff2coverage.py < <features.gff>

Or the equivalent command (it is called "genomic" method):
$ python gff2coverage.py -m genomic < <features.gff>

If you provide a genome sequence in fasta format, the coverage per feature is also computed:
$ python gff2coverage.py -m genomic -g <genome.fasta> < <features.gff>

Histogram method
++++++++++++++++

On the other hand, it is possible to create a cumulative histogram per feature (option: -f <feature-list>):
$ python gff2coverage.py -m histogram -b <int> -g <genome.fasta> -f <feature_t1,feature_t2,feature_t3> < <features.gff>

With the above options, you are using a sequence file (option: -g <genome.fasta>) and a histogram with <int> bins will be created per feature. Otherwise you may also provide a window size (option: -w <int>):
$ python gff2coverage.py -m histogram -g <genome.fasta> -w <int> -f <feature_t1,feature_t2,feature_t3> < <features.gff>

By default the histogram data is saved to a file called "<contig>.hist" but you can also modify this (option: -o <"%s.name">):
$ python gff2coverage.py -m histogram -g <genome.fasta> -f feature_t1,feature_t2,feature_t3 -o <"%s.name">  < <features.gff>

where "%s" will be replaced by the contig name.


Examples
--------

Output summary statistics for all features in file hg19.chr19.gtf:
$ python gff2coverage.py -g hg19.chr19.fasta < hg19.chr19.gtf

Create an histogram for all the features in the file hg19.chr19.gtf and modify the output file name to "19.all"
$ python scripts/gff2coverage.py -m histogram -g hg19.chr19.fasta -o "%s.all" -f exon,CDS,start_codon,stop_codon < hg19.chr19.gtf


Type
----

   python gff2coverage.py --help

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
import CGAT.GTF as GTF



def printValues( contig, max_size, window_size, values, options ):
    """
    Prints the histogram values in a file for each feature.
    """

    outfile = open( options.output_pattern % contig, "w" )

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
    

def checkHistogramConfig(histogram_length, window_size, num_bins):
    """ 
    Check whether the window size and the number of bins are appropriate to create an histogram.
    If no errors, it returns a list with two values: [window_size, num_bins]
    """

    result = []
    ws = 0
    nb = 0

    if window_size and num_bins:
        if window_size < 1 or window_size > histogram_length:
           raise ValueError( " Selected window_size (%i) is not valid for this histogram (length=%i). Select an appropriate value with option -w " % (window_size, histogram_length) )
        elif num_bins < 1 or num_bins > histogram_length:
           raise ValueError( " Selected num_bins (%i) is not valid for this histogram (length=%i). Select an appropriate value with option -b " % (num_bins, histogram_length) )
        elif num_bins != histogram_length / window_size:
           raise ValueError( " Selected num_bins (%i) is not compatible with window_size (%i) for this histogram (length=%i). Select appropriate values with options -b and -w respectively." % (num_bins, window_size, histogram_length) )
        else:
           ws = window_size
           nb = num_bins
    elif window_size:
        if window_size < 1 or window_size > histogram_length:
           raise ValueError( " Selected window_size (%i) is not valid for this histogram (length=%i). Select an appropriate value with option -w " % (window_size, histogram_length) )
        else:
           ws = window_size
           nb = int(math.ceil( (float(histogram_length) / window_size)))
    elif num_bins:
        if num_bins < 1 or num_bins > histogram_length:
           raise ValueError( " Selected num_bins (%i) is not valid for this histogram (length=%i). Select an appropriate value with option -b " % (num_bins, histogram_length) )
        else:
           ws = int(math.ceil( float(histogram_length) / num_bins))
           nb = num_bins
    else:
        raise ValueError( " Error: you have to select either a window_size (option: -w) or a number of bins (-b) or both of them. " )

    result.append(ws)
    result.append(nb)
    return result


def nonOverlappingIntervals(i1, i2):
    """
    Returns true if two intervals are non-overlapping.
    """
    # IMPORTANT: it has been checked that given a feature called x, x.start points to the previous coordinate where the feature really starts
    # For example, if a feature is defined within range <1000,1200>, then x.start will have a value equal to 999
    # That is why I have added +1 when comparing different features
    return i1.end < i2.start + 1 or i1.start + 1 > i2.end


def processChunk( contig, chunk, options, fasta = None ):
    """
    Computes how the features are distributed along the X axis.
    This function requires segments to be non-overlapping.
    """
    
    if len(chunk) == 0: return

    # check whether there are overlapping features or not
    checked = []
    for feature in chunk:
        checked.append(feature)
	others = [x for x in chunk if x not in checked]
	for otherFeature in others:
	   if nonOverlappingIntervals(feature, otherFeature) == False:
		#print feature
		#print otherFeature
		raise ValueError( " Histogram could not be created since the file contains overlapping features! \n%s\n%s  " % (feature, otherFeature) )
    # clear auxiliary list
    del checked[:]

    # compute min_coordinate, max_coordinate for the histogram
    min_coordinate = 1
    max_coordinate = max( map( lambda x: x.end, chunk) )
    if fasta:
	contig_length = fasta.getLength( contig )
        assert max_coordinate <= contig_length, "maximum coordinate (%i) larger than contig size (%i) for contig %s" % (max_coordinate, contig_length, contig )
        max_coordinate = contig_length
    else:
        min_coordinate = min( map( lambda x: x.start, chunk) )

    histogram_length = max_coordinate - min_coordinate + 1

    [window_size, num_bins] = checkHistogramConfig(histogram_length, options.window_size, options.num_bins)

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
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gff2coverage.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome [default = None]."  )

    parser.add_option( "-o", "--output-pattern", dest="output_pattern", type="string",
                       help="output pattern for the histogram filename. %s is substituted with the contig name [default = \"%s.hist\"]."  )

    parser.add_option( "-f", "--features", dest="features", type="string", action="append",
                       help="a comma separated list of features to collect [default = None]. IMPORTANT: the features in the gff file must be non-overlapping"  )

    parser.add_option( "-m", "--method", dest="method", type="string",
                       help="method to use: genomic | histogram [default = genomic]." )

    parser.add_option( "-w", "--window", dest="window_size", type="int", 
                       help="window size to create the histogram [default = None]." )

    parser.add_option( "-b", "--bins", dest="num_bins", type="int", 
                       help="number of bins to create the histogram [default = None]. IMPORTANT: you need to specify either the number of bins or the window size of both!" )


    parser.set_defaults(
        genome_file = None,
        output_pattern = "%s.hist",
        window_size = None,
        num_bins = None,
        value_format = "%6.4f",
        features = [],
        method = "genomic",
        )

    (options, args) = E.Start( parser )

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
        for entry in GTF.iterator( sys.stdin ):
            intervals[ (entry.contig, entry.source, entry.feature) ] += 1
            bases[ (entry.contig, entry.source, entry.feature) ] += entry.end - entry.start

        options.stdout.write( "contig\tsource\tfeature\tintervals\tbases" )
        if fasta:
            options.stdout.write( "\tcontig_percent_coverage\tgenome_percent_coverage\n" )
            total_genome_size = sum( fasta.getContigSizes( with_synonyms = False ).values() )
        else:
            options.stdout.write( "\n" )

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
                                     
    else:
	raise ValueError( "ERROR in the selected method: " + options.method + "\nValid methods are: \"genomic\" or \"histogram\"\n" )
    

    E.Stop()

