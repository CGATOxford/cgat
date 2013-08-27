################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
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
gff2table.py - compute features for intersection of two gff files
=================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

collect intervals from two gff files and compute features based 
on their intersection. The script is intended to compute properties
for a set of non-overlapping windows.

Transforms:
   * none:        no transform
   * overlap:     overlap between set1 and set2
   * complement:  part of set1 that is not covered by set2
   * third_codon: only takes every third position. Needs frame information
                  in the gff file.

Decorators:
   * GC:            G+C content of intervals
   * count:         number of windows
   * mean-length:   mean length of intervals overlapping with window

Usage
-----

Example::

   python gff2table.py --help

Type::

   python gff2table.py --help

for command line help.

Documentation
-------------

Code
----

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

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome

import CGAT.Intervals as Intervals
import CGAT.Stats as Stats
import CGAT.GTF as GTF

#########################################################
#########################################################
#########################################################
## Decorators: decorate a window through intervals
#########################################################
def decorator_counts( intervals, start, end, contig, fasta ):
    """compute length distribution."""
    d = Stats.DistributionalParameters( map( lambda x: x[1] - x[0], intervals ))
    return d['nval'], str(d)

def decorator_percent_coverage( intervals, start, end, contig, fasta ):
    """compute length of intervals."""
    d = Stats.DistributionalParameters( map( lambda x: x[1] - x[0], intervals ))
    return 100.0 * float(d['sum']) / (end - start), str(d)

def decorator_mean_length( intervals, start, end, contig, fasta ):
    """compute length distribution."""
    d = Stats.DistributionalParameters( map( lambda x: x[1] - x[0], intervals ))
    return d['mean'], str(d)

def decorator_median_length( intervals, start, end, contig, fasta ):
    """compute length distribution."""
    d = Stats.DistributionalParameters( map( lambda x: x[1] - x[0], intervals ))
    return d['median'], str(d)

def decorator_percent_gc( intervals, start, end, contig, fasta ):
    """compute G+C content in intervals.
    """
    l, ngc = 0, 0

    # fetch sequence of the complete window first
    sequence = fasta.getSequence( contig, "+", start, end )
    
    for istart, iend in intervals:
        ngc += len( filter( lambda x: x in "GCgc", sequence[istart-start:iend-start]) )
        l += iend - istart
    
    return 100.0 * ngc / l, None

def decorator_median_score( values, start, end, contig ):
    """compute median of values."""
    d = Stats.DistributionalParameters( values )
    return d['median'], str(d)

def decorator_mean_score( values, start, end, contig ):
    """compute mean of values."""
    d = Stats.DistributionalParameters( values )
    return d['mean'], str(d)

def decorator_stddev_score( values, start, end, contig ):
    """compute stddev of values."""
    d = Stats.DistributionalParameters( values )
    return d['stddev'], str(d)

def decorator_min_score( values, start, end, contig ):
    """compute minumum of values."""
    d = Stats.DistributionalParameters( values )
    return d['min'], str(d)

def decorator_max_score( values, start, end, contig ):
    """compute minumum of values."""
    d = Stats.DistributionalParameters( values )
    return d['max'], str(d)

#########################################################
#########################################################
#########################################################
## Transformers: manipulate a set of intervals
#########################################################

def transform_overlap( start, end, intervals_with_gff ):
    """transform: overlap of intervals in x with y."""
    y = Intervals.combineIntervals( map( lambda x: (x[0], x[1]), intervals_with_gff ) )
    return Intervals.pruneIntervals( y, start, end )

def transform_complement( start, end, intervals_with_gff ):
    y = Intervals.combineIntervals( map( lambda x: (x[0], x[1]), intervals_with_gff ) )
    return Intervals.complementIntervals( y, start, end )

def transform_third_codon( start, end, intervals_with_gff):
    """transform: only return nucleotide positions in window (start, end) 
    that are in third codon position.
    """
    intervals = []
    for istart, iend, gff in intervals_with_gff: 

        if gff.frame == ".":
            raise ValueError("need a frame for third codon positions.")
        
        # frame = nucleotides from start to next codon
        frame = int(gff.frame)

        # to make life easier, convert to 0-based coordinates, 
        # with zero starting at first position in window
        # re-arrange positions on negative strand
        if Genomics.IsNegativeStrand( gff.strand ):
            # convert to negative strand coordinates counting from 0
            coordinate_offset = end
            reverse = True
            istart, iend = end - iend, end - istart
        else: 
            istart, iend = istart - start, iend - start
            reverse = False
            coordinate_offset = start
            
        # make sure that you start on a second codon position and within window
        if istart < 0:
            frame = (frame + istart) % 3
            istart = 0
        if frame != 0: istart -= (3-frame)
        istart += 2

        iend = min( iend, end - start)

        for x in range( istart, iend, 3 ):

            if reverse:
                c = coordinate_offset - x - 1
            else:
                c = coordinate_offset + x
            intervals.append( (c, c+1) )

    return Intervals.combineIntervals( intervals )

def test_transform_third_codon():

    def test_entry( frame, strand, xfrom, xto, start, end, ref ):

        entry = GTF.Entry()
        entry.frame = frame
        entry.strand = strand
        entry.start = xfrom
        entry.end = xto
        
        intervals = transform_third_codon( start, end, [ (xfrom, xto, entry) ] )
        if ref != intervals:
            print "failed:", ref != intervals
        
    test_entry( 0, "+", 1, 7, 0, 6, [(3,4)] )
    test_entry( 0, "-", 1, 7, 0, 6, [(1,2), (4,5)])
    test_entry( 1, "+", 1, 7, 0, 6, [(1,2), (4,5)])
    test_entry( 2, "+", 1, 7, 0, 6, [(2,3), (5,6)])
    test_entry( 1, "-", 1, 7, 0, 6, [(3,4)] )
    test_entry( 2, "-", 1, 7, 0, 6, [(2,3), (5,6)] )

    sys.exit(0)

def annotateWindows( contig, windows, gff_data, fasta, options ):
    """annotate windows."""

    index = IndexedGenome.IndexedGenome()
    for g in gff_data:
        index.add( g.contig, g.start, g.end, g )
    
    is_gtf = options.is_gtf

    if options.transform == "none":
        transform = lambda x,y,z: map( lambda x: (x[0], x[1]), z)
    elif options.transform == "overlap":
        transform = transform_overlap
    elif options.transform == "complement":
        transform = transform_complement
    elif options.transform == "third_codon":
        transform = transform_third_codon
    else:
        raise ValueError("unknown transform %s" % options.transform)

    work_on_intervals = True
    if options.decorator == "counts":
        decorator = decorator_counts
    elif options.decorator == "mean-length":
        decorator = decorator_mean_length
    elif options.decorator == "median-length":
        decorator = decorator_median_length
    elif options.decorator == "percent-coverage":
        decorator = decorator_percent_coverage
    elif options.decorator == "gc":
        decorator = decorator_percent_gc
    elif options.decorator == "median-score":
        decorator = decorator_median_score
        work_on_intervals = False
    elif options.decorator == "mean-score":
        decorator = decorator_mean_score
        work_on_intervals = False
    elif options.decorator == "stddev-score":
        decorator = decorator_stddev_score
        work_on_intervals = False
    elif options.decorator == "min-score":
        decorator = decorator_min_score
        work_on_intervals = False
    elif options.decorator == "max-score":
        decorator = decorator_max_score
        work_on_intervals = False
    else:
        raise ValueError("unknown decorator %s" % options.decorator)

    for start,end in windows:
        
        # counts/length before/after transformation
        n1, l1, n2, l2 = 0, 0, 0, 0

        values, intervals_with_gff, genes, transcripts = [], [], set(), set()
        
        try:
            for istart, iend, value in index.get( contig, start, end ):
                n1 += 1
                l1 += iend - istart
                intervals_with_gff.append( (istart, iend, value) )
                values.append( value.score )
                if is_gtf:
                    genes.add( value.gene_id )
                    transcripts.add( value.transcript_id )
        except KeyError:
            pass
        
        if n1 == 0 and options.skip_empty: continue

        if work_on_intervals:

            if options.loglevel >= 3:
                options.stdlog.write("# intervals in window %i:%i before transformation: %s\n" % (start, end, str(intervals)))

            intervals = transform( start, end, intervals_with_gff )

            for xstart, xend in intervals:
                n2 += 1
                l2 += xend - xstart

            if options.loglevel >= 3:
                options.stdlog.write("# intervals in window %i:%i after transformation: %s\n" % (start,end, str(intervals)))

            score, extra_info = decorator( intervals, start, end, contig, fasta )
            
        else:
            if len(values) > 0:
                values = map(float, values)
                score, extra_info = decorator( values, start, end, contig )
            else:
                score, extra_info = 0, None
                
            l2 = 0
            n2 = 0

        if is_gtf:
            ngenes, ntranscripts = len(genes), len(transcripts)
        else:
            ngenes, ntranscripts = 0, 0

        if extra_info:
            extra_info = re.sub("\t", ";", extra_info )
        options.stdout.write( "\t".join( \
            map(str, (contig, start, end,
                      ngenes, ntranscripts,
                      n1, l1,
                      n2, l2,
                      score,
                      extra_info) ) ) + "\n" )
                      
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gff2table.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"])

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome (indexed)."  )

    parser.add_option( "-w", "--filename-windows", dest="filename_windows", type="string",
                       help="gff file with windows to use."  )

    parser.add_option( "-d", "--filename-data=", dest="filename_data", type="string",
                       help="gff file with data to use." )

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                       help="filename-data is gtf file [default=%default." )

    parser.add_option( "-f", "--features=", dest="features", type="choice", action="append",
                       choices=("GC", ),
                       help="features to compute." )

    parser.add_option( "-c", "--decorator=", dest="decorator", type="choice", 
                       choices=( "counts", "gc", "gc3", "mean-length", "median-length", "percent-coverage",
                                 "median-score", "mean-score", "stddev-score", "min-score", "max-score" ),
                       help="decorators to use." )

    parser.add_option( "-e", "--skip-empty", dest="skip_empty", action="store_true",
                       help="skip empty windows." )

    parser.add_option( "-t", "--transform=", dest="transform", type="choice", 
                       choices=("none", "overlap", "complement", "third_codon"),
                       help="transform to use when mapping overlapping regions onto window." )


    parser.set_defaults(
        genome_file = None,
        filename_windows = None,
        filename_data = None,
        features = [],
        skip_empty = False,
        decorator = "counts",
        transform = "none",
        is_gtf = False,
        )

    (options, args) = E.Start( parser )

    #    test_transform_third_codon()
    
    if not options.filename_windows:
        raise ValueError("please supply a gff file with window information." )
    
    if options.loglevel >= 1:
        options.stdlog.write("# reading windows..." )
        options.stdlog.flush()
        
    windows = GTF.readAsIntervals( GFF.iterator( IOTools.openFile(options.filename_windows, "r" ) ) )

    if options.loglevel >= 1:
        options.stdlog.write("done\n" )
        options.stdlog.flush()

    if options.filename_data:
        if options.loglevel >= 1:
            options.stdlog.write("# reading data..." )
            options.stdlog.flush()

        if options.is_gtf:
            gff_data = GTF.readFromFile( IOTools.openFile( options.filename_data, "r" ) )
        else: 
            gff_data = GTF.readFromFile( IOTOols.openFile( options.filename_data, "r" ) )

        if options.loglevel >= 1:
            options.stdlog.write("done\n" )
            options.stdlog.flush()
        
        data_ranges = GTF.SortPerContig( gff_data )
    else:
        ## use windows to compute properties
        ## by supplying no data and asking for the complement = original window
        gff_data = None
        data_ranges = None
        options.transform = "complement"

    map_contig2size = {}

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        map_contig2size = fasta.getContigSizes()
    else:
        for contig, values in windows.items():
            map_contig2size[contig] = max( lambda x:x[1], values )
        fasta = None

    contigs = map_contig2size.keys()
    contigs.sort()

    ## proceed contig wise
    noutput_contigs, ncontigs_skipped_windows, ncontigs_skipped_data = 0, 0, 0

    options.stdout.write( "\t".join( \
        map(str, ("contig", "start", "end",
                  "ngenes", "ntranscripts",
                  "n1", "l1",
                  "n2", "l2",
                  "score",
                  "extra_info") ) ) + "\n" )
    
    for contig in contigs:
        
        skip = False
        if contig not in windows: 
            ncontigs_skipped_windows += 1
            skip = True

        if data_ranges and contig not in data_ranges: 
            ncontigs_skipped_data += 1
            skip = True
            
        if skip: continue
        
        noutput_contigs += 1 
        if data_ranges:
            annotateWindows( contig,
                             windows[contig],
                             gff_data[data_ranges[contig][0]:data_ranges[contig][1]],
                             fasta,
                             options )
        else:
            annotateWindows( contig,
                             windows[contig],
                             [],
                             fasta,
                             options )
            

    if options.loglevel >= 1:
        options.stdout.write( "# ninput_windows=%i, noutput_contigs=%i, ninput_contigs=%i, nskipped_windows=%i, nskipped_data=%i\n" %\
                                  ( len(windows), noutput_contigs, len(contigs), ncontigs_skipped_windows, ncontigs_skipped_data ) )

    E.Stop()


