################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
bam2peakshape.py - compute peak shape features from a bam-file
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a :term:`bed` formatted file with regions of interest,
for example binding intervals from a ChIP-Seq experiment. Using a collection
of aligned reads is a :term:`bam` formatted file, the script outputs a collection
of features describing the peak shape.

Among the features output are:

peak-height - number of reads at peak
peak-median - median coverage compared to peak height within interval
interval-width - width of interval
peak-max-width - width of peak. Distance to peak of furthest half-height position.
peak-min-width - width of peak. Distance to peak of closest half-height position.

histogram - a histogram of read depth within the interval. 

.. todo::
   
   paired-endedness is not fully implemented.

-c/--control-file
   if given, two peakshapers are computed, one for the primary and one for the control
   :term:`bam` file. The control file is centered around the same
   base as the primary file and output in the same sort order as the primary profile.

Usage
-----

Example::

   python script_template.py --help

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

For read counts to be correct the NH flag to be set correctly.

Code
----

'''

import os, sys, re, optparse, collections
import Experiment as E
import IOTools
import pysam
import Bed
import numpy

import pyximport
pyximport.install(build_in_temp=False)
import _bam2peakshape

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-o", "--only-interval", dest="only_interval", action="store_true",
                       help = "only count tags strictly in interval. Otherwise, use window around peak center "
                              " [%default]" )

    parser.add_option( "-w", "--window-size", dest="window_size", type = "int",
                       help = "window size to use"
                              "[%default]" )

    parser.add_option( "-b", "--bin-size", dest="bin_size", type = "int",
                       help = "bin-size for histogram of read depth. "
                              "[%default]" )

    parser.add_option( "-s", "--sort", dest="sort", type = "choice", action = "append",
                       choices = ("peak-height", "peak-width", "unsorted", "interval-width", "interval-score" ),
                       help = "output sort order for matrices. "
                              "[%default]" )

    parser.add_option( "-i", "--shift", dest="shift", type = "int",
                       help = "shift for reads. Reads will be shifted upstream/downstream by this amount. "
                              "[%default]" )

    parser.add_option( "-c", "--control-file", dest="control_file", type = "string",
                       help = "control file "
                              "[%default]" )

    parser.add_option( "-r", "--random-shift", dest="random_shift", action="store_true",
                       help = "shift intervals in random direction directly up/downstream of interval "
                              "[%default]" )

    parser.add_option( "-e", "--centring-method", dest="centring_method", type = "choice",
                       choices = ("reads", "middle"),
                       help = "centring method (reads=use reads to determine peak, middle=use middle of interval" 
                              "[%default]" )
                              
    parser.add_option( "-n", "--normalization", dest="normalization", type = "choice",
                       choices = ("none", "sum" ),
                       help = "normalisation to perform. "
                              "[%default]" )                                 

    parser.set_defaults(
        remove_rna = False,
        ignore_pairs = False,
        input_reads = 0,
        force_output = False,
        bin_size = 10,
        shift = 0,
        window_size = 1000,
        sort = [],
        centring_method = "reads",
        control_file = None,
        random_shift = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if len(args) != 2:
        raise ValueError("please specify a bam and bed file" )

    bamfile, bedfile = args

    pysam_in = pysam.Samfile( bamfile, "rb" )
    shift = options.shift

    if options.control_file:
        E.info("using control file %s" % options.control_file)
        pysam_control = pysam.Samfile( options.control_file, "rb" )
    else:
        pysam_control = None

    options.stdout.write( "\t".join( ("contig", 
                                      "start",
                                      "end",
                                      "name",
                                      "\t".join(_bam2peakshape.PeakShapeResult._fields) ) ) + "\n" )
    
    if options.window_size:
        # bins are centered at peak-center and then stretching outwards.
        bins = numpy.arange( -options.window_size + options.bin_size // 2, 
                              +options.window_size, 
                              options.bin_size )        
        
    contigs = set(pysam_in.references)

    result =[]
    c = E.Counter()
    for bed in Bed.iterator( IOTools.openFile( bedfile ) ):
        c.input += 1

        if bed.contig not in contigs: 
            c.skipped += 1
            continue
        
        features = _bam2peakshape.count( pysam_in, bed.contig, bed.start, bed.end, 
                                         shift = shift,
                                         window_size = options.window_size,
                                         bins = bins,
                                         only_interval = options.only_interval,
                                         centring_method = options.centring_method )

        if pysam_control:
            control = _bam2peakshape.countAroundPos(pysam_control, 
                                                    bed.contig,
                                                    features.peak_center,
                                                    shift = shift,
                                                    bins = features.bins )
        else:
            control = None

        if options.random_shift:
            direction = numpy.random.randint( 0, 2 )
            if direction: pos = features.peak_center + 2 * bins[0]
            else: pos = features.peak_center + 2 * bins[-1]
            shifted = _bam2peakshape.countAroundPos(pysam_in, 
                                                    bed.contig,
                                                    pos,
                                                    shift = shift,
                                                    bins = features.bins )
        else:
            shifted = None

        result.append( (features, bed, control, shifted) )
        c.added += 1

    E.info( "interval processing: %s" % c )

    # center bins
    out_bins = bins[:-1] + options.bin_size

    def writeMatrix( result, sort ):

        outfile_matrix = E.openOutputFile( "matrix_%s.gz" % re.sub("-", "_", sort ) )            
        outfile_matrix.write( "name\t%s\n" % "\t".join( map(str, out_bins )))
        
        if result[0][2] != None:
            outfile_control = E.openOutputFile( "control_%s.gz" % re.sub("-", "_", sort ) )            
            outfile_control.write( "name\t%s\n" % "\t".join( map(str, out_bins )))

        if result[0][3] != None:
            outfile_shift = E.openOutputFile( "shift_%s.gz" % re.sub("-", "_", sort ) )            
            outfile_shift.write( "name\t%s\n" % "\t".join( map(str, out_bins )))
            
        n = 0
        for features, bed, control, shifted in result:
            n += 1
            if "name" in bed: name = bed.name
            else:name = str(n)
            bins, counts = features[-2], features[-1]
            outfile_matrix.write( "%s\t%s\n" % (name, "\t".join(map(str,counts))) )
            if control:
                outfile_control.write( "%s\t%s\n" % (name, "\t".join(map(str,control.counts))) )
            if shifted:
                outfile_shift.write( "%s\t%s\n" % (name, "\t".join(map(str,shifted.counts))) )

        outfile_matrix.close()

    n = 0
    for features, bed, control, shifted in result:
        n += 1
        if "name" in bed: name = bed.name
        else: name = str(n)
        options.stdout.write( "%s\t%i\t%i\t%s\t" % (bed.contig, bed.start, bed.end, name ) )

        options.stdout.write( "\t".join( map(str, features[:-2] ) ) )
        bins, counts = features[-2], features[-1]
        options.stdout.write( "\t%s" % ",".join( map(str, bins )))
        options.stdout.write( "\t%s" % ",".join( map(str, counts)))
        options.stdout.write( "\n" )

    norm_result = []
    if options.normalization == "sum":
        E.info("Starting sum normalization")
        # get total counts across all intervals
        norm = 0.0
        for features, bed, control, shifted in result:
            counts = features[-1]
            norm += sum(counts)
        norm /= 1000000
        E.info("norm = %i" % norm)
        
        # normalise
        for features, bed, control, shifted in result:
            counts = features[-1]
            norm_counts = []
            for c in counts:
                norm_counts.append( c/(norm) )
            new_features = features._replace( counts=norm_counts )
            norm_result.append( (new_features, bed, control, shifted) )
    else:
        E.info("No normalization performed")
        norm_result = result
        
    # output sorted matrices
    if not options.sort: writeMatrix( norm_result, "unsorted" )

    for sort in options.sort: 

        if sort == "peak-height":
            norm_result.sort( key = lambda x: x[0].peak_height )
            
        elif sort == "peak-width":
            norm_result.sort( key = lambda x: x[0].peak_width )

        elif sort == "interval-width":
            norm_result.sort( key = lambda x: x[1].end - x[1].start )

        elif sort == "interval-score":
            try:
                norm_result.sort( key = lambda x: x[1].score )
            except IndexError:
                E.warn("score field not present - no output" )

        writeMatrix( norm_result, sort )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


