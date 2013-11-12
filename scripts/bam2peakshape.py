'''
bam2peakshape.py - compute peak shape features from a bam-file
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Intervals BAM BED Summary

Purpose
-------

This script takes a :term:`bed` formatted file with regions of interest,
for example binding intervals from a ChIP-Seq experiment. Using a collection
of aligned reads is a :term:`bam` formatted file, the script outputs a collection
of features describing the peak shape.

This script is designed with a slight emphasis on ChIP-Seq datasets. For example,
for the purpose of visualizing binding profiles of transcription factors
ChIP-Seq purely based on the peaks regions defined by MACS, and without the need to use any genomic
annotations (e.g. ENSEMBL, refseq). The reason is that, bam2peakshape.py is able to
center the counting window at the summit of every individual peak. bam2peakshape.py is also able
to: (1) plot the control RNA-Seq library to enable side-by-side comparison of treatment vs control;
(2) randomly shift the set of input resions to generate a background set of regions in order to
provide a peaks profile that can be used as the control.

Alternatively,  you may also consider using :doc:`bam2geneprofile`, which is
designed with a slight emphasis on RNA-Seq datasets, which takes care of spliced reads,
by using the CIGAR string in the BAM file to accurately define the covered bases (when
the --base-accurate-off is not specified, currently it is not specified by default).

The script outputs a tab-separated table on stdout containing features for each interval.
A peak is defined as the location of the highest density in an interval. The width of 
the peak (peak_width) is defined as the region around the peak in which the density does
not drop below a threshold of peak_heigt * 90%.

Among the features output are:

+---------------------+----------------------------------------------------------+
|*Column*             |*Content*                                                 |
+---------------------+----------------------------------------------------------+
|peak_height          |number of reads at peak                                   |
+---------------------+----------------------------------------------------------+
|peak_median          |median coverage compared to peak height                   | 
+---------------------+----------------------------------------------------------+
|interval_width       |width of interval                                         |
+---------------------+----------------------------------------------------------+
|peak_width           |width of peak                                             |
+---------------------+----------------------------------------------------------+
|bins                 |bins for a histogram of densities within the interval.    |
+---------------------+----------------------------------------------------------+
|npeaks               |number of density peaks in interval.                      |
+---------------------+----------------------------------------------------------+
|peak_center          |point of highest density in interval                      |
+---------------------+----------------------------------------------------------+
|peak_relative_pos    |point of highest density in interval coordinates          |
+---------------------+----------------------------------------------------------+
|counts               |counts for a histogram of densities within the interval   |   
+---------------------+----------------------------------------------------------+
|furthest_half_height |Distance of peak center to furthest half-height position  |
+---------------------+----------------------------------------------------------+
|closest_half_height  |Distance of peak center to closest half-height position   |
+---------------------+----------------------------------------------------------+

Additionally, the script outputs a set of matrixes with densities over intervals that
can be used for plotting. The default filenames are ``(matrix|control)_<sortorder>.tsv.gz``,
The names can be controlled with the ``--output-filename-pattern`` option.

.. todo::
   
   paired-endedness is not fully implemented.

Usage
-----

Example::

   python bam2peakshape.py bamfile bedfile

Type::

   python bam2peakshape.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.Bed as Bed
import numpy
import bx.bbi.bigwig_file

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2peakshape
except ImportError:
    import CGAT._bam2peakshape as _bam2peakshape

def buildOptionParser( argv ):
    
    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-f", "--format", dest="format", type = "choice",
                       choices = ("bam", "bigwig" ),
                       help = "format of genomic input files for densities "
                              "[%default]" )

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
                       choices = ("peak-height", "peak-width", "unsorted", 
                                  "interval-width", "interval-score" ),
                       help = "output sort order for matrices. "
                              "[%default]" )

    parser.add_option( "-c", "--control-file", dest="control_file", type = "string",
                       help = "control file. If given, two peakshapes are computed, "
                       "one for the primary data and one for the control data. "
                       "The control file is centered around the same "
                       "base as the primary file and output in the same "
                       "sort order as the primary profile to all side-by-side. "
                       "comparisons. "
                       "[%default]" )

    parser.add_option( "-r", "--random-shift", dest="random_shift", action="store_true",
                       help = "shift intervals in random direction up/downstream of interval "
                              "[%default]" )
 
    parser.add_option( "-e", "--centring-method", dest="centring_method", type = "choice",
                       choices = ("reads", "middle"),
                       help = "centring method. Available are: "
                              "reads=use density to determine peak, "
                              "middle=use middle of interval " 
                              "[%default]" )
                              
    parser.add_option( "-n", "--normalization", dest="normalization", type = "choice",
                       choices = ("none", "sum" ),
                       help = "normalisation to perform. "
                              "[%default]" )                                 

    parser.add_option( "--use-strand", dest="strand_specific", action="store_true",
                       help = "use strand information in intervals. Intervals on the "
                       "negative strand are flipped "
                       "[%default]" )
 
    parser.add_option( "-i", "--shift", dest="shift", type = "int",
                       help = "shift for reads. When processing bam files, "
                       "reads will be shifted upstream/downstream by this amount. "
                       "[%default]" )

    parser.set_defaults(
        bin_size = 10,
        shift = 0,
        window_size = 1000,
        sort = [],
        centring_method = "reads",
        control_file = None,
        random_shift = False,
        strand_specific = False,
        format = "bam",
        report_step = 100,
        )

    return parser

def outputResults( result, bins, options ):
    '''ouput results from density profiles.'''

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
                norm_result.sort( key = lambda x: float(x[1].score) )
            except IndexError:
                E.warn("score field not present - no output" )

        writeMatrix( norm_result, sort )


def buildResults( bedfile, fg_file, control_file, counter, options ):
    '''compute densities and peakshape parameters.'''

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
        
    #contigs = set(pysam_in.references)
 
    strand_specific = options.strand_specific

    result =[]
    c = E.Counter()
    c.input = 0

    for bed in Bed.iterator( IOTools.openFile( bedfile ) ):
        c.input += 1

        #if bed.contig not in contigs: 
        #    c.skipped += 1
        #    continue

        if c.input % options.report_step == 0:
            E.info( "iteration: %i" % c.input )
        
        features = counter.countInInterval( fg_file, bed.contig, bed.start, bed.end, 
                                            window_size = options.window_size,
                                            bins = bins,
                                            only_interval = options.only_interval,
                                            centring_method = options.centring_method )
     

        if control_file:
            control = counter.countAroundPos(control_file, 
                                             bed.contig,
                                             features.peak_center,
                                             bins = features.bins )

        else:
            control = None

        if options.random_shift:
            direction = numpy.random.randint( 0, 2 )
            if direction: pos = features.peak_center + 2 * bins[0]
            else: pos = features.peak_center + 2 * bins[-1]
            shifted = counter.countAroundPos(fg_file, 
                                             bed.contig,
                                             pos,
                                             bins = features.bins )
        else:
            shifted = None

        if strand_specific and bed.strand == "-":
            features._replace( hist=hist[::-1] )
            if control: control._replace( hist=hist[::-1] )
            if shifted: shift._replace( hist=hist[::-1] )

        result.append( (features, bed, control, shifted) )
        c.added += 1

    E.info( "interval processing: %s" % c )

    return result, bins

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    parser = buildOptionParser( argv )
    
    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if len(args) != 2:
        raise ValueError("please specify a bam and bed file" )

    infile, bedfile = args
    control_file = None
    if options.control_file:
        E.info("using control file %s" % options.control_file)

    if options.format == "bigwig":
        fg_file = bx.bbi.bigwig_file.BigWigFile( open(infile) )
        if options.control_file:
            control_file = bx.bbi.bigwig_file.BigWigFile( open(options.control_file) )
        counter = _bam2peakshape.CounterBigwig()

    elif options.format == "bam":
        fg_file = pysam.Samfile( infile, "rb" )
        if options.control_file:
            control_file = pysam.Samfile( options.control_file, "rb" )
        counter = _bam2peakshape.CounterBam( shift = options.shift)

    result, bins = buildResults( bedfile,
                                 fg_file, 
                                 control_file, 
                                 counter,
                                 options )

    if len(result) == 0:
        E.warn( "no data - no output" )
        E.Stop()
        return

    outputResults( result, bins, options )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


