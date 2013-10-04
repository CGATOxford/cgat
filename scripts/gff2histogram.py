'''
gff2histogram.py - compute histograms from intervals in gff or bed format
=========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals GFF Summary

Purpose
-------

This script computes distributions of interval sizes, intersegmental distances
and interval overlap from a list of intervals in :term:`gff` or :term:`bed` 
format. 

The output will be written into separate files. Filenames are given by 
``--ouput-filename-pattern``.

Available methods are:

hist
    Output a histogram of interval sizes and distances between intervals
    in nucleotides.

stats
    Output summary statistics of interval sizes and distances between
    intervals

values
    Output distances, sizes, and overlap values to separate files.

all
    all of the above.

Usage
-----

Example::

   python gff2histogram.py < in.gff

Type::

   python gff2histogram.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import re
import optparse

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.Histogram as Histogram
import CGAT.Stats as Stats

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gff2histogram.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-b", "--bin-size", dest="bin_size", type="string",
                      help="bin size."  )
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    parser.add_option("--no-empty-bins", dest="no_empty_bins", action="store_true",
                      help="do not display empty bins.")
    parser.add_option("--with-empty-bins", dest="no_empty_bins", action="store_false",
                      help="display empty bins.")
    parser.add_option("--ignore-out-of-range", dest="ignore_out_of_range", action="store_true",
                      help="ignore values that are out of range (as opposed to truncating them to range border.")
    parser.add_option("--missing", dest="missing_value", type="string",
                      help="entry for missing values [%default]." )
    parser.add_option("--dynamic-bins", dest="dynamic_bins", action="store_true",
                      help="each value constitutes its own bin." )
    parser.add_option( "--format", dest="format", type="choice", 
                       choices=( "gff", "gtf", "bed" ),
                       help="input file format [%default].")
    parser.add_option( "--method", dest="methods", type="choice", action="append",
                       choices=( "all", "hist", "stats", "overlaps", "values" ),
                       help="methods to apply [%default].")
    parser.add_option( "--data", dest="data", type="choice",
                       choices=( "all", "size", "distance" ),
                       help="data to compute [%default].")

    parser.set_defaults(
        no_empty_bins = True,
        bin_size = None,
        dynamic_bins = False,
        ignore_out_of_range = False,
        min_value = None,
        max_value = None,
        nonull = None,
        missing_value = "na",
        output_filename_pattern="%s",
        methods = [],
        data = "all",
        format = "gff",
        )

    (options, args) = E.Start( parser, add_output_options = True )

    if "all" in options.methods:
        options.methods = ("hist", "stats", "overlaps")
        if not options.output_filename_pattern: options.output_filename_pattern = "%s"

    if len(options.methods) == 0:
        raise ValueError( "please provide counting method using --method option" )

    if options.format in ( "gff", "gtf" ):
        gffs = GTF.iterator( options.stdin )
    elif options.format == "bed":
        gffs = Bed.iterator( options.stdin )

    values_between = []
    values_within = []
    values_overlaps = []

    if "overlaps" in options.methods:
        if not options.output_filename_pattern: 
            options.output_filename_pattern = "%s"
        outfile_overlaps = E.openOutputFile( "overlaps" )
    else:
        outfile_overlaps = None

    last = None
    ninput, noverlaps = 0,0
    for this in gffs:
        ninput += 1
        values_within.append( this.end - this.start )

        if last and last.contig == this.contig:
            if this.start < last.end:
                noverlaps += 1
                if outfile_overlaps:
                    outfile_overlaps.write( "%s\t%s\n" % (str(last), str(this)) )
                values_overlaps.append( min(this.end, last.end) - max(last.start, this.start) )
                if this.end > last.end:
                    last = this
                continue
            else:
                values_between.append( this.start - last.end )
                # if this.start - last.end < 10: 
                #     print str(last)
                #     print str(this)
                #     print "=="
                values_overlaps.append( 0 )

        last = this

    if "hist" in options.methods:
        outfile = E.openOutputFile( "hist" )
        h_within = Histogram.Calculate( values_within,
                                        no_empty_bins = options.no_empty_bins,
                                        increment = options.bin_size,
                                        min_value = options.min_value,
                                        max_value = options.max_value,
                                        dynamic_bins = options.dynamic_bins,
                                        ignore_out_of_range = options.ignore_out_of_range )

        h_between = Histogram.Calculate( values_between,
                                         no_empty_bins = options.no_empty_bins,
                                         increment = options.bin_size,
                                         min_value = options.min_value,
                                         max_value = options.max_value,
                                         dynamic_bins = options.dynamic_bins,
                                         ignore_out_of_range = options.ignore_out_of_range )

        if "all" == options.data:
            outfile.write( "residues\tsize\tdistance\n" )
            combined_histogram = Histogram.Combine( [h_within, h_between], missing_value = options.missing_value )
            Histogram.Write( outfile, combined_histogram, nonull = options.nonull )        
        elif options.data == "size":
            outfile.write( "residues\tsize\n" )
            Histogram.Write( outfile, h_within, nonull = options.nonull )        
        elif options.data == "distance":
            outfile.write( "residues\tdistance\n" )
            Histogram.Write( outfile, h_between, nonull = options.nonull )        

        outfile.close()

    if "stats" in options.methods:
        outfile = E.openOutputFile( "stats" )
        outfile.write( "data\t%s\n" % Stats.Summary().getHeader() )
        if options.data in ("size", "all"):
            outfile.write( "size\t%s\n" % str(Stats.Summary(values_within)) )
        if options.data in ("distance", "all"):
            outfile.write( "distance\t%s\n" % str(Stats.Summary(values_between)) )
        outfile.close()

    if "values" in options.methods:
        outfile = E.openOutputFile( "distances" )
        outfile.write( "distance\n%s\n" % "\n".join( map(str, values_between) ) )
        outfile.close()
        outfile = E.openOutputFile( "sizes" )
        outfile.write( "size\n%s\n" % "\n".join( map(str, values_within) ) )
        outfile.close()
        outfile = E.openOutputFile( "overlaps" )
        outfile.write( "overlap\n%s\n" % "\n".join( map(str, values_overlaps) ) )
        outfile.close()

    E.info( "ninput=%i, ndistance=%i, nsize=%i, noverlap=%i" % (ninput, 
                                                                len(values_between),
                                                                len(values_within),
                                                                noverlaps) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
