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
bam2geneprofile.py - build coverage profile for a set of transcripts/genes
===========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a :term:`gtf` formatted file and computes density profiles
over various annotations derived from the :term:`gtf` file. 

The densities can be computed from :term:`bam` or :term:`bed` formatted files.
:term:`bam` files need to be sorted by coordinate and indexed. If a :term:`bed` 
formatted file is supplied, it must be compressed with and indexed with :file:`tabix`.

.. todo::
   
   * paired-endedness is ignored. Both ends of a paired-ended read are treated 
   individually.

   * use CIGAR string to accurately define aligned regions. Currently only the
   full extension of the alignment is used. Thus, results from RNASEQ experiments
   might be mis-leading.
  
-n/--normalization
   normalize counts before aggregation. The normalization options are:
   * sum: sum of counts within a region
   * max: maximum count within a region
   * total-sum: sum of counts across all regions
   * total-max: maximum count in all regions
   * none: no normalization

-N/--normalize-profiles
   normalize profiles when outputting. The normalization methods are:
   * none: no normalization
   * area: normalize such that the area under the profile is 1.
   * counts: normalize by number of features (genes,tss) that have been counted

   several options can be combined.

-m/--method
   utrprofile - aggregate over gene models with UTR. The areas are
                upstream - UTR5 - CDS - UTR3 - downstream
   geneprofile - aggregate over exonic gene models. The areas are
                upstream - EXON - downstream
   tssprofile  - aggregate over transcription start/stop
                upstream - TSS/TTS - downstream
   intervalprofile - aggregate over interval
                upstream - EXON - downstream
   midpointprofile - aggregate over midpoint of gene model. The areas are:
                upstream - downstream

-a/--merge-pairs
   merge pairs in paired-ended data.
   
Usage
-----

Example::

   python bam2geneprofile.py reads.bam genes.gtf

Type::

   python bam2geneprofile.py --help

for command line help.

Documentation
-------------


Code
----

'''

import os
import sys
import re
import optparse
import collections
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.GTF as GTF
import numpy

import pyximport
pyximport.install(build_in_temp=False)
import _bam2geneprofile
from bx.bbi.bigwig_file import BigWigFile

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--method", dest="methods", type = "choice", action = "append",
                       choices = ("geneprofile", "tssprofile", "utrprofile", 
                                  "intervalprofile", "midpointprofile",
                                  "),
                       help = "counters to use. Counters describe the meta-gene structure to use"
                              " for counting. Several counters can be chosen."
                              " [%default]" )

    parser.add_option( "-b", "--bamfile", "--bedfile", "--bigwigfile", dest="infiles", 
                       type = "string", action = "append",
                       help = "BAM/bed/bigwig files to use. Do not mix different types"
                              "[%default]" )

    parser.add_option( "-g", "--gtffile", dest="gtffile", type = "string",
                       help = "GTF file to use. "
                              "[%default]" )

    parser.add_option( "-n", "--normalization", dest="normalization", type = "choice",
                       choices = ("none", "max", "sum", "total-max", "total-sum"),
                       help = "normalization to use. Normalize the counts before aggregating "
                              "[%default]" )

    parser.add_option( "-p", "--normalize-profile", dest="profile_normalizations", type = "choice", action="append",
                       choices = ("none", "area", "counts"),
                       help = "profile normalization to use. "
                              "[%default]" )

    parser.add_option( "-r", "--reporter", dest="reporter", type = "choice",
                       choices = ("gene", "transcript"  ),
                       help = "report results for genes or transcripts."
                              " When 'genes` is chosen, exons across all transcripts for"
                              " a gene are merged. When 'transcript' is chosen, counts are"
                              " computed for each transcript separately with each transcript"
                              " contributing equally to the meta-gene profile."
                              " [%default]" )

    parser.add_option( "-i", "--shift", dest="shifts", type = "int", action = "append",
                       help = "shift reads in :term:`bam` formatted file before computing densities (ChIP-Seq). "
                              "[%default]" )

    parser.add_option( "-a", "--merge-pairs", dest="merge_pairs", action = "store_true",
                       help = "merge pairs in :term:`bam` formatted file before computing"
                              " densities (ChIP-Seq)."
                              "[%default]" )

    parser.add_option( "-u", "--base-accuracy", dest="base_accuracy", action = "store_true",
                       help = "compute densities with base accuracy. The default is to"
                              " only use the start and end of the aligned region (RNA-Seq)"
                              " [%default]" )

    parser.add_option( "-e", "--extend", dest="extends", type = "int", action = "append",
                       help = "extend reads in :term:`bam` formatted file (ChIP-Seq). "
                              "[%default]" )

    parser.add_option("--resolution-upstream", dest="resolution_upstream", type = "int",
                       help = "resolution of upstream region in bp "
                              "[%default]" )

    parser.add_option("--resolution-downstream", dest="resolution_downstream", type = "int",
                       help = "resolution of downstream region in bp "
                              "[%default]" )

    parser.add_option("--resolution-upstream-utr", dest="resolution_upstream_utr", type = "int",
                       help = "resolution of upstream UTR region in bp "
                              "[%default]" )

    parser.add_option("--resolution-downstream-utr", dest="resolution_downstream_utr", type = "int",
                       help = "resolution of downstream UTR region in bp "
                              "[%default]" )

    parser.add_option("--resolution-cds", dest="resolution_cds", type = "int",
                       help = "resolution of cds region in bp "
                              "[%default]" )

    parser.add_option("--extension_upstream", dest = "extension_upstream", type = "int",
                       help = "extension upstream from the first exon in bp"
                              "[%default]" )

    parser.add_option("--extension_downstream", dest = "extension_downstream", type = "int",
                       help = "extension downstream from the last exon in bp"
                              "[%default]" )

    parser.add_option("--extension_inward", dest="extension_inward", type = "int",
                       help = "extension inward from a TSS start site in bp"
                              "[%default]" )

    parser.add_option("--extension_outward", dest="extension_outward", type = "int",
                       help = "extension outward from a TSS start site in bp"
                              "[%default]" )
                       
    parser.add_option("--scale_flank_length", dest="scale_flanks", type = "int",
                       help = "scale flanks to (integer multiples of) gene length"
                              "[%default]" )

    parser.set_defaults(
        remove_rna = False,
        ignore_pairs = False,
        input_reads = 0,
        force_output = False,
        bin_size = 10,
        extends = [],
        shifts = [],
        sort = [],
        reporter = "transcript",
        resolution_cds = 1000,
        resolution_upstream_utr = 1000,
        resolution_downstream_utr = 1000,
        resolution_upstream = 1000,
        resolution_downstream = 1000,
        # mean length of transcripts: about 2.5 kb
        extension_upstream = 2500,
        extension_downstream = 2500,
        extension_inward = 3000,
        extension_outward = 3000,
        plot = True,
        methods = [],
        infiles = [],
        gtffile = None,
        profile_normalizations = [],
        normalization = None,
        scale_flanks = 0,
        merge_pairs = False,
        min_insert_size = 0,
        max_insert_size = 1000,
        base_accuracy = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    # Keep for backwards compatability
    if len(args) == 2:
        infile, gtf = args
        options.infiles.append(infile)
        options.gtffile = gtf

    if not options.gtffile:
        raise ValueError("no GTF file specified" )

    if len(options.infiles) == 0:
        raise ValueError("no bam/wig/bed files specified" )

    if options.reporter == "gene":
        gtf_iterator = GTF.flat_gene_iterator( GTF.iterator( IOTools.openFile( options.gtffile ) ) )
    elif options.reporter == "transcript":
        gtf_iterator = GTF.transcript_iterator( GTF.iterator( IOTools.openFile( options.gtffile ) ) )
        
    # Select rangecounter based on file type
    if len(options.infiles) > 0:
        if options.infiles[0].endswith( ".bam" ):
            bamfiles = [ pysam.Samfile( x, "rb" ) for x in options.infiles ]
            format = "bam"
            if options.merge_pairs:
                range_counter = _bam2geneprofile.RangeCounterBAM( bamfiles, 
                                                                  shifts = options.shifts, 
                                                                  extends = options.extends,
                                                                  merge_pairs = options.merge_pairs,
                                                                  min_insert_size = options.min_insert_size,
                                                                  max_insert_size = options.max_insert_size )
            elif options.shifts or options.extends:
                range_counter = _bam2geneprofile.RangeCounterBAM( bamfiles, 
                                                                  shifts = options.shifts, 
                                                                  extends = options.extends )
            elif options.base_accuracy:
                range_counter = _bam2geneprofile.RangeCounterBAMBaseAccuracy( bamfiles )
            else:
                range_counter = _bam2geneprofile.RangeCounterBAM( bamfiles )
            
                                                              
        elif options.infiles[0].endswith( ".bed.gz" ):
            bedfiles = [ pysam.Tabixfile( x ) for x in options.infiles ]
            format = "bed"
            range_counter = _bam2geneprofile.RangeCounterBed( bedfiles )

        elif options.infiles[0].endswith( ".bw" ):
            wigfiles = [ BigWigFile(file=open(x)) for x in options.infiles ]
            format = "bigwig"
            range_counter = _bam2geneprofile.RangeCounterBigWig( wigfiles )

        else:
            raise NotImplementedError( "can't determine file type for %s" % bamfile )

    counters = []
    for method in options.methods:
        if method == "utrprofile":
            counters.append( _bam2geneprofile.UTRCounter( range_counter, 
                                                          options.resolution_upstream,
                                                          options.resolution_upstream_utr,
                                                          options.resolution_cds,
                                                          options.resolution_downstream_utr,
                                                          options.resolution_downstream,
                                                          options.extension_upstream,
                                                          options.extension_downstream ) )
        elif method == "geneprofile":
            counters.append( _bam2geneprofile.GeneCounter( range_counter, 
                                                           options.resolution_upstream,
                                                           options.resolution_cds,
                                                           options.resolution_downstream,
                                                           options.extension_upstream,
                                                           options.extension_downstream,
                                                           options.scale_flanks ) )

        elif method == "tssprofile":
            counters.append( _bam2geneprofile.TSSCounter( range_counter, 
                                                           options.extension_outward,
                                                           options.extension_inward ) )

        elif method == "intervalprofile":
            counters.append( _bam2geneprofile.RegionCounter( range_counter, 
                                                             options.resolution_upstream,
                                                             options.resolution_cds,
                                                             options.resolution_downstream,
                                                             options.extension_upstream,
                                                             options.extension_downstream ) )

        elif method == "midpointprofile":
            counters.append( _bam2geneprofile.MidpointCounter( range_counter, 
                                                               options.resolution_upstream,
                                                               options.resolution_downstream,
                                                               options.extension_upstream,
                                                               options.extension_downstream ) )

    # set normalization
    for c in counters:
        c.setNormalization( options.normalization )

    E.info( "starting counting with %i counters" % len(counters) )

    _bam2geneprofile.count( counters, gtf_iterator )

    for method, counter in zip(options.methods, counters):
        if not options.profile_normalizations:
            options.profile_normalizations.append( "none" )

        for norm in options.profile_normalizations:
            outfile = IOTools.openFile( E.getOutputFile( counter.name ) + ".%s.tsv.gz" % norm, "w")
            counter.writeMatrix( outfile, normalize=norm )
            outfile.close()
            
        outfile = IOTools.openFile( E.getOutputFile( counter.name ) + ".lengths.tsv.gz", "w")
        counter.writeLengthStats( outfile )
        outfile.close()

    if options.plot:

        import matplotlib
        # avoid Tk or any X
        matplotlib.use( "Agg" )
        import matplotlib.pyplot as plt
        
        for method, counter in zip(options.methods, counters):

            if method in ("geneprofile", "utrprofile", "intervalprofile" ):

                plt.figure()
                plt.subplots_adjust( wspace = 0.05)
                max_scale = max( [max(x) for x in counter.aggregate_counts ] )

                for x, counts in enumerate( counter.aggregate_counts ):
                    plt.subplot( 5, 1, x+1)
                    plt.plot( range(len(counts)), counts )
                    plt.title( counter.fields[x] )
                    plt.ylim( 0, max_scale )

                figname = counter.name + ".full"
                
                fn = E.getOutputFile( figname ) + ".png"
                plt.savefig( os.path.expanduser(fn) )

                plt.figure()

                points = []
                cuts = []
                for x, counts in enumerate( counter.aggregate_counts ):
                    points.extend( counts )
                    cuts.append( len( counts ) )
                                 
                plt.plot( range(len(points)), points )
                xx,xxx = 0, []
                for x in cuts:
                    xxx.append( xx + x // 2 )
                    xx += x
                    plt.axvline( xx, color = "r", ls = "--" )

                plt.xticks( xxx, counter.fields )

                figname = counter.name + ".detail"
                
                fn = E.getOutputFile( figname ) + ".png"
                plt.savefig( os.path.expanduser(fn) )

            elif method == "tssprofile":

                plt.figure()
                plt.subplot( 1, 3, 1)
                plt.plot( range(-options.extension_outward, options.extension_inward), counter.aggregate_counts[0] )
                plt.title( counter.fields[0] )
                plt.subplot( 1, 3, 2)
                plt.plot( range(-options.extension_inward, options.extension_outward), counter.aggregate_counts[1] )
                plt.title( counter.fields[1] )
                plt.subplot( 1, 3, 3)
                plt.title( "combined" )
                plt.plot( range(-options.extension_outward, options.extension_inward), counter.aggregate_counts[0] )
                plt.plot( range(-options.extension_inward, options.extension_outward), counter.aggregate_counts[1] )
                plt.legend( counter.fields[:2] )

                fn = E.getOutputFile( counter.name ) + ".png"
                plt.savefig( os.path.expanduser(fn) )

            elif method == "midpointprofile":

                plt.figure()
                plt.plot( numpy.arange(-options.resolution_upstream, 0), counter.aggregate_counts[0] )
                plt.plot( numpy.arange(0, options.resolution_downstream), counter.aggregate_counts[1] )

                fn = E.getOutputFile( counter.name ) + ".png"
                plt.savefig( os.path.expanduser(fn) )
        
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

    
