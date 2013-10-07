'''
bam2geneprofile.py - build meta-gene profile for a set of transcripts/genes
===========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Genesets Intervals GTF BAM Summary

Purpose
-------

This script takes a :term:`gtf` formatted file and computes meta-gene profiles
over various annotations derived from the :term:`gtf` file. 

A meta-gene profile is an abstract genomic entity over which reads stored in a 
:term:`bam` formatted file have been counted. A meta-gene might be an idealized
eukaryotic gene (upstream, exonic sequence, downstream) or any other genomic landmark 
of interest such as transcription start sites.

The script can be used to visualize binding profiles of a chromatin mark in gene 
bodies, binding of transcription factors in promotors or 3' bias in RNASeq data.

This script is designed with a slight emphasis on RNA-Seq datasets. For example, it
takes care of spliced reads, by using the CIGAR string in the BAM file to accurately
define aligned bases (when the --base-accurate-off is not specified, currently it is not specified by default).

Alternatively, for the purpose of visualizing binding profiles of transcription factors
ChIP-Seq without the need to use any genomic annotations (ENSEMBL, or refseq), you may
also consider using :doc:`bam2peakshape`, which is
designed with a slight emphasis on Chip-Seq datasets, for example, :doc:`bam2peakshape` is able to
center the counting window to the summit of every individual peak. :doc:`bam2peakshape` is also able
to: (1) plot the control RNA-Seq library to enable side-by-side comparison; (2) randomize
the given regions to provide a semi-control.

Background
---------- 

The :file:`bam2geneprofile.py` script reads in a set of transcripts from a :term:`gtf` 
formatted file. For each transcript, overlapping reads from the provided :term:`bam` file 
are collected. The counts within the transcript are then mapped onto the meta-gene
structure and counts are aggregated over all transcripts in the :term:`gtf` file.

:term:`Bam` files need to be sorted by coordinate and indexed. 

A meta-gene structure has two components - regions of variable size, such as exons, introns, etc,
which nevertheless have a fixed start and end coordinate in a transcript. The other component
are regions of fixed width, such a regions of a certain size upstream or downstream of a landmark such as
a transcription start site.

The size of the former class, regions of variable size, can be varied with ``--resolution`` options.
For example, the option ``--resolution-upstream-utr=1000`` will create a meta-gene with 
a 1000bp upstream UTR region. UTRs that are larger will be compressed, but UTRs that are smaller, will
be stretched to fit the 1000bp meta-gene UTR region.

The size of fixed-width regions can be set with ``--extension`` options. For example,
the options ``--extension-upstream`` will set the size of the uptsream extension region
to 1000bp. Note that no scaling is required when counting reads towards the meta-gene profile.

Options
-------

The behaviour of the script can be modified by several options.

Profiles
++++++++

Different profiles are accessible through the ``--counter`` option. Multiple
counters can be applied at the same time. While ``upstream`` and ``downstream``
typically have a fixed size, the other regions such as ``CDS``, ``UTR`` will be
scaled to a common size.

utrprofile
    UPSTREAM - UTR5 - CDS - UTR3 - DOWNSTREAM
    gene models with UTR. Separate the coding section from the non-coding part.

geneprofile 
    UPSTREAM - EXON - DOWNSTREAM
    simple exonic gene models

geneprofilewithintrons 
    UPSTREAM - EXON - INTRON - DOWNSTREAM
    gene models containing also intronic sequence, only correct if used with --base-accuracy
     
geneprofileabsolutedistancefromthreeprimeend
    UPSTREAM - EXON(AbsoluteFromThreePrimeEnd) - INTRON(AbsoluteFromThreePrimeEnd) - DOWNSTREAM
    For EXON region, the script counts over the mRNA transcript only, skipping introns. 
    Designed to visualize the 3 prime bias in RNASeq data, only correct if used with --base-accuracy
    
tssprofile
    UPSTREAM - DOWNSTREAM
    transcription start/stop sites

intervalprofile 
    UPSTREAM - INTERVAL - DOWNSTREAM
    Similar to geneprofile, but count over the complete span of the gene
    (including introns).

midpointprofile
    UPSTREAM  - DOWNSTREAM
    aggregate over midpoint of gene model

Genes versus transcripts
++++++++++++++++++++++++

The default is to collect reads on a per-transcript level. Alternatively, the script can
merge all transcripts of a gene into a single virtual transcript. Note that this virtual transcript
might not be a biologically plausible transcript. It is usually better to provide :file:`bam2geneprofile.py`
with a set of representative transcripts per gene in order to avoid up-weighting genes with
multiple transcripts.

Meta-gene structures
++++++++++++++++++++

The script provides a variety of different meta-gene structures (``--methods``).

Normalization
+++++++++++++

Normalization can be applied in two stages of the computation.

Before adding counts to the meta-gene profile, the profile for the individual
transcript can be normalized. Without normalization, highly expressed genes
will contribute more to the meta-gene profile than lowly expressed genes.
With normalization, each gene contributes an equal amount (``-normalize``).

When outputting the meta-gene profile, the meta-gene profile itself can be normalized.
Normalization a profile can help comparing the shapes of profiles between different
experiments independent of the number of reads or transcripts used in the construction
of the meta-gene profile (``--normalize-profile``).

Bed and wiggle files
++++++++++++++++++++

The densities can be computed from :term:`bed` or :term:`wiggle` formatted files.
If a :term:`bed` formatted file is supplied, it must be compressed with and indexed with :file:`tabix`.

.. note::
   
   Paired-endedness is ignored. Both ends of a paired-ended read are treated 
   individually.

Usage
-----

The following command will generate the gene profile from RNA-Seq data::

   python bam2geneprofile.py --bamfile=rnaseq.bam -g geneset.gtf.gz 
                     --method=geneprofilewithintrons  --reporter=gene 
                     --resolution-cds=1400 --resolution-introns=2000

The output will contain read coverage over genes, merging all transcripts within
a gene into a single chain of exons. The profile will contain four separate
segments:

1. the upstream region of a gene (default = 1000bp)

2. the transcribed region of a gene. The transcribed region of every gene will be 
   scaled to 1400 bp (``--resolution-cds=1400``), 
   shrinking longer transcripts and expanding shorter transcripts.

3. the intronic regions of a gene. These will be scaled to 2kb (``-resolution-introns=2000``).

4. the downstream region of a gene (default = 1000bp).

Type::

   python bam2geneprofile.py --help

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
import CGAT.GTF as GTF
import numpy

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2geneprofile
except ImportError:
    import CGAT._bam2geneprofile as _bam2geneprofile


from bx.bbi.bigwig_file import BigWigFile

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--method", dest="methods", type = "choice", action = "append",
                       choices = ("geneprofile", "tssprofile", "utrprofile", 
                                  "intervalprofile", "midpointprofile",
                                  "geneprofilewithintrons", 
                                  "geneprofileabsolutedistancefromthreeprimeend",
                                  ),
                       help = 'counters to use. Counters describe the meta-gene structure to use '
                       '[%default]. \n Note using geneprofilewithintrons, or geneprofileabsolutedistancefromthreeprimeend will automatically turn on the --base-accuracy option' )

    parser.add_option( "-b", "--bamfile", "--bedfile", "--bigwigfile", dest="infiles", 
                       metavar = "BAM",
                       type = "string", action = "append",
                       help = "BAM/bed/bigwig files to use. Do not mix different types"
                              "[%default]" )

    parser.add_option( "-g", "--gtffile", dest="gtffile", type = "string",
                       metavar = "GTF",
                       help = "GTF file to use. "
                              "[%default]" )

    parser.add_option( "-n", "--normalization", dest="normalization", type = "choice",
                       choices = ("none", "max", "sum", "total-max", "total-sum"),
                       help = """normalization to apply on each transcript profile before adding to meta-gene profile. 

The options are:

* none: no normalization
* sum: sum of counts within a region
* max: maximum count within a region
* total-sum: sum of counts across all regions
* total-max: maximum count in all regions
[%default]""" )

    parser.add_option( "-p", "--normalize-profile", dest="profile_normalizations", type = "choice", action="append",
                       choices = ("none", "area", "counts"),
                       help = """normalization to apply on meta-gene profile normalization. 

The options are:
* none: no normalization
* area: normalize such that the area under the meta-gene profile is 1.
* counts: normalize by number of features (genes,tss) that have been counted
[%default]""" )

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

    parser.add_option("--resolution-introns", dest="resolution_introns", type = "int",
                       help = "resolution of introns region in bp "
                              "[%default]" )
                              
    parser.add_option("--resolution-exons-absolute-distance-topolya", dest="resolution_exons_absolute_distance_topolya", type = "int",
                       help = "resolution of exons absolute distance topolya in bp "
                              "[%default]" )                    
                       
    parser.add_option("--resolution-introns-absolute-distance-topolya", dest="resolution_introns_absolute_distance_topolya", type = "int",
                       help = "resolution of introns absolute distance topolya in bp "
                              "[%default]" )
                       
    parser.add_option("--extension-exons-absolute-distance-topolya", dest="extension_exons_absolute_distance_topolya", type = "int",
                       help = "extension for exons from the absolute distance from the topolya in bp"
                              "[%default]" )

    parser.add_option("--extension-introns-absolute-distance-topolya", dest="extension_introns_absolute_distance_topolya", type = "int",
                       help = "extension for introns from the absolute distance from the topolya in bp"
                              "[%default]" )                
                              
    parser.add_option("--extension-upstream", dest = "extension_upstream", type = "int",
                       help = "extension upstream from the first exon in bp"
                              "[%default]" )

    parser.add_option("--extension-downstream", dest = "extension_downstream", type = "int",
                       help = "extension downstream from the last exon in bp"
                              "[%default]" )

    parser.add_option("--extension-inward", dest="extension_inward", type = "int",
                       help = "extension inward from a TSS start site in bp"
                              "[%default]" )

    parser.add_option("--extension-outward", dest="extension_outward", type = "int",
                       help = "extension outward from a TSS start site in bp"
                              "[%default]" )
                       
    parser.add_option("--scale-flank-length", dest="scale_flanks", type = "int",
                       help = "scale flanks to (integer multiples of) gene length"
                              "[%default]" )

    parser.add_option( "--matrix-format", dest="matrix_format", type = "choice", 
                       choices = ("multiple", "single" ),
                       help = "matrix output format, either 'multiple' files or a 'single' file "
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
        resolution_introns = 1000,
        resolution_exons_absolute_distance_topolya = 3000,    #3kb is a good balance of seeing long enough 3 prime bias and not omit too many genes. Tim 31th Aug 2013
        resolution_introns_absolute_distance_topolya = 500,   #introns is only for assess the noise level, thus do ont need a long region, a long region has the side effect of omit more genes. Tim 31th Aug 2013
        extension_exons_absolute_distance_topolya = 3000,     #3kb is a good balance of seeing long enough 3 prime bias and not omit too many genes. Tim 31th Aug 2013
        extension_introns_absolute_distance_topolya = 500,    #introns is only for assess the noise level, thus do ont need a long region, a long region has the side effect of omit more genes. Tim 31th Aug 2013
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
        matrix_format = "single",
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
    
    for methodsRequiresBaseAccuracy in ["geneprofilewithintrons", 
                                        "geneprofileabsolutedistancefromthreeprimeend",
                                       ]:
        # If you implemented any methods that you do not want the spliced out introns 
        # or exons appear to be covered by non-existent reads, it is better you let those
        # methods imply --base-accurarcy by add them here.
        if methodsRequiresBaseAccuracy in options.methods:
            options.base_accuracy=True

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
        elif method == "geneprofilewithintrons":
            counters.append( _bam2geneprofile.GeneCounterWithIntrons( range_counter, 
                                                           options.resolution_upstream,
                                                           options.resolution_cds,
                                                           options.resolution_introns,
                                                           options.resolution_downstream,
                                                           options.extension_upstream,
                                                           options.extension_downstream,
                                                           options.scale_flanks ) )

        elif method == "geneprofileabsolutedistancefromthreeprimeend":
            counters.append( _bam2geneprofile.GeneCounterAbsoluteDistanceFromThreePrimeEnd( range_counter, 
                                                           options.resolution_upstream,                                                           
                                                           options.resolution_downstream,
                                                           options.resolution_exons_absolute_distance_topolya,
                                                           options.resolution_introns_absolute_distance_topolya,
                                                           # options.resolution_exons_absolute_distance_tostartsite,   
                                                           # options.resolution_introns_absolute_distance_tostartsite,             
                                                           # Tim 31th Aug 2013: a possible feature for future,  if five prime bias is of your interest. 
                                                           #(you need to create another class). It is not very difficult to derive from this class, but is not implemented yet
                                                           # This future feature is slightly different the TSS profile already implemented, because in this future feature introns are skipped, 
                                                           options.extension_upstream,
                                                           options.extension_downstream,
                                                           options.extension_exons_absolute_distance_topolya,
                                                           options.extension_introns_absolute_distance_topolya,
                                                           # options.extension_exons_absolute_distance_tostartsite,   
                                                           # options.extension_introns_absolute_distance_tostartsite,                                           
                                                           # Tim 31th Aug 2013: a possible feature for future,  if five prime bias is of your interest. 
                                                           #(you need to create another class). It is not very difficult to derive from this class, but is not implemented yet
                                                           # This future feature is slightly different the TSS profile already implemented, because in this future feature introns are skipped, 
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

    # output matrices
    if not options.profile_normalizations:
        options.profile_normalizations.append( "none" )
    elif "all" in options.profile_normalizations:
        options.profile_normalizations = ["none", "area", "counts" ]

    for method, counter in zip(options.methods, counters):
        if options.matrix_format == "multiple":
            for norm in options.profile_normalizations:
                with IOTools.openFile( E.getOutputFile( counter.name ) + ".%s.tsv.gz" % norm, "w") as outfile:
                    counter.writeMatrix( outfile, normalize=norm )

        elif options.matrix_format == "single":
            # build a single output
            matrices = []
            for norm in options.profile_normalizations:
                matrix = counter.buildMatrix( normalize = norm )
                nrows, ncols = matrix.shape
                matrix.shape = (nrows * ncols, 1 )
                matrices.append( matrix )

            for x in range(1,len(matrices)):
                assert matrices[0].shape == matrices[x].shape
                
            # build a single matrix
            matrix = numpy.hstack( matrices )
            nrows, ncols = matrix.shape
            with IOTools.openFile( E.getOutputFile( counter.name ) + ".matrix.tsv.gz", "w" ) as outfile:
                outfile.write( "bin\tregion\tregion_bin\t%s\n" % "\t".join( options.profile_normalizations) )
                fields = []
                bins = []
                for field, nbins in zip( counter.fields, counter.nbins ):
                    fields.extend( [field] * nbins )
                    bins.extend( list(range(nbins)) )

                for row, cols in enumerate(zip( fields, bins, matrix)):
                    outfile.write( "%i\t%s\t" % (row, "\t".join( [ str(x) for x in cols[:-1]  ]) ))
                    outfile.write( "%s\n" % ("\t".join( [ str(x) for x in cols[-1]  ]) ))

        with IOTools.openFile( E.getOutputFile( counter.name ) + ".lengths.tsv.gz", "w") as outfile:
            counter.writeLengthStats( outfile )

    if options.plot:

        import matplotlib
        # avoid Tk or any X
        matplotlib.use( "Agg" )
        import matplotlib.pyplot as plt
        
        for method, counter in zip(options.methods, counters):

            if method in ("geneprofile", "geneprofilewithintrons", "geneprofileabsolutedistancefromthreeprimeend", "utrprofile", "intervalprofile" ):

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

    
