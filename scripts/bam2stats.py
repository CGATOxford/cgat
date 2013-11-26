'''
bam2stats.py - compute stats from a bam-file
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Summary BAM

Purpose
-------

This script takes a bam file as input and computes a few metrics by
iterating over the file. The metrics output are:

+------------------------+------------------------------------------+
|*Category*              |*Content*                                 |
+------------------------+------------------------------------------+
|total                   |total number of alignments in bam file    |
+------------------------+------------------------------------------+
|alignments_mapped       |alignments mapped to a chromosome (bam    |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|alignments_unmapped     |alignments unmapped (bam flag)            |
+------------------------+------------------------------------------+
|qc_fail                 |alignments failing QC (bam flag)          |
+------------------------+------------------------------------------+
|mate_unmapped           |alignments in which the mate is unmapped  |
|                        |(bam flag)                                |
+------------------------+------------------------------------------+
|reverse                 |alignments in which read maps to reverse  |
|                        |strand (bam flag)                         |
+------------------------+------------------------------------------+
|mate_reverse            |alignments in which mate maps to reverse  |
|                        |strand (bam flag)                         |
+------------------------+------------------------------------------+
|proper_pair             |alignments in which both pairs have been  |
|                        |mapped properly (according to the mapper) |
|                        |(bam flag)                                |
+------------------------+------------------------------------------+
|read1                   |alignments for 1st read of pair (bam flag)|
+------------------------+------------------------------------------+
|paired                  |alignments of reads that are paired (bam  |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|duplicate               |read is PCR or optical duplicate (bam     |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|read2                   |alignment is for 2nd read of pair (bam    |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|secondary               |alignment is not primary alignment        |
+------------------------+------------------------------------------+
|alignments_rna          |alignments mapping to regions specified in|
|                        |a .gff file                               |
+------------------------+------------------------------------------+
|alignments_no_rna       |alignments not mapping to regions in a    |
|                        |.gff file (if --remove-rna has been set,  |
|                        |otherwise equal to mapped)                |
+------------------------+------------------------------------------+
|alignments_duplicates   |number of alignments mapping to the same  |
|                        |location                                  |
+------------------------+------------------------------------------+
|alignments_unique       |number of alignments mapping to unique    |
|                        |locations                                 |
+------------------------+------------------------------------------+
|reads_total             |number of reads in file. Either given via |
|                        |--input-reads or deduc ed as the sum of   |
|                        |mappend and unmapped reads                |
+------------------------+------------------------------------------+
|reads_mapped            |number of reads mapping in file. Derived  |
|                        |from the total number o f alignments and  |
|                        |removing counts for multiple              |
|                        |matches. Requires the NH flag to be set   |
|                        |correctly.                                |
+------------------------+------------------------------------------+
|reads_unmapped          |number of reads unmapped in file. Assumes |
|                        |that there is only one                    |
|                        | entry per unmapped read.                 |
+------------------------+------------------------------------------+
|reads_missing           |number of reads missing, if number of     |
|                        |reads given by --input-rea ds. Otherwise  |
|                        |0.                                        |
+------------------------+------------------------------------------+
|reads_norna             |reads not mapping to repetetive RNA       |
|                        |regions.                                  |
+------------------------+------------------------------------------+
|pairs_total             |number of total pairs - this is the number|
|                        |of reads_total divided by two. If there   |
|                        |were no pairs, pairs_total will be 0.     |
+------------------------+------------------------------------------+
|pairs_mapped            |number of mapped pairs - this is the same |
|                        |as the number of proper pairs.            |
+------------------------+------------------------------------------+

Additionally, the script outputs histograms for the following tags and
scores.

* NM: number of mismatches in alignments.
* NH: number of hits of reads.
* mapq: mapping quality of alignments.

If a fastq file is supplied (``--filename-fastq``), the script will
compute some additional summary statistics. However, as it builds a dictionary
of all sequences, it will also require a good  amount of memory. The additional
metrics output are:

+-------------------------+----------------------------------------+
|*Category*               |*Content*                               |
+-------------------------+----------------------------------------+
|pairs_total              |total number of pairs in input data     |
+-------------------------+----------------------------------------+
|pairs_mapped             |pairs in which both reads map           |
+-------------------------+----------------------------------------+
|pairs_unmapped           |pairs in which neither read maps        |
+-------------------------+----------------------------------------+
|pairs_proper_unique      |pairs which are proper and map uniquely.|
+-------------------------+----------------------------------------+
|pairs_incomplete         |pairs in which one of the reads maps    |
|                         |uniquely, but the other does not map.   |
+-------------------------+----------------------------------------+
|pairs_proper_duplicate   |pairs which are proper and unique, but  |
|                         |marked as duplicates.                   |
+-------------------------+----------------------------------------+
|pairs_proper_multimapping|pairs which are proper, but map to      |
|                         |multiple locations.                     |
+-------------------------+----------------------------------------+
|pairs_not_proper_unique  |pairs mapping uniquely, but not flagged |
|                         |as proper                               |
+-------------------------+----------------------------------------+
|pairs_other              |pairs not in any of the above categories|
+-------------------------+----------------------------------------+

Usage
-----

Example::

   cat in.bam | python bam2stats.py -

This command will generate various statistics based on the supplied BAM file, such as
percentage reads mapped and percentage reads mapped in pairs.

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

Reads are not counted via read name, but making use of NH and HI flags
when present.  To recap, NH is the number of reported alignments that
contain the query in the current record, while HI is the hit index and
ranges from 0 to NH-1.

Unfortunately, not all aligners follow this convention. For example,
gsnap seems to set NH to the number of reportable alignments, while
the actual number of reported alignments in the file is less. Thus, if
the HI flag is present, the maximum HI is used to correct the NH
flag. The assumption is, that the same reporting threshold has been
used for all alignments.

If no NH flag is present, it is assumed that all reads have only been reported once.

Multi-matching counts after filtering are really guesswork. Basically,
the assumption is that filtering is consistent and will tend to remove
all alignments of a query.

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
import CGAT.GTF as GTF
import pysam

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2stats
except ImportError:
    import CGAT._bam2stats as _bam2stats

FLAGS = {
    1: 'paired',
    2: 'proper_pair',
    4: 'unmapped',
    8: 'mate_unmapped',
    16: 'reverse',
    32: 'mate_reverse',
    64: 'read1',
    128: 'read2',
    256: 'secondary',
    512: 'qc_fail',
    1024: 'duplicate'  }

def computeMappedReadsFromAlignments( total_alignments, nh, max_hi ):
    '''compute number of reads alignment from total number of alignments'''
    nreads_mapped = total_alignments
    if len(nh) > 0:
        max_nh = max(nh.keys())
        if max_hi > 0:
            for x in xrange( 2, min( max_nh + 1, max_hi )): 
                nreads_mapped -= (nh[x] / x) * (x-1)
            for x in xrange( max_hi, max_nh + 1): 
                nreads_mapped -= (nh[x] / max_hi) * (max_hi-1)
        else:
            for x in xrange( 2, max(nh.keys() ) + 1 ): nreads_mapped -= (nh[x] / x) * (x-1)

    return nreads_mapped

def writeNH( outfile, nh, max_hi ):
    '''output nh array, correcting for max_hi if less than nh'''

    # need to remove double counting
    # one read matching to 2 positions is only 2

    max_nh = max(nh.keys())
    if max_hi > 0:
        for x in xrange( 1, min( max_nh + 1, max_hi) ): 
            if nh[x] == 0: continue
            outfile.write("%i\t%i\n" % (x, nh[x] / x))
        for x in xrange( max_hi, max_nh + 1): 
            if nh[x] == 0: continue
            outfile.write("%i\t%i\n" % (x, nh[x] / max_hi))
    else:
        for x in xrange( 1, max_nh + 1):
            if nh[x] == 0: continue
            outfile.write("%i\t%i\n" % (x, nh[x] / x))


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-r", "--filename-rna", dest="filename_rna", type="string", metavar='GFF',
                       help = "gff formatted file with rna locations. Note that the computation currently "
                              "does not take into account indels, so it is an approximate count only [%default]" )
    parser.add_option( "-f", "--remove-rna", dest="remove_rna", action="store_true",
                       help = "as well as counting, also remove rna reads for duplicate and other counts [%default]" )
    parser.add_option( "-i", "--input-reads", dest="input_reads", type="int",
                       help = "the number of reads - if given, used to provide percentages [%default]" )
    parser.add_option( "--force-output", dest="force_output", action="store_true",
                       help = "output nh/nm stats even if there is only a single count [%default]" )
    parser.add_option( "-d", "--output-details", dest="output_details", action="store_true",
                       help = "output per-read details [%default]" )
    parser.add_option( "-q", "--filename-fastq", dest = "filename_fastq",
                       help = "filename with sequences and quality scores. This file is only "
                       "used to collect sequence identifiers. Thus, for paired end data a "
                       "single file is sufficient [%default]" )

    parser.set_defaults(
        filename_rna = None,
        remove_rna = False,
        input_reads = 0,
        force_output = False,
        filename_fastq = None,
        output_details = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if options.filename_rna:
        rna = GTF.readAndIndex( GTF.iterator( IOTools.openFile( options.filename_rna ) ) )
    else:
        rna = None

    pysam_in = pysam.Samfile( "-", "rb" )

    if options.output_details:
        outfile_details = E.openOutputFile( "details", "w")
    else:
        outfile_details = None

    if options.filename_fastq and not os.path.exists( options.filename_fastq ):
        raise IOError("file %s does not exist" % options.filename_fastq)

    counter, flags_counts, nh_filtered, nh_all, nm_filtered, nm_all, mapq, mapq_all, max_hi = \
        _bam2stats.count( pysam_in, 
                          options.remove_rna, 
                          rna,
                          filename_fastq = options.filename_fastq,
                          outfile_details = outfile_details )

    if max_hi > 0 and max_hi != max( nh_all.keys() ):
        E.warn( "max_hi(%i) is inconsistent with max_nh (%i) - counts will be corrected" \
                    % (max_hi, max(nh_all.keys())))
    flags = sorted(flags_counts.keys())

    outs = options.stdout
    outs.write( "category\tcounts\tpercent\tof\n" )
    outs.write( "alignments_total\t%i\t%5.2f\ttotal\n" % (counter.input, 100.0 ) )
    if counter.input == 0: 
        E.warn( "no input - skipped" )
        E.Stop()
        return

    nalignments_unmapped = flags_counts["unmapped"]
    nalignments_mapped = counter.input - nalignments_unmapped
    outs.write( "alignments_mapped\t%i\t%5.2f\ttotal\n" % \
                    (nalignments_mapped, 100.0 * nalignments_mapped / counter.input ) )
    outs.write( "alignments_unmapped\t%i\t%5.2f\ttotal\n" % \
                    ( nalignments_unmapped, 100.0 * nalignments_unmapped / counter.input ) )

    if nalignments_mapped == 0: 
        E.warn( "no alignments - skipped" )
        E.Stop()
        return

    for flag, counts in flags_counts.iteritems():
        if flag == "unmapped": continue
        outs.write( "%s\t%i\t%5.2f\talignments_mapped\n" % ( flag, counts, 100.0 * counts / nalignments_mapped ) )

    if options.filename_rna:
        outs.write( "alignments_rna\t%i\t%5.2f\talignments_mapped\n" % (counter.rna, 100.0 * counter.rna / nalignments_mapped ) )
        outs.write( "alignments_no_rna\t%i\t%5.2f\talignments_mapped\n" % (counter.no_rna, 100.0 * counter.no_rna / nalignments_mapped ) )

    outs.write( "alignments_filtered\t%i\t%5.2f\talignments_mapped\n" % (counter.filtered, 100.0 * counter.filtered / nalignments_mapped ) )
    if counter.filtered == nalignments_mapped:
        normby = "alignments_mapped"
    else:
        normby = "alignments_filtered"

    if counter.filtered > 0:
        outs.write( "alignments_duplicates\t%i\t%5.2f\t%s\n" % (counter.duplicates, 100.0* counter.duplicates / counter.filtered, 
                                                                normby))
        outs.write( "alignments_unique\t%i\t%5.2f\t%s\n" % (counter.filtered - counter.duplicates,
                                                            100.0*(counter.filtered - counter.duplicates)/counter.filtered,
                                                            normby) )

    # derive the number of mapped reads in file from alignment counts
    nreads_unmapped = flags_counts["unmapped"]
    nreads_mapped = computeMappedReadsFromAlignments( nalignments_mapped, nh_all, max_hi )
    
    nreads_missing = 0
    if options.input_reads:
        nreads_total = options.input_reads
        # unmapped reads in bam file?
        if nreads_unmapped: 
            nreads_missing = nreads_total - nreads_unmapped - nreads_mapped
        else: 
            nreads_unmapped = nreads_total - nreads_mapped

    elif nreads_unmapped:
        # if unmapped reads are in bam file, take those
        nreads_total = nreads_mapped + nreads_unmapped
    else:
        # otherwise normalize by mapped reads
        nreads_unmapped = 0
        nreads_total = nreads_mapped

    outs.write( "reads_total\t%i\t%5.2f\treads_total\n" % (nreads_total, 100.0 ) )
    outs.write( "reads_mapped\t%i\t%5.2f\treads_total\n" % (nreads_mapped, 100.0 * nreads_mapped / nreads_total ) )
    outs.write( "reads_unmapped\t%i\t%5.2f\treads_total\n" % (nreads_unmapped, 100.0 * nreads_unmapped / nreads_total ) )
    outs.write( "reads_missing\t%i\t%5.2f\treads_total\n" % (nreads_missing, 100.0 * nreads_missing / nreads_total ) )

    if len(nh_all) > 1:
        outs.write( "reads_unique\t%i\t%5.2f\treads_mapped\n" % (nh_all[1], 100.0 * nh_all[1] / nreads_mapped ) )

    # compute after filtering
    # not that these are rough guesses
    if options.filename_rna:
        nreads_norna = computeMappedReadsFromAlignments( counter.filtered, nh_filtered, max_hi )

        outs.write( "reads_norna\t%i\t%5.2f\treads_mapped\n" % (nreads_norna, 100.0 * nreads_norna / nreads_mapped ) )

        if len(nh_filtered) > 1:
            outs.write( "reads_norna_unique\t%i\t%5.2f\treads_norna\n" % (nh_filtered[1], 100.0 * nh_filtered[1] / nreads_norna ) )

    pysam_in.close()

    # output paired end data 
    if flags_counts["read2"] > 0:
        if options.filename_fastq:
            pairs_mapped = counter.total_pairs - counter.total_pair_is_unmapped
            outs.write( "pairs_total\t%i\t%5.2f\tpairs_total\n" % \
                            (counter.total_pairs, 100.0 * counter.total_pairs / counter.total_pairs ) )
            outs.write( "pairs_mapped\t%i\t%5.2f\tpairs_total\n" % \
                            (pairs_mapped, 100.0 * pairs_mapped / counter.total_pairs))
            outs.write( "pairs_unmapped\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_unmapped, 100.0 * counter.total_pair_is_unmapped / counter.total_pairs ) )
            outs.write( "pairs_proper_unique\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_proper_uniq, 100.0 * counter.total_pair_is_proper_uniq / counter.total_pairs ) )
            outs.write( "pairs_incomplete\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_incomplete, 100.0 * counter.total_pair_is_incomplete / counter.total_pairs ) )
            outs.write( "pairs_proper_duplicate\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_proper_duplicate, 100.0 * counter.total_pair_is_proper_duplicate / counter.total_pairs ) )
            outs.write( "pairs_proper_multimapping\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_proper_mmap, 100.0 * counter.total_pair_is_proper_mmap / counter.total_pairs ) )
            outs.write( "pairs_not_proper_unique\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_not_proper_uniq, 100.0 * counter.total_pair_not_proper_uniq / counter.total_pairs ) )
            outs.write( "pairs_other\t%i\t%5.2f\tpairs_total\n" % \
                            ( counter.total_pair_is_other, 100.0 * counter.total_pair_is_other / counter.total_pairs ) )
        else:
            # approximate counts
            pairs_total = nreads_total // 2
            pairs_mapped = flags_counts["proper_pair"] // 2
            outs.write( "pairs_total\t%i\t%5.2f\tpairs_total\n" % \
                            (pairs_total, 100.0))
            outs.write( "pairs_mapped\t%i\t%5.2f\tpairs_total\n" % \
                            (pairs_mapped, 100.0 * pairs_mapped / pairs_total))
    else:
        pairs_total = pairs_mapped = 0
        outs.write( "pairs_total\t%i\t%5.2f\tpairs_total\n" % (pairs_total,0.0))
        outs.write( "pairs_mapped\t%i\t%5.2f\tpairs_total\n" % (pairs_mapped, 0.0))

    if options.force_output or len(nm_filtered) > 0:
        outfile = E.openOutputFile( "nm", "w" )
        outfile.write( "NM\talignments\n" )
        if len(nm_filtered) > 0:
            for x in xrange( 0, max( nm_filtered.keys() ) + 1 ): outfile.write("%i\t%i\n" % (x, nm_filtered[x]))
        else:
            outfile.write( "0\t%i\n" % (counter.filtered) )
        outfile.close()

    if options.force_output or len(nh_all) > 1:
        outfile = E.openOutputFile( "nh_all", "w")
        outfile.write( "NH\treads\n" )
        if len(nh_all) > 0:
            writeNH( outfile, nh_all, max_hi )
        else:
            # assume all are unique if NH flag not set
            outfile.write( "1\t%i\n" % (counter.mapped_reads) )
        outfile.close()

    if options.force_output or len(nh_filtered) > 1:
        outfile = E.openOutputFile( "nh", "w")
        outfile.write( "NH\treads\n" )
        if len(nh_filtered) > 0:
            writeNH( outfile, nh_filtered, max_hi )
        else:
            # assume all are unique if NH flag not set
            outfile.write( "1\t%i\n" % (counter.filtered) )
        outfile.close()

    if options.force_output or len(mapq_all) > 1:
        outfile = E.openOutputFile( "mapq", "w")
        outfile.write( "mapq\tall_reads\tfiltered_reads\n" )
        for x in xrange( 0, max( mapq_all.keys() ) + 1 ):       
            outfile.write("%i\t%i\t%i\n" % (x, mapq_all[x], mapq[x]))
        outfile.close()
        

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

    
