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
bam2stats.py - compute stats from a bam-file
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   paired-read libraries need to be treated correctly. Currently
   the second read is ignored to avoid over-counting.

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
import GFF

import pyximport
pyximport.install(build_in_temp=False)
import _bam2stats

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

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-r", "--filename-rna", dest="filename_rna", type="string",
                       help = "gff formatted file with rna locations. Note that the computation currently does not take"
                              " into account indels, so it is an approximate count only [%default]" )
    parser.add_option( "-f", "--remove-rna", dest="remove_rna", action="store_true",
                       help = "remove rna reads for duplicate and other counts [%default]" )
    parser.add_option( "-i", "--input-reads", dest="input_reads", type="int",
                       help = "the number of reads - if given, used to provide percentages [%default]" )
    parser.add_option( "--force-output", dest="force_output", type="int",
                       help = "output nh/nm stats even if there is only a single count [%default]" )

#    parser.add_option( "-p", "--ignore-pairs", dest="ignore_pairs", action="store_true",
#                       help = "if set, pairs will be counted individually. The default is to count a pair as one [%default]" )
                       
    parser.set_defaults(
        filename_rna = None,
        remove_rna = False,
        ignore_pairs = False,
        input_reads = 0,
        force_output = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if options.filename_rna:
        rna = GFF.readAndIndex( GFF.iterator( IOTools.openFile( options.filename_rna ) ) )
    else:
        rna = None

    pysam_in = pysam.Samfile( "-", "rb" )

    c, flags_counts, nh, nh_all, nm, nm_all = _bam2stats.count( pysam_in, options.remove_rna, rna )

    flags = sorted(flags_counts.keys())

    outs = options.stdout
    outs.write( "category\tcounts\tpercent\tof\n" )
    outs.write( "total\t%i\t%5.2f\ttotal\n" % (c.input, 100.0 ) )
    if c.input == 0: 
        E.warn( "no input - skipped" )
        E.Stop()
        return

    nmapped = c.input - flags_counts["unmapped"]
    outs.write( "mapped\t%i\t%5.2f\ttotal\n" % (nmapped, 100.0 * nmapped / c.input ) )
    if nmapped == 0: 
        E.warn( "no mapped reads - skipped" )
        E.Stop()
        return

    for flag, counts in flags_counts.iteritems():
        outs.write( "%s\t%i\t%5.2f\tmapped\n" % ( flag, counts, 100.0 * counts / c.input ) )

    if options.filename_rna:
        outs.write( "rna\t%i\t%5.2f\tmapped\n" % (c.rna, 100.0 * c.rna / nmapped ) )
        outs.write( "no_rna\t%i\t%5.2f\tmapped\n" % (c.filtered, 100.0 * c.filtered / nmapped ) )
        normby = "norna"
    else:
        normby = "mapped"
    
    if c.filtered > 0:
        outs.write( "duplicates\t%i\t%5.2f\t%s\n" % (c.duplicates, 100.0* c.duplicates / c.filtered, normby))
        outs.write( "unique\t%i\t%5.2f\t%s\n" % (c.filtered - c.duplicates,
                                                     100.0*(c.filtered - c.duplicates)/c.filtered,
                                                       normby) )

    # derive the number of mapped reads in file from alignment counts
    nreads_unmapped = flags_counts["unmapped"]
    nreads_mapped = nmapped
    if len(nh_all) > 1:
        for x in xrange( 2, max(nh_all.keys() ) + 1 ): nreads_mapped -= (nh_all[x] / x) * (x-1)

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

    if options.filename_rna:
        nreads_norna = c.filtered
        if len(nh) > 1:
            for x in xrange( 2, max(nh.keys() ) + 1 ): nreads_norna -= (nh[x] / x) * (x-1)

        outs.write( "reads_norna\t%i\t%5.2f\treads_mapped\n" % (nreads_norna, 100.0 * nreads_norna / nreads_mapped ) )

        if len(nh) > 1:
            outs.write( "reads_norna_unique\t%i\t%5.2f\treads_norna\n" % (nh[1], 100.0 * nh[1] / nreads_norna ) )

    pysam_in.close()

    if options.force_output or len(nm) > 0:
        outfile = E.openOutputFile( "nm", "w" )
        outfile.write( "NM\talignments\n" )
        for x in xrange( 0, max( nm.keys() ) + 1 ): outfile.write("%i\t%i\n" % (x, nm[x]))
        outfile.close()

    if options.force_output or len(nh) > 1:
        # need to remove double counting
        # one read matching to 2 positions is only 2
        outfile = E.openOutputFile( "nh", "w")
        outfile.write( "NH\treads\n" )
        for x in xrange( 1, max( nh.keys() ) + 1 ): 
            outfile.write("%i\t%i\n" % (x, nh[x] / x))
        outfile.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

    
