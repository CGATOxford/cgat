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
   
   describe purpose of the script.

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
                       
    parser.set_defaults(
        output_duplicates = False,
        filename_rna = None,
        remove_rna = False
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if options.filename_rna:
        rna = GFF.readAndIndex( GFF.iterator( IOTools.openFile( options.filename_rna ) ) )
    else:
        rna = None

    pysam_in = pysam.Samfile( "-", "rb" )

    if options.output_duplicates:
        outs_dupl = open( outfile + ".duplicates", "w" )
        outs_dupl.write( "contig\tpos\tcounts\n" )
        outs_hist = open( outfile + ".histogram", "w" )
        outs_hist.write( "duplicates\tcounts\tcumul\tfreq\tcumul_freq\n" )
    else:
        outs_dupl = False
        outs_hist = False

    last_contig, last_pos = None, None
    ninput, nduplicates = 0, 0

    duplicates = collections.defaultdict( int )
    counts = collections.defaultdict( int )
    count = 0

    # count nh, nm tags
    nh, nm = collections.defaultdict( int ), collections.defaultdict( int )
    nh_all, nm_all = collections.defaultdict( int ), collections.defaultdict( int )
    flags_counts = collections.defaultdict( int )
    flags = sorted(FLAGS.keys())
    nrna, nfiltered = 0, 0

    contig, remove_rna = None, options.remove_rna

    for read in pysam_in:

        ninput += 1

        f = read.flag
        for x in flags:
            if f & x: flags_counts[x] += 1

        try:
            nh_all[read.opt("NH")] += 1
            nm_all[read.opt("NM")] += 1
        except KeyError:
            # ignore missing tags
            pass

        # skip unmapped reads
        if read.is_unmapped: continue

        if read.rname != last_contig:
            contig = pysam_in.getrname( read.rname )

        # note: does not take into account gaps within reads
        if rna and rna.contains( contig, read.pos, read.pos + read.qlen ):
            nrna += 1
            if remove_rna: continue
        
        nfiltered += 1

        try:
            nh[read.opt("NH")] += 1
            nm[read.opt("NM")] += 1
        except KeyError:
            pass

        if read.rname == last_contig and read.pos == last_pos:
            count += 1
            nduplicates += 1
            continue

        if count > 1:
            counts[count] += 1
            if outs_dupl:
                outs_dupl.write("%s\t%i\t%i\n" % (last_contig, last_pos, count) )

        count = 1
        last_contig, last_pos = read.rname, read.pos

    outs = options.stdout
    outs.write( "category\tcounts\tpercent\tof\n" )
    outs.write( "total\t%i\t%5.2f\ttotal\n" % (ninput, 100.0 ) )
    nmapped = ninput - flags_counts[4]
    outs.write( "mapped\t%i\t%5.2f\ttotal\n" % (nmapped, 100.0 * nmapped / ninput ) )

    for x in flags:
        outs.write( "%s\t%i\t%5.2f\tmapped\n" % ( FLAGS[x], flags_counts[x], 100.0 * flags_counts[x] / ninput ) )

    outs.write( "rna\t%i\t%5.2f\tmapped\n" % (nrna, 100.0 * nrna / nmapped ) )
    outs.write( "no_rna\t%i\t%5.2f\tmapped\n" % (nfiltered, 100.0 * nfiltered / nmapped ) )
    outs.write( "duplicates\t%i\t%5.2f\tno_rna\n" % (nduplicates, 100.0* nduplicates / nfiltered))
    outs.write( "unique\t%i\t%5.2f\tno_rna\n" % (nfiltered - nduplicates,
                                                 100.0*(nfiltered - nduplicates)/nfiltered))

    nreads = nmapped
    outs.write( "reads_total\t%i\t%5.2f\treads_total\n" % (nreads, 100.0 ) )
    nreads_mapped = nreads - flags_counts[4]
    outs.write( "reads_mapped\t%i\t%5.2f\treads_total\n" % (nreads_mapped, 100.0 * nreads_mapped / nreads ) )

    if len(nh_all) > 1:
        for x in xrange( 2, max(nh_all.keys() ) + 1 ): nreads -= (nh_all[x] / x) * (x-1)
        outs.write( "reads_unique\t%i\t%5.2f\treads_mapped\n" % (nh_all[1], 100.0 * nh_all[1] / nreads_mapped ) )

    nreads_norna = nfiltered
    outs.write( "reads_norna\t%i\t%5.2f\treads_mapped\n" % (nreads_norna, 100.0 * nreads_norna / nreads_mapped ) )

    if len(nh) > 1:
        for x in xrange( 2, max(nh.keys() ) + 1 ): nreads_norna -= (nh[x] / x) * (x-1)
        outs.write( "reads_norna_unique\t%i\t%5.2f\treads_norna\n" % (nh[1], 100.0 * nh[1] / nreads_norna ) )

    pysam_in.close()

    if outs_dupl:
        keys = counts.keys()
        # count per position (not the same as nduplicates, which is # of reads)
        c = 0
        total = sum( counts.values() )
        for k in sorted(keys):
            c += counts[k]
            outs_hist.write("%i\t%i\t%i\t%f\t%f\n" % (k, counts[k], c, 
                                                      100.0 * counts[k] / total,
                                                      100.0 * c / total) )
        outs_dupl.close()
        outs_hist.close()

    if len(nm) > 0:
        outfile = E.openOutputFile( "nm", "w" )
        outfile.write( "NM\talignments\n" )
        for x in xrange( 0, max( nm.keys() ) + 1 ): outfile.write("%i\t%i\n" % (x, nm[x]))
        outfile.close()

    if len(nh) > 1:
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

    
