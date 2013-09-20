################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
bam2bam.py - modify bam files
=============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a :term:`bam` formatted file from stdin and outputs a
modified :term:`bam` formatted file on stdout.

.. note::
   You need to redirect logging information to a file or turn it off
   via -v 0 in order to get a valid sam/bam file.

Usage
-----

Example::

   python bam2bam.py < in.bam > out.bam

Type::

   python bam2bam.py --help

for command line help.

Documentation
-------------

--set-nh: set the NH flag. Some tools (bowtie, bwa) do not set the NH flag. 
   If set, this option will set the NH flag (for mapped reads).
   This option requires the bam/sam file to be sorted by read name.

--unset-unmapped_mapq: some tools set the mapping quality of unmapped reads. This
   causes a violation in the Picard tools.

--remove-better: remove alignments that are worse than alignments in the :term:`bam`
   formatted file ``--filter-bam``. 

Code
----

'''

import os
import sys
import re
import optparse
import collections
import itertools

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2bam
except ImportError:
    import CGAT._bam2bam as _bam2bam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "--set-nh", dest="set_nh", action="store_true",
                       help = "sets the NH flag. The file needs to be sorted by readname [%default]" )

    parser.add_option( "--unset-unmapped-mapq", dest="unset_unmapped_mapq", action="store_true",
                       help = "sets the mapping quality of unmapped reads to 0 [%default]" )

    parser.add_option( "--strip", dest="strip", type="choice",
                       choices = ("sequence", "quality" ),
                       help = "remove parts of the bam-file. Note that stripping the sequence will "
                       "also strip the quality values [%default]" )

    parser.add_option( "--filter", dest="filter", action="append", type="choice",
                       choices=('NM', 'CM', 'mapped', 'unique' ),
                       help = "filter bam file. The option denotes the property that is  "
                       "used to determine better match [%default]" )

    parser.add_option( "--reference-bam", dest="reference_bam", type="string",
                       help = "bam-file to filter with [%default]" )
    
    parser.add_option( "--sam", dest="output_sam", action="store_true",
                       help = "output in sam format [%default]" )

    parser.set_defaults(
        filter = [],
        set_nh = False,
        unset_unmapped_mapq = False,
        output_sam = False,
        reference_bam = None,
        strip = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # reading bam from stdin does not work with only "r" tag
    pysam_in = pysam.Samfile( "-", "rb" )

    if options.output_sam:
        pysam_out = pysam.Samfile( "-", "wh", template = pysam_in )
    else:
        pysam_out = pysam.Samfile( "-", "wb", template = pysam_in )

    if options.filter:

        if "NM" in options.filter:
            remove_mismatches = True
            colour_mismatches = False

        elif "CM" in options.filter:
            remove_mismatches = True
            colour_mismatches = True

        if remove_mismatches:
            if not options.reference_bam:
                raise ValueError( "requiring reference bam file for removing by mismatches")

            pysam_ref = pysam.Samfile( options.reference_bam, "rb" )
        else:
            pysam_ref = None

        c = _bam2bam.filter_bam( pysam_in, pysam_out, pysam_ref,
                                 remove_nonunique = "unique" in options.filter,
                                 remove_contigs = None,
                                 remove_unmapped = "mapped" in options.filter,
                                 remove_mismatches = remove_mismatches,
                                 colour_mismatches = colour_mismatches )
        
        options.stdlog.write( "category\tcounts\n%s\n" % c.asTable() )

    else:
        # set up the modifying iterators
        it = pysam_in.fetch( until_eof = True )

        if options.unset_unmapped_mapq:
            def unset_unmapped_mapq( i ):
                for read in i:
                    if read.is_unmapped:
                        read.mapq = 0
                    yield read
            it = unset_unmapped_mapq( it )

        if options.set_nh and False:
            def set_nh( i ):

                for key, reads in itertools.groupby( i, lambda x: x.qname ):
                    l = list(reads)
                    nh = len(l)
                    for read in l:
                        if not read.is_unmapped:
                            t = dict( read.tags )
                            t['NH'] = nh
                            read.tags = list(t.iteritems())
                        yield read
            it = set_nh( it )

        if options.strip != None:
            def strip_sequence( i ):
                for read in i:
                    read.seq = None
                    yield read

            def strip_quality( read ):
                for read in i:
                    read.qual = None
                    yield read

            if options.strip == "sequence":
                it = strip_sequence( it )
            elif options.strip == "quality":
                it = strip_quality( it )

        if options.set_nh:
            it = _bam2bam.SetNH( it )

        # read and output
        for read in it: 
            pysam_out.write( read )

        pysam_in.close()
        pysam_out.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

