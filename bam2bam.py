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
bam2bam.py - modify bam files
=============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script modifies bam files by going through the whole
file.

.. note::
   You need to redirect logging information to a file or turn it off
   via -v 0 in order to get a valid sam/bam file.

Usage
-----

Example::

   python script_template.py --help

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

--set-nh: set the NH flag. Some tools (bowtie, bwa) do not set the NH flag. 
   If set, this option will set the NH flag (for mapped reads).
   This option requires the bam/sam file to be sorted by read name.

--unset-unmapped_mapq: some tools set the mapping quality of unmapped reads. This
   causes a violation in the Picard tools.

Code
----

'''

import os, sys, re, optparse, collections, itertools

import Experiment as E
import IOTools
import pysam
import GFF

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "--set-nh", dest="set_nh", action="store_true",
                       help = "sets the NH flag. The file needs to be sorted by readname [%default]" )

    parser.add_option( "--unset-unmapped-mapq", dest="unset_unmapped_mapq", action="store_true",
                       help = "sets the mapping quality of unmapped to 0 [%default]" )

    parser.add_option( "--sam", dest="output_sam", action="store_true",
                       help = "output in sam format [%default]" )

    parser.set_defaults(
        set_nh = False,
        unset_unmapped_mapq = False,
        output_sam = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    pysam_in = pysam.Samfile( "-" )
    if options.output_sam:
        pysam_out = pysam.Samfile( "-", "wh", template = pysam_in )
    else:
        pysam_out = pysam.Samfile( "-", "wb", template = pysam_in )


    # set up the modifying iterators
    it = pysam_in

    if options.unset_unmapped_mapq:
        def unset_unmapped_mapq( i ):
            for read in i:
                if read.is_unmapped:
                    read.mapq = 0
                yield read
        it = unset_unmapped_mapq( it )

    if options.set_nh:
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

    # read and output
    for read in it: pysam_out.write( read )

    pysam_in.close()
    pysam_out.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

