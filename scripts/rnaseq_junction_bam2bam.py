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
rnaseq_bams2bam.py - convert mappings against junctions to genomic coordinates
================================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes as input a BAM file resulting from reads mapped against
a junction database and outputs a :term:`bam` formatted file in genomic
coordinates.

The contigs should be of the format 
<chromosome>|<start>|<exon-end>-<exon-start>|<end>|<splice>|<strand>.

<start> - 0-based coordinate of first base
<exon-end> - 0-based coordinate of last base in exon
<exon-start> - 0-based coordinate of first base in exon
<end> - 0-based coordinate of base after last base

Strand can be either ``fwd`` or ``rev``, though sequences in the database
and coordinates are all on the forward strand.

For example ``chr1|1244933|1244982-1245060|1245110|GTAG|fwd`` tranlates to the
intron ``chr1:1244983-1245060`` in python coordinates.

The input bam-file is supposed to be sorted by read. Only the best matches
are output for each read, were best is defined both in terms of number 
of mismatches and number of colour mismatches.

Usage
-----

Example::

   cat input.bam | python rnaseq_junction_bam2bam.py - --log=log > output.bam

Type::

   python rnaseq_junction_bam2bam.py --help

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
import time
import itertools

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import pysam

import pyximport
pyximport.install(build_in_temp=False)
import _rnaseq_bams2bam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-t", "--template-bam", dest="filename_genome_bam", type="string",
                       help = "input bam file for header information [%default]" )

    parser.add_option( "-s", "--contig-sizes", dest="filename_contigs", type="string",
                       help = "filename with contig sizes [%default]" )

    parser.add_option( "-o", "--colour", dest="colour_mismatches", action="store_true",
                       help = "mismatches will use colour differences (CM tag) [%default]" )

    parser.add_option( "-i", "--ignore-mismatches", dest="ignore_mismatches", action="store_true",
                       help = "ignore mismatches [%default]" )

    parser.add_option( "-c", "--remove-contigs", dest="remove_contigs", type="string",
                       help = "','-separated list of contigs to remove [%default]" )

    parser.add_option( "-f", "--force", dest="force", action = "store_true",
                       help = "force overwriting of existing files [%default]" )

    parser.add_option( "-u", "--unique", dest="unique", action = "store_true",
                       help = "remove reads not matching uniquely [%default]" )

    parser.set_defaults(
        filename_genome_bam = None,
        filename_gtf = None,
        filename_mismapped = None,
        remove_contigs = None,
        force = False,
        unique = False,
        colour_mismatches = False,
        ignore_mismatches = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    genomefile, referencenames, referencelengths = None, None, None

    if options.filename_genome_bam:
        genomefile = pysam.Samfile( options.filename_genome_bam, "rb" )
    elif options.filename_contigs:
        contigs = IOTools.ReadMap( IOTools.openFile( options.filename_contigs ) )
        data = zip( *list( contigs.iteritems() ) )
        referencenames, referencelengths = data[0], map(int, data[1])
    else:
        raise ValueError("please provide either --template-bam or --contig-sizes" )

    infile = pysam.Samfile( "-", "rb" )
    outfile = pysam.Samfile( "-", "wb", template = genomefile,
                             referencenames = referencenames,
                             referencelengths = referencelengths )
    
    if options.colour_mismatches: 
        tag = "CM"
    else:
        tag = "NM"

    nambiguous = 0
    ninput = 0
    nunmapped = 0
    ncigar = 0
    nfull = 0
    noutput = 0

    contig2tid = dict( [ (y,x) for x,y in enumerate( outfile.references ) ] )

    for qname, readgroup in itertools.groupby( infile, lambda x: x.qname ):
        ninput += 1
        reads = list(readgroup)
        if reads[0].is_unmapped:
            nunmapped += 1
            continue

        # filter for best match
        best = min( [ x.opt(tag) for x in reads ] )
        reads = [ x for x in reads if x.opt(tag) == best ]
        if len(reads) > 1:
            nambiguous += 1
            continue
        
        read = reads[0]

        # reject complicated matches (indels, etc)
        # to simplify calculations below.
        if len(read.cigar) > 1:
            ncigar += 1
            continue 
        
        # set NH flag to latest count
        t = dict( read.tags )
        t['NH'] = 1
        read.tags = list(t.iteritems())

        sname = infile.getrname( read.tid )

        contig, first_exon_start, middle, last_exon_end, splice, strand = sname.split("|")
        first_exon_end, last_exon_start = middle.split("-")
        first_exon_start, first_exon_end, last_exon_start, last_exon_end = map(int, (\
                first_exon_start, first_exon_end, last_exon_start, last_exon_end ) )
        first_exon_end += 1

        total = first_exon_end - first_exon_start + last_exon_end - last_exon_start
        first_exon_length = first_exon_end - first_exon_start

        match1 = first_exon_length - read.pos
        intron_length = last_exon_start - first_exon_end
        match2 = read.qlen - match1

        # match lies fully in one exon - ignore
        if match1 <= 0 or match2 <= 0:
            nfull += 1
            continue

        # increment pos
        read.pos = first_exon_start + read.pos
        read.tid = contig2tid[contig]
        # 3 = BAM_CREF_SKIP
        read.cigar = [( 0, match1), (3, intron_length), (0, match2) ]

        outfile.write( read )
      
        noutput += 1

    outfile.close()
    if genomefile:
        genomefile.close()

    c = E.Counter()
    c.input = ninput
    c.output = noutput
    c.full = nfull
    c.cigar = ncigar
    c.ambiguous = nambiguous
    c.unmapped = nunmapped

    E.info( "%s" % str(c))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
