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
rnaseq_bams2bam.py - merge genomic and transcriptome mapped bamfiles
====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes as input two BAM files from an RNASeq experiment. 
The first bam file (:file:`bamG`) should contain reads mapped against the genome
using a mapper permitting splicing (e.g. tophat). The second bam file
(:file:`bamT`) should contain reads mapped against known transcripts. This script
will write a new bam file that removes reads from :term:`bamG` that 
map to regions that are conflicting with those in :term:`bamT`.

.. note:: 
   Note that if junctions are supplied, the resultant bam files will not
   be sorted by position.

.. glossary::

   bamG
      :term:`bam` formatted file with reads mapped against the genome

   bamT
      :term:`bam` formatted file with reads mapped against transcripts

Usage
-----

Example::

   python cgat_script_template.py bamT.bam bamG.bam

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

The script needs to look-up reads via their names. It thus builds
an index of reads mapping 

This script requires the NM attributes to be set. If it is not 
set, you will need to set a policy.

Code
----

'''

import os
import sys
import re
import optparse
import time

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import pysam
import CGAT.IndexedGenome as IndexedGenome

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _rnaseq_bams2bam
except ImportError:
    import CGAT._rnaseq_bams2bam as _rnaseq_bams2bam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-g", "--filename-gtf", dest="filename_gtf", type="string",
                       help = "filename with gene models in gtf format [%default]" )

    parser.add_option( "-m", "--filename-mismapped", dest="filename_mismapped", type="string",
                       help = "output bam file for mismapped reads [%default]" )

    parser.add_option( "-j", "--filename-junctions", dest="filename_junctions", type="string",
                       help = "bam file with reads mapped across junctions [%default]" )

    parser.add_option( "-r", "--filename-regions", dest="filename_regions", type="string",
                       help = "filename with regions to remove in bed format [%default]" )

    parser.add_option( "-t", "--filename-transcriptome", dest="filename_transcriptome", type="string",
                       help = "bam file with reads mapped against transcripts [%default]" )

    parser.add_option( "-p", "--filename-map", dest="filename_map", type="string",
                       help = "filename mapping transcript numbers (used by --filename-transciptome) to transcript names (used by --filename-gtf) [%default]" )

    parser.add_option( "-s", "--filename-stats", dest="filename_stats", type="string",
                       help = "filename to output stats to [%default]" )

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

    parser.add_option( "--sam", dest="output_sam", action="store_true",
                       help = "output in sam format [%default]" )

    parser.set_defaults(
        filename_gtf = None,
        filename_mismapped = None,
        filename_junctions = None,
        filename_transcriptome = None,
        filename_map = None,
        remove_contigs = None,
        force = False,
        unique = False,
        colour_mismatches = False,
        ignore_mismatches = False,
        output_sam = False,
        filename_table = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 1:
        raise ValueError( "please supply one bam file" )

    bamfile_genome = args[0]
    genome_samfile = pysam.Samfile( bamfile_genome, "rb" )

    if options.remove_contigs:
        options.remove_contigs = options.remove_contigs.split(",")

    if options.filename_map:
        E.info( "reading map" )
        id_map = IOTools.readMap( IOTools.openFile( options.filename_map), has_header = True )
	id_map = dict( [ (y,x) for x,y in id_map.iteritems() ])
    else:
        id_map = None
        
    transcripts = {}
    if options.filename_gtf:
        E.info( "indexing geneset" )
        mapped, missed = 0, 0
        for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile(options.filename_gtf) )):
            gtf.sort( key = lambda x: x.start )
            transcript_id = gtf[0].transcript_id
            if id_map:
                try:
                    transcript_id = id_map[transcript_id]
                    mapped += 1
                except KeyError:
                    missed += 1
                    continue
            transcripts[transcript_id] = gtf

        E.info( "read %i transcripts from geneset (%i mapped, %i missed)" % (len(transcripts), mapped, missed) )

    regions_to_remove = None
    if options.filename_regions:
        E.info( "indexing regions" )
        regions_to_remove = IndexedGenome.Simple()
        for bed in Bed.iterator( IOTools.openFile( options.filename_regions )):
            regions_to_remove.add( bed.contig, bed.start, bed.end )
        E.info( "read %i regions" % len(regions_to_remove) )

    if options.filename_transcriptome:
        transcripts_samfile = pysam.Samfile( options.filename_transcriptome, 
                                             "rb" )
    else:
        transcripts_samfile = None

    if not options.force and os.path.exists( bamfile_output ):
        raise IOError( "output file %s already exists" % bamfile_output )

    if options.output_sam:
        output_samfile = pysam.Samfile( "-", "wh", template = genome_samfile )
    else:
        output_samfile = pysam.Samfile( "-", "wb", template = genome_samfile )

    if options.filename_mismapped:
        if not options.force and os.path.exists( options.filename_mismapped ):
            raise IOError( "output file %s already exists" % options.filename_mismapped )
        output_mismapped = pysam.Samfile( options.filename_mismapped, 
                                          "wb", 
                                          template = genome_samfile )
    else:
        output_mismapped = None

    if options.filename_junctions:
        junctions_samfile = pysam.Samfile( options.filename_junctions, 
                                           "rb" )
    else:
        junctions_samfile = None

    c = _rnaseq_bams2bam.filter( genome_samfile, 
                                 output_samfile, 
                                 output_mismapped,
                                 transcripts_samfile,
                                 junctions_samfile,
                                 transcripts,
                                 regions = regions_to_remove,
                                 unique = options.unique,
                                 remove_contigs = options.remove_contigs,
                                 colour_mismatches = options.colour_mismatches,
                                 ignore_mismatches = options.ignore_mismatches,
                                 ignore_transcripts = transcripts_samfile == None,
                                 ignore_junctions = junctions_samfile == None)

    if options.filename_stats:
        outf = IOTools.openFile( options.filename_stats, "w" )
        outf.write( "category\tcounts\n%s\n" % c.asTable() )
        outf.close()

    if options.filename_transcriptome:
        transcripts_samfile.close()
    
    genome_samfile.close()
    output_samfile.close()
    if output_mismapped:
        output_mismapped.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

