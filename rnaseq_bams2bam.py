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

Usage
-----

Example::

   python script_template.py bamT.bam bamG.bam

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

The script needs to look-up reads via their names. It thus builds
an index of reads mapping 

This script requires the NM attributes to be set.


Code
----

'''

import os, sys, re, optparse, time

import Experiment as E
import GTF, IOTools
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
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-g", "--filename-gtf", dest="filename_gtf", type="string",
                       help = "filename with gene models [%default]" )

    parser.add_option( "-m", "--filename-mismapped", dest="filename_mismapped", type="string",
                       help = "output bam file for mismapped reads [%default]" )

    parser.add_option( "-c", "--remove-contigs", dest="remove_contigs", type="string",
                       help = "','-separated list of contigs to remove [%default]" )

    parser.add_option( "-f", "--force", dest="force", action = "store_true",
                       help = "force overwriting of existing files [%default]" )

    parser.add_option( "-u", "--unique", dest="unique", action = "store_true",
                       help = "remove reads not matching uniquely [%default]" )

    parser.set_defaults(
        filename_gtf = None,
        filename_mismapped = None,
        remove_contigs = None,
        force = False,
        unique = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 3:
        raise ValueError( "please supply three bam files" )
    
    bamfile_transcripts, bamfile_genome, bamfile_output = args

    if options.remove_contigs:
        options.remove_contigs = options.remove_contigs.split(",")

    if options.filename_gtf == None:
        raise ValueError( "please supply a file with the gene set." )

    E.info( "indexing geneset" )

    transcripts = {}
    for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile(options.filename_gtf) )):
        transcripts[gtf[0].transcript_id] = gtf

    E.info( "read %i transcripts from geneset" % len(transcripts) )

    transcripts_samfile = pysam.Samfile( bamfile_transcripts, "rb" )
    genome_samfile = pysam.Samfile( bamfile_genome, "rb" )

    if not options.force and os.path.exists( bamfile_output ):
        raise IOError( "output file %s already exists" % bamfile_output )
    output_samfile = pysam.Samfile( bamfile_output, "wb", template = genome_samfile )

    if options.filename_mismapped:
        if not options.force and os.path.exists( options.filename_mismapped ):
            raise IOError( "output file %s already exists" % options.filename_mismapped )
        output_mismapped = pysam.Samfile( options.filename_mismapped, 
                                          "wb", 
                                          template = genome_samfile )
    else:
        output_mismapped = None

    c = _rnaseq_bams2bam.filter( genome_samfile, transcripts_samfile,
                                 output_samfile, output_mismapped,
                                 transcripts,
                                 unique = options.unique,
                                 remove_contigs = options.remove_contigs )
    
    options.stdout.write( "category\tcounts\n%s\n" % c.asTable() )

    transcripts_samfile.close()
    genome_samfile.close()
    output_samfile.close()
    if output_mismapped:
        output_mismapped.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

