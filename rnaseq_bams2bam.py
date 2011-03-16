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

    parser.add_option( "-f", "--force", dest="force", action = "store_true",
                       help = "force overwriting of existing files [%default]" )

    parser.set_defaults(
        filename_gtf = None,
        filename_mismapped = None,
        force = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 3:
        raise ValueError( "please supply three bam files" )
    
    bamfile_transcripts, bamfile_genome, bamfile_output = args

    if options.filename_gtf == None:
        raise ValueError( "please supply a file with the gene set." )

    E.info( "indexing geneset" )

    transcripts = {}
    for gtf in GTF.transcript_iterator( GTF.iterator( IOTools.openFile(options.filename_gtf) )):
        transcripts[gtf[0].transcript_id] = gtf

    E.info( "read %i transcripts from geneset" % len(transcripts) )

    E.info( "building index of transcriptome reads" )
    
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

    index = pysam.IndexedReads( transcripts_samfile )
    index.build()
    
    E.info( "finished building index" )

    c = E.Counter()

    E.info( "starting filtering" )

    for read in genome_samfile:
        c.input += 1

        if read.is_unmapped:
            c.unmapped_genome += 1
            c.output += 1
            output_samfile.write( read )
            continue

        # get transcripts that read matches to
        try:
            matches = list( index.find( read.qname ) )
        except KeyError:
            c.notfound += 1
            c.output += 1
            output_samfile.write( read )
            continue
        
        g_contig = genome_samfile.getrname( read.tid )

        # set mapped = True, if read is mapped to transcripts
        #
        # set location_ok = True, if read matches in expected location
        # according to transcripts
        location_ok = False
        mapped = False
        read_mismatches = read.opt("NM")

        for match in matches:

            # ignore reads that are not mapped to transcripts
            if match.is_unmapped: continue
            # ignore reads that are mapped to transcripts with
            # more mismatches than the genomic location
            if match.opt("NM") > read_mismatches: continue

            mapped = True

            # find transcript
            transcript_id = transcripts_samfile.getrname( match.tid )
            gtfs = transcripts[transcript_id]
            t_contig, t_start, t_end = gtfs[0].contig, gtfs[0].start, gtfs[-1].end
            
            # does read map to genomic transcript location?
            if g_contig == t_contig and t_start <= read.pos <= t_end:
                location_ok = True
                break

        if location_ok:
            c.location_ok += 1
            c.output += 1
            output_samfile.write( read )
        elif mapped:
            c.mismapped += 1
            if output_mismapped:
                output_mismapped.write( read )
        else:
            c.unmapped_transcript += 1
            c.output += 1
            output_samfile.write( read )
            
        
    E.info( "finished filtering" )

    E.info( "%s" % str(c))

    transcripts_samfile.close()
    genome_samfile.close()
    output_samfile.close()
    if output_mismapped:
        output_mismapped.close()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

