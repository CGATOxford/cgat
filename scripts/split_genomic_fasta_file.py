################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
split_genomic_fasta_file.py - 
======================================================

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

   python split_genomic_fasta_file.py --help

Type::

   python split_genomic_fasta_file.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse

import CGAT.Experiment as E

USAGE="""python %s [OPTIONS] [genomic_sequence] [ < genomic sequence]

Version: $Id: split_genomic_fasta_file.py 1225 2007-04-10 15:13:11Z andreas $

If the file contains multiple sequence entries, it is split
at entry boundaries. Currently, individual entries are not
split.

""" % sys.argv[0]

global_written2file = 0
global_outfile = None
global_file_id = 0

def maskSequence( sequence, regions, mask_char = "N" ):
    """mask sequence with regions."""
    nmasked = 0
    for start,end in regions:
        sequence[start:end] = mask_char * (end-start)
        nmasked += end-start
        
    return sequence, nmasked

def processSequence( key, description_rest, sequence, options, mask_regions = None ):

    global global_outfile
    global global_file_id
    global global_written2file

    written = 0

    lsequence = len(sequence)

    if mask_regions:
        sequence, nmasked = maskSequence( sequence, mask_regions[key] )
    else:
        nmasked = 0

    while written < lsequence:

        ## size of current chunk
        chunk_size = min( options.chunk_size - global_written2file, lsequence - written)

        first_res = max(0, written )
        last_res  = min(written + chunk_size + options.extend, lsequence )

        if options.loglevel >= 2:
            options.stdlog.write( "# file_id=%i, chunk_size=%i, global_written2file=%i, written=%i, lsequence=%i, first_res=%i, last_res=%i\n" % \
                                  (global_file_id, chunk_size, global_written2file, written, lsequence, first_res, last_res ) )
            

        ## check if we need to open a new sequence file:
        if global_outfile == None:
            global_file_id += 1
            global_outfile = open( options.output_pattern % global_file_id, "w")
                
        offset_positive_strand = first_res
        offset_negative_strand = lsequence - last_res

        if options.format == "fasta":        
            global_outfile.write(">%s_%i_%i %s\n%s\n" % (key,
                                                         offset_positive_strand, offset_negative_strand,                                                     
                                                         description_rest, sequence[first_res:last_res]))
        elif options.format == "ranges":
            global_outfile.write("%s\t%i\t%i\n" % (key, first_res, last_res ) )

        global_written2file += chunk_size
        
        ## if we have reached the chunksize: close file
        if global_written2file >= options.chunk_size:
            global_outfile.close()
            global_outfile = None
            global_written2file = 0

        written += chunk_size

    return len(sequence), nmasked

if __name__ == '__main__':


    parser = E.OptionParser( version = "%prog version: $Id: split_genomic_fasta_file.py 1225 2007-04-10 15:13:11Z andreas $")

    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend sequences by # residues."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for the genomic files. Should contain a %i."  )

    parser.add_option("-s", "--sequence", dest="filename_sequence", type="string",
                      help="sequence filename."  )

    parser.add_option("-c", "--chunk-size", dest="chunk_size", type="int",
                      help="chunk size in nucleotides."  )

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("fasta", "ranges"),
                      help="output format options."  )

    parser.add_option("-m", "--mask-regions", dest="filename_mask_regions", type="string",
                      help="mask regions - parameter is a filename of regions with forward strand coordinates in the format 'contig <tab> start <tab> end'."  )

    parser.set_defaults(
        extend = 0,
        output_pattern = "genome_%5i.fasta",
        format = "fasta",
        filename_sequence = None,
        filename_mask_regions = None,
        mask_char = "N" )

    (options, args) = E.Start( parser, add_pipe_options = True )

    ## read segments to mask
    
    if options.filename_mask_regions:
        mask_regions = {}
        infile = open(options.filename_mask_regions, "r")
        nregions = 0
        for line in infile:
            if line[0] == "#": continue
            contig, start, end = line[:-1].split("\t")
            if contig not in mask_regions: mask_regions[contig] = []
            mask_regions[contig].append( (int(start),int(end) ) )
            nregions += 1

        if options.loglevel >= 1:
            options.stdlog.write("# masking %i regions on %i contigs\n" % (nregions, len(mask_regions)))
            
        for contig, regions in mask_regions.items():
            regions.sort()
            
        infile.close()
    else:
        mask_regions = None
        
    if options.filename_sequence:
        infile = open(options.filename_sequence, "r")
    else:
        infile = sys.stdin

    ## For large genomic sequences, reading the whole data is
    ## a stretch for memory resources. Thus proceed sequence by
    ## sequence

    # read complete genomic sequences
    fragments = []
    sequences = []
    key = None
    description_rest = None
    total_keys, total_written, total_masked = 0, 0, 0
    
    for line in infile:
        if line[0] == "#": continue
        if line[0] == ">":
            if fragments:
                nwritten, nmasked = processSequence( key, description_rest, string.join(fragments, ""), options, mask_regions)
                total_written += nwritten
                total_masked += nmasked
                total_keys += 1
                
            fragments = []
            
            x = re.split( "\s", line[1:-1])
            key = x[0]
            description_rest = string.join(x[1:], " ")
            continue
        fragments.append( line[:-1] )
        
    nwritten, nmasked = processSequence( key, description_rest, string.join(fragments, ""), options, mask_regions)        

    total_written += nwritten
    total_masked += nmasked
    total_keys += 1
    
    if options.filename_sequence: infile.close()

    if options.loglevel >= 1:
        options.stdlog.write( "# nkeys=%i, nwritten=%i, nmasked=%i\n" % (total_keys, total_written, total_masked))

    E.Stop()
