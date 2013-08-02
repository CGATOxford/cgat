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
mask_fasta.py - mask fasta formatted sequences
==============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

If the file contains multiple sequence entries, it is split
at entry boundaries. Currently, individual entries are not
split.

Usage
-----

Example::

   python mask_fasta.py --help

Type::

   python mask_fasta.py --help

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
import array

import CGAT.Experiment as E

def maskSequence( sequence, regions, mask_char = "N" ):
    """mask sequence with regions."""
    nmasked = 0
    v = array.array( 'c' )
    v.fromstring( sequence)
    errors = []
    for start,end in regions:

        if start < 0 or end > len(sequence):
            errors.append( "out of bounds: %i-%i" % (start, end) )
            continue
        
        for x in range(start,end):
            if v[x] != mask_char:
                v[x] = mask_char
                nmasked += 1
        
    return v.tostring(), nmasked, errors


if __name__ == '__main__':


    parser = E.OptionParser( version = "%prog version: $Id: mask_fasta.py 2782 2009-09-10 11:40:29Z andreas $")

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
    total_keys, total_written, total_masked, nerrors = 0, 0, 0, 0
    outfile = options.stdout
    
    for line in infile:
        if line[0] == "#": continue
        if line[0] == ">":
            if fragments:

                if options.loglevel >= 1:
                    options.stdlog.write( "# processing sequence %s\n" % key )
                    options.stdlog.flush()
                    
                if key in mask_regions:
                    sequence, nmasked, errors = maskSequence( string.join(fragments, ""), mask_regions[key], options.mask_char)
                    if errors:
                        nerrors += len(errors)
                        options.stdlog.write("# %i errors while masking\n" % len(errors))
                        for e in errors:
                            options.stdlog.write("# %s\n" % e )
                else:
                    sequence = string.join(fragments, "")
                    
                if description_rest:
                    outfile.write(">%s %s\n%s\n" % (key, description_rest, sequence ) )
                else:
                    outfile.write(">%s\n%s\n" % (key, sequence ) )                    
                total_masked += nmasked
                total_keys += 1
                total_written += len(sequence)
                
            fragments = []
            
            x = re.split( "\s", line[1:-1])
            key = x[0]
            description_rest = string.join(x[1:], " ")
            continue
        fragments.append( line[:-1] )

    if fragments:
        if key in mask_regions:
            sequence, nmasked, errors = maskSequence( string.join(fragments, ""), mask_regions[key], options.mask_char)
            if errors:
                nerrors += len(errors)                
                options.stdlog.write("# %i errors while masking\n" % len(errors))
                for e in errors:
                    options.stdlog.write("# %s\n" % e )
        else:
            sequence = string.join(fragments, "")

    total_masked += nmasked
    total_keys += 1
    total_written += len(sequence)

    if description_rest:
        outfile.write(">%s %s\n%s\n" % (key, description_rest, sequence ) )
    else:
        outfile.write(">%s\n%s\n" % (key, sequence ) )                    
    
    if options.filename_sequence: infile.close()

    if options.loglevel >= 1:
        options.stdlog.write( "# nkeys=%i, nwritten=%i, nmasked=%i, nerrors=%i\n" % (total_keys, total_written, total_masked, nerrors) )

    E.Stop()
