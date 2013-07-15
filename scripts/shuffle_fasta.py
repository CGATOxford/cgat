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
shuffle_fasta.py - shuffle sequences
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

shuffle a set of fasta sequences.

Usage
-----

Example::

   python shuffle_fasta.py --help

Type::

   python shuffle_fasta.py --help

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
import getopt
import time
import optparse
import math
import tempfile
import random
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import Bio
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: shuffle_fasta.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--codons", dest="codons", action="store_true",
                      help="make sure that shuffled sequences only contain valid codons."  )

    parser.add_option("-a", "--conserve-aminos", dest="conserve_aminos", action="store_true",
                      help="conserve amino acids."  )

    parser.add_option("-b", "--bias", dest="bias", type="float",
                      help="introduce bias into codon usage choice. Complete bias is 1.0, while no bias is 0.0."  )

    parser.add_option("-i", "--biased-codon-usage", dest="filename_biased_codon_usage", type="string",
                      help="Filename with reference codon usage table for biased codon usage."  )

    parser.add_option("-u", "--bulk-codon-usage", dest="filename_bulk_codon_usage", type="string",
                      help="Filename with reference codon usage table for unbiased codon usage."  )

    parser.set_defaults(
        codons = False,
        conserve_aminos = False,
        bias = 0.0,
        filename_biased_codon_usage= None,
        filename_bulk_codon_usage= None,        
        stop_codons = ("TAG", "TAA", "TGA"),
        precision = 10000,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    iterator = FastaIterator.FastaIterator( sys.stdin )

    ## get map of amino acids to codons
    map_aa2codons = Genomics.GetMapAA2Codons()

    ## for codon based shuffling: build ranges based on strength of bias and on reference codon usage
    ## Bias switches from completely biased to unbiased. Unbiased is uniform usage.
    if options.filename_biased_codon_usage:
        
        map_codon2frequency = IOTools.ReadMap( open(options.filename_biased_codon_usage, "r"), map_functions = (str,float), has_header = True )

        if options.filename_bulk_codon_usage:
            map_codon2frequency_bulk = IOTools.ReadMap( open(options.filename_bulk_codon_usage, "r"), map_functions = (str,float), has_header = True )
            
        codon_ranges = {}
        for aa in map_aa2codons.keys():
            c = []
            x = 0
            for codon in map_aa2codons[aa]:

                if options.filename_bulk_codon_usage:
                    u = map_codon2frequency_bulk[codon]
                else:
                    # uniform usage                    
                    u = 1.0 / len(map_aa2codons[aa])
                
                g = map_codon2frequency[codon]
                f = g + (u - g) * (1.0 - options.bias)
                x += f * options.precision
                c.append( x )
            codon_ranges[aa] = c
            
    while 1:
        cur_record = iterator.next()
        
        if cur_record is None: break
        
        sequence = re.sub( " ", "", cur_record.sequence)
        l = len(sequence)

        if options.conserve_aminos:
            n = []
            for codon in [ sequence[x:x+3] for x in range(0, l, 3) ]:
                aa = Genomics.MapCodon2AA( codon )
                if aa not in map_aa2codons: continue
                if options.bias or options.filename_biased_codon_usage:
                    ## get random number from 0 to precision 
                    v = random.randint( 0, options.precision )
                    ## find the corresponding intervall:
                    l = len(map_aa2codons[aa])
                    x = 0
                    while x < l - 1:
                        if v < codon_ranges[aa][x]:
                            break
                        x += 1
                else:
                    x = random.randint( 0, len(map_aa2codons[aa])-1)
                n.append( map_aa2codons[aa][x] )
            sequence = "".join( n )
        else:
            sequence = list(sequence)            
            if options.codons:
                while 1:
                    random.shuffle(sequence)
                    for codon in [ sequence[x:x+3] for x in range(0, l, 3) ]:
                        if codon in options.stop_codons:
                            redo = True
                            break
                    else:
                        break
            else:
                random.shuffle( sequence )
            sequence = "".join(sequence)
        options.stdout.write( ">%s\n%s\n" % (cur_record.title, "".join(sequence)))
        
    E.Stop()
