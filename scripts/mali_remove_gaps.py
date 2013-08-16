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
mali_remove_gaps.py - 
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

   python mali_remove_gaps.py --help

Type::

   python mali_remove_gaps.py --help

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
import math

USAGE="""python %s [OPTIONS] < exonerate_output > filtered

Prune a nucelotide multiple alignment according to a master sequence.

1. Go in codon steps through the multiple alignment according
to the master sequence.

2. Remove all columns in other sequences, that
        1. fall out of frame
        2. are incomplete codons

Version = $Id: mali_remove_gaps.py 2782 2009-09-10 11:40:29Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-o, --file-output               output
""" % sys.argv[0]

param_long_options=["verbose=", "help", "file-output=", "version" ]

param_short_options="v:hm:e:p:c"

param_loglevel = 1

param_gap_char = "-"
param_mask_char = "x"

param_filename_output = None

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.MaliIO as MaliIO
import CGAT.Exons as Exons

##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-o", "--file-output"):
            param_filename_output = a
            
    ## 1. read multiple alignment in fasta format
    mali, identifiers = MaliIO.readFasta( sys.stdin )

    if param_loglevel >= 1:
        print "# read mali with %i entries." % len(identifiers)

    print E.GetHeader()
    print E.GetParams()

    ## 1. remove gaps in multiple alignment

    mali = MaliIO.removeGaps( mali )
    
    if param_master:
        frame_columns = GetFrameColumns( mali, param_master )
    elif param_master_pattern:
        columns = []
        for id in identifiers:
            if re.search( param_master_pattern, id ):
                columns += GetFrameColumns( mali, id )

        if len(columns) == 0:
            columns += GetFrameColumns( mali, identifiers[0] )
                
        # sort all columns by tuple. The "shortest" codon will be first (1,2,3) before (1,2,100)
        columns.sort()

        # select codons
        frame_columns = []
        last_codon = columns[0]
        for codon in columns[1:]:
            # skip identical codons
            if codon == last_codon: continue

            # take first (shortest) codon in case of identical first residue
            if codon[0] == last_codon[0]: continue
                
            # if not overlapping, keep
            if codon[0] > last_codon[2]:
                frame_columns.append( last_codon )

            # if overlapping, but out of register: skip
            last_codon = codon
            
        frame_columns.append( last_codon )        

    ## translate characters to upper/lower case according to exon info.
    if exons:
        for id in mali:
            if id in exons:
                mali[id] = AddExonInformation( mali[id], exons[id], mask_char = param_mask_char )

    if param_loglevel >= 1:
        print "# found %i columns" % (len(frame_columns))

    mask_chars = ( string.upper( param_mask_char ), string.lower( param_mask_char ) )

    
    for id in mali.keys():
        sequence = mali[id]
        fragments = []
        nstops, ncodons, naligned = 0, 0, 0
        for a,b,c in frame_columns:
            codon = sequence[a] + sequence[b] + sequence[c]

            codon_is_aligned = False
            codon_is_ok = True
            for x in codon:
                ## a codon will be masked, if it either
                ## 1. contains a gap character
                ## 2. is an unaligned character, i.e.,
                ##     exons and masked, or no exons and lowerwase
                residue_is_unaligned = (x == param_gap_char) or \
                                       (not exons and x in string.lowercase) or \
                                       (exons and x in mask_chars)
                codon_is_aligned = codon_is_aligned or not residue_is_unaligned
                codon_is_ok = codon_is_ok and not residue_is_unaligned

            if codon_is_aligned: naligned += 1
            
            if codon_is_ok:
                ncodons += 1
                if string.upper(codon) in ("TAG", "TAA", "TGA"):
                    if param_remove_stops:
                        fragments.append( param_gap_char * 3 )
                    else:
                        fragments.append( codon )                        
                    nstops += 1
                else:
                    fragments.append( codon )                                            
            else:
                fragments.append( param_gap_char * 3 )
                
        mali[id] = string.join(fragments, "")
        if param_loglevel >= 1:
            print "# sequence: %s\tpositions: %i\taligned:%i\tcodons: %i\t stops: %i" % (id, len(fragments), naligned, ncodons, nstops)
            sys.stdout.flush()
            
    for id in mali.keys():
        if param_mark_codons:
            a = mali[id]
            f = lambda x: a[x:x+3]
            s = string.join( [ f(x) for x in range(0,len(a),3)], " ")
        else:
            s = mali[id]
        print ">%s\n%s" % (id, s )

    if param_filename_translation:
        outfile = open( param_filename_translation, "w")
        for id in mali.keys():
            outfile.write( ">%s\n%s\n" % (id, Genomics.TranslateDNA2Protein(mali[id])))
        outfile.close()

    print E.GetFooter()
    
    
