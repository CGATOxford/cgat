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
mali2predictions.py - 
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

   python mali2predictions.py --help

Type::

   python mali2predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
USAGE="""convert a multiply aligned genomic segment to predictions.

lower case characters are exons, upper-case characters are introns.
Introns of a size smaller than 6 are frameshifts.
"""

import sys
import os
import string
import re
import optparse

import CGAT.Experiment as E
import CGAT.Prediction as Prediction
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.Mali as Mali
import CGAT.IOTools as IOTools

def processCodon( codon ):
    """process a codon.

    A codon is a list of residues belonging to a codon.
    """
    alignment = []

    max_d = 0
    max_x = 0
    for x in range(1, len(codon)):
        d = codon[x][0] - codon[x-1][0]
        if d > max_d:
            max_d = d
            max_x = x

    if max_d > 6:
        # process introns
        frame = max_x % 3
        if frame:
            alignment.append( ["S", 0, frame] )
        alignment.append( ["5", 0, 2] )
        alignment.append( ["I", 0, max_d - 5] )
        alignment.append( ["3", 0, 2] )        
        if frame:
            alignment.append( ["S", 1, 3-frame] )

    elif max_d == 1 and len(codon) == 3:
        alignment.append( ["M", 1, 3] )

    elif len(codon) < 3:
        alignment.append( ["F", 1, len(codon) ] )

    elif len(codon) == 3:
        d = codon[-1][0] - codon[0][0] + 1
        alignment.append( ["F", 1, d] )
    else:
        raise "untreated"        
    
    # print codon, alignment

    return alignment

##------------------------------------------------------------------------
##------------------------------------------------------------------------
##------------------------------------------------------------------------
def getAlignedColumns( mali, exons ):

    # columns that are aligned, maps to the aligned residue in master
    aligned_columns = {}
    # residues in master sequence that are aligned and exonic, maps to the aligned exonic residue
    aligned_exons = {}
    
    if options.master:
        
        try:
            s = mali[options.master]
        except KeyError:
            raise KeyError("master sequence %s not found in mali." % s )

        n = 0
        r = 0
        for x in range(len(s)):
            c = s[x]
            if c in options.gap_chars:
                continue

            r += 1
            aligned_columns[ x ] = r
            
            if c in string.uppercase:
                continue
            
            aligned_exons[ x ] = n
            n += 1
    else:
        raise "please specify a master to define frames."
    
    return aligned_columns, aligned_exons

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: mali2predictions.py 2781 2009-09-10 11:33:14Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome."  )

    parser.add_option( "-l", "--filename-locations", dest="filename_locations", type="string",
                       help="filename with locations"  )

    parser.add_option( "-m", "--master", dest="master", type="string",
                       help="the master determines the frame."  )

    parser.set_defaults(
        filename_locations = None,
        gap_chars = "-.",
        master = None
        )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    mali = Mali.Mali()

    mali.readFromFile( sys.stdin )

    identifiers = mali.getIdentifiers()

    aligned_columns, aligned_exons = getAlignedColumns( mali, options )

    map_id2location = {}

    if options.filename_locations:
        map_id2location = IOTools.ReadMap( open(options.filename_locations, "r") )
        
    options.stdout.write( Prediction.Prediction().getHeader() + "\n" )

    nid = 1
    
    for identifier in identifiers:

        if options.loglevel >= 2:
            options.stdlog.write("# processing %s\n" % (identifier) )

        entry = mali.getEntry( identifier )
        
        sequence = entry.mString
        if sequence[0] not in string.lowercase:
            raise "all sequences should start with an exon."

        was_exon = True
        d = 0
        alignment = []
        carry_over = 0

        last_codon = []
        codon = []
        nchars_in_codon = 0
        n = 0

        last_master_residue = 0
        master_residue = 0
        for column in range(len(sequence)):

            c = sequence[column]
            is_gap = c in options.gap_chars
            is_aligned = column in aligned_columns
            is_exon = column in aligned_exons

            if is_gap:
                continue

            if is_exon:
                master_residue = aligned_exons[column]                
                codon.append( (n, master_residue) )

            n += 1

            # check if we have a complete codon
            if is_exon:
                # A codon is complete, if it ends at frame 2 or
                # it spans more than one codons in the master.
                # Gaps in the master that are a multiple of 3 are ignored
                d = master_residue - last_master_residue - 1
                
                if master_residue % 3 == 2 or (d % 3 != 0 and d > 0):
                    
                    if last_codon:
                        d = codon[0][0] - last_codon[-1][0] - 1
                        if d > 0:
                            # add in-frame introns
                            if d > 10:
                                alignment.append( ["5", 0, 2] )
                                alignment.append( ["I", 0, d - 4] )
                                alignment.append( ["3", 0, 2] )                                
                            else:
                                raise "untreated case"
                        
                    alignment += processCodon( codon )
                    last_codon = codon
                    codon = []
                    
            last_master_residue = master_residue
            
        last = alignment[0]
        new_alignment = []
        for this in alignment[1:]:
            if this[0] == last[0]:
                last[1] += this[1]
                last[2] += this[2]
                continue
            
            new_alignment.append( last )
            last = this
            
        new_alignment.append(last)

        if options.loglevel >= 4:
            options.stdlog.write("# output=%s\n" % (str(new_alignment)))

        assert( new_alignment[-1][2] % 3 == 0)

        lalignment = sum( map( lambda x: x[2], new_alignment ) )
        
        prediction = Prediction.Prediction()
        
        prediction.mQueryToken = identifier
        
        genomic_sequence = re.sub("[%s]" % options.gap_chars, "", mali[identifier])

        prediction.mPredictionId = nid
        nid += 1
        
        if identifier in map_id2location:

            prediction.mSbjctToken, prediction.mSbjctStrand, sfrom, sto = map_id2location[identifier].split(":")[:4]

            prediction.mSbjctGenomeFrom = int(sfrom) + entry.mFrom
            prediction.mSbjctGenomeTo = int(sto)            

        else:
            prediction.mSbjctToken = "unk"
            prediction.mSbjctStrand = "+"
            prediction.mSbjctGenomeFrom = 0

        prediction.mQueryCoverage = 100
        prediction.mPercentIdentity = 100
        prediction.mPercentSimilarity = 100

        prediction.mQueryLength = prediction.mQueryTo
        
        prediction.mSbjctGenomeTo = prediction.mSbjctGenomeFrom + lalignment
            
        prediction.mMapPeptide2Genome = new_alignment
        prediction.mAlignmentString = string.join( map( \
                                      lambda x: string.join(map(str, x), " "),
                                      prediction.mMapPeptide2Genome), " ")

        prediction.mMapPeptide2Translation, prediction.mTranslation = Genomics.Alignment2PeptideAlignment( \
               prediction.mMapPeptide2Genome, 0, 0, genomic_sequence )

        (prediction.mNIntrons, prediction.mNFrameShifts, prediction.mNGaps, prediction.mNSplits, prediction.mNStopCodons, disruptions) = \
                          Genomics.CountGeneFeatures( 0,
                                                      prediction.mMapPeptide2Genome,
                                                      genomic_sequence )
        
        options.stdout.write( str(prediction) + "\n" )
    
    E.Stop()
