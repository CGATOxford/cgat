################################################################################
#   Gene prediction pipeline 
#
#   $Id: Prediction.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
"""
Prediction.py - Gene prediction with exonerate/genewise
========================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

"""
import os, sys, string, re, getopt, tempfile, copy

import Genomics
import alignlib_lite

## number of nucleotides to ignore for counting stop-codons
## This used to be three, but is now set to 0 in order to keep
## accounting for stop codons consistent with the disruptions
## table.
GlobalBorderStopCodons = 0
GlobalStopCodons = ("TAG", "TAA", "TGA")

class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class AlignmentError(Error):
    """Exception raised by inconsistent alignments.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class Prediction:
    
    def __init__(self, expand = 1):

        self.mExpand = expand
        
        self.mPredictionId = 0
        self.mQueryToken = 0
        self.mQueryFrom = 0
        self.mQueryTo = 0
        self.mSbjctToken = 0
        self.mSbjctStrand = 0
        self.mSbjctFrom = 0
        self.mSbjctTo = 0
        self.mRank = 0
        self.score = 0
        self.mQueryLength = 0
        self.mQueryCoverage = 0
        self.mNGaps = 0
        self.mNFrameShifts = 0
        self.mNIntrons = 0
        self.mNSplits = 0
        self.mNStopCodons = 0
        self.mPercentIdentity = 0
        self.mPercentSimilarity = 0
        self.mTranslation = ""
        self.mSbjctGenomeFrom = 0
        self.mSbjctGenomeTo = 0
        self.mAlignmentString = ""
        self.mQueryAli = ""
        self.mSbjctAli = ""
        
        if self.mExpand:
            self.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()
            self.mMapPeptide2Genome = []
        else:
            self.mMapPeptide2Translation = None
            self.mMapPeptide2Genome = None
        self.mNAssembled = 0

    def __del__(self):
        if self.mExpand:
            self.mMapPeptide2Translation.clear()
            self.mMapPeptide2Genome = []
        self.mTranslation = ""
        self.mAlignmentString = ""
        
    def __eq__(self, other ):

        return self.mQueryToken == other.mQueryToken and \
               self.mSbjctToken == other.mSbjctToken and \
               self.mSbjctStrand == other.mSbjctStrand and \
               self.score == other.score and \
               self.mQueryFrom == other.mQueryFrom and \
               self.mQueryTo == other.mQueryTo and \
               self.mSbjctGenomeFrom == other.mSbjctGenomeFrom and \
               self.mSbjctGenomeTo == other.mSbjctGenomeTo and \
               self.mAlignmentString == other.mAlignmentString 

    def cleanAlignmentInfo( self ):
        """remove alignment information (to save space)."""
        self.mTranslation, self.mAlignmentString, self.mQueryAli, self.mSbjctAli = None, None, None, None
        
    def getCopy( self ):
        """return a new copy.
        """

        new_entry = Prediction()

        new_entry.mExpand = self.mExpand 
        
        new_entry.mPredictionId = self.mPredictionId 
        new_entry.mQueryToken = self.mQueryToken 
        new_entry.mQueryFrom = self.mQueryFrom 
        new_entry.mQueryTo = self.mQueryTo 
        new_entry.mSbjctToken = self.mSbjctToken 
        new_entry.mSbjctStrand = self.mSbjctStrand 
        new_entry.mSbjctFrom = self.mSbjctFrom 
        new_entry.mSbjctTo = self.mSbjctTo 
        new_entry.mRank = self.mRank 
        new_entry.score = self.score 
        new_entry.mQueryLength = self.mQueryLength 
        new_entry.mQueryCoverage = self.mQueryCoverage 
        new_entry.mNGaps = self.mNGaps 
        new_entry.mNFrameShifts = self.mNFrameShifts 
        new_entry.mNIntrons = self.mNIntrons 
        new_entry.mNSplits = self.mNSplits 
        new_entry.mNStopCodons = self.mNStopCodons 
        new_entry.mPercentIdentity = self.mPercentIdentity 
        new_entry.mPercentSimilarity = self.mPercentSimilarity 
        new_entry.mTranslation = self.mTranslation 
        new_entry.mSbjctGenomeFrom = self.mSbjctGenomeFrom 
        new_entry.mSbjctGenomeTo = self.mSbjctGenomeTo 
        new_entry.mAlignmentString = self.mAlignmentString 
        new_entry.mQueryAli = self.mQueryAli 
        new_entry.mSbjctAli = self.mSbjctAli 

        if self.mExpand:
            new_entry.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()
            alignlib_lite.py_copyAlignment( new_entry.mMapPeptide2Translation, self.mMapPeptide2Translation)
            new_entry.mMapPeptide2Genome = Genomics.String2Alignment( new_entry.mAlignmentString) 
        else:
            new_entry.mMapPeptide2Translation = self.mMapPeptide2Translation = None
            new_entry.mMapPeptide2Genome = self.mMapPeptide2Genome = None

        return new_entry

    def Expand( self ):
        self.mExpand = True
        
        if self.mMapPeptide2Translation.getLength() > 0:
            f = alignlib_lite.py_AlignmentFormatEmissions( self.mMapPeptide2Translation )
            self.mQueryAli, self.mSbjctAli = f.mRowAlignment, f.mColAlignment
            self.mQueryFrom = self.mMapPeptide2Translation.getRowFrom()
            self.mQueryTo = self.mMapPeptide2Translation.getRowTo()
            self.mSbjctFrom = self.mMapPeptide2Translation.getColFrom()
            self.mSbjctTo = self.mMapPeptide2Translation.getColTo()

        self.mMapPeptide2Genome = Genomics.String2Alignment( self.mAlignmentString )

    def Add( self, const_other,
             combine_contig = False,
             allow_overlap = False,
             contig_size = 0,
             combine_queries = False,
             as_intron = False ):
        """add one entry to another.

        This procedure allows to add
        
        - predictions on different contigs if combine_contig = True
        - overlapping predictions on the same query if allow_overlap = True
        - results from different queries if combine_queries = True

        - if as_intron is set to true, the new fragment is added as an intron.
        
        """

        ## create working copies of each prediction
        other = const_other.getCopy()
        this  = self.getCopy()

        other.Expand()
        this.Expand()

        if as_intron:
            code = "I"
        else:
            code = "P"

        ## check for query overlaps
        if this.mQueryToken == other.mQueryToken:

            query_overlap = max( 0, min(this.mQueryTo, other.mQueryTo) -\
                                 max(this.mQueryFrom, other.mQueryFrom) + 1)

            if query_overlap > 0:

                if allow_overlap:
                    overlap = query_overlap
                    ## if queries overlap, truncate this before adding the other
                    this.mMapPeptide2Translation.removeRowRegion( this.mQueryTo - overlap + 1, this.mQueryTo )
                    other.mMapPeptide2Translation.moveAlignment( 0, -overlap )
                    this.mQueryTo -= overlap
                    this.mTranslation = this.mTranslation[:-overlap]

                    ## remove aligned residues from the back
                    for x in range(len(this.mMapPeptide2Genome) - 1, 0, -1):
                        if this.mMapPeptide2Genome[x][1] <= overlap:
                            overlap -= this.mMapPeptide2Genome[x][1]
                            del this.mMapPeptide2Genome[x]
                        else:
                            break
                    this.mMapPeptide2Genome[-1] = (this.mMapPeptide2Genome[-1][0],
                                                   this.mMapPeptide2Genome[-1][1] - overlap,
                                                   this.mMapPeptide2Genome[-1][2] - overlap * 3)
                else:
                    raise ValueError, "refusing to add overlapping entries: overlap = %i, queries:\n%s\n%s\n, set allow_overlap = True " % (query_overlap, str(this), str(other))


        else:
            if not combine_queries:
                raise ValueError, "refusing to add different queries - set combine_queries = True."

        if this.mSbjctToken != other.mSbjctToken or \
               this.mSbjctStrand != other.mSbjctStrand :
            if combine_contig:
                this.mSbjctToken += "-" + other.mSbjctToken
                this.mSbjctStrand += other.mSbjctStrand
            else:
                raise ValueError, "can not add different sbjct."                

        sbjct_overlap = max(0, min(this.mSbjctGenomeTo, other.mSbjctGenomeTo) -\
                            max(this.mSbjctGenomeFrom, other.mSbjctGenomeFrom), 0)

        if sbjct_overlap > 0:
            if not combine_contig:
                raise ValueError, "refusing to add overlapping entries: overlap = %i, sbjct:\n%s\n%s\n" % (sbjct_overlap, str(this), str(other))

        if this.mSbjctToken == other.mSbjctToken:

            ## set precedence
            if this.mSbjctGenomeFrom < other.mSbjctGenomeFrom:
                first = this
                second = other
            else:
                first = other
                second = this

            ## get length of gap
            d_na = second.mSbjctGenomeFrom - first.mSbjctGenomeTo

            if this.mQueryToken != other.mQueryToken:
                d_aa = first.mQueryLength - first.mQueryTo                 
                # create a new virtual query by concatenating
                # the two queries
                this.mQueryToken += "-" + other.mQueryToken

                # sort out the alignment
                second.mMapPeptide2Translation.moveAlignment( first.mQueryLength, 0 )

                this.mQueryLength = first.mQueryLength + second.mQueryLength

            else:
                d_aa = second.mQueryFrom - first.mQueryTo - 1
            
            this.mSbjctGenomeFrom = min(this.mSbjctGenomeFrom, other.mSbjctGenomeFrom )
            this.mSbjctGenomeTo = max(this.mSbjctGenomeTo, other.mSbjctGenomeTo )        

            this.mMapPeptide2Genome = first.mMapPeptide2Genome + [(code, d_aa, d_na)] + second.mMapPeptide2Genome
            this.mTranslation = first.mTranslation + second.mTranslation

            second.mMapPeptide2Translation.moveAlignment( 0, first.mSbjctTo - 1 )
            
        else:
            ## join on different contigs
            d_na = contig_size - this.mSbjctGenomeTo + other.mSbjctGenomeFrom + query_overlap * 3
            d_aa = other.mQueryFrom - this.mQueryTo - 1
            this.mMapPeptide2Genome += [(code, d_aa, d_na),] + other.mMapPeptide2Genome
            this.mTranslation += other.mTranslation 
            other.mMapPeptide2Translation.moveAlignment( 0, this.mSbjctTo - 1 )

            this.mSbjctGenomeFrom = this.mSbjctGenomeFrom
            this.mSbjctGenomeTo = contig_size + other.mSbjctGenomeTo

        ## now fill self from first and this
        self.mQueryToken = first.mQueryToken
        self.mQueryLength = this.mQueryLength
        
        nthis  = this.mMapPeptide2Translation.getLength() - this.mMapPeptide2Translation.getNumGaps()
        nother = other.mMapPeptide2Translation.getLength() - other.mMapPeptide2Translation.getNumGaps()

        self.mMapPeptide2Genome = first.mMapPeptide2Genome
        self.mSbjctGenomeFrom = this.mSbjctGenomeFrom
        self.mSbjctGenomeTo= this.mSbjctGenomeTo
        
        ## there might be some reference counting issues, thus
        ## do it the explicit way.
        alignlib_lite.py_addAlignment2Alignment( this.mMapPeptide2Translation, other.mMapPeptide2Translation)
        self.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()
        alignlib_lite.py_addAlignment2Alignment( self.mMapPeptide2Translation, this.mMapPeptide2Translation )
        
        self.mTranslation = this.mTranslation
        
        self.mQueryFrom = self.mMapPeptide2Translation.getRowFrom()
        self.mQueryTo = self.mMapPeptide2Translation.getRowTo()
        self.mSbjctFrom = self.mMapPeptide2Translation.getColFrom()
        self.mSbjctTo = self.mMapPeptide2Translation.getColTo()
        
        self.mQueryCoverage = 100.0 * (self.mQueryTo - self.mQueryFrom + 1) / float(self.mQueryLength)

        self.mAlignmentString = string.join( map( \
                                      lambda x: string.join(map(str, x), " "),
                                      self.mMapPeptide2Genome), " ")

        f = alignlib_lite.py_AlignmentFormatEmssions( self.mMapPeptide2Translation )
        self.mQueryAli, self.mSbjctAli = f.mRowAlignment, f.mColAlignment

        ## summary parameters
        self.mRank = max( this.mRank, other.mRank)
        self.score += other.score
        self.mNGaps += other.mNGaps
        self.mNFrameShifts += other.mNFrameShifts
        self.mNIntrons += other.mNIntrons + 1
        self.mNStopCodons += other.mNStopCodons
        
        nnew = self.mMapPeptide2Translation.getLength() - self.mMapPeptide2Translation.getNumGaps()
        
        self.mPercentIdentity = min( 100.0, (self.mPercentIdentity * nthis + other.mPercentIdentity * nother) / nnew )
        self.mPercentSimilarity = min( 100.0, (self.mPercentSimilarity * nthis + other.mPercentSimilarity * nother) / nnew )

        self.mNAssembled += 1 + other.mNAssembled

    def InvertGenomicCoordinates( self, lgenome ):
        """invert genomic alignment on sequence.

        Negative strand is calculated from the other end.
        """
        if self.mSbjctStrand == "-":
            x = self.mSbjctGenomeFrom 
            self.mSbjctGenomeFrom = lgenome - self.mSbjctGenomeTo
            self.mSbjctGenomeTo   = lgenome - x
        
    def shiftGenomicRegion( self, forward_offset, reverse_offset = 0 ):
        """shift aligned genomic region by offset.

        if reverse offset is given, negative strand predictions
        are shifted by this amount.
        """
        if reverse_offset:
            if self.mSbjctStrand == "+":
                self.mSbjctGenomeFrom += forward_offset
                self.mSbjctGenomeTo += forward_offset
            else:
                self.mSbjctGenomeFrom += reverse_offset                
                self.mSbjctGenomeTo += reverse_offset
        else:
            self.mSbjctGenomeFrom += forward_offset
            self.mSbjctGenomeTo += forward_offset

    def shiftQueryRegion( self, offset ):
        """shift aligned genomic region by offset."""
        self.mQueryFrom += offset
        self.mQueryTo += offset
        if self.mExpand:
            self.mMapPeptide2Translation.moveAlignment( offset, 0)

    def setStrand( self, strand ):
        """set strand."""
        self.mSbjctStrand = strand 

    def getHeader(self):
        """return header."""
        return "\t".join( ("id",
                          "query_token",
                          "sbjct_token",
                          "sbjct_strand",
                          "rank",
                          "score",
                          "query_from",
                          "query_to",
                          "query_ali",
                          "sbjct_from",
                          "sbjct_to",
                          "sbjct_ali",
                          "query_length",
                          "query_coverage",
                          "ngaps",
                          "nframeshifts",
                          "nintrons",
                          "nsplits",
                          "nstops",
                          "pidentity",
                          "psimilarity",
                          "translation",
                          "sbjct_genome_from",
                          "sbjct_genome_to",
                          "map_peptide2genome",
                          "nassembled" ) )
    def __str__( self ):
        """get a string representation of results."""

            
        if self.mExpand:
            if self.mMapPeptide2Translation.getLength() > 0:
                f = alignlib_lite.py_AlignmentFormatEmissions( self.mMapPeptide2Translation )
                row_ali, col_ali = f.mRowAlignment, f.mColAlignment
                self.mQueryFrom = self.mMapPeptide2Translation.getRowFrom()
                self.mQueryTo = self.mMapPeptide2Translation.getRowTo()
                self.mSbjctFrom = self.mMapPeptide2Translation.getColFrom()
                self.mSbjctTo = self.mMapPeptide2Translation.getColTo()
            else:
                row_ali, col_ali = "", ""
        else:
            row_ali = self.mQueryAli
            col_ali = self.mSbjctAli
            
        if self.mPredictionId:
            return string.join( map(str, (\
                self.mPredictionId,
                self.mQueryToken,
                self.mSbjctToken,
                self.mSbjctStrand,
                self.mRank,
                self.score,
                self.mQueryFrom,
                self.mQueryTo,
                row_ali,
                self.mSbjctFrom,
                self.mSbjctTo,
                col_ali,
                self.mQueryLength,
                self.mQueryCoverage, 
                self.mNGaps,
                self.mNFrameShifts,
                self.mNIntrons,
                self.mNSplits,
                self.mNStopCodons,
                "%5.2f" % self.mPercentIdentity,
                "%5.2f" % self.mPercentSimilarity,
                self.mTranslation,
                self.mSbjctGenomeFrom,
                self.mSbjctGenomeTo,
                self.mAlignmentString,
                self.mNAssembled,
                )), "\t")
        else:
            return string.join( map(str, (\
              self.mQueryToken,
              self.mSbjctToken,
              self.mSbjctStrand,
              self.mRank,
              self.score,
              self.mQueryFrom,
              self.mQueryTo,
              row_ali,
              self.mSbjctFrom,
              self.mSbjctTo,
              col_ali,
              self.mQueryLength,
              self.mQueryCoverage, 
              self.mNGaps,
              self.mNFrameShifts,
              self.mNIntrons,
              self.mNSplits,
              self.mNStopCodons,
                "%5.2f" % self.mPercentIdentity,
                "%5.2f" % self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom,
              self.mSbjctGenomeTo,
              self.mAlignmentString,
              self.mNAssembled,
              )), "\t")

    def fillFromTable( self, table_row ):

        if len(table_row) == 25:
            ( self.mPredictionId,
              self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString) = table_row
        elif len(table_row) == 26:
            ( self.mPredictionId,
              self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString,
              self.mNAssembled) = table_row[:26]
        elif len(table_row) > 26:
            ( self.mPredictionId,
              self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString,
              self.mNAssembled) = table_row[:26]
        else:
            raise ValueError, "unknown format: %i fields" % len(data)
            sys.exit(0)
            
        if self.mExpand:
            self.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()

            if self.mQueryAli != "" and self.mSbjctAli != "":
                alignlib_lite.py_AlignmentFormatEmissions( self.mQueryFrom, self.mQueryAli,
                                                   self.mSbjctFrom, self.mSbjctAli ).copy( self.mMapPeptide2Translation )

            self.mMapPeptide2Genome = Genomics.String2Alignment( self.mAlignmentString )

    def setTranslation( self, genomic_sequence ):
        """set translation from genomic sequence."""
        self.mMapPeptide2Translation, self.mTranslation = Genomics.Alignment2PeptideAlignment( \
                self.mMapPeptide2Genome, self.mQueryFrom, self.mSbjctGenomeFrom, genomic_sequence )
        
    def getExonBoundariesOnGenome( self ):
        """return a list of exons for this prediction.
        """
        return map( lambda x: (x[3], x[4]),
                    Genomics.Alignment2ExonBoundaries( self.mMapPeptide2Genome,
                                                       self.mQueryFrom,
                                                       self.mSbjctGenomeFrom ))
        
    def write( self ):
        """write a tab separated list of results.
        """
        print str(self)
        
    def read( self, line ):

        data = string.split( line[:-1], "\t")
        
        if len(data) == 26:
            ( self.mPredictionId,
              self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString, self.mNAssembled,
              ) = data
        elif len(data) == 25:
            ( self.mPredictionId,
              self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString,
              ) = data
        elif len(data) == 24:            
            ( self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              self.mAlignmentString,
              ) = data
        elif len(data) == 23:
            ( self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
              self.mRank, self.score,
              self.mQueryFrom, self.mQueryTo, self.mQueryAli,
              self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
              self.mQueryLength, self.mQueryCoverage, 
              self.mNGaps, self.mNFrameShifts, self.mNIntrons,
              self.mNSplits, self.mNStopCodons,
              self.mPercentIdentity, self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
              ) = data
            self.mAlignmentString = ""
        else:
            raise ValueError, "unknown format: %i fields in line %s" % (len(data), line[:-1])

        (self.score, self.mQueryCoverage, self.mPercentIdentity, self.mPercentSimilarity) = map (\
            float, (self.score, self.mQueryCoverage, self.mPercentIdentity, self.mPercentSimilarity))

        (self.mPredictionId, 
         self.mQueryFrom, self.mQueryTo, self.mQueryLength,
         self.mSbjctFrom, self.mSbjctTo,
         self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
         self.mNGaps, self.mNIntrons, self.mNSplits, self.mNStopCodons,
         self.mNFrameShifts, self.mNAssembled) = map (\
            int, ( self.mPredictionId,
                   self.mQueryFrom, self.mQueryTo, self.mQueryLength,
                   self.mSbjctFrom, self.mSbjctTo,
                   self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
                   self.mNGaps, self.mNIntrons, self.mNSplits, self.mNStopCodons,
                   self.mNFrameShifts, self.mNAssembled))

        if self.mExpand:        
            self.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()

            if self.mQueryAli != "" and self.mSbjctAli != "":
                
                alignlib_lite.py_AlignmentFormatExplicit(
                    self.mQueryFrom, self.mQueryAli,
                    self.mSbjctFrom, self.mSbjctAli).copy( self.mMapPeptide2Translation )

            self.mMapPeptide2Genome = Genomics.String2Alignment( self.mAlignmentString )

##--------------------------------------------------------------------------
class Predictions :
    mLogLevel = 0

    def __init__( self ):
        self.mMatches = []

    def __eq__(self, other):
        if not other: return 0        

        if not other.mMatches: return 0
        if not self.mMatches: return 0
        return self.mMatches[0] == other.mMatches[0]

    def __getitem__(self, key ):
        return self.mMatches[key]

    def __setitem__(self, key, value ):
        self.mMatches[key] = value

    def __len__(self):
        return len(self.mMatches)

    def __str__(self):
        return string.join( map(str, self.mMatches), "\n")

    def append( self, entry ):
        return self.mMatches.append( entry )
    
    def clear( self ):
        for x in self.mMatches:
            x.__del__()
        self.mMatches = []
    
    def combine(self, other):
        self.mMatches += other.mMatches

    def add( self, prediction ):
        """prediction to predictions."""
        self.mMatches.append(prediction)

    def reRank( self, offset = 1 ):
        """rerank predictions based on current sorting order."""
        
        for x in range(len(self.mMatches)):
            self.mMatches[x].mRank = offset + x

    def getBestMatch( self ):
        self.sort( lambda x,y: cmp(-x.score, -y.score)  )
        return self.mMatches[0]
        
    def getNumMatches( self ):
        """return number of matches."""
        return len(self.mMatches)

    def shiftGenomicRegion( self, forward_offset, reverse_offset = 0):
        """shift aligned genomic region by offset for all entries.

        if reverse offset is given, negative strand predictions
        are shifted by this amount.
        """
        for x in self.mMatches:
            x.shiftGenomicRegion( forward_offset, reverse_offset)

    def shiftQueryRegion( self, offset ):
        """shift aligned genomic region by offset for all entries."""
        for x in self.mMatches:
            x.shiftQueryRegion( offset )

    def setStrand( self, strand ):
        """set strand."""
        for x in self.mMatches:
            x.setStrand( strand )
    
    def write( self ):
        """write a tab separated list of results.
        """
        for x in self.mMatches:
            x.write()

    def sort( self, cmp_function = None ):
        """sort predictions.
        """
        self.mMatches.sort( cmp_function )

    def getRange( self ):
        """return a prediction covering all predictions."""

        p = self.mMatches[0].getCopy()

        p.mExpand = False
        
        if len(self.mMatches) == 1: return p
        
        for pp in self.mMatches[1:]:

            p.mQueryFrom       = min(p.mQueryFrom      , pp.mQueryFrom )
            p.mQueryTo         = max(p.mQueryTo        , pp.mQueryTo )
            p.mSbjctGenomeFrom = min(p.mSbjctGenomeFrom, pp.mSbjctGenomeFrom )
            p.mSbjctGenomeTo   = max(p.mSbjctGenomeTo  , pp.mSbjctGenomeTo )
            p.score += pp.score

        return p


def iterate( infile ):
    '''iterate over predicions in infile.'''

    for line in infile:
        if line.startswith("#"): continue
        if line.startswith("id\tquery_token"): continue
        p = Prediction()
        p.read( line )
        yield p
        

    
