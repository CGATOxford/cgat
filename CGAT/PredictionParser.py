################################################################################
#   Gene prediction pipeline 
#
#   $Id: PredictionParser.py 2764 2009-09-04 16:55:07Z andreas $
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
PredictionParser.py - Parser for exonerate/genewise output
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

"""
import os, sys, string, re, getopt, tempfile, copy

import Genomics, GTF
try: import alignlib_lite
except ImportError: pass

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


class PredictionParserEntry :

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

    def CleanAlignmentInfo( self ):
        """remove alignment information (to save space)."""
        self.mTranslation, self.mAlignmentString, self.mQueryAli, self.mSbjctAli = None, None, None, None
        
    def GetCopy( self ):
        """return a new copy.
        """

        new_entry = PredictionParserEntry()

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

    def Add( self, const_other,
             combine_contig = False,
             allow_overlap = False,
             contig_size = 0):
        """add one entry to another.

        This procedure allows to add
        - predictions on different contigs if combine_contig = True
        - overlapping predictions on the same query if allow_overlap = True
        """

        other = const_other.GetCopy()

        if not self.mExpand or not other.mExpand:
            raise ValueError, "both object must be expanded."
        
        if self.mQueryToken != other.mQueryToken:
            raise ValueError, "can not add different query."

        if self.mSbjctToken != other.mSbjctToken or \
               self.mSbjctStrand != other.mSbjctStrand :
            if combine_contig:
                self.mSbjctToken += "-" + other.mSbjctToken
                self.mSbjctStrand += other.mSbjctStrand
            else:
                raise ValueError, "can not add different sbjct."                

        query_overlap = max( 0, min(self.mQueryTo, other.mQueryTo) -\
                             max(self.mQueryFrom, other.mQueryFrom) + 1)
        
        if query_overlap > 0:
            
            if allow_overlap:
                overlap = query_overlap
                ## if queries overlap, truncate this before adding the other
                self.mMapPeptide2Translation.removeRowRegion( self.mQueryTo - overlap + 1, self.mQueryTo )
                other.mMapPeptide2Translation.moveAlignment( 0, -overlap )
                self.mQueryTo -= overlap
                self.mTranslation = self.mTranslation[:-overlap]
                
                ## remove aligned residues from the back
                for x in range(len(self.mMapPeptide2Genome) - 1, 0, -1):
                    if self.mMapPeptide2Genome[x][1] <= overlap:
                        overlap -= self.mMapPeptide2Genome[x][1]
                        del self.mMapPeptide2Genome[x]
                    else:
                        break
                self.mMapPeptide2Genome[-1] = (self.mMapPeptide2Genome[-1][0],
                                               self.mMapPeptide2Genome[-1][1] - overlap,
                                               self.mMapPeptide2Genome[-1][2] - overlap * 3)
            else:
                raise ValueError, "can not add overlapping entries: overlap = %i, queries:\n%s\n%s\n" % (query_overlap, str(self), str(other))

        sbjct_overlap = max(0, min(self.mSbjctGenomeTo, other.mSbjctGenomeTo) -\
                            max(self.mSbjctGenomeFrom, other.mSbjctGenomeFrom), 0)

        if sbjct_overlap > 0:
            if not combine_contig:
                raise ValueError, "can not add overlapping entries: overlap = %i, sbjct:\n%s\n%s\n" % (sbjct_overlap, str(self), str(other))
        
        if self.mSbjctToken == other.mSbjctToken:
            ## join on the same contig
            if self.mSbjctGenomeFrom < other.mSbjctGenomeFrom:
                ## get length of gap
                d_na = other.mSbjctGenomeFrom - self.mSbjctGenomeTo
                d_aa = other.mQueryFrom       - self.mQueryTo - 1
                self.mMapPeptide2Genome += [("P", d_aa, d_na)] + other.mMapPeptide2Genome
                self.mTranslation += other.mTranslation
                other.mMapPeptide2Translation.moveAlignment( 0, self.mSbjctTo - 1 )
            else:
                d_na = self.mSbjctGenomeFrom - other.mSbjctGenomeTo
                d_aa = self.mQueryFrom - other.mQueryTo - 1
                self.mTranslation = other.mTranslation + self.mTranslation
                self.mMapPeptide2Genome = other.mMapPeptide2Genome + [("P", d_aa, d_na)] + self.mMapPeptide2Genome
                self.mMapPeptide2Translation.moveAlignment( 0, other.mSbjctTo - 1 )
                
            self.mSbjctGenomeFrom = min(self.mSbjctGenomeFrom, other.mSbjctGenomeFrom )
            self.mSbjctGenomeTo = max(self.mSbjctGenomeTo, other.mSbjctGenomeTo )        
                
        else:
            ## join on different contigs
            d_na = contig_size - self.mSbjctGenomeTo + other.mSbjctGenomeFrom + query_overlap * 3
            d_aa = other.mQueryFrom - self.mQueryTo - 1
            self.mMapPeptide2Genome += [("P", d_aa, d_na),] + other.mMapPeptide2Genome
            self.mTranslation += other.mTranslation 
            other.mMapPeptide2Translation.moveAlignment( 0, self.mSbjctTo - 1 )

            self.mSbjctGenomeFrom = self.mSbjctGenomeFrom
            self.mSbjctGenomeTo = contig_size + other.mSbjctGenomeTo

        self.mRank = max( self.mRank, other.mRank)
        self.score += other.score
        self.mQueryCoverage += other.mQueryCoverage
        self.mNGaps += other.mNGaps
        self.mNFrameShifts += other.mNFrameShifts
        self.mNIntrons += other.mNIntrons + 1
        self.mNStopCodons += other.mNStopCodons

        nself  = self.mMapPeptide2Translation.getLength() - self.mMapPeptide2Translation.getNumGaps()
        nother = other.mMapPeptide2Translation.getLength() - other.mMapPeptide2Translation.getNumGaps()
        alignlib_lite.py_addAlignment2Alignment( self.mMapPeptide2Translation, other.mMapPeptide2Translation )

        self.mQueryFrom = self.mMapPeptide2Translation.getRowFrom()
        self.mQueryTo = self.mMapPeptide2Translation.getRowTo()
        self.mSbjctFrom = self.mMapPeptide2Translation.getColFrom()
        self.mSbjctTo = self.mMapPeptide2Translation.getColTo()
        
        nnew = self.mMapPeptide2Translation.getLength() - self.mMapPeptide2Translation.getNumGaps()
        
        self.mPercentIdentity = (self.mPercentIdentity * nself + other.mPercentIdentity * nother) / nnew
        self.mPercentSimilarity = (self.mPercentSimilarity * nself + other.mPercentSimilarity * nother) / nnew

        self.mAlignmentString = string.join( map( \
                                      lambda x: string.join(map(str, x), " "),
                                      self.mMapPeptide2Genome), " ")

        self.mNAssembled += 1 + other.mNAssembled

    def InvertGenomicCoordinates( self, lgenome ):
        """invert genomic alignment on sequence.

        Negative strand is calculated from the other end.
        """
        if self.mSbjctStrand == "-":
            x = self.mSbjctGenomeFrom 
            self.mSbjctGenomeFrom = lgenome - self.mSbjctGenomeTo
            self.mSbjctGenomeTo   = lgenome - x
        
    def ShiftGenomicRegion( self, offset ):
        """shift aligned genomic region by offset."""
        self.mSbjctGenomeFrom += offset
        self.mSbjctGenomeTo += offset

    def ShiftQueryRegion( self, offset ):
        """shift aligned genomic region by offset."""
        self.mQueryFrom += offset
        self.mQueryTo += offset
        if self.mExpand:
            self.mMapPeptide2Translation.moveAlignment( offset, 0)

    def SetStrand( self, strand ):
        """set strand."""
        self.mSbjctStrand = strand 

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
                self.mPercentIdentity,
                self.mPercentSimilarity,
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
              self.mPercentIdentity,
              self.mPercentSimilarity,
              self.mTranslation,
              self.mSbjctGenomeFrom,
              self.mSbjctGenomeTo,
              self.mAlignmentString,
              self.mNAssembled,
              )), "\t")

    def FillFromTable( self, table_row ):

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
                f = alignlib_lite.py_AlignmentFormatEmissions()
                f.mRowFrom = self.mQueryFrom
                f.mColFrom = self.SbjctFrom
                f.mRowAlignment = self.mQueryAli
                f.mColAlignment = self.mSbjctAli

                f.copy( self.mMapPeptide2Translation )

            self.mMapPeptide2Genome = Genomics.String2Alignment( self.mAlignmentString )

    def SetTranslation( self, genomic_sequence ):
        """set translation from genomic sequence."""
        self.mMapPeptide2Translation, self.mTranslation = Genomics.Alignment2PeptideAlignment( \
                self.mMapPeptide2Genome, self.mQueryFrom, self.mSbjctGenomeFrom, genomic_sequence )
        
    def GetExonBoundariesOnGenome( self ):
        """return a list of exons for this prediction.
        """
        return map( lambda x: (x[3], x[4]),
                    Genomics.Alignment2ExonBoundaries( self.mMapPeptide2Genome,
                                                       self.mQueryFrom,
                                                       self.mSbjctGenomeFrom ))
        
    def Write( self ):
        """write a tab separated list of results.
        """
        print str(self)
        
    def Read( self, line ):

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

        (self.mQueryFrom, self.mQueryTo, self.mQueryLength,
         self.mSbjctFrom, self.mSbjctTo,
         self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
         self.mNGaps, self.mNIntrons, self.mNSplits, self.mNStopCodons,
         self.mNFrameShifts, self.mNAssembled) = map (\
            int, ( self.mQueryFrom, self.mQueryTo, self.mQueryLength,
                   self.mSbjctFrom, self.mSbjctTo,
                   self.mSbjctGenomeFrom, self.mSbjctGenomeTo,
                   self.mNGaps, self.mNIntrons, self.mNSplits, self.mNStopCodons,
                   self.mNFrameShifts, self.mNAssembled))

        if self.mExpand:        
            self.mMapPeptide2Translation = alignlib_lite.py_makeAlignmentVector()

            if self.mQueryAli != "" and self.mSbjctAli != "":
                f = alignlib_lite.py_AlignmentFormatEmissions()
                f.mRowFrom = self.mQueryFrom
                f.mColFrom = self.mSbjctFrom
                f.mRowAlignment = self.mQueryAli
                f.mColAlignment = self.mSbjctAli

                f.copy( self.mMapPeptide2Translation )

            self.mMapPeptide2Genome = Genomics.String2Alignment( self.mAlignmentString )



##-----------------------------------------------------------------------------------------------
class PredictionParserEntryCaleb( PredictionParserEntry ):
    """a parser to convert Calebs predictions.

    The format is:
    ENSG00000000005 ENST00000003616 688.97  96.85   100     1       318     317     17 61 108 142 193 249   NULL    NULL    951     100     X       +       77495591        77511098      MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGGKHFWPETPKKTYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENEEITTTFFEQSVIWVPAEKPIENRDFLKNSKILEICDNVTMYWINPTLIAVSELQDFEEDGEDLHFPTNEKKGIEQNEQWVVPQVKMEKIRHTRQASEEELPINDYTENGLEFDPMLDERGYCCIYCRRGNRYCRRVCEPLLGYYPYPYCYQGGRVICRVIMPCNWWVARMLGRV atggcaaagaatcctccagagaactgtgaggactgtcacattctaaatgcagaagcttttaaatccaagaagatatgtaaatcacttaagatttgtggattggtgtttggtatcctggctctaactctaattgtcctgttttgggggggcaagcacttctggccggagacacccaaaaaaacatatgacatggagcacactttctacagcaatggagagaagaagaagatttacatggaaattgatcctgtgaccagaactgaaatattcagaagtggaaatggcactgatgaaacattggaagtacatgactttaaaaatggttacactggcatctactttgtaggtcttcaaaaatgcttcatcaaaactcagattaaagtgattcctgaattttctgaaccagaggaggaaatagatgagaatgaagaaatcaccacaactttttttgaacagtctgtgatttgggtcccagcagaaaagcctattgaaaatcgagactttctgaaaaattccaaaattctggaaatttgtgataatgtgaccatgtattggatcaatcccactctaatagcagtttcagagttacaagactttgaggaagatggtgaagatcttcactttcctaccaatgaaaaaaaaggtattgaacaaaacgagcagtgggtggtcccccaagtgaagatggaaaagatccgccacaccagacaagcaagtgaggaagaacttccaataaatgactatactgaaaatggactagaatttgatcccatgctggatgagagaggttattgttgtatttactgccgtcgaggcaaccgctactgccgccgcgtctgtgaacctttactaggttactacccgtatccatactgctatcaaggagggcgggttatctgtcgtgtcatcatgccttgcaactggtgggtagcccgcatgctggggagggtc     MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGSKHFWPEVPKKAYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENEEITTTFFEQSVIWVPAEKPIENRDFLKNSKILEICDNVTMYWINPTLISVSELQDFEEEGEDLHFPANEKKGIEQNEQWVVPQVKVEKTRHARQASEEELPINDYTENGIEFDPMLDERGYCCIYCRRGNRYCRRVCEPLLGYYPYPYCYQGGRVICRVIMPCNWWVARMLGRV       MAKNPPENCEDCHILNAEAFKSKKICKSLKICGLVFGILALTLIVLFWGGKHFWPETPKKTYDMEHTFYSNGEKKKIYMEIDPVTRTEIFRSGNGTDETLEVHDFKNGYTGIYFVGLQKCFIKTQIKVIPEFSEPEEEIDENEEITTTFFEQSVIWVPAEKPIENRDFLKNSKILEICDNVTMYWINPTLIAVSELQDFEEDGEDLHFPTNEKKGIEQNEQWVVPQVKMEKIRHTRQASEEELPINDYTENGLEFDPMLDERGYCCIYCRRGNRYCRRVCEPLLGYYPYPYCYQGGRVICRVIMPCNWWVARMLGRV       0       3904_results    0       0.0167  0.2019  0.0825  678.5   272.5   0.00517 0.03308 4.8054  -1484.5 0.2092  NULL    0       317     317     1       0    317      317     317     NULL    NULL    794.08

    Note: Caleb's protein positions sometimes start 0.
    
    """
    def __init__(self):
        
        PredictionParserEntry.__init__(self)
        
    def Read( self, line ):

        data = string.split( line[:-1], "\t")

        (self.mQueryGene, self.mQueryToken, self.score,
         self.mPercentIdentity, self.mQueryCoverage,
         self.mQueryFrom, self.mQueryTo, self.mQueryLength,
         introns, self.mNFrameShifts, self.mNStopCodons,
         psl_score, psl_identity, self.mSbjctToken, self.mSbjctStrand, self.mSbjctGenomeFrom,
         self.mSbjctGenomeTo, self.mTranslation, cdna, query_aligned_string, sbjct_aligned_string, errors,
         f_name, self.mRank) = data[:24]
        
        self.score = float(self.score)
        self.mPercentIdentity = float(self.mPercentIdentity)
        self.mQueryCoverage = float(self.mQueryCoverage)

        self.mSbjctToken = "chr%s" % self.mSbjctToken
        
        (self.mQueryFrom, self.mQueryTo, self.mSbjctGenomeFrom, self.mSbjctGenomeTo) = map(int, \
            (self.mQueryFrom, self.mQueryTo, self.mSbjctGenomeFrom, self.mSbjctGenomeTo))

        if self.mNFrameShifts == "NULL":
            self.mNFrameShifts = 0
        else:
            self.mNFrameShifts = len(string.split(self.mNFrameShifts, " "))

        if self.mNStopCodons == "NULL":
            self.mNStopCodons = 0
        else:
            self.mNStopCodons = len(string.split(self.mNStopCodons, " "))
            
        if introns == "NULL":
            self.mNIntrons == 0
        else:
            self.mNIntrons = len(string.split(introns, " "))

        ## correct bug in parsing: entries like A:F[ttt] or -:F[ttt] in strings
        self.mTranslation = re.sub( "\s*\S{0,1}:(\S)\[\S{3}\]\s*", lambda x: x.groups()[0], self.mTranslation)

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
    
    def Clear( self ):
        for x in self.mMatches:
            x.__del__()
        self.mMatches = []
    
    def Combine(self, other):
        self.mMatches += other.mMatches

    def Add( self, prediction ):
        """prediction to predictions."""
        self.mMatches.append(prediction)

    def ReRank( self, offset = 1 ):
        """rerank predictions based on current sorting order."""
        
        for x in range(len(self.mMatches)):
            self.mMatches[x].mRank = offset + x

    def GetBestMatch( self ):
        return self.mMatches[0]
        
    def GetNumMatches( self ):
        """return number of matches."""
        return len(self.mMatches)

    def ShiftGenomicRegion( self, offset ):
        """shift aligned genomic region by offset for all entries."""
        for x in self.mMatches:
            x.ShiftGenomicRegion( offset )

    def ShiftQueryRegion( self, offset ):
        """shift aligned genomic region by offset for all entries."""
        for x in self.mMatches:
            x.ShiftQueryRegion( offset )

    def SetStrand( self, strand ):
        """set strand."""
        for x in self.mMatches:
            x.SetStrand( strand )
    
    def Write( self ):
        """write a tab separated list of results.
        """
        for x in self.mMatches:
            x.Write()

    def Sort( self, cmp_function = None ):
        """sort predictions.
        """
        self.mMatches.sort( cmp_function )

    def GetRange( self ):
        """return a prediction covering all predictions."""

        p = self.mMatches[0].GetCopy()

        p.mExpand = False
        
        if len(self.mMatches) == 1: return p
        
        for pp in self.mMatches[1:]:

            p.mQueryFrom       = min(p.mQueryFrom      , pp.mQueryFrom )
            p.mQueryTo         = max(p.mQueryTo        , pp.mQueryTo )
            p.mSbjctGenomeFrom = min(p.mSbjctGenomeFrom, pp.mSbjctGenomeFrom )
            p.mSbjctGenomeTo   = max(p.mSbjctGenomeTo  , pp.mSbjctGenomeTo )
            p.score += pp.score

        return p
##--------------------------------------------------------------------------
class PredictionParser :
    
    mLogLevel = 0
    
    def __init__( self ):
        self.mMatches = []
        self.mErrors = []
        self.mNumInput = 0
        self.mNumOutput = 0
        
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

    def GetNumErrors( self ):
        """return number of errors."""
        return len(self.mErrors)

    def GetNumInput( self ):
        """return number of input."""
        return self.mNumInput

    def GetNumOutput( self ):
        """return number of output."""
        return self.mNumOutput

    def Parse( self, lines, peptide_sequence = None, genomic_sequence = None ):
        """parse lines.
        """
        if self.mLogLevel >= 3:
            print "# received"
            print string.join(lines)

        return self.ExtractData( lines, peptide_sequence = peptide_sequence, genomic_sequence = genomic_sequence)
    
    def Write( self ):
        """write a tab separated list of results.
        """
        for x in self.mMatches:
            x.Write()

    def Sort( self, cmp_function = None ):
        """sort predictions.
        """
        self.mMatches.sort( cmp_function )

def iterator_predictions( infile ):
    """a simple iterator over all entries in a file."""
    while 1:
        line = infile.readline()
        if line.startswith("id\t"): continue
        if not line: return
        if line[0] == "#": continue
        e = PredictionParserEntry()
        e.Read( line )
        yield e

##--------------------------------------------------------------------------
class PredictionParserGenewise( PredictionParser ):
    """Parse output from genewise.
    """
    
    def __init__(self):
        PredictionParser.__init__( self )

        self.mParsers = [ self.ParseSummary, self.ParseTranslation, self.ParsePeptide, self.ParseGFF, self.ParseAlb ]

    def GetNew( self ):
        """return object of the same type."""
        return PredictionParserGenewise()

    def ExtractData( self, lines, peptide_sequence, genomic_sequence ):
        """call the various parsers for a chunk of data.
        """

        if len(lines) == 0:
            return 0

        matches = Predictions()
        self.mMatches = matches
        
        # get info about chromosomal location of chunk submitted to
        # genewise.
        if lines[0][:2] == "#>":
            (self.mQueryToken, self.mSbjctToken, self.mSbjctStrand,
             self.mSbjctFrom, self.mSbjctTo) = re.split( "\s+", lines[0][2:-1])[:5]
            self.mSbjctFrom, self.mSbjctTo = map(int, (self.mSbjctFrom, self.mSbjctTo))
            first_index = 1
        else:
            self.mQueryToken, self.mSbjctToken, self.mSbjctStrand, self.mSbjctFrom, self.mSbjctTo = "", "", "", 0, 0
            first_index = 0

        if len(lines) == 1:
            print "## ERROR: no data for %s vs %s" % (self.mQueryToken, self.mSbjctToken)
            return

        
        entry = PredictionParserEntry()
        entry.mRank = 1
    
        parser_counter = 0

        ## call the various parsers
        for index in range(first_index, len(lines)):

            if lines[index][:2] == "//":
                self.mParsers[parser_counter](lines[first_index:index])

                parser_counter += 1
                first_index = index + 1

        if first_index < index:
            entry.mParsers[parser_counter](lines[first_index:index]) 

        ## set the fields according to parsed output
        # fill in summary info for each entry:
        for entry in self.mMatches:
            entry.mQueryToken = self.mQueryToken
            entry.mSbjctToken = self.mSbjctToken            
            entry.mSbjctStrand = self.mSbjctStrand

            if peptide_sequence:
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                alignlib_lite.py_rescoreAlignment( entry.mMapPeptide2Translation, row_seq, col_seq )
        
                entry.mQueryLength = len(peptide_sequence)            
                entry.mPercentIdentity = alignlib_lite.py_calculatePercentIdentity( entry.mMapPeptide2Translation, row_seq, col_seq ) * 100
                entry.mPercentSimilarity = alignlib_lite.py_calculatePercentSimilarity( entry.mMapPeptide2Translation ) * 100
                entry.mQueryCoverage = ( entry.mMapPeptide2Translation.getRowTo() - \
                                         entry.mMapPeptide2Translation.getRowFrom() + 1 ) * 100 /\
                                         entry.mQueryLength

            entry.mAlignmentString = string.join( map( \
                                      lambda x: string.join(map(str, x), " "),
                                      entry.mMapPeptide2Genome), " ")

            if len(entry.mMapPeptide2Genome) == 0:
                print "### PARSING ERROR: empty alignment"
                print str(entry)
                print string.join(lines, "")
                sys.exit(1)
            
            # get regions matching
            # use summary information, as it is correct for genomic data.
            # when -u is set (in contrast to -alb format)
            # genewise starts counting at 1, thus subtract 1 from first position
            # Note: When position starts with a gap, genewise says alignments starts
            # from position 1 while the first aligned residue is position 2
            if entry.mMapPeptide2Genome[0][0] == "G": entry.mQueryFrom += entry.mMapPeptide2Genome[0][1]
            if entry.mMapPeptide2Genome[-1][0] == "G": entry.mQueryTo -= entry.mMapPeptide2Genome[0][1]

            # a sanity check
            if entry.mQueryFrom != entry.mMapPeptide2Translation.getRowFrom():
                print "## PARSING ERROR: wrong start at %s vs %s: %i <> %i" %\
                      (entry.mQueryToken, entry.mSbjctToken,
                       entry.mQueryFrom,
                       entry.mMapPeptide2Translation.getRowFrom() )
                print str(entry)
                print string.join(lines,"")
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                print alignlib_lite.py_AlignmentFormatExplicit( entry.mMapPeptide2Translation, row_seq, col_seq )
                sys.exit(1)
            if entry.mQueryTo != entry.mMapPeptide2Translation.getRowTo() and \
                   entry.mMapPeptide2Genome[-1][0] != "G":
                print "## PARSING ERROR: wrong end at %s vs %s: %i <> %i" %\
                  (entry.mQueryToken, entry.mSbjctToken,
                   entry.mQueryTo,
                   entry.mMapPeptide2Translation.getRowTo() )
                print str(entry)                
                print string.join(lines,"")
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                print alignlib_lite.py_AlignmentFormatExplicit( entry.mMapPeptide2Translation, row_seq, col_seq )
                sys.exit(1)

            (entry.mNIntrons, entry.mNFrameShifts, entry.mNGaps, entry.mNSplits, entry.mNStopCodons, disruptions) = \
                              Genomics.CountGeneFeatures( entry.mSbjctGenomeFrom,
                                                          entry.mMapPeptide2Genome,
                                                          genomic_sequence,
                                                          GlobalBorderStopCodons,
                                                          GlobalStopCodons )

        matches.Sort( lambda x,y : cmp( -x.score, -y.score))
        return matches
    
    def ParseGFF( self, lines ):
        return "gff", None

    def ParseSummary( self, lines ):
        """parse summary entry of genewise."""

        for line in lines[1:]:
            data = re.split( "\s+", line[:-1])
            entry = PredictionParserEntry()
            entry.score = float(data[0])
            entry.mQueryFrom = int(data[2])
            entry.mQueryTo = int(data[3])
            entry.mSbjctGenomeFrom = int(data[5]) - 1
            entry.mSbjctGenomeTo = int(data[6])
            entry.mIndels = int(data[7])
            if data[8] == "":
                entry.mIntrons = 0
            else:
                entry.mIntrons = int(data[8])
            self.mMatches.append( entry )

    def ParseTranslation( self, lines ):
        """read a translation entry from genewise.

        return a tuple (description, sequence)
        """

        if len(lines) < 1:
            raise "not enough lines"
        
        n = 0
        sequence = ""
        for line in lines:
            if line[0] in (">", "#"):
                if sequence:
                    self.mMatches[n].mTranslation = sequence
                    n += 1
                if line[0] != "#":
                    self.mMatches[n].mDescription = line[1:-1]
                else:
                    self.mMatches[n].mDescription = "pseudogene"                    
                sequence = ""
                continue
            elif line[:6] == "Making":
                continue
            sequence += line[:-1]

        self.mMatches[n].mTranslation = sequence
        
    def ParsePeptide( self, lines ):
        """read a peptide entry from genewise.

        return a tuple (description, sequence, modifications)
        """

        self.ParseTranslation(lines)
        
    def ParseAlb( self, lines ):
        """Parse genewise logical block alignment and return aligment.

        :returns: tuple with ``(code, score, ngaps, nframeshifts,
          nintrons, nphase0, nphase1, nphase2,
          alignment of peptide to genome,
          alignment of peptide to translation)``

        """

        map_genewise2code = \
                          { ( "MATCH_STATE", "CODON") : "M",
                            ( "INTRON_STATE", "5SS_PHASE_0"): "5",
                            ( "INTRON_STATE", "5SS_PHASE_1"): "5",                        
                            ( "INTRON_STATE", "5SS_PHASE_2"): "5",
                            ( "INTRON_STATE", "CENTRAL_INTRON"): "I",
                            ( "MATCH_STATE", "3SS_PHASE_0"): "3",
                            ( "MATCH_STATE", "3SS_PHASE_1"): "3",
                            ( "MATCH_STATE", "3SS_PHASE_2"): "3",                        
                            ( "INSERT_STATE", "3SS_PHASE_0") : "3",
                            ( "INSERT_STATE", "3SS_PHASE_1") : "3",
                            ( "INSERT_STATE", "3SS_PHASE_2") : "3",                                                
                            ( "DELETE_STATE", "INSERT") : "G",
                            ( "INSERT_STATE", "CODON") : "G",
                            ( "INSERT_STATE", "SEQUENCE_DELETION") : "F", # a codon with less than 3 nucleotides at an insert
                            ( "MATCH_STATE", "SEQUENCE_INSERTION") : "F", # a codon with more than 3 nucleotides
                            ( "MATCH_STATE", "SEQUENCE_DELETION")  : "F", # a codon with less than 3 nucleotides
                            ( "LOOP_STATE", "RANDOM_SEQUENCE") : None,
                            ( "END", "END") : None,
                            }
        
        raise "might be broken with new alignlib_lite.py_- check"

        ## note: all the numbers including the ranges can be negative (-1 for first position) !
        rx = re.compile( "^(\S+)\s+\[(\S+):(\S+)\s+\"(\S+)\"\s+(\S+)\],\[(\S+):(\S+)\s+\"(\S+)\"\s+(\S+)\]" )
        
        query_last_state = None
        query_first_from = 0

        sbjct_last_state = None
        sbjct_genomic_residue_first_from = 0
        sbjct_peptide_residue_first_from = 0

        states = []
        total_bin_score = 0.0

        nintrons = 0
        nframeshifts = 0
        nphases = [0,0,0]
        ngaps = 0

        ## parsing errors
        nerrors = 0

        for line in lines:

            try:
                (bin_score,
                 query_from, query_to, query_state, query_score,
                 sbjct_from, sbjct_to, sbjct_state, sbjct_score ) =  rx.match( line[:-1] ).groups()
            except AttributeError:
                raise "parsing error in line %s" % line

            query_from, query_to, sbjct_from, sbjct_to = map(int, (query_from, query_to, sbjct_from, sbjct_to))

            total_bin_score += string.atof( bin_score )

            if query_last_state != query_state or sbjct_last_state != sbjct_state:

                if query_last_state:
                    states.append( ( query_last_state, query_first_from, query_from,
                                     sbjct_last_state, sbjct_first_from, sbjct_from ) )

                
                query_first_from = query_from
                sbjct_first_from = sbjct_from
                query_last_state = query_state
                sbjct_last_state = sbjct_state

        states.append( ( query_last_state, query_first_from, query_to,
                         sbjct_last_state, sbjct_first_from, sbjct_to) )

        ## translate genewise alignment to exonerate alignment codes
        n = 0
        
        ## residue positions on the query/sbjct for peptide sequences
        ## note: genewise would start the first state at -1:0. Thus, add
        ## 1 to get first aligned position in peptide
        query_peptide_from = states[0][1] + 1
        sbjct_peptide_from = 0

        for x in range( len(states)):

            query_state, query_from, query_to, sbjct_state, sbjct_from, sbjct_to = states[x]

            query_len = query_to - query_from
            sbjct_len = sbjct_to - sbjct_from
            
            ## build alignment between peptide sequence and genomic sequence
            try:
                code = map_genewise2code[ (query_state, sbjct_state)]
            except KeyError:
                print "## ERROR: unknown genewise code: %s %s" % (query_state, sbjct_state)
                nerrors += 1
                code = "?"

            if query_state == "LOOP_STATE" and x > 0:
                n += 1
                query_peptide_from = states[x+1][1] + 1
                sbjct_peptide_from = 0
                continue

            ## add code for split codons
            ## move matches to second split codon, if
            ## a split codon is on a gap: add a Gap
            if sbjct_state == "5SS_PHASE_1":
                self.mMatches[n].mMapPeptide2Genome.append( ("S", 0, 1) )
                sbjct_len -= 1
            elif sbjct_state == "5SS_PHASE_2":
                self.mMatches[n].mMapPeptide2Genome.append( ("S", 0, 2) )
                sbjct_len -= 2
            elif sbjct_state == "3SS_PHASE_1":
                query_len = 0
                sbjct_len -= 2
            elif sbjct_state == "3SS_PHASE_2":
                query_len = 0                
                sbjct_len -= 1

            if code:
                if query_from != query_peptide_from - 1:
                    print "Mismatch", query_state, query_from, query_to, sbjct_state,
                    sbjct_from, sbjct_to, query_peptide_from
                
                self.mMatches[n].mMapPeptide2Genome.append( (code, query_len, sbjct_len) )

            if sbjct_state == "3SS_PHASE_1":
                self.mMatches[n].mMapPeptide2Genome.append( ("S", 1, 2) )
                sbjct_len = 2
                if query_state == "MATCH_STATE":
                    query_len = 1
                else:
                    query_len = 0
                    
            elif sbjct_state == "3SS_PHASE_2":
                self.mMatches[n].mMapPeptide2Genome.append( ("S", 1, 1) )                
                sbjct_len = 1
                if query_state == "MATCH_STATE":
                    query_len = 1
                else:
                    query_len = 0
                
            ## build alignment between peptide sequences
            query_increment = 0
            sbjct_increment = 0

            if query_state == "MATCH_STATE":
                 query_increment = query_len

            elif query_state == "DELETE_STATE":
                query_increment = query_len

            elif query_state == "INSERT_STATE":
                query_increment = query_len

            elif query_state == "UNKNOWN_LABEL":
                query_peptide_from += query_len
                ## there are three X's in the sbjct
                sbjct_peptide_from += 3             
                continue
            elif query_state == "LOOP_STATE":
                query_increment += query_len

            if sbjct_state == "CODON":
                sbjct_increment = sbjct_len 
            elif sbjct_state in ("3SS_PHASE_1", "3SS_PHASE_2"):
                sbjct_increment = 3
            elif sbjct_state == "INSERT":
                sbjct_increment = sbjct_len
            elif sbjct_state == "RANDOM_SEQUENCE":
                sbjct_increment = 0

            ## for a frameshift: genewise outputs two characters in the
            ## translated sequence? I do not know. 3 makes more sense.
            elif sbjct_state == "SEQUENCE_DELETION":
                sbjct_increment = 3

            ## for a frameshift: if there is a sequence insertion, no
            ## character is output.
            elif sbjct_state == "SEQUENCE_INSERTION":
                if sbjct_len > 3:
                    sbjct_increment = 3
                else:
                    sbjct_increment = 0

            if code and query_increment and sbjct_increment:
                alignlib_lite.py_addDiagonal2Alignment(  self.mMatches[n].mMapPeptide2Translation,
                                                 query_peptide_from, query_peptide_from + query_increment,
                                                 sbjct_peptide_from - query_peptide_from
                                                 )
                
            query_peptide_from += query_increment
            sbjct_peptide_from += sbjct_increment / 3

    
            
##--------------------------------------------------------------------------
class PredictionParserExonerate (PredictionParser):
    """Parse output from exonerate.
    """
    
    def __init__(self):
        PredictionParser.__init__( self )

        self.mRestrictToStrand = None

    def SetRestrictToStrand( self, v ):
        self.mRestrictToStrand = v
            
    def ExtractData( self, lines, peptide_sequence = None, genomic_sequence = None):

        matches = Predictions()
        
        nmatches = 0
        for line in lines:
            
            if not re.match("diy", line): continue

            nmatches += 1
            data = string.split( line[:-1], "\t")

            entry = PredictionParserEntry()
            entry.mRank = nmatches 

            (sugar, entry.mQueryLength, rank, percent_identity, percent_similarity, vulgar ) = data[1:]

            (entry.mQueryToken, entry.mQueryFrom, entry.mQueryTo, entry.mQueryStrand, 
             entry.mSbjctToken, entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
             entry.mSbjctStrand, entry.score) = string.split(sugar, " " )

            if self.mRestrictToStrand and entry.mSbjctStrand != self.mRestrictToStrand:
                continue

            try:
                (entry.mQueryFrom, entry.mQueryTo,
                 entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                 entry.score, entry.mQueryLength) = map( \
                    int , (entry.mQueryFrom, entry.mQueryTo,
                           entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                           entry.score, entry.mQueryLength))
                
            except ValueError:
                raise InputError, "error while parsing ints: %s" % " ".join( (entry.mQueryFrom, entry.mQueryTo,
                                                                            entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                                                                            entry.score, entry.mQueryLength) )
            
            try:            
                entry.mPercentIdentity, entry.mPercentSimilarity = map( float,
                                                                        (percent_identity, percent_similarity ))
            except ValueError:
                raise InputError, "error while parsing floats: %s, %s" % (percent_identity, percent_similarity)
            
            entry.mAlignmentString = vulgar
            try:
                entry.mMapPeptide2Genome = Genomics.String2Alignment( vulgar )
            except ValueError:
                raise InputError, "error while parsing alignment %s" % vulgar
            
            entry.mQueryCoverage = (entry.mQueryTo - entry.mQueryFrom) * 100 / entry.mQueryLength

            entry.mMapPeptide2Translation, entry.mTranslation = Genomics.Alignment2PeptideAlignment( \
                entry.mMapPeptide2Genome, entry.mQueryFrom, entry.mSbjctGenomeFrom, genomic_sequence )

            ## add 1 for peptide alignments
            entry.mQueryFrom += 1

            (entry.mNIntrons, entry.mNFrameShifts, entry.mNGaps, entry.mNSplits, entry.mNStopCodons, disruptions) = \
                              Genomics.CountGeneFeatures( entry.mSbjctGenomeFrom,
                                                          entry.mMapPeptide2Genome,
                                                          genomic_sequence,
                                                          GlobalBorderStopCodons,
                                                          GlobalStopCodons )

            matches.append(entry)
            
        return matches

##--------------------------------------------------------------------------
def Blocks2Alignment( query_block_starts,
                      sbjct_block_starts,
                      block_sizes ):
    """convert a list of blocks to a peptid2genome alignment.
    """

    alignment = []
    sbjct_pos = 0
    sbjct_offset = sbjct_block_starts[0]
    current_phase = 0

    map_peptide2translation = alignlib_lite.py_makeAlignmentVector()

    query_peptide_pos = query_block_starts[0]
    sbjct_peptide_pos = 0

    for x in range( len(block_sizes) - 1):

        lquery = block_sizes[x]
        lsbjct = block_sizes[x]
        lgap = query_block_starts[x+1] - ( query_block_starts[x] + lquery )
        lintron = sbjct_block_starts[x+1] - ( sbjct_block_starts[x] + lsbjct * 3 )
        
        ## do peptide2translation alignment
        alignlib_lite.py_addDiagonal2Alignment( map_peptide2translation,
                                        query_peptide_pos,
                                        query_peptide_pos + lquery,
                                        sbjct_peptide_pos - query_peptide_pos )
        query_peptide_pos += lquery + lgap
        sbjct_peptide_pos += lsbjct 
        
        
        ## do peptide2genome alignment
        lsbjct *= 3
        sbjct_pos += lsbjct
        phase = sbjct_pos % 3
        if phase: lquery -= 1
        
        alignment.append( ("M", lquery, lsbjct) )
        if phase: alignment.append( ("S", 0, phase) )
        alignment.append( ("5", 0, 0) )
        alignment.append( ("I", 0, lintron) )
        alignment.append( ("3", 0, 0) )
        if phase: alignment.append( ("S", 1, 3-phase) )
        if lgap: alignment.append( ("G", lgap, 0) )

    lquery = block_sizes[-1]
    lsbjct = block_sizes[-1]
    alignlib_lite.py_addDiagonal2Alignment( map_peptide2translation,
                                    query_peptide_pos,
                                    query_peptide_pos + lquery,
                                    sbjct_peptide_pos - query_peptide_pos )


    lsbjct *= 3
    alignment.append( ("M", lquery, lsbjct) )
    return alignment, map_peptide2translation

##--------------------------------------------------------------------------
def Blocks2AlignmentCDNA( query_block_starts,
                          sbjct_block_starts,
                          block_sizes ):
    """convert a list of blocks to a peptid2genome alignment.

    Both blocks are in nucleotides.

    Just get the exon structure correct.
    
    """
        
##     for x in range( len(query_block_starts)):
##         print query_block_starts[x], sbjct_block_starts[x], block_sizes[x]

    alignment = []
    sbjct_pos = 0
    sbjct_offset = sbjct_block_starts[0]
    query_pos = 0
    current_phase = 0

    map_peptide2translation = alignlib_lite.py_makeAlignmentVector()

    query_peptide_pos = query_block_starts[0]/3
    sbjct_peptide_pos = 0

    query_overhang = 0
    ## phase of last block
    phase = 0
    matches = 0
    splits = 0
    intron = 0
    for x in range( len(block_sizes) - 1):

##         print "x=", x, query_block_starts[x], sbjct_block_starts[x], block_sizes[x],\
##               query_block_starts[x+1] - ( query_block_starts[x] + block_sizes[x] ), \
##               sbjct_block_starts[x+1] - ( sbjct_block_starts[x] + block_sizes[x] )

        # count in nucleotides
        nbases  = block_sizes[x] 
        lgap    = query_block_starts[x+1] - ( query_block_starts[x] + nbases )
        lintron = sbjct_block_starts[x+1] - ( sbjct_block_starts[x] + nbases )

        ## add slipped residues at the ends
        ## add to query, make gap in sbjct smaller.
        ## gaps are thus artificially reduced to multiples of 3
##         nbases  += lgap
##         lsbjct  += lgap
##         lintron -= lgap
##         lgap    -= lgap 

        naminos = (nbases - phase) / 3

        sbjct_pos += nbases
        query_pos += nbases 

        phase = query_pos % 3        

##         print "lgap=", lgap, "lintron=",lintron, "phase=", phase,\
##               "qp=", query_peptide_pos, "sp=", sbjct_peptide_pos, "phase=", phase,\
##               "nbases=",nbases, "naminos=",naminos

        alignment.append( ("M", naminos, naminos * 3) )
        matches += naminos * 3

        # treat gaps in query
        # make sure that gaps only occur at codon positions!
        # if there is overhang (phase > 0) at one more gap.
        # and move position to gap in sbjct
        # new_phase is pseudo-phase, makes sure, that extra
        # nucleotides are taken off from next region.
        new_phase = 0
        if lgap:
            if phase:
                query_pos -= phase
                lgap += phase
                lintron += phase
            # increasing lgap to next codon end
            # remembering to take off residues
            if lgap % 3:
                new_phase = 3 - lgap % 3
                lgap += new_phase
                lintron += new_phase 

            query_pos += lgap
            alignment.append( ("G", lgap / 3, 0) )
            phase = query_pos % 3
            if phase != 0: raise "what?"
            
        if lintron:

            if phase == 0 and lintron % 3 == 0 and lintron < 50:
                alignment.append( ("G", 0, lintron ) )
                intron += lintron
            else:
                if phase: alignment.append( ("S", 0, phase) )
                alignment.append( ("5", 0, 0) )
                alignment.append( ("I", 0, lintron) )
                intron += lintron
                alignment.append( ("3", 0, 0) )

                if phase: alignment.append( ("S", 1, 3-phase) )
                if phase: phase = 3 - phase
                if phase: splits += 1

        if new_phase:
            phase = new_phase
        
        # do peptide2translation alignment
        alignlib_lite.py_addDiagonal2Alignment( map_peptide2translation,
                                        query_peptide_pos,
                                        query_peptide_pos + naminos,
                                        sbjct_peptide_pos - query_peptide_pos)
        
        query_peptide_pos += naminos + lgap / 3
        sbjct_peptide_pos += naminos

        if phase:
            map_peptide2translation.add( query_peptide_pos,
                                         sbjct_peptide_pos,
                                         0 )
            query_peptide_pos += 1
            sbjct_peptide_pos += 1
        
    nbases = block_sizes[-1]
    
    naminos = (nbases - phase) / 3

    alignlib_lite.py_addDiagonal2Alignment( map_peptide2translation,
                                    query_peptide_pos,
                                    query_peptide_pos + naminos,
                                    sbjct_peptide_pos - query_peptide_pos )


    sbjct_pos += nbases
    query_pos += nbases 

    naminos = (nbases-phase)/ 3
    phase = query_pos % 3

    alignment.append( ("M", naminos, naminos * 3 ) )
    matches += naminos * 3

##     print "intron=", intron    
##     print "final phase=" , phase
##     print "splits=", splits
##     print "matches=", matches, matches / 3
##     print "query_pos=", query_pos
##     print "sum=", splits * 3 + intron + matches
    return alignment, map_peptide2translation

##--------------------------------------------------------------------------
class PredictionParserBlatCDNA (PredictionParser):
    """Parse output from psl.

    Sample line::

      psLayout version 3

      match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes
              qStarts  tStarts
              match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
      :: ---------------------------------------------------------------------------------------------------------------------------------------------------------------
      885     0       0       0       0       0       3       85      +       ENSCAFT00000000002      885     0       885     chr1.fa 124897793       3228799 3229769 4       41,577,103,164,       0,41,618,721,   3228799,3228845,3229427,3229605,        aactactttggattcacccacagtggggcggcggcagcagc,ggctgcggcccagtatagccagcatccagcttctggtgtagcctactctcatccaactacagttgctagctacgctatccatcaggctccagtagctgctcacacagttactgcggcccatgcaccagcaggcaccacagttgcagttgccaggcctgccccagtagctgttgcagctgctgcaatagctgctgcttatggaggctaccccactggacatacagcaagtgactatggctctatccagagacaacaagaagcaccaccaccagcacccccagctactacacagaactaccaggactcatactcatatgtaagatccactgctcctgctgtagcatatgatagtaagcagctgcctgcccagcctcagccttctgttgctgaaacctactatcagattgcccccaaagcaggttatagccaaggtgcaactcagtatacacaagcccagcaaactcgacaagtgacagccataaaaccagccacaccaagtccagctaccactactttctccatttatcccgtatcttccatcgttcagtcagtagcagttgcagctactgtggtgcca,tatacccagagtgctacttatagtaccatggcagttacttattctggtacatcttattcaggttatgaagcagcagtatattcagctgcatccttctactacc,cagctgcctggacagggaccatctttactaaaaaagcaccattccaaaataaacaactgaaaccaaaacagcctcccaaaccggcccagatagattatcgtgatgtttgtaagattagctgtgctgaaccacagacttataaagaacatttagaaggacaaaaa,     aactactttggattcacccacagtggggcggcggcagcagc,ggctgcggcccagtatagccagcatccagcttctggtgtagcctactctcatccaactacagttgctagctacgctatccatcaggctccagtagctgctcacacagttactgcggcccatgcaccagcaggcaccacagttgcagttgccaggcctgccccagtagctgttgcagctgctgcaatagctgctgcttatggaggctaccccactggacatacagcaagtgactatggctctatccagagacaacaagaagcaccaccaccagcacccccagctactacacagaactaccaggactcatactcatatgtaagatccactgctcctgctgtagcatatgatagtaagcagctgcctgcccagcctcagccttctgttgctgaaacctactatcagattgcccccaaagcaggttatagccaaggtgcaactcagtatacacaagcccagcaaactcgacaagtgacagccataaaaccagccacaccaagtccagctaccactactttctccatttatcccgtatcttccatcgttcagtcagtagcagttgcagctactgtggtgcca,tatacccagagtgctacttatagtaccatggcagttacttattctggtacatcttattcaggttatgaagcagcagtatattcagctgcatccttctactacc,cagctgcctggacagggaccatctttactaaaaaagcaccattccaaaataaacaactgaaaccaaaacagcctcccaaaccggcccagatagattatcgtgatgtttgtaagattagctgtgctgaaccacagacttataaagaacatttagaaggacaaaaa,

    Run is cDNA versus genome.
    """
    
    def __init__(self):
        PredictionParser.__init__( self )
            
    def ExtractData( self, lines, peptide_sequence, genomic_sequence):
        """extract data from psl, tab-formatted output.
        """

        matches = Predictions()
        
        nmatches = 0
        for line in lines:

            if line[0] == "#": continue
            if not re.match("^[0-9]", line): continue

            entry = PredictionParserEntry()
            entry.mRank = nmatches 

            data = string.split(line[:-1], "\t")
            
            ( nmatches, nmismatches, nrepmatches, nns,
              query_ngaps_counts, query_ngaps_bases,
              sbjct_ngaps_counts, sbjct_ngaps_bases,
              sbjct_strand,
              entry.mQueryToken, query_length, query_from, query_to, 
              entry.mSbjctToken, sbjct_length, sbjct_from, sbjct_to, 
              nblocks, block_sizes,
              query_block_starts, sbjct_block_starts,
              query_block_seqs, sbjct_block_seqs) = data

            ( nmatches, nmismatches, nrepmatches, nns,
              query_ngaps_counts, query_ngaps_bases,
              sbjct_ngaps_counts, sbjct_ngaps_bases,
              query_length, query_from, query_to, 
              sbjct_length, sbjct_from, sbjct_to,
              nblocks ) = map(int, \
                              ( nmatches, nmismatches, nrepmatches, nns,
                                query_ngaps_counts, query_ngaps_bases,
                                sbjct_ngaps_counts, sbjct_ngaps_bases,
                                query_length, query_from, query_to, 
                                sbjct_length, sbjct_from, sbjct_to,
                                nblocks ))


            ## compile peptide2genome alignment
            query_block_starts = map(int, string.split( query_block_starts[:-1], ","))
            sbjct_block_starts = map(int, string.split( sbjct_block_starts[:-1], ","))
            query_block_seqs = string.split( query_block_seqs[:-1], ",")
            sbjct_block_seqs = string.split( sbjct_block_seqs[:-1], ",")            
            block_sizes = map(int, string.split( block_sizes[:-1], ","))
            
            if sbjct_strand == "-":
                q = []
                s = []
                for x in range(len(query_block_starts)):
                    q.append( query_length - query_block_starts[x] - block_sizes[x])
                    s.append( sbjct_length - sbjct_block_starts[x] - block_sizes[x])
                q.reverse()
                s.reverse()
                query_block_starts = q
                sbjct_block_starts = s
                block_sizes.reverse()
                sbjct_block_seqs.reverse()
                query_block_seqs.reverse()

            # correct for starts not at codon
            if query_from % 3:
                correct_query = 3 - query_from % 3
                query_block_starts[0] += correct_query
                sbjct_block_starts[0] += correct_query
                block_sizes[0] -= correct_query
                query_from += correct_query
                sbjct_from += correct_query
                print "## correcte from to %i" % query_from
            # correct for stop not at codon
            if query_to % 3:
                correct_query = query_to % 3
                block_sizes[-1] -= correct_query
                query_to -= correct_query
                sbjct_to -= correct_query
                print "## correcte to to %i" % query_to
            
            sbjct_ali_length = sbjct_to - sbjct_from
            query_ali_length = (query_to - query_from) / 3
                            
            query_from /= 3 + 1
            query_to /= 3
            query_length /= 3

            entry.mNIntrons = nblocks - 1
            entry.mQueryLength = query_length 
            entry.mQueryCoverage = 100 * (query_to - query_from) / query_length
            entry.mSbjctStrand = sbjct_strand
            entry.mNGaps = query_ngaps_counts + sbjct_ngaps_counts
            entry.mQueryFrom = query_from
            entry.mQueryTo = query_to
            
            if sbjct_strand == "+":
                entry.mSbjctGenomeFrom = sbjct_from
                entry.mSbjctGenomeTo = sbjct_to
            else:
                entry.mSbjctGenomeFrom = sbjct_length - sbjct_to
                entry.mSbjctGenomeTo = sbjct_length - sbjct_from
            
            entry.mMapPeptide2Genome, entry.mMapPeptide2Translation = \
                                      Blocks2AlignmentCDNA( query_block_starts, sbjct_block_starts, block_sizes )
            entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )

            entry.mTranslation = string.join( sbjct_block_seqs, "")
            entry.mSbjctFrom = 1
            entry.mSbjctTo = len(entry.mTranslation)
            
            lquery, lsbjct = Genomics.GetAlignmentLength( entry.mMapPeptide2Genome )
            
            if 0:
                print str(entry)
                if 0 and (lquery != (entry.mQueryTo - entry.mQueryFrom + 1) or lquery != query_ali_length):
                    raise AlignmentError( "alignment length discrepancy in query: %i vs %i vs %i" %\
                                      ( lquery, entry.mQueryTo - entry.mQueryFrom + 1, query_ali_length))

                elif lsbjct != (entry.mSbjctGenomeTo - entry.mSbjctGenomeFrom) or lsbjct != sbjct_ali_length:
                    raise AlignmentError( "# alignment length discrepancy in sbjct: %i vs %i vs %i" %\
                                          ( lsbjct, entry.mSbjctGenomeTo - entry.mSbjctGenomeFrom, sbjct_ali_length))

                elif 0 and (entry.mQueryFrom != entry.mMapPeptide2Translation.getRowFrom() or\
                         entry.mQueryTo != entry.mMapPeptide2Translation.getRowTo()):
                    raise AlignmentError( "# alignment length discrepancy in query: %i-%i vs %i-%i" % \
                                          ( entry.mQueryFrom, entry.mQueryTo,
                                            entry.mMapPeptide2Translation.getRowFrom(),
                                            entry.mMapPeptide2Translation.getRowTo()) )

                elif entry.mSbjctFrom != entry.mMapPeptide2Translation.getColFrom() or\
                         entry.mSbjctTo != entry.mMapPeptide2Translation.getColTo():
                    raise AlignmentError( "# alignment length discrepancy in sbjct: %i-%i vs %i-%i" % \
                                          ( entry.mSbjctFrom, entry.mSbjctTo,
                                            entry.mMapPeptide2Translation.getColFrom(),
                                            entry.mMapPeptide2Translation.getColTo()))

            entry.mPercentIdentity = 100 * nmatches / (nmatches + nmismatches)
            entry.mPercentSimilarity = 100 * nmatches / (nmatches + nmismatches)                

##             peptide_sequence = ["X" * (query_block_starts[0]) + query_block_seqs[0]]
##             p = query_block_starts[0] + block_sizes[0]
            
##             for x in range(1,len(query_block_starts)):
##                 peptide_sequence.append("X" * (query_block_starts[x] - p) + query_block_seqs[x])
##                 p = query_block_starts[x] + block_sizes[x]
                
##             peptide_sequence = string.join( peptide_sequence, "")
##             row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
##             col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
##             alignlib_lite.py_rescoreAlignment( entry.mMapPeptide2Translation, row_seq, col_seq )
        
##             entry.mPercentIdentity = alignlib_lite.py_calculatePercentIdentity( entry.mMapPeptide2Translation, row_seq, col_seq ) * 100
##             entry.mPercentSimilarity = alignlib_lite.py_calculatePercentSimilarity( entry.mMapPeptide2Translation ) * 100
        
            matches.append( entry )

        return matches

    
class PredictionParserBlatTrans (PredictionParser):
    """Parse output from psl.

    Query is protein, database is 6frame translation
    
    Input looks like this::

        psLayout version 3

        match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes
                qStarts  tStarts
                match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
        :: ---------------------------------------------------------------------------------------------------------------------------------------------------------------
        98      0       0       0       0       0       2       3113    +-      185061  98      0       98      chr19.fa        56914383        6068186 6071593 3       28,27,43,        0,28,55,        50842790,50845539,50846068,     IFSCRQNCVEFYPIFLVTLWMAGWYFNQ,VFATCLGLVYIYARHQYFWGYSEAAKK,RITGFRLSLGCLALLTVLGALGIANSFLDEYLDLNVIKKLRHF,   IFSCRQNCVEFYPIFLVTLWMAGWYFNQ,VFATCLGLVYIYARHQYFWGYSEAAKK,RITGFRLSLGCLALLTVLGALGIANSFLDEYLDLNVIKKLRHF,
    """
    
    def __init__(self):
        PredictionParser.__init__( self )
            
    def ExtractData( self, lines, peptide_sequence, genomic_sequence):
        """extract data from psl, tab-formatted output.
        """

        matches = Predictions()
        
        nmatches = 0
        for line in lines:

            if line[0] == "#": continue
            if not re.match("^[0-9]", line): continue
            
            entry = PredictionParserEntry()
            entry.mRank = nmatches 

            data = string.split(line[:-1], "\t")
            
            ( nmatches, nmismatches, nrepmatches, nns,
              query_ngaps_counts, query_ngaps_bases,
              sbjct_ngaps_counts, sbjct_ngaps_bases,
              sbjct_strand,
              entry.mQueryToken, query_length, query_from, query_to, 
              entry.mSbjctToken, sbjct_length, sbjct_from, sbjct_to, 
              nblocks, block_sizes,
              query_block_starts, sbjct_block_starts,
              query_block_seqs, sbjct_block_seqs) = data

            ( nmatches, nmismatches, nrepmatches, nns,
              query_ngaps_counts, query_ngaps_bases,
              sbjct_ngaps_counts, sbjct_ngaps_bases,
              query_length, query_from, query_to, 
              sbjct_length, sbjct_from, sbjct_to,
              nblocks ) = map(int, \
                              ( nmatches, nmismatches, nrepmatches, nns,
                                query_ngaps_counts, query_ngaps_bases,
                                sbjct_ngaps_counts, sbjct_ngaps_bases,
                                query_length, query_from, query_to, 
                                sbjct_length, sbjct_from, sbjct_to,
                                nblocks ))

            entry.score = nmatches
            entry.mNIntrons = nblocks - 1
            entry.mQueryLength = query_length 
            entry.mQueryCoverage = 100 * (query_to - query_from) / query_length
            ## take sbjct strand
            entry.mSbjctStrand = sbjct_strand[1]
            entry.mNGaps = query_ngaps_counts + sbjct_ngaps_counts
            entry.mQueryFrom = query_from + 1
            entry.mQueryTo = query_to

            ## reverse strand according to our coordinates
            if entry.mSbjctStrand == "+":
                entry.mSbjctGenomeFrom = sbjct_from
                entry.mSbjctGenomeTo = sbjct_to
            else:
                entry.mSbjctGenomeFrom = sbjct_length - sbjct_to
                entry.mSbjctGenomeTo = sbjct_length - sbjct_from


            ## compile peptide2genome alignment
            query_block_starts = map(int, string.split( query_block_starts[:-1], ","))
            sbjct_block_starts = map(int, string.split( sbjct_block_starts[:-1], ","))
            query_block_seqs = string.split( query_block_seqs[:-1], ",")
            sbjct_block_seqs = string.split( sbjct_block_seqs[:-1], ",")            
            block_sizes = map(int, string.split( block_sizes[:-1], ","))

            entry.mMapPeptide2Genome, entry.mMapPeptide2Translation = \
                                      Blocks2Alignment( query_block_starts, sbjct_block_starts, block_sizes )
            entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )

            entry.mTranslation = string.join( sbjct_block_seqs, "")
            entry.mSbjctFrom = 1
            entry.mSbjctTo = len(entry.mTranslation)

            lquery, lsbjct = Genomics.GetAlignmentLength( entry.mMapPeptide2Genome )
            error = None
            if lquery != (entry.mQueryTo - entry.mQueryFrom + 1):
                raise AlignmentError( "alignment length discrepancy in query: %i vs %i vs %i" %\
                                  ( query_length,  lquery, entry.mQueryTo - entry.mQueryFrom + 1))

            elif lsbjct != (entry.mSbjctGenomeTo - entry.mSbjctGenomeFrom):
                raise AlignmentError( "# alignment length discrepancy in sbjct: %i vs %i vs %i" %\
                                      ( sbjct_length, lsbjct, entry.mSbjctGenomeTo - entry.mSbjctGenomeFrom))
                                                                                                 
                                                                                                 
            elif entry.mQueryFrom != entry.mMapPeptide2Translation.getRowFrom() or\
                     entry.mQueryTo != entry.mMapPeptide2Translation.getRowTo():
                raise AlignmentError( "# alignment length discrepancy in query: %i-%i vs %i-%i" % \
                                      ( entry.mQueryFrom, entry.mQueryTo,
                                        entry.mMapPeptide2Translation.getRowFrom(),
                                        entry.mMapPeptide2Translation.getRowTo()) )
            elif entry.mSbjctFrom != entry.mMapPeptide2Translation.getColFrom() or\
                     entry.mSbjctTo != entry.mMapPeptide2Translation.getColTo():
                raise AlignmentError( "# alignment length discrepancy in sbjct: %i-%i vs %i-%i" % \
                                      ( entry.mSbjctFrom, entry.mSbjctTo,
                                        entry.mMapPeptide2Translation.getColFrom(),
                                        entry.mMapPeptide2Translation.getColTo()))

            peptide_sequence = ["X" * (query_block_starts[0]) + query_block_seqs[0]]
            p = query_block_starts[0] + block_sizes[0]
            
            for x in range(1,len(query_block_starts)):
                peptide_sequence.append("X" * (query_block_starts[x] - p) + query_block_seqs[x])
                p = query_block_starts[x] + block_sizes[x]
                
            peptide_sequence = string.join( peptide_sequence, "")
            row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
            col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
            alignlib_lite.py_rescoreAlignment( entry.mMapPeptide2Translation, row_seq, col_seq )
        
            entry.mPercentIdentity = alignlib_lite.py_calculatePercentIdentity( entry.mMapPeptide2Translation, row_seq, col_seq ) * 100
            entry.mPercentSimilarity = alignlib_lite.py_calculatePercentSimilarity( entry.mMapPeptide2Translation ) * 100
            
            matches.append( entry )

        return matches

class PredictionParserExons (PredictionParser):
    """Parse output from exons file.
    
    Input looks like this::

       CG9926-RA       chr3R   -1      0       1       0       696     9776126 9776819
       CG9927-RA       chr3R   -1      0       1       0       398     9771125 9771523
       CG9927-RA       chr3R   -1      1       2       398     1026    9770442 9771067
       CG9928-RA       chr2L   -1      0       1       0       321     13131417        13131735
       CG9929-RA       chr3R   -1      0       1       0       117     9769026 9769143
       CG9929-RA       chr3R   -1      0       2       117     347     9768742 9768972
       CG9929-RA       chr3R   -1      1       3       347     942     9768091 9768683
       CG9930-RA       chr3R   -1      0       1       0       880     9699040 9699920
       CG9930-RA       chr3R   -1      2       2       880     1038    9694339 9694497
       CG9930-RA       chr3R   -1      0       3       1038    1575    9693241 9693775

    Negative strand exons are given as positive strand coordinates, but are still sorted
    in increasing order.

    Coordinates are zero-based in open/closed notation.

    The score is set to the length of the transcript.
    """
    
    def __init__(self, contig_sizes = {}):
        PredictionParser.__init__( self )
        self.contigSizes = contig_sizes
        
    def ExtractData( self,
                     lines, peptide_sequence,
                     genomic_sequence ):
        """extract data from exons file, tab-formatted output.
        """

        matches = Predictions()
        
        nmatches = 0
        
        entry = None

        id = 0 
        last_query = None
        exons = []
        for line in lines:

            data = line[:-1].split()

            (query, sbjct_token, sbjct_strand,
             phase, exon_id, peptide_from, peptide_to,
             sbjct_genome_from, sbjct_genome_to ) = data

            ( sbjct_genome_from, sbjct_genome_to,
              peptide_from, peptide_to, phase ) = \
              map(int, ( sbjct_genome_from, sbjct_genome_to,
                         peptide_from, peptide_to, phase ) )

            ## sanity check
            if sbjct_genome_to <= sbjct_genome_from:
                ## Ensembl had small exons of wrong orientation.
                ## Those have been removed in later releases.
                ## Simply ignore them here.
                continue
            
            if sbjct_strand in ("-1", "0", "-"):
                sbjct_strand = "-"
            elif sbjct_strand in ( "1", "+", "+1") :
                sbjct_strand = "+"                

            if sbjct_token in self.contigSizes:
                sbjct_length = self.contigSizes[sbjct_token]
                if sbjct_strand == "-":
                    sbjct_genome_from, sbjct_genome_to = sbjct_length - sbjct_genome_to, sbjct_length - sbjct_genome_from
                    
            if last_query != query:
                
                if entry:
                    self.mNumInput += 1
                    try:
                        entry.mMapPeptide2Genome = Genomics.Exons2Alignment( exons )
                    except ValueError, msg:
                        self.mErrors.append( (last_query, msg) )
                    else:
                        self.mNumOutput += 1
                        entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )
                        matches.append( entry )
                    exons = []

                    
                id += 1
                entry = PredictionParserEntry()
                
                entry.mRank = 0 
                entry.mQueryToken = query
                entry.mPredictionId = query

                entry.mSbjctToken = sbjct_token
                entry.mSbjctStrand = sbjct_strand
                entry.mSbjctGenomeFrom = sbjct_genome_from
                entry.mSbjctGenomeTo = sbjct_genome_to

                entry.mPercentIdentity = 100
                entry.mPercentSimilarity = 100                
                
                entry.mQueryCoverage = 100                

                last_query = query
                
            else:
                entry.mSbjctGenomeFrom = min(entry.mSbjctGenomeFrom, sbjct_genome_from)
                entry.mSbjctGenomeTo = max(entry.mSbjctGenomeTo, sbjct_genome_to)
                
            exons.append( (sbjct_genome_from, sbjct_genome_to) )
                    
        if entry:
            self.mNumInput += 1
            try:
                entry.mMapPeptide2Genome = Genomics.Exons2Alignment( exons )
            except ValueError, msg:
                self.mErrors.append( (last_query, msg) )
            else:
                entry.score = entry.mSbjctTo - entry.mSbjctFrom                
                entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )
                matches.append( entry )
                self.mNumOutput += 1
                
        return matches

class PredictionParserGFF(PredictionParser):
    """Parse GFF files.
    """
    
    def __init__(self, contig_sizes = {}):
        PredictionParser.__init__( self )
        self.contigSizes = contig_sizes
        
    def ExtractData( self,
                     lines, peptide_sequence,
                     genomic_sequence ):
        """extract data from psl, tab-formatted output.
        """
        raise "this module is incomplete"
    
        matches = Predictions()

        for line in lines:
            if line[0] == "#": continue
        
            gff_entry = GTF.Entry()
            gff_entry.read( line )

            entry = PredictionParserEntry()
            entry.mSbjctToken = gff_entry.name
            entry.mSbjctGenomeFrom = gff_entry.start
            entry.mSbjctGenomeTo = gff_entry.end
            
            (query, sbjct_token, sbjct_strand,
             phase, exon_id, peptide_from, peptide_to,
             sbjct_genome_from, sbjct_genome_to ) = data

            ( sbjct_genome_from, sbjct_genome_to,
              peptide_from, peptide_to, phase ) = \
            map(int, ( sbjct_genome_from, sbjct_genome_to,
                       peptide_from, peptide_to, phase ) )
            
            if sbjct_strand in ("-1", "0"):
                sbjct_strand = "-"
            elif sbjct_strand == "1":
                sbjct_strand = "+"                

            if sbjct_token in self.contigSizes:
                sbjct_length = self.contigSizes[sbjct_token]
                if sbjct_strand == "-":
                    sbjct_genome_from, sbjct_genome_to = sbjct_length - sbjct_genome_to, sbjct_length - sbjct_genome_from

            if last_query != query:

                if entry:
                    entry.mMapPeptide2Genome = Genomics.Exons2Alignment( exons )
                    entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )

                    matches.append( entry )
                    exons = []
                    
                id += 1
                entry = PredictionParserEntry()
                
                entry.mRank = 0 
                entry.mQueryToken = query
                entry.mPredictionId = id

                entry.mSbjctToken = sbjct_token
                entry.mSbjctStrand = sbjct_strand
                entry.mSbjctGenomeFrom = sbjct_genome_from
                entry.mSbjctGenomeTo = sbjct_genome_to

                entry.mPercentIdentity = 100
                entry.mPercentSimilarity = 100                
                
                entry.mQueryCoverage = 100                
                
                
                last_query = query
                
            else:
                entry.mSbjctGenomeFrom = min(entry.mSbjctGenomeFrom, sbjct_genome_from)
                entry.mSbjctGenomeTo = max(entry.mSbjctGenomeTo, sbjct_genome_to)

            exons.append( (sbjct_genome_from, sbjct_genome_to) )
                    
        if entry:
            entry.mMapPeptide2Genome = Genomics.Exons2Alignment( exons )
            entry.mAlignmentString = Genomics.Alignment2String( entry.mMapPeptide2Genome )

            matches.append( entry )
            
        return matches


if __name__ == "__main__":

    p = PredictionParserGenewise()
    p.Parse( sys.stdin.readlines() )
    print str(p)
                                                      
