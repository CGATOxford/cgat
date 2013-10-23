################################################################################
#   Gene prediction pipeline 
#
#   $Id: AlignedPairs.py 2781 2009-09-10 11:33:14Z andreas $
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
AlignedPairs.py - 
==================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

"""
import sys, os, string
try: import alignlib_lite
except ImportError: pass

import WrapperDialign
import WrapperDBA
import WrapperBlastZ
import Genomics

# Wrapper for clustal needs to be updated
# import WrapperClustal


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

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
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class UnalignedPair:

    def __init__(self, other = None, sequence1 = "", sequence2 = "", token1 = 0, token2 = 0):

        if other:
            (self.mCategory, self.mMethod,
             self.mToken1, self.mId1, self.mNum1, self.mLen1,
             self.mToken2, self.mId2, self.mNum2, self.mLen2,
             self.mSequence1, self.mSequence2,
             self.mFrom1, self.mTo1,
             self.mFrom2, self.mTo2 ) = \
             (other.mCategory, other.mMethod,
              other.mToken1, other.mId1, other.mNum1, other.mLen1,
              other.mToken2, other.mId2, other.mNum2, other.mLen2,
              other.mSequence1, other.mSequence2,
              other.mFrom1, other.mTo1,
              other.mFrom2, other.mTo2 )
        else:
            self.mCategory = ""
            self.mMethod = "unaligned"
            self.mToken1 = token1
            self.mId1 = 0
            self.mNum1 = 0
            self.mLen1 = 0
            self.mToken2 = token2
            self.mId2 = 0
            self.mNum2 = 0
            self.mLen2 = 0
            self.mSequence1 = ""
            self.mSequence2 = ""        
            self.mFrom1 = 0
            self.mTo1 = 0
            self.mFrom2 = 0
            self.mTo2 = 0

        if sequence1:
            self.mSequence1 = sequence1
            self.mFrom1 = 0
            self.mTo1 = len(sequence1)
            
        if sequence2:
            self.mSequence2 = sequence2        
            self.mFrom2 = 0
            self.mTo2 = len(sequence2)
        
    def Read( self, line ):
        """read data from tab-separated line."""
        
        data = string.split( line[:-1], "\t")
        
        (self.mCategory, self.mMethod,
         self.mToken1, self.mId1, self.mNum1, self.mLen1,
         self.mToken2, self.mId2, self.mNum2, self.mLen2,
         self.mFrom1, self.mTo1, self.mSequence1, 
         self.mFrom2, self.mTo2, self.mSequence2, ) = data[:16]

        ( self.mId1, self.mNum1, self.mLen1,
          self.mId2, self.mNum2, self.mLen2,
          self.mFrom1, self.mTo1,
          self.mFrom2, self.mTo2) = \
          map(int, ( self.mId1, self.mNum1, self.mLen1,
                     self.mId2, self.mNum2, self.mLen2,
                     self.mFrom1, self.mTo1,
                     self.mFrom2, self.mTo2, ))
        
    def __str__( self ):
        return "\t".join( map(str,\
                              ( self.mCategory, self.mMethod,
                                self.mToken1, self.mId1, self.mNum1, self.mLen1,
                                self.mToken2, self.mId2, self.mNum2, self.mLen2,
                                self.mFrom1, self.mTo1, self.mSequence1,
                                self.mFrom2, self.mTo2, self.mSequence2 ) ) )
                              

    def GetHeader( self ):
        return "\t".join( \
            ( "CATEGORY", "METHOD",
              "TOKEN1", "ID1", "TOTAL1", "LEN1",
              "TOKEN2", "ID2", "TOTAL2", "LEN2",               
              "FROM1", "TO1", "SEQ1",
              "FROM2", "TO2", "SEQ2" ) )


    def GetLegend( self ):
        return """# CATEGORY:       category
# METHOD:         alignment method
# TOKEN1:         name
# ID1:            id of segment
# TOTAL1:         number of segments
# LEN1:           length of segment
# TOKEN2:         name
# ID2:            id of segment
# TOTAL2:         number of segments
# LEN2:           length of segment
# FROM1:          start of segment1
# TO1:            end of segment1
# SEQ1:           sequence of segment1
# FROM2:          start of segment1
# TO2:            end of segment1
# SEQ2:           sequence of segment1"""

        
class AlignedPair( UnalignedPair ):

    mOptionsBlastZ = "C=2 B=0 T=0 W=6 K=2200"
    mOptionsDialign = "-n"
    mOptionsDialignLGS = "-n -it -thr 2 -lmax 30 -smin 8"

    def __init__(self, other = None, sequence1="", sequence2=""):
        
        UnalignedPair.__init__( self, other, sequence1, sequence2 )

        if other and isinstance(other, AlignedPair ):
            ( self.mIdentical,
              self.mTransitions,
              self.mTransversions,
              self.mNumGaps,
              self.mLength,
              self.mAligned,
              self.mUnaligned,
              self.mAlignmentFrom1,
              self.mAlignmentTo1,
              self.mAlignment1,
              self.mAlignmentFrom2,
              self.mAlignmentTo2,
              self.mAlignment2,
              self.mUnaligned,
              self.mBlockSizes,
              self.mGapSizes1,
              self.mGapSizes2 ) = \
              ( other.mIdentical,
                other.mTransitions,
                other.mTransversions,
                other.mNumGaps,
                other.mLength,
                other.mAligned,
                other.mUnaligned,
                other.mAlignmentFrom1,
                other.mAlignmentTo1,
                other.mAlignment1,
                other.mAlignmentFrom2,
                other.mAlignmentTo2,
                other.mAlignment2,
                other.mUnaligned,
                other.mBlockSizes,
                other.mGapSizes1,
                other.mGapSizes2 ) 
        else:
            self.mIdentical = 0
            self.mTransitions = 0
            self.mTransversions = 0
            self.mNumGaps = 0
            self.mLength = 0
            self.mAligned = 0        
            self.mUnaligned = 1
            self.mAlignmentFrom1 = 0
            self.mAlignmentTo1 = 0
            self.mAlignment1 = ""
            self.mAlignmentFrom2 = 0
            self.mAlignmentTo2 = 0
            self.mAlignment2 = ""        
            self.mBlockSizes = []
            self.mGapSizes1 = []
            self.mGapSizes2 = []


    ##-----------------------------------------------------------------------------------
    def __str__(self):
        """string representation of object."""
        return UnalignedPair.__str__(self) + "\t" +\
               "\t".join( map( str, \
                               ( self.mIdentical,
                                 self.mTransitions,
                                 self.mTransversions,
                                 self.mNumGaps,
                                 self.mLength,
                                 self.mAligned,
                                 self.mUnaligned,
                                 self.mAlignmentFrom1,
                                 self.mAlignmentTo1,
                                 self.mAlignment1,
                                 self.mAlignmentFrom2,
                                 self.mAlignmentTo2,
                                 self.mAlignment2,
                                 ",".join(map(str,self.mBlockSizes)),
                                 ",".join(map(str,self.mGapSizes1)),
                                 ",".join(map(str,self.mGapSizes2)),
                                 ) ) )

    ##-----------------------------------------------------------------------------------
    def Read( self, line ):
        """read data from tab-separated line."""
        UnalignedPair.Read( self, line )
        
        data = string.split( line[:-1], "\t")

        ( self.mIdentical,
          self.mTransitions,
          self.mTransversions,
          self.mNumGaps,
          self.mLength,
          self.mAligned,
          self.mUnaligned,
          self.mAlignmentFrom1,
          self.mAlignmentTo1,
          self.mAlignment1,
          self.mAlignmentFrom2,
          self.mAlignmentTo2,
          self.mAlignment2,
          self.mBlockSizes,
          self.mGapSizes1,
          self.mGapSizes2) = data[16:]

        ( self.mIdentical,
          self.mTransitions,
          self.mTransversions,
          self.mNumGaps,
          self.mLength,
          self.mAligned,
          self.mUnaligned,
          self.mAlignmentFrom1,
          self.mAlignmentTo1,
          self.mAlignmentFrom2,
          self.mAlignmentTo2 ) = map(int, \
                                     ( self.mIdentical,
                                       self.mTransitions,
                                       self.mTransversions,
                                       self.mNumGaps,
                                       self.mLength,
                                       self.mAligned,
                                       self.mUnaligned,
                                       self.mAlignmentFrom1,
                                       self.mAlignmentTo1,
                                       self.mAlignmentFrom2,
                                       self.mAlignmentTo2 ))

        self.mBlockSizes = map(int, self.mBlockSizes.split(","))
        if self.mGapSizes1:
            self.mGapSizes1 = map(int, self.mGapSizes1.split(","))
        else:
            self.mGapSizes1 = []
        if self.mGapSizes2:
            self.mGapSizes2 = map(int, self.mGapSizes2.split(","))
        else:
            self.mGapSizes2 = []
        
    ##-----------------------------------------------------------------------------------
    def GetHeader( self ):
        return UnalignedPair.GetHeader(self) + "\t" + \
               "\t".join( ( \
            "IDENTICAL",
            "TRANSITIONS",
            "TRANSVERSIONS",
            "NUMGAPS",
            "LENGTH",
            "ALIGNED",
            "UNALIGNED",
            "ALIFROM1", "ALITO1", "ALIGNED1",
            "ALIFROM2", "ALITO2", "ALIGNED2",
            "BLOCKSIZES",
            "GAPSIZES1", "GAPSIZES2" ) )
    
    ##-----------------------------------------------------------------------------------
    def GetLegend( self ):
        return UnalignedPair.GetLegend(self) + "\n" +\
               """# IDENTICAL:       number of identical positions
# TRANSITIONS:    number of transitions
# TRANSVERSIONS:  number of transversion
# NUMGAPS:        number of gaps
# LENGTH:         length of alignment
# ALIGNED:        number of aligned positions
# UNALIGNED:      number of unaligned positions
# ALIFROM1:       start of alignment in sequence1
# ALIEND1:        end of alignment in sequence1
# ALIGNED1:       alignment of sequence1
# ALIFROM2:       start of alignment in sequence2
# ALIEND2:        end of alignment in sequence2
# ALIGNED2:       alignment of sequence2
# BLOCKSIZES:     alignment, length of blocks
# GAPSIZES1:      length of gap sizes in sequence1
# GAPSIZES2:      length of gap sizes in sequence2"""

    ##-----------------------------------------------------------------------------------
    def GetMap( self ):
        """return map between the two segments."""
        if self.mAlignmentFrom1 and self.mAlignmentFrom2:
            map_a2b = alignlib_lite.py_makeAlignmentVector()
            alignlib_lite.py_AlignmentFormatEmissions( 
                self.mAlignmentFrom1, self.mAlignment1,
                self.mAlignmentFrom2, self.mAlignment2 ).copy( map_a2b )
            return map_a2b
        else:
            return None

    ##-----------------------------------------------------------------------------------
    def Align( self, method, anchor = 0, loglevel = 1 ):
        """align a pair of sequences.
        get rid of this and use a method class instead in the future
        """
        
        map_a2b = alignlib_lite.py_makeAlignmentVector()
        s1 = "A" * anchor + self.mSequence1 + "A" * anchor
        s2 = "A" * anchor + self.mSequence2 + "A" * anchor    

        self.strand = "+"

        if method == "dialign":
            dialign = WrapperDialign.Dialign( self.mOptionsDialign )
            dialign.Align( s1, s2, map_a2b )
        elif method == "blastz":
            blastz = WrapperBlastZ.BlastZ( self.mOptionsBlastZ )
            blastz.Align( s1, s2, map_a2b )
            if blastz.isReverseComplement():
                self.strand = "-"
                self.mSequence2 = Genomics.complement( self.mSequence2 )

        elif method == "dialignlgs":
            dialignlgs = WrapperDialign.Dialign( self.mOptionsDialignLGS )
            dialignlgs.Align( s1, s2, map_a2b ) 
        elif method == "dba":
            dba = WrapperDBA.DBA()
            dba.Align( s1, s2, map_a2b )
        elif method == "clustal":
            raise NotImplementedError( "clustal wrapper needs to be updated")
            clustal = WrapperClustal.Clustal()
            clustal.Align( s1, s2, map_a2b )
        elif method == "nw":
            seq1 = alignlib_lite.py_makeSequence( s1 )
            seq2 = alignlib_lite.py_makeSequence( s2 )
            alignator = alignlib_lite.py_makeAlignatorDPFull( alignlib_lite.py_ALIGNMENT_GLOBAL,
                                                      gop=-12.0,
                                                      gep=-2.0 )
            alignator.align( map_a2b, seq1, seq2 )
        elif method == "sw":                        
            seq1 = alignlib_lite.py_makeSequence( s1 )
            seq2 = alignlib_lite.py_makeSequence( s2 )
            alignlib_lite.py_performIterativeAlignment( map_a2b, seq1, seq2, alignator_sw, min_score_sw )
        else:
            ## use callback function
            method(s1, s2, map_a2b)

        if map_a2b.getLength() == 0:
            raise AlignmentError("empty alignment")

        if anchor:
            map_a2b.removeRowRegion( anchor + len(self.mSequence1) + 1, map_a2b.getRowTo() )
            map_a2b.removeRowRegion( 1, anchor)        
            map_a2b.removeColRegion( anchor + len(self.mSequence2) + 1, map_a2b.getColTo() )        
            map_a2b.removeColRegion( 1, anchor)
            map_a2b.moveAlignment( -anchor, -anchor )

        f = alignlib_lite.py_AlignmentFormatExplicit( map_a2b, 
                                              alignlib_lite.py_makeSequence( self.mSequence1),
                                              alignlib_lite.py_makeSequence( self.mSequence2) )

        self.mMethod = method
        self.mAlignment = map_a2b
        self.mAlignedSequence1, self.mAlignedSequence2 = f.mRowAlignment, f.mColAlignment
        f = alignlib_lite.py_AlignmentFormatEmissions( map_a2b )
        self.mAlignment1, self.mAlignment2 = f.mRowAlignment, f.mColAlignment
        self.mAlignmentFrom1 = map_a2b.getRowFrom()
        self.mAlignmentTo1 = map_a2b.getRowTo()        
        self.mAlignmentFrom2 = map_a2b.getColFrom()
        self.mAlignmentTo2 = map_a2b.getColTo()        
        self.mNumGaps, self.mLength = map_a2b.getNumGaps(), map_a2b.getLength()
        self.mAligned = self.mLength - self.mNumGaps

        self.SetPercentIdentity()
        self.SetBlockSizes()

    ##-----------------------------------------------------------------------------------
    def SetBlockSizes( self, gap_char = "-" ):
        """return the block-sizes in the alignment.
        """
        block_sizes = []
        left_gap_sizes = []
        right_gap_sizes = []    

        was_left_gap = self.mAlignedSequence1[0] == gap_char
        was_block = self.mAlignedSequence1[0] != gap_char and self.mAlignedSequence2[0] != gap_char
        size = 1

        for x in range(1,min(len(self.mAlignedSequence1), len(self.mAlignedSequence2))):
            is_left_gap  = self.mAlignedSequence1[x] == gap_char
            is_right_gap = self.mAlignedSequence2[x] == gap_char
            is_block = not is_left_gap and not is_right_gap
            if is_block and not was_block:
                if was_left_gap:
                    left_gap_sizes.append(size)
                    size = 0
                else:
                    right_gap_sizes.append(size)
                    size = 0
            elif not is_block and was_block:
                block_sizes.append(size)
                size = 0
            elif is_left_gap and not was_left_gap:
                right_gap_sizes.append(size)
                size = 0
            elif is_right_gap and was_left_gap:
                left_gap_sizes.append(size)
                size = 0
            elif is_left_gap and is_right_gap:
                raise "double gap"

            was_block = is_block
            was_left_gap = is_left_gap

            size += 1

        if was_block:
            block_sizes.append(size)
        else:
            if was_left_gap:
                left_gap_sizes.append(size)
            else:
                right_gap_sizes.append(size)

        self.mBlockSizes, self.mGapSizes1, self.mGapSizes2 = block_sizes, left_gap_sizes, right_gap_sizes

    ##-----------------------------------------------------------------------------------
    def SetPercentIdentity( self, gap_char = "-" ):
        """return number of idential and transitions/transversions substitutions
        in the alignment.
        """
        transitions   = ("AG", "GA", "CT", "TC")
        transversions = ("AT", "TA", "GT", "TG", "GC", "CG", "AC", "CA" )    

        nidentical = 0
        naligned = 0
        ndifferent = 0
        ntransitions = 0
        ntransversions = 0
        nunaligned = (self.mAlignmentFrom1 - self.mFrom1) +\
                     (self.mTo1 - self.mAlignmentTo1 ) +\
                     (self.mAlignmentFrom2 - self.mFrom2) +\
                     (self.mTo2 - self.mAlignmentTo2 )                      
        
        for x in range(min(len(self.mAlignedSequence1), len(self.mAlignedSequence2))):
            if self.mAlignedSequence1[x] != gap_char and \
               self.mAlignedSequence2[x] != gap_char:
                naligned += 1
                if self.mAlignedSequence1[x] == self.mAlignedSequence2[x]:
                    nidentical += 1
                else:
                    ndifferent += 1

                if (self.mAlignedSequence1[x] + self.mAlignedSequence2[x]) in transitions:
                    ntransitions += 1
                if (self.mAlignedSequence1[x] + self.mAlignedSequence2[x]) in transversions:
                    ntransversions += 1
            else:
                nunaligned += 1

        self.mIdentical, self.mTransitions, self.mTransversions, self.mAligned, self.mUnaligned = \
                         nidentical, ntransitions, ntransversions, naligned, nunaligned





