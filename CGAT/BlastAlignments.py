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
BlastAlignments.py - tools for working with alignments
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import re, string, sys

try: import alignlib_lite
except ImportError: pass

class Map:
    '''a blast alignment.'''

    def __init__(self ):
        (self.mQueryToken, self.mSbjctToken, self.mEvalue,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli ) =\
         ("", "", 0, 0, 0, "", 0, 0, "" )
        self.mIsExpanded = False
        self.mMapQuery2Sbjct = None

    def Read( self, line ):
        
        (self.mQueryToken, self.mSbjctToken, self.mEvalue,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.mEvalue = float(self.mEvalue)

    def __str__( self ):

        return string.join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.mEvalue,
            self.mQueryFrom, self.mQueryTo, self.mQueryAli,
            self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli)), "\t")

    def MapRange( self, query_token, query_from, query_to ):
        """map something."""
        map_query2sbjct = alignlib_lite.py_makeAlignataVector()
        alignlib_lite.py_fillAlignataCompressed( map_query2sbjct,
                                         self.mQueryFrom, self.mQueryAli,
                                         self.mSbjctFrom, self.mSbjctAli)
        
        new_from = 0

        if query_from <= map_query2sbjct.getRowTo():
            x = max(query_from, self.mQueryFrom)
            while map_query2sbjct.mapRowToCol(x) == 0:
                x += 1
                if x > map_query2sbjct.getRowTo(): break
            else:
               new_from = map_query2sbjct.mapRowToCol(x)

        new_to = 0

        if query_to >= map_query2sbjct.getRowFrom():
            x = min(query_to, self.mQueryTo)
            while map_query2sbjct.mapRowToCol(x) == 0:
                x -= 1
                if x < map_query2sbjct.getRowFrom(): break
            else:
               new_to = map_query2sbjct.mapRowToCol(x)
        
        return self.mSbjctToken, new_from, new_to


    def Expand( self ):
        if not self.mIsExpanded:
            self.mMapQuery2Sbjct = alignlib_lite.py_makeAlignataVector()
            alignlib_lite.py_fillAlignataCompressed( self.mMapQuery2Sbjct,
                                             self.mQueryFrom, self.mQueryAli,
                                             self.mSbjctFrom, self.mSbjctAli )
        self.mIsExpanded = True
        
    def Clear( self ):
        if self.mIsExpanded:
            self.mMapQuery2Sbjct = None
            self.mIsExpanded = False
        
    def Contract( self ):
        if self.mIsExpanded:
            self.mQueryAli, self.mSbjctAli = alignlib_lite.py_writeAlignataCompressed( self.mMapQuery2Sbjct )
            self.mQueryFrom = self.mMapQuery2Sbjct.getRowFrom()
            self.mQueryTo = self.mMapQuery2Sbjct.getRowTo()        
            self.mSbjctFrom = self.mMapQuery2Sbjct.getColFrom()
            self.mSbjctTo = self.mMapQuery2Sbjct.getColTo()
            self.mIsExpanded = False

    def GetClone( self ):
        """get copy of self.
        """

        m = Map()
        
        (m.mQueryToken, m.mSbjctToken, m.mEvalue,
         m.mQueryFrom, m.mQueryTo, m.mQueryAli,
         m.mSbjctFrom, m.mSbjctTo, m.mSbjctAli,
         m.mIsExpanded) = \
         (self.mQueryToken, self.mSbjctToken, self.mEvalue,
          self.mQueryFrom, self.mQueryTo, self.mQueryAli,
          self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
          self.mIsExpanded)

        if self.mIsExpanded:
            m.mMapQuery2Sbjct = alignlib_lite.py_makeAlignataVector()
            alignlib_lite.py_copyAlignata( m.mMapQuery2Sbjct, self.mMapQuery2Sbjct )
            
        return m
        
    def MapAlignment( self, map_query = None, map_sbjct = None ):

        self.Expand()
        
        tmp = alignlib_lite.py_makeAlignataVector()

        if map_query:
            tmp.Clear()
            map_query.Expand()
            alignlib_lite.py_combineAlignata( tmp,
                                      map_query.mMapQuery2Sbjct,
                                      self.mMapQuery2Sbjct,
                                      alignlib_lite.py_RR )
            map_query.Clear()
            alignlib_lite.py_copyAlignata( self.mMapQuery2Sbjct, tmp )
            self.mQueryToken = map_query.mSbjctToken
            
        if map_sbjct:
            tmp.Clear()
            map_sbjct.Expand()
            alignlib_lite.py_combineAlignata( tmp,
                                      self.mMapQuery2Sbjct,
                                      map_sbjct.mMapQuery2Sbjct,
                                      alignlib_lite.py_CR )
            map_sbjct.Clear()
            alignlib_lite.py_copyAlignata( self.mMapQuery2Sbjct, tmp )
            self.mSbjctToken = map_sbjct.mSbjctToken
            
        self.Contract()
        
        return self.IsOk()

    def Reverse( self ):
        (self.mQueryToken, self.mSbjctToken, 
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli ) = \
        ( self.mSbjctToken, self.mQueryToken,
          self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli, 
          self.mQueryFrom, self.mQueryTo, self.mQueryAli )
        
    def IsOk( self ):
        return self.mMapQuery2Sbjct.getRowFrom() != 0 and \
               self.mMapQuery2Sbjct.getRowTo() != 0 and \
               self.mMapQuery2Sbjct.getColFrom() != 0 and \
               self.mMapQuery2Sbjct.getColTo() != 0

        
class SelfMap( Map ):
    def __init__(self):
        Map.__init__(self)
    def __str__(self):
        return ""
    def Read(line):
        pass
    def MapRange( self, query_token, query_from, query_to):
        return query_token, query_from, query_to
        
def ReadMap( file, multiple = False ):
    """read a map from a file.

    If multiple is true, return a list of mappings.
    """
    m = {}
    for line in file:
        if line[0] == "#": continue
        mm = Map()
        mm.Read(line)
        if multiple:
            if mm.mQueryToken not in m: m[mm.mQueryToken] = []
            m[mm.mQueryToken].append(mm)
        else:
            m[mm.mQueryToken] = mm
    return m
        

class Link( Map ):

    mFormat = 0
    
    def __init__(self):
        pass

    def Read( self, line ):

        data = line[:-1].split("\t")

        try:
            (self.mQueryToken, self.mSbjctToken, self.mEvalue,
             self.mQueryFrom, self.mQueryTo, self.mQueryAli,
             self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
             self.score, self.mPercentIdentity, 
             self.mQueryLength, self.mSbjctLength) = data[:13]

            (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo,
             self.score, self.mPercentIdentity, self.mQueryLength, self.mSbjctLength) = map(
                int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo,
                      self.score, self.mPercentIdentity, self.mQueryLength, self.mSbjctLength))
            self.mQueryLength = max( self.mQueryTo, self.mQueryLength )
            self.mSbjctLength = max( self.mSbjctTo, self.mSbjctLength )         
            self.mEvalue = float(self.mEvalue)
        except ValueError:
            raise ValueError, line

        self.mBitScore = 0.0
        self.mFormat = 1
        
        if (len(data) == 14):
            self.mFormat = 2
            self.mBitScore = float(data[13])

    def Reverse( self ):
        """reverse alignment."""
        Map.Reverse( self )
        self.mQueryLength, self.mSbjctLength = self.mSbjctLength, self.mQueryLength
        
    def __str__( self ):

        if self.mFormat == 1:
            return string.join( map(str, (
                self.mQueryToken, self.mSbjctToken, self.mEvalue,
                self.mQueryFrom, self.mQueryTo, self.mQueryAli,
                self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
                self.score, self.mPercentIdentity, 
                self.mQueryLength, self.mSbjctLength)), "\t")
        elif self.mFormat == 2:
            return string.join( map(str, (
                self.mQueryToken, self.mSbjctToken, self.mEvalue,
                self.mQueryFrom, self.mQueryTo, self.mQueryAli,
                self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli,
                self.score, self.mPercentIdentity, 
                self.mQueryLength, self.mSbjctLength, self.mBitScore)), "\t")
            
def iterator_links( infile ):
    """a simple iterator over all entries in a file."""
    while 1:
        line = infile.readline()
        if not line: return
        if line[0] == "#": continue
        v = Link()
        v.Read( line )
        yield v

    
