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
ScopTester.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys
import re
import string
import os
import time 

from Pairsdb import *

import alignlib
import pairsdblib

from MessagePairsdb import MessagePairsdb
from TableDomainsScopTest import TableDomainsScopTest
from TablePairsdbNeighbours import TablePairsdbNeighbours
from Pairsdb import *
import Tools

#-------------------------------------------
# Class:	       ScopTest
# Superclasses:  Message
# Subclasses:    
# Function:      update ScopTest-database
#
# Author:	       Andreas Heger
#-------------------------------------------

class ScopTester:
    
    def __init__ (self, dbhandle, alignator, table_scop_test ):

        self.mLogLevel = 1
        
        self.dbhandle = dbhandle 
        # C++ connection to pairsdb
        self.mDatabaseNamePairsdb = dbhandle.GetDatabase()

        self.mConnectionPairsdb  = pairsdblib.Connection( dbhandle.GetHost(),
                                                          dbhandle.GetUser(),
                                                          dbhandle.GetPassword(),
                                                          dbhandle.GetPort())
        
        self.mConnectionPairsdb.Connect( self.mDatabaseNamePairsdb )

        self.mTableScopTest = table_scop_test
        
        self.mAlignanda = []
        self.mInformation = []
        
        self.mAlignator = alignator

        self.startAt = 0
        
    ##--------------------------------------------------------------------------------------
    def CalculateMatches( self ):
        """calculate all-vs-all alignments.
        """

        if not self.mAlignanda:
            self.GetAlignanda()
            
        if self.mLogLevel >= 1:
            print "# --> calculating alignments for %i entries" % len(self.mAlignanda)
            print "# --> starting at:", Tools.GetTimeStamp()

        nalignanda = len(self.mAlignanda)
        
        for a1 in range(self.startAt, nalignanda-1):
            if self.mLogLevel >= 1:
                print "# %5i/%5i at %s" % (a1, nalignanda, Tools.GetTimeStamp())
                sys.stdout.flush()
                
            for a2 in range(a1+1,nalignanda):
                
                if self.mLogLevel >= 3:
                    print "#    aligning to %i" % (a2), self.mInformation[a2]
                    sys.stdout.flush()

                result = alignlib.makeAlignataVector()

                self.mAlignator.Align( self.mAlignanda[a1], self.mAlignanda[a2], result )

                info = self.mAlignator.CheckResult( result, self.mInformation[a1], self.mInformation[a2] )
                
                if info:
                    r = tuple(self.mInformation[a1]) + tuple(self.mInformation[a2]) + tuple(info)
                    print string.join(r, "\t" )
                
                sys.stdout.flush()
                
            self.mAlignanda[a1].Release()
            self.mAlignanda[a1] = None
            
        if self.mLogLevel >= 1:
            print "# --> finished at:", Tools.GetTimeStamp()

    ##--------------------------------------------------------------------------------------            
    def SetAlignator( self, alignator ):
        self.mAlignator = alignator

    ##--------------------------------------------------------------------------------------
    def GetAlignanda( self ):
        """retrieve alignandum-objects.
        """

        ## do not get additional info so that non-redundant table can be used.
        domains = self.mTableScopTest.GetAllDomains( all = 0 )

        if self.mLogLevel >= 1:
            print "# --> retrieving %i entries at %s" % (len(domains), Tools.GetTimeStamp())
            sys.stdout.flush()
            
        for domain in domains:

            if self.mLogLevel >= 2:
                print "# retrieving", domain
                sys.stdout.flush()

            (nid, nrdb_from, nrdb_to, scop_class) = domain
            
            ## dummy values
            pdb_id = "test"
            region_nr = 0
            
            if scop_class[0:3] not in ("00a", "00b", "00c", "00d"):
                if self.mLogLevel >= 2:
                    print "# skipped because not in first four classes"
                    sys.stdout.flush()
                continue
            
            # if nid not in (47268, 74355): continue
            # if nid not in (36388, 148361): continue
            # if nid not in (3503, 115681): continue
            # if nid not in (17, 1011060): continue                        
            alignandum, info = self.GetAlignandum( nid, nrdb_from, nrdb_to )

            if alignandum:
                if info:
                    self.mInformation.append( ( "%i_%i_%i" % (nid, nrdb_from, nrdb_to),
                                                scop_class, pdb_id, str(region_nr)) + tuple(info) )
                    
                else:
                    self.mInformation.append( ( "%i_%i_%i" % (nid, nrdb_from, nrdb_to),
                                                scop_class, pdb_id, str(region_nr)) )
                    
                self.mAlignanda.append(alignandum)
            else:
                if self.mLogLevel >= 2:
                    print "# skipped because no alignandum found"
                    sys.stdout.flush()

        if self.mLogLevel >= 1:
            print "# --> retrieved %i entries at %s" % (len(self.mAlignanda), Tools.GetTimeStamp())
            sys.stdout.flush()
                    
##--------------------------------------------------------------------------------------                                
class ScopTesterSequences( ScopTester ):
    
    def GetAlignandum( self, nid, nrdb_from, nrdb_to ):

        alignandum = pairsdblib.makeSequenceFromPairsdb( self.mConnectionPairsdb, nid )
        alignandum.useSegment( nrdb_from, nrdb_to )

        return alignandum, None

##--------------------------------------------------------------------------------------                                
class ScopTesterProfiles( ScopTester ):

    def __init__(self,
                 dbhandle,
                 alignator,
                 table_scop_test,
                 min_profile_size = 20,
                 min_level = 30,
                 max_level = 90,
                 neighbours = "pairsdb_90x90"):
        
        self.mMinProfileSize = min_profile_size
        self.mMinLevel = min_level
        self.mMaxLevel = max_level
        self.mTableNameNeighbours = neighbours
        
        ScopTester.__init__( self, dbhandle, alignator, table_scop_test )

        self.mTableNeighbours = TablePairsdbNeighbours( self.dbhandle )
        self.mTableNeighbours.SetName( self.mTableNameNeighbours)

        self.mBlastL             = 0.3                  # lambda    
        self.mLogOddorScaleFactor = self.mBlastL
        self.mLogOddor    = alignlib.makeLogOddorDirichlet( self.mLogOddorScaleFactor )
        self.mMaxLinesMali = 1000
        self.mRegularizor = alignlib.makeRegularizorDirichletPrecomputed()  

    def GetAlignandum( self, nid, nrdb_from, nrdb_to ):

        n = self.mTableNeighbours.GetNumNeighbours( nid )
        
        if n >= self.mMinProfileSize:
            
            profile = alignlib.makeEmptyProfile( self.mRegularizor, self.mLogOddor )

            pairsdblib.fillProfileNeighbours( profile,
                                              self.mConnectionPairsdb,
                                              nid,
                                              self.mTableNameNeighbours,
                                              self.mMaxLinesMali,
                                              0,
                                              self.mMaxLevel,
                                              self.mMinLevel )
            if self.mLogLevel >= 3:
                print "# ------> using profile for rep %i" % nid
                
        else:
            profile = pairsdblib.makeSequenceFromPicasso( self.mConnectionPairsdb, nid )
            if self.mLogLevel >= 3:
                print "# ------> using sequence for rep %i" % nid
            
        profile.useSegment( nrdb_from, nrdb_to )
                                          
        return profile, (str(n),)

##--------------------------------------------------------------------------------------                                
class ScopTesterFullProfiles( ScopTesterProfiles ):
    """use full length profiles.
    beware of multidomain-proteins, use iterative multiple alignment
    method.
    """
    def __init__(self,
                 dbhandle,
                 alignator,
                 table_scop_test,
                 min_profile_size = 20,
                 min_level = 30,
                 max_level = 90,
                 neighbours = "pairsdb_90x90"):

        ScopTesterProfiles.__init__( self,dbhandle, alignator, table_scop_test,
                                     min_profile_size, min_level, max_level, neighbours)

        self.mAddLength = 500
        self.mMaxLength = 2000
        
    def GetAlignandum( self, nid, nrdb_from, nrdb_to ):

        profile, x = ScopTesterProfiles.GetAlignandum( self, nid, nrdb_from, nrdb_to)

        profile.useFullLength()

        # add some context around
        xfrom = max( 1, nrdb_from - self.mAddLength)
        xto = min( profile.getLength(), nrdb_to + self.mAddLength)

        if self.mLogLevel >= 3:
            print "using segment %i-%i" % (xfrom, xto)
            sys.stdout.flush()
        
        profile.useSegment( xfrom, xto)

        return profile, x

##--------------------------------------------------------------------------------------
class Alignator:
    """
    aligns two sequences and returns result.
    """
    
    def __init__( self, alignator ):
        self.mAlignator = alignator

    ##--------------------------------------------------------------------------------------
    def Align( self, a1, a2, result ):
        self.mAlignator.Align( a1, a2, result )
        
    ##--------------------------------------------------------------------------------------
    def CheckResult( self,
                     result,
                     info1 = None,
                     info2 = None):
        """check if result is ok. The function below returns everything.
        return tuple of strings as result.
        """

        if (result.getLength() > 0):
            row_ali, col_ali = alignlib.writeAlignataCompressed( result )
            return map(str, (result.getScore(),
                             result.getLength(),
                             result.getNumGaps(),
                             alignlib.calculatePercentSimilarity( result ),
                             result.getRowFrom(), result.getRowTo(), row_ali,
                             result.getColFrom(), result.getColTo(), col_ali ) )
        else:
            return ("0",) * 12
        
##--------------------------------------------------------------------------------------
class AlignatorIterative(Alignator):
    """
    aligns two sequences iteratively, checks if alignment regions are overlapping with
    domain regions and returns result only for those overlapping. This is useful if you have
    several domains in a sequence, but you need only compare to one.
    """

    def __init__( self, alignator, min_score, min_overlap, gop, gep ):
        self.mAlignator = alignator
        self.mMinScore = float(min_score)
        self.mMinOverlap = min_overlap
        self.mGop = gop
        self.mGep = gep

    ##--------------------------------------------------------------------------------------
    def Align( self, a1, a2, result ):
        """align repetetively. Take highest scoring alignment, that overlaps with domains
        and put it in result. Note: performIterativeAlignment does not work, as it is linear.
        It requires domains to be in the same order.

        Result is empty, fragments are saved in object.
        """
        
        fragmentor = alignlib.makeFragmentorRepetitive( self.mAlignator, self.mMinScore )

        ## align iteratively and convert fragments to Alignata-objects
        val = fragmentor.Fragment( a1, a2, result)
        
        self.mFragments = map( lambda x: alignlib.AlignataPtr(x), val)

        for fragment in self.mFragments:
            fragment.thisown = 1
                    
##         alignlib.performIterativeAlignmentNonConst( result,
##                                                     a1, a2,
##                                                     self.mAlignator,
##                                                     self.mMinScore )
                                                    
    ##--------------------------------------------------------------------------------------
    def CheckResult( self,
                     result,
                     info1, info2):
        """check if result is ok. Check for each fragment, if it overlaps
        with the domains to be tested and dump if ok. This simulates
        psiblast.
        """

        row_from, row_to = map(string.atoi, info1[1:3])
        col_from, col_to = map(string.atoi, info2[1:3])

        ## check for overlap
        for fragment in self.mFragments:
            
            # print alignlib.writeAlignataTable( fragment, 8, 1)
            
            xcol_from = Tools.MapRight(fragment, row_from )
            xcol_to   = Tools.MapLeft(fragment, row_to )

            overlap = min(col_to, xcol_to) - max(col_from, xcol_from)

            # print self.mMinOverlap, overlap, xcol_from, xcol_to, col_from, col_to

            if overlap > self.mMinOverlap:
            
                return map(str, (fragment.getScore(),
                                 fragment.getLength(),
                                 fragment.getNumGaps(),
                                 alignlib.calculatePercentSimilarity( fragment ),
                                 fragment.getRowFrom(), fragment.getRowTo(),
                                 fragment.getColFrom(), fragment.getColTo(),
                                 overlap, xcol_from, xcol_to,
                                 (xcol_to - xcol_from) - (col_to - col_from)) )
            
        return ("0",) * 12
        
##--------------------------------------------------------------------------------------
                
if __name__ == '__main__':

    dbhandle = Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed"
	sys.exit(1)
    
    a = alignlib.makeFullDP( -10.0, -2.0 )
    alignator = Alignator( a )

    x = ScopTesterSequences( dbhandle, alignator )
    
    x.Process()

    if param_alignator == 0:
        a = alignlib.makeFullDP( param_gop, param_gep)

        alignator = Alignator( a )
    if param_entities == 0:
        tester = ScopTesterSequences( dbhandle, alignator )

        tester.mLogLevel = param_loglevel

        matches = a.CalculateMatches()
                

    



