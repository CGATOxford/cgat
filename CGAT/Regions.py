################################################################################
#   Gene prediction pipeline 
#
#   $Id: Regions.py 2781 2009-09-10 11:33:14Z andreas $
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
Regions.py - helper functions for working with genomic segments
===============================================================

Version: $Id: Regions.py 2781 2009-09-10 11:33:14Z andreas $
"""

import os, string, re, bisect

class RegionFilter:
    """Filter class based on regions."""

    def __init__(self):
        
        self.mRegions = {}
    
    def readFromFile( self, infile, ignore_strand = False ):
        """read regions from a file."""
    
        self.mForwardRegions = {}
        self.mReverseRegions = {}
        self.mRegions = []
        self.mIgnoreStrand = ignore_strand
        n = 0
        for line in infile:
            if line[0] == "#": continue

            token, sbjct_token, sbjct_strand, sbjct_from, sbjct_to = line[:-1].split("\t")[:5]

            if ignore_strand:
                key = sbjct_token
            else:
                key = "%s-%s" % (sbjct_token, sbjct_strand)
                
            if key not in self.mForwardRegions:
                self.mForwardRegions[key] = []
                self.mReverseRegions[key] = []
                
            self.mForwardRegions[key].append( (int(sbjct_from), n) )
            self.mReverseRegions[key].append( (int(sbjct_to), n) )            
            self.mRegions.append( (token, sbjct_from, sbjct_to) )
            n += 1

        for k, v in self.mForwardRegions.items():
            v.sort()
            self.mForwardRegions[k] = ( map(lambda x: x[0], v),
                                        map(lambda x: x[1], v) )
                                       
        for k, v in self.mReverseRegions.items():
            v.sort()
            self.mReverseRegions[k] = ( map(lambda x: x[0], v),
                                        map(lambda x: x[1], v) )
                                       

    def getOverlaps( self, sbjct_token, sbjct_strand, sbjct_from, sbjct_to):
        """return overlapping regions with region."""

        if self.mIgnoreStrand:
            key = sbjct_token
        else:
            key = "%s-%s" % (sbjct_token, sbjct_strand)
        
        if key not in self.mForwardRegions: return []

        ## index of intervalls starting after sbjct_to
        index_from = bisect.bisect_right( self.mForwardRegions[key][0], sbjct_to - 1)
        ## index of intervalls stopping before sbjct_from
        index_to   = bisect.bisect_left( self.mReverseRegions[key][0], sbjct_from + 1)
        
        set1 = set(self.mForwardRegions[key][1][:index_from])
        set2 = set(self.mReverseRegions[key][1][index_to:])

        overlaps = set1.intersection(set2)
        
#         for o in overlaps:
#             print "overlaps", overlaps
            
#             print "set1: regions not starting after %i" % sbjct_to
#             for x in set1:
#                 print self.mRegions[x]
                
#             print "set1: regions not stopping before %i" % sbjct_from
#             for x in set2:
#                 print self.mRegions[x]                

#             for x in overlaps:
#                 print sbjct_token, sbjct_strand, sbjct_from, sbjct_to, "overlaps with", self.mRegions[x]

        return overlaps
