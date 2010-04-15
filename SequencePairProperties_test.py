################################################################################
#   Gene prediction pipeline 
#
#   $Id: SequencePairProperties_test.py 2792 2009-09-16 15:27:03Z andreas $
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
"""unit testing module for the Tree.py class."""
import random, re


from SequencePairProperties import *
import unittest

class SequencePairPropertiesNATest(unittest.TestCase):
    
    def testCounts1(self):
        seq1 = "A" * 5 + "C" * 5 + "-" * 5 + "G" * 5 + "T"*5
        seq2 = seq1
        
        c = SequencePairPropertiesCountsNa()
        c.loadPair( seq1, seq2 )
        
        self.assertEqual( c.mNAligned, 20 )
        self.assertEqual( c.mPercentGC, "50.00" )

    def testCounts2(self):
        seq1 = "A" * 5 + "C" * 5 + "-" * 5 + "G" * 5 + "T"*5
        seq2 = "T" * 5 + "G" * 5 + "-" * 5 + "C" * 5 + "A"*5
        
        c = SequencePairPropertiesCountsNa()
        c.loadPair( seq1, seq2 )
        
        self.assertEqual( c.mNAligned, 20 )
        self.assertEqual( c.mPercentGC, "50.00" )

    def testCounts3(self):
        seq1 = "A" * 20
        seq2 = "T" * 20
        
        c = SequencePairPropertiesCountsNa()
        c.loadPair( seq1, seq2 )
        
        self.assertEqual( c.mNAligned, 20 )
        self.assertEqual( c.mPercentGC, " 0.00" )

    def testCounts4(self):
        seq1 = "G" * 20
        seq2 = "C" * 20
        
        c = SequencePairPropertiesCountsNa()
        c.loadPair( seq1, seq2 )
        
        self.assertEqual( c.mNAligned, 20 )
        self.assertEqual( c.mPercentGC, "100.00" )

    def testCounts5(self):
        seq1 = "G" * 20
        seq2 = "A" * 20
        
        c = SequencePairPropertiesCountsNa()
        c.loadPair( seq1, seq2 )
        
        self.assertEqual( c.mNAligned, 20 )
        self.assertEqual( c.mPercentGC, "50.00" )

    def testCounts5(self):

        a, b = [], []
        for x in range(0,100):
            a.append( random.choice( "ACGT" ) ) 
            b.append( random.choice( "ACGT" ) )
        a,b = "".join(a), "".join(b)

        gc = 100.0 * float(len(re.sub( "[^GC]", "", a)) + len(re.sub( "[^GC]", "", b))) / 100.0 / 2.0

        c = SequencePairPropertiesCountsNa()
        c.loadPair( a, b )
        
        self.assertEqual( c.mNAligned, 100 )
        self.assertEqual( c.mPercentGC, "%5.2f" %  gc )
        
        


            

if __name__ == "__main__":
    unittest.main()
