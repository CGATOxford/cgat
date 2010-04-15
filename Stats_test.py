################################################################################
#   Gene prediction pipeline 
#
#   $Id: Stats_test.py 2784 2009-09-10 11:41:14Z andreas $
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

from Stats import *
import numpy
import scipy.stats
import unittest
from rpy import r as R

class StatsCheck(unittest.TestCase):

    mNumSamples = 100
    mNumReplicates = 1000
    mSignificance = 0.05

    mPvalues = [0.033500000000000002, 0.035099999999999999, 0.055, 0.058500000000000003, 0.039199999999999999, 0.045199999999999997, 0.017600000000000001, 0.019099999999999999, 
                0.0001, 0.0001, 0.97299999999999998, 0.99429999999999996, 0.99709999999999999, 0.0001, 0.0001, 0.98829999999999996, 0.99529999999999996, 0.0001, 
                0.0001, 0.0001, 0.00050000000000000001, 0.00059999999999999995, 0.00050000000000000001, 0.00069999999999999999, 0.40250000000000002, 0.43680000000000002, 
                0.35560000000000003, 0.33410000000000001, 0.41260000000000002, 0.3644, 0.12039999999999999, 0.0001, 0.0001, 0.0465, 0.0001, 0.0001, 0.0001, 0.0001, 0.0086, 
                0.046300000000000001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0137, 0.0137, 0.0843, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 
                0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 
                0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]

    def testLRT(self):
        """test that the false positive rate is in the same order as mSignificance.

        Sample from a normal distribution and compare two models:

        1. mean estimated = complex model (1 df)
        2. mean given     = simple model  (0 df)

        Likelihood = P(model | data)
        """
        simple_np = 0
        complex_np = 1

        npassed = 0

        for replicate in range(0,self.mNumReplicates):
            sample = scipy.stats.norm.rvs(size=self.mNumSamples, loc = 0.0, scale = 1.0)
            mean = scipy.mean( sample )

            complex_ll = numpy.sum( numpy.log( scipy.stats.norm.pdf( sample, loc = mean, scale = 1.0 ) ) )
            simple_ll = numpy.sum( numpy.log( scipy.stats.norm.pdf( sample, loc = 0.0, scale = 1.0 ) ) )
            
            a =  doLogLikelihoodTest( complex_ll, complex_np,
                                      simple_ll, simple_np,
                                      significance_threshold = self.mSignificance )
            
            if a.mPassed: npassed += 1

        r = float(npassed) / self.mNumReplicates

        self.assertAlmostEqual( self.mSignificance, r, 2, "number of passed tests does not reflect chi-square distn" )

    def checkFDR( self, pi0_method ):

        result = doFDR( self.mPvalues, fdr_level = 0.05, pi0_method=pi0_method )
        R("""require ('qvalue')""")
        qvalues = R.qvalue( self.mPvalues, fdr_level = 0.05, pi0_method=pi0_method )
        self.assertEqual( len(result.mQValues), len(qvalues["qvalues"]) )
        self.assertEqual( len(result.mLambda), len(qvalues["lambda"]) )
        self.assertEqual( result.mPi0, qvalues["pi0"] )
        for a,b in zip( result.mQValues, qvalues["qvalues"] ):
            self.assertAlmostEqual( a, b, 2, "unequal: %f != %f" % (a,b) )

        for a,b in zip( result.mPassed, qvalues["significant"] ):
            self.assertEqual( a, b, "threshold-passed flag not equal: %s != %s" % (a,b) )

    def testFDRSmoother( self ):
        self.checkFDR( "smoother" )

    def testFDRBootstrap( self ):
        self.checkFDR( "bootstrap" )

if __name__ == "__main__":
    unittest.main()
