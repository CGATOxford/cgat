##########################################################################
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
##########################################################################
"""unit testing module for the Tree.py class."""


import numpy
import scipy.stats
import unittest

#import Stats_old as Stats
#from rpy import r as R

import CGAT.Stats as Stats
from rpy2.robjects import r as R
import rpy2.robjects as ro


class TestStats(unittest.TestCase):

    mNumSamples = 100
    mNumReplicates = 1000
    mSignificance = 0.05
    nplaces = 1

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

        for replicate in range(0, self.mNumReplicates):
            sample = scipy.stats.norm.rvs(
                size=self.mNumSamples, loc=0.0, scale=1.0)
            mean = scipy.mean(sample)

            complex_ll = numpy.sum(
                numpy.log(scipy.stats.norm.pdf(sample, loc=mean, scale=1.0)))
            simple_ll = numpy.sum(
                numpy.log(scipy.stats.norm.pdf(sample, loc=0.0, scale=1.0)))

            a = Stats.doLogLikelihoodTest(complex_ll, complex_np,
                                          simple_ll, simple_np,
                                          significance_threshold=self.mSignificance)

            if a.mPassed:
                npassed += 1

        r = float(npassed) / self.mNumReplicates

        self.assertAlmostEqual(self.mSignificance, r, places=self.nplaces)


class TestFDRRAgainstR(unittest.TestCase):

    '''test python against qvalue implementation.
    '''

    mPvalues = [0.033500000000000002, 0.035099999999999999, 0.055, 0.058500000000000003, 0.039199999999999999, 0.045199999999999997, 0.017600000000000001, 0.019099999999999999,
                0.0001, 0.0001, 0.97299999999999998, 0.99429999999999996, 0.99709999999999999, 0.0001, 0.0001, 0.98829999999999996, 0.99529999999999996, 0.0001,
                0.0001, 0.0001, 0.00050000000000000001, 0.00059999999999999995, 0.00050000000000000001, 0.00069999999999999999, 0.40250000000000002, 0.43680000000000002,
                0.35560000000000003, 0.33410000000000001, 0.41260000000000002, 0.3644, 0.12039999999999999, 0.0001, 0.0001, 0.0465, 0.0001, 0.0001, 0.0001, 0.0001, 0.0086,
                0.046300000000000001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0137, 0.0137, 0.0843, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
                0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
                0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]

    def checkFDR(self, pi0_method):

        result = Stats.doFDR(
            self.mPvalues, fdr_level=0.05, pi0_method=pi0_method)
        R("""require ('qvalue')""")
        qvalues = R.qvalue(ro.FloatVector(self.mPvalues),
                           fdr_level=0.05,
                           pi0_method=pi0_method)

        assert qvalues.names[1] == "pi0"
        assert qvalues.names[2] == "qvalues"
        assert qvalues.names[5] == "significant"
        assert qvalues.names[6] == "lambda"

        r_qvalues = qvalues[2]
        r_pi0 = qvalues[1][0]

        self.assertEqual(len(result.mQValues), len(qvalues[2]))
        self.assertEqual(len(result.mLambda), len(qvalues[6]))
        self.assertEqual(result.mPi0, r_pi0)
        for a, b in zip(result.mQValues, r_qvalues):
            self.assertAlmostEqual(a, b, 2, "unequal: %f != %f" % (a, b))

        for a, b in zip(result.mPassed, qvalues[5]):
            self.assertEqual(
                a, b, "threshold-passed flag not equal: %s != %s" % (a, b))

    def testFDRSmoother(self):
        self.checkFDR("smoother")

    def testFDRBootstrap(self):
        self.checkFDR("bootstrap")


def getRelativeError(a, b):
    return abs(a - b) / a


class TestFDRPythonAgainstRDataset1(unittest.TestCase):

    '''test full python implemenation against R.'''

    # there are differences in the smoothing
    # procedure leading to different values of
    # pi0 - hence the resultant qvalues can
    # be slightly different.
    nplaces = 1

    # maximum 10% error
    max_error = 0.10

    def setUp(self):

        R('''library(qvalue)''')
        R('''data(hedenfalk)''')

        self.pvalues = R('''hedenfalk''')
        # print sorted(self.pvalues)
        # print self.pvalues[-1]
        # self.pvalues = numpy.array(range(0,100) + range(0,100),
        # dtype=float) / 100.0
        # print self.pvalues

    def testAgainstQValue(self):

        R.assign("pvalues", self.pvalues)
        qvalue = R('''qvalue( pvalues )''')
        r_qvalues = qvalue[2]
        r_pi0 = qvalue[1][0]

        new = Stats.doFDRPython(self.pvalues)
        self.assertTrue(getRelativeError(r_pi0, new.mPi0) < self.max_error)

        for a, b in zip(r_qvalues, new.mQValues):
            self.assertAlmostEqual(a, b, places=self.nplaces)

    def checkFDR(self, **kwargs):

        old = Stats.doFDR(self.pvalues, **kwargs)
        # print old.mQValues[:10]
        # print old.mPi0
        new = Stats.doFDRPython(self.pvalues, **kwargs)
        # print new.mQValues[:10]
        # print new.mPi0
        # self.assertAlmostEqual( old.mPi0, new.mPi0, places=3)
        self.assertTrue(getRelativeError(old.mPi0, new.mPi0) < self.max_error)

        for pvalue, a, b in zip(self.pvalues, old.mQValues, new.mQValues):
            self.assertTrue(getRelativeError(a, b) < self.max_error,
                            "qvalues: relative error %f > %f (pvalue=%f, %f, %f)" %
                            (getRelativeError(a, b),
                             self.max_error,
                             pvalue, a, b))

    def testFDR(self):
        self.checkFDR()

    def testFDRSmoother(self):
        self.checkFDR(pi0_method="smoother")

    def testFDRBootstrap(self):
        self.checkFDR(pi0_method="bootstrap")

    def testFDRLambda(self):
        self.checkFDR(vlambda=(0.5,))


class TestPValueAdust(unittest.TestCase):

    def setUp(self):
        self.pvalues = ro.FloatVector(
            [3.5703e-01, 1.8901e-01,
             3.0227e-01, 6.8011e-01, 7.1580e-01, 5.0043e-01,
             7.8114e-01, 2.1019e-01, 7.9668e-01, 8.5777e-01,
             8.1345e-01, 5.3087e-01, 7.5629e-01, 7.5666e-01,
             5.2840e-01, 6.6289e-01, 7.7629e-01, 6.4036e-01,
             5.3431e-01, 5.5584e-01, 2.3265e-01, 8.1606e-01,
             6.2403e-01, 6.1609e-01, 7.9606e-01, 6.3858e-01,
             7.2817e-01, 4.1527e-01, 5.0852e-01, 7.8213e-01,
             8.1430e-01, 6.7900e-03, 1.7310e-01, 6.7662e-01,
             7.6978e-01, 2.2915e-01, 7.7845e-01, 8.1132e-01,
             6.5324e-01, 7.0884e-01, 2.7733e-01, 2.2420e-01,
             2.1176e-01, 6.0016e-01, 1.9747e-01, 1.4250e-01,
             7.5438e-01, 5.0737e-01, 7.7629e-01, 8.5746e-01,
             6.7486e-01, 6.8187e-01, 4.9900e-01, 5.5313e-01,
             5.0102e-01, 7.0573e-01, 7.8015e-01, 4.2552e-01,
             5.1993e-01, 7.8006e-01, 2.2966e-01, 7.1619e-01,
             4.7140e-02, 3.2127e-01, 5.8572e-01, 4.7865e-01,
             4.5794e-01, 7.0762e-01, 8.5290e-01, 7.7477e-01,
             4.9756e-01, 4.3011e-01, 3.8470e-02, 2.5932e-01,
             4.0840e-01, 5.8163e-01, 7.7615e-01, 1.5576e-01,
             5.9640e-02, 7.1280e-01, 6.9584e-01, 7.8302e-01,
             9.3500e-02, 3.0228e-01, 1.6430e-02, 5.0919e-01,
             4.6662e-01, 7.9593e-01, 2.2595e-01, 3.4016e-01,
             7.1660e-01, 7.4866e-01, 2.4574e-01, 3.5032e-01,
             5.9100e-01, 6.2572e-01, 7.5720e-01, 7.0475e-01,
             3.0123e-01, 3.5238e-01, 7.8066e-01, 8.5794e-01,
             3.2544e-01, 8.1513e-01, 6.8518e-01, 6.3880e-01,
             5.3906e-01, 4.8648e-01, 7.1588e-01, 2.1637e-01,
             7.9461e-01, 7.1245e-01, 8.8084e-01, 2.7758e-01,
             7.6641e-01, 2.0526e-01, 5.8600e-01, 6.6236e-01,
             6.8191e-01, 7.0945e-01, 5.2711e-01, 6.2196e-01,
             7.0475e-01, 6.7714e-01, 1.2524e-01, 4.2886e-01,
             8.0450e-02, 3.0297e-01, 5.6986e-01, 6.3553e-01,
             6.4007e-01, 6.1447e-01, 6.4272e-01, 8.1268e-01,
             1.9693e-01, 2.8661e-01, 4.5727e-01, 6.0581e-01,
             2.9422e-01, 3.5555e-01, 2.8332e-01, 5.6513e-01,
             4.7400e-01, 6.4527e-01, 7.0818e-01, 6.1335e-01,
             1.9530e-01, 6.1503e-01, 7.0788e-01, 1.0805e-01,
             6.8783e-01, 5.6482e-01, 5.9258e-01, 5.6684e-01,
             6.0419e-01, 5.9690e-01, 1.0411e-01, 2.1848e-01,
             4.8744e-01, 4.4486e-01, 5.9832e-01, 6.2815e-01,
             7.7340e-02, 3.3220e-02, 6.0470e-01, 9.9230e-02,
             6.4443e-01, 5.5593e-01, 3.7584e-01, 6.1314e-01,
             6.4730e-01, 6.8665e-01, 6.1246e-01, 6.4015e-01,
             4.7346e-01, 5.9839e-01, 8.0669e-01, 6.0338e-01,
             6.4697e-01, 4.1835e-01, 3.0939e-01, 4.2020e-01,
             2.8626e-01, 4.1615e-01, 1.4539e-01, 2.5152e-01,
             6.1258e-01, 1.9615e-01, 2.0814e-01, 3.5945e-01,
             2.7513e-01, 1.6680e-01, 7.5770e-02, 6.8493e-01,
             6.3498e-01, 7.4494e-01, 4.2371e-01, 5.6501e-01,
             2.3194e-01, 8.1081e-01, 8.7850e-02, 2.5838e-01,
             1.8457e-01, 5.9697e-01, 5.8886e-01, 6.3387e-01,
             6.4380e-01, 5.9856e-01, 1.5285e-01, 7.0722e-01,
             2.8330e-01, 8.1095e-01, 4.9184e-01, 5.2554e-01,
             5.6463e-01, 4.0671e-01, 2.8000e-04, 5.8700e-03,
             2.1250e-02, 2.5720e-02, 2.6050e-02, 6.8560e-02,
             3.9890e-02, 8.4750e-02, 8.7820e-02, 1.0286e-01,
             7.4510e-02, 1.6138e-01, 1.8390e-01, 3.2850e-02,
             1.1312e-01, 1.1556e-01, 8.5350e-02, 1.2327e-01,
             1.0807e-01, 1.2647e-01, 6.3700e-03, 1.3564e-01,
             1.0702e-01, 2.4620e-02, 1.3361e-01, 2.5921e-01,
             1.3132e-01, 2.7540e-02, 3.9260e-02, 2.2444e-01,
             9.6680e-02, 4.0250e-02, 1.8486e-01, 4.9320e-02,
             1.1898e-01, 2.5090e-02, 7.4260e-02, 1.0508e-01,
             3.4140e-02, 9.8970e-02, 1.0011e-01, 1.2910e-01,
             1.6164e-01, 7.5470e-02, 1.3989e-01, 2.5257e-01,
             2.9892e-01, 2.8947e-01, 1.9050e-02, 2.6359e-01,
             1.8922e-01, 2.4522e-01, 1.0101e-01, 4.0510e-02,
             3.2285e-01, 2.6995e-01, 2.7198e-01, 2.9808e-01,
             3.2103e-01, 1.8112e-01, 2.9776e-01, 3.5214e-01,
             2.5325e-01, 3.7140e-02, 1.8400e-01, 7.5500e-02,
             3.4361e-01, 4.2253e-01, 3.7713e-01, 2.3674e-01,
             2.6463e-01, 3.5779e-01, 3.1044e-01, 4.0422e-01,
             2.9973e-01, 3.9765e-01, 4.1503e-01, 2.4038e-01,
             3.5544e-01, 4.1207e-01, 1.0202e-01, 3.0463e-01,
             2.1991e-01, 4.1872e-01, 2.5065e-01, 2.6679e-01,
             3.9336e-01, 3.3450e-01, 3.0440e-01, 2.6198e-01,
             3.0682e-01, 2.5341e-01, 4.2416e-01, 5.0201e-01,
             3.8728e-01, 4.2616e-01, 2.5298e-01, 1.9127e-01,
             4.3884e-01, 2.5816e-01, 4.8537e-01, 4.3609e-01,
             3.4732e-01, 3.7736e-01, 4.3901e-01, 4.5836e-01,
             4.0918e-01, 4.8454e-01, 2.5898e-01, 4.7588e-01,
             4.2705e-01, 4.6887e-01, 9.1640e-02, 4.6272e-01,
             2.4527e-01, 3.8905e-01, 3.2614e-01, 5.1200e-01,
             2.8509e-01, 4.7901e-01, 3.2064e-01, 3.0500e-01,
             4.0306e-01, 5.0636e-01, 3.0535e-01, 3.9964e-01,
             2.4321e-01, 4.0180e-01, 5.5109e-01, 3.8770e-01,
             5.5091e-01, 5.0608e-01, 3.9073e-01, 4.7774e-01,
             4.4110e-01, 3.8644e-01, 4.3418e-01, 4.1289e-01,
             4.5496e-01, 5.4656e-01, 4.7470e-01, 6.1053e-01,
             5.5539e-01, 3.4375e-01, 3.9771e-01, 5.5518e-01,
             4.7669e-01, 4.3686e-01, 4.8281e-01, 4.4194e-01,
             4.5973e-01, 4.7424e-01, 4.9175e-01, 4.8020e-01,
             4.1258e-01, 5.2228e-01, 4.8238e-01, 4.9247e-01,
             4.9724e-01, 5.0410e-01, 5.1969e-01, 4.9581e-01,
             5.5244e-01, 5.1085e-01, 4.4022e-01, 3.7369e-01,
             4.4847e-01, 4.2010e-01, 4.4931e-01, 4.1517e-01,
             4.1648e-01, 4.3190e-01, 4.1962e-01, 4.1750e-01,
             3.9514e-01, 4.1189e-01, 4.0184e-01, 4.1238e-01,
             4.0556e-01, 3.9327e-01, 3.6691e-01, 3.6744e-01,
             3.9382e-01, 3.8568e-01, 3.1613e-01, 3.7000e-01,
             3.0003e-01, 3.4318e-01, 3.5342e-01, 3.3962e-01,
             2.8939e-01, 3.6290e-01, 3.6277e-01, 3.6685e-01,
             3.3982e-01, 2.7396e-01, 3.5878e-01, 2.9111e-01,
             3.3638e-01, 3.5901e-01, 3.1641e-01, 2.2482e-01,
             3.0998e-01, 3.6248e-01, 3.1717e-01, 3.0584e-01,
             1.6863e-01, 3.0002e-01, 1.2076e-01, 2.7596e-01,
             1.7485e-01, 2.8420e-01, 2.8503e-01, 2.5287e-01,
             3.1915e-01, 3.0764e-01, 2.2573e-01, 3.3076e-01,
             3.0768e-01, 2.6496e-01, 1.0692e-01, 2.4507e-01,
             3.0225e-01, 3.0225e-01, 3.3532e-01, 2.2439e-01,
             1.7248e-01, 2.7667e-01, 2.5522e-01, 2.7345e-01,
             2.6334e-01, 2.3398e-01, 2.4848e-01, 1.3712e-01,
             2.2280e-02, 1.2419e-01, 2.4924e-01, 6.9550e-02,
             3.2560e-02, 1.4077e-01, 2.5225e-01, 1.8327e-01,
             2.9285e-01, 3.2305e-01, 1.3596e-01, 1.8016e-01,
             1.4832e-01, 2.1049e-01, 2.1347e-01, 2.5555e-01,
             2.0492e-01, 2.5225e-01, 1.9864e-01, 1.8680e-01,
             9.0980e-02, 2.2762e-01, 8.1700e-02, 1.1379e-01,
             9.9880e-02, 1.3229e-01, 3.9850e-02, 9.9840e-02,
             1.0044e-01, 2.7016e-01, 2.4448e-01, 6.0010e-02,
             2.4254e-01, 5.4230e-02, 4.9920e-02, 1.8194e-01,
             2.2501e-01, 2.1155e-01, 1.0778e-01, 1.2855e-01,
             1.5890e-02, 1.6487e-01, 1.8918e-01, 2.5900e-02,
             1.0000e-05, 3.8720e-02, 9.2480e-02, 9.7000e-04,
             7.6400e-02, 1.0210e-01, 2.8290e-01, 5.5830e-02,
             2.2050e-02, 8.2920e-02, 1.2863e-01, 8.5360e-02,
             9.4570e-02, 1.0503e-01, 5.3360e-02, 1.4112e-01,
             8.1950e-02, 1.6769e-01, 1.0428e-01, 1.2118e-01,
             7.0000e-04, 2.0936e-01, 1.6847e-01, 3.7100e-03,
             1.3900e-03, 7.2370e-02, 1.8379e-01, 7.0540e-02,
             5.8000e-03, 2.3900e-02, 1.2930e-01, 2.5640e-02,
             7.9580e-02, 1.6998e-01, 9.8000e-03, 1.0145e-01,
             8.7760e-02, 2.8000e-04, 1.4651e-01, 1.1680e-02,
             7.3870e-02, 4.5500e-03, 6.7930e-02, 6.5300e-03,
             1.7852e-01, 1.1650e-02, 3.8330e-02, 3.7310e-02,
             9.9200e-03, 1.3350e-02, 5.6200e-03, 5.7980e-02,
             5.0000e-05, 1.8160e-01, 5.4170e-02, 1.7340e-02,
             1.1130e-01, 1.4980e-02, 4.0000e-04, 7.2970e-02,
             4.9500e-02, 3.9100e-03, 3.7200e-02, 8.7000e-03,
             7.7300e-03, 1.5073e-01, 8.2870e-02, 5.8910e-02,
             3.6800e-03, 3.0700e-03, 1.7353e-01, 1.4780e-01,
             4.0280e-02, 6.1780e-02, 1.5082e-01, 1.2130e-02,
             1.6710e-02, 1.4060e-02, 2.1920e-02, 1.7668e-01,
             5.9400e-02, 2.7370e-02, 1.0050e-01, 1.6980e-02,
             1.1070e-02, 1.0000e-05, 1.6590e-02, 1.2208e-01,
             5.5260e-02, 1.0000e-05, 8.0300e-03, 2.9380e-02,
             4.6300e-03, 4.6700e-02, 1.5253e-01, 8.2330e-02,
             1.3758e-01, 2.7940e-02, 7.1300e-03, 2.3010e-02,
             4.5520e-02, 4.0650e-02, 2.8700e-03, 1.5419e-01,
             7.1000e-04, 2.6500e-03, 4.6600e-03, 8.8840e-02,
             9.8790e-02, 1.5520e-02, 4.4000e-03, 3.8060e-02,
             1.8500e-02, 6.0900e-03, 3.2000e-03, 1.0000e-05,
             9.0900e-03, 3.9170e-02, 9.7300e-03, 3.6000e-04,
             4.8390e-02, 1.1080e-02, 5.3800e-02, 2.5650e-02,
             1.0000e-05, 1.0000e-05, 6.9420e-02, 8.9200e-03,
             7.7000e-03, 1.4830e-02, 4.6280e-02, 4.7000e-04,
             1.9910e-02, 9.0140e-02, 3.6700e-02, 1.8690e-02,
             1.0000e-05, 1.8300e-03, 4.9400e-03, 9.0280e-02,
             3.5920e-02, 2.3400e-03, 8.9980e-02, 1.1668e-01,
             3.2300e-02, 1.5540e-02, 6.6130e-02, 1.4090e-02,
             1.5170e-02, 3.0700e-03, 3.5500e-03, 8.2000e-03,
             5.2000e-04, 1.0560e-02, 6.9420e-02, 1.0000e-05,
             2.9000e-04, 9.7590e-02, 5.1190e-02, 4.5440e-02,
             9.3480e-02, 1.7000e-04, 1.0000e-05, 1.0000e-05,
             1.4600e-03, 3.0400e-03, 1.1000e-04, 8.1900e-03,
             2.8200e-03, 5.6000e-04, 4.0400e-03, 8.6880e-02,
             8.6250e-02, 1.0000e-05, 1.0000e-05, 6.0400e-03,
             1.1000e-04, 3.0500e-03, 1.0000e-05, 3.0000e-05,
             1.2400e-03, 4.7600e-03, 8.5000e-03, 1.0000e-05,
             8.2200e-03, 2.0000e-05, 5.2200e-03, 1.3000e-04,
             1.0000e-05, 1.2500e-02, 2.5500e-03, 2.6100e-02,
             9.0000e-04, 2.0000e-04, 2.4200e-03, 3.2190e-02,
             3.2170e-02, 3.2070e-02, 1.8100e-03, 3.8500e-03,
             1.7700e-03, 1.0000e-05, 1.1000e-04, 1.0000e-05,
             1.0000e-05, 1.1700e-03, 1.0000e-05, 1.0000e-05,
             3.2300e-03, 2.1000e-04, 3.3700e-03, 3.1400e-03,
             1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
             8.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
             2.1500e-03, 1.0000e-05, 1.0000e-05, 1.0000e-05,
             1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
             1.0000e-05, 1.0000e-05])

    def check(self, method):
        '''check for length equality and elementwise equality.'''
        a = R['p.adjust'](self.pvalues, method=method)
        b = Stats.adjustPValues(self.pvalues, method=method)
        self.assertEqual(len(a), len(b))
        for x, y in zip(a, b):
            self.assertAlmostEqual(x, y)

    def testBonferroni(self):
        self.check("bonferroni")

    def testHolm(self):
        self.check("holm")

    def testHommel(self):
        pass
        # self.check("hommel")

    def testHochberg(self):

        # code for checking
        R.assign("p", self.pvalues)
        R('''
        lp = length(p)
        n = length(p)
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin((n - i + 1L) * p[o]))[ro]
        ''')

        self.check("hochberg")

    def testBH(self):
        self.check("BH")

    def testBY(self):
        self.check("BY")

    def testNone(self):
        self.check("none")

if __name__ == "__main__":
    unittest.main()
