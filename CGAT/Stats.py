'''
Stats.py - statistical utility functions
========================================

:Tags: Python

Code
----

'''
import math
import numpy
import scipy
import scipy.stats
import scipy.interpolate
import collections
from rpy2.robjects import r as R
import rpy2.robjects as ro
from functools import reduce


def getSignificance(pvalue, thresholds=[0.05, 0.01, 0.001]):
    """return cartoon of significance of a p-Value."""
    n = 0
    for x in thresholds:
        if pvalue > x:
            return "*" * n
        n += 1
    return "*" * n


class Result(object):

    '''allow both member and dictionary access.'''
    slots = ("_data")

    def __init__(self):
        object.__setattr__(self, "_data", dict())

    def fromR(self, take, r_result):
        '''convert from an *r_result* dictionary using map *take*.

        *take* is a list of tuples mapping a field to the corresponding
        field in *r_result*.
        '''
        for x, y in take:
            if y:
                self._data[x] = r_result.rx(y)[0][0]
            else:
                self._data[x] = r_result.rx(x)[0][0]

            # if y:
            #     self._data[x] = r_result[y]
            # else:
            #     self._data[x] = r_result[x]

        return self

    def __getattr__(self, key):
        if not key.startswith("_"):
            try:
                return object.__getattribute__(self, "_data")[key]
            except KeyError:
                pass
        return getattr(self._data, key)

    def keys(self):
        return list(self._data.keys())

    def values(self):
        return list(self._data.values())

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return self._data.__len__()

    def __str__(self):
        return str(self._data)

    def __contains__(self, key):
        return key in self._data

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __setattr__(self, key, value):
        if not key.startswith("_"):
            self._data[key] = value
        else:
            object.__setattr__(self, key, value)

    def __getstate__(self):
        # required for correct pickling/unpickling
        return object.__getattribute__(self, "_data")

    def __setstate__(self, d):
        # required for correct unpickling, otherwise
        # maximum recursion threshold will be reached
        object.__setattr__(self, "_data", d)


class LogLikelihoodTest:

    def __init__(self):
        pass


def doLogLikelihoodTest(complex_ll, complex_np,
                        simple_ll, simple_np,
                        significance_threshold=0.05):
    """perform log-likelihood test between model1 and model2.
    """

    assert complex_ll >= simple_ll, "log likelihood of complex model smaller than for simple model: %f > %f" % (
        complex_ll, simple_ll)

    chi = 2 * (complex_ll - simple_ll)
    df = complex_np - simple_np

    if df <= 0:
        raise ValueError("difference of degrees of freedom not larger than 0")

    p = scipy.stats.chisqprob(chi, df)

    l = LogLikelihoodTest()

    l.mComplexLogLikelihood = complex_ll
    l.mSimpleLogLikelihood = simple_ll
    l.mComplexNumParameters = complex_np
    l.mSimpleNumParameters = simple_np
    l.mSignificanceThreshold = significance_threshold
    l.mProbability = p
    l.mChiSquaredValue = chi
    l.mDegreesFreedom = df

    if p < significance_threshold:
        l.mPassed = True
    else:
        l.mPassed = False

    return l


class BinomialTest:

    def __init__(self):
        pass


def doBinomialTest(p, sample_size, observed, significance_threshold=0.05):
    """perform a binomial test.

    Given are p: the probability of the NULL hypothesis, the sample_size
    and the number of observed counts.
    """


class ChiSquaredTest:

    def __init__(self):
        pass


def doChiSquaredTest(matrix, significance_threshold=0.05):
    '''perform chi-squared test on a matrix.

    The observed/expected values are in rows, the categories are in
    columns, for example:

    +---------+--------------+--------+----------+
    |set      |protein_coding|intronic|intergenic|
    +---------+--------------+--------+----------+
    |observed |92            |90      |194       |
    +---------+--------------+--------+----------+
    |expected |91            |10      |15        |
    +---------+--------------+--------+----------+

    If there are only two categories (one degrees of freedom) the
    Yates correction is applied.  For each entry (observed-expected),
    the value 0.5 is subtracted ignoring the sign of the difference.

    The test throws an exception if

    1. one or more expected categories are less than 1 (it does not
    matter what the observed values are)

    2. more than one-fifth of expected categories are less than 5

    '''

    nrows, ncols = matrix.shape

    if nrows != 2:
        raise NotImplementedError(
            "chi-square currently only implemented for 2xn tables.")

    n = 0
    for x in range(ncols):
        if matrix[1][x] < 1:
            raise ValueError("matrix contains expected counts < 1")
        if matrix[1][x] < 5:
            n += 1
    if 100.0 * n / ncols > 20.0:
        raise ValueError(
            "more than 20% of expected categories are less than 5")

    row_sums = [sum(matrix[x, :]) for x in range(nrows)]
    col_sums = [sum(matrix[:, x]) for x in range(ncols)]
    sample_size = float(sum(row_sums))
    chi = 0.0

    df = (nrows - 1) * (ncols - 1)
    # Yates correction applies for a 2x2 table only (df==1)
    if df == 1:
        correction = 0.5 * 0.5
    else:
        correction = 0

    for x in range(nrows):
        for y in range(ncols):
            expected = row_sums[x] * col_sums[y] / sample_size
            # compute difference and apply Yates correction
            d = abs(matrix[x, y] - expected) - correction
            chi += (d * d) / expected

    result = ChiSquaredTest()

    result.mProbability = scipy.stats.chisqprob(chi, df)
    result.mDegreesFreedom = df
    result.mChiSquaredValue = chi
    result.mPassed = result.mProbability < significance_threshold
    result.mSignificance = getSignificance(result.mProbability)
    result.mSampleSize = sample_size
    result.mPhi = math.sqrt(result.mChiSquaredValue / result.mSampleSize)
    return result


def doPearsonChiSquaredTest(p, sample_size, observed,
                            significance_threshold=0.05):
    """perform a pearson chi squared test.

    Given are p: the probability of the NULL hypothesis, the sample_size
    and the number of observed counts.

    For large sample sizes, this test is a continuous approximation to
    the binomial test.
    """
    e = float(p) * sample_size
    d = float(observed) - e
    chi = d * d / e
    df = 1

    result = ChiSquaredTest()

    result.mProbability = scipy.stats.chisqprob(chi, df)
    result.mDegreesFreedom = df
    result.mChiSquaredValue = chi
    result.mPassed = result.mProbability < significance_threshold
    result.mSignificance = getSignificance(result.mProbability)
    result.mSampleSize = sample_size
    result.mPhi = math.sqrt(result.mChiSquaredValue / result.mSampleSize)
    result.mObserved = observed
    result.mExpected = e
    return result


class DistributionalParameters:

    """a collection of distributional parameters. Available properties
    are:

    mMean, mMedian, mMin, mMax, mSampleStd, mSum, mCounts

    This method is deprecated - use :class:`Summary` instead.
    """

    def __init__(self, values=None, format="%6.4f", mode="float"):

        self.mMean, self.mMedian, self.mMin, self.mMax, self.mSampleStd, self.mSum, self.mCounts, self.mQ1, self.mQ3 = \
            (0, 0, 0, 0, 0, 0, 0, 0, 0)

        if values is not None and len(values) > 0:
            self.updateProperties(values)
        self.mFormat = format
        self.mMode = mode
        self.mNErrors = 0

    def updateProperties(self, values):
        """update properties.

        If values is an vector of strings, each entry will be converted
        to float. Entries that can not be converted are ignored.
        """
        values = [x for x in values if x is not None]

        if len(values) == 0:
            raise ValueError("no data for statistics")

        # convert
        self.mNErrors = 0
        if type(values[0]) not in (int, float):
            n = []
            for x in values:
                try:
                    n.append(float(x))
                except ValueError:
                    self.mNErrors += 1
        else:
            n = values

        if len(n) == 0:
            raise ValueError("no data for statistics")

        # use a non-sort algorithm later.
        n.sort()
        self.mQ1 = n[len(n) / 4]
        self.mQ3 = n[len(n) * 3 / 4]

        self.mCounts = len(n)
        self.mMin = min(n)
        self.mMax = max(n)
        self.mMean = scipy.mean(n)
        self.mMedian = scipy.median(n)
        self.mSampleStd = scipy.std(n)
        self.mSum = reduce(lambda x, y: x + y, n)

    def getZScore(self, value):
        """return zscore for value."""
        if self.mSampleStd > 0:
            return (value - self.mMean) / self.mSampleStd
        else:
            return 0

    def setFormat(self, format):
        """set number format."""
        self.mFormat = format

    def getHeaders(self):
        """returns header of column separated values."""
        return ("nval", "min", "max", "mean", "median",
                "stddev", "sum", "q1", "q3")

    def getHeader(self):
        """returns header of column separated values."""
        return "\t".join(self.getHeaders())

    def items(self):
        return [(x, self.__getitem__(x)) for x in self.getHeaders()]

    def __getitem__(self, key):

        if key == "nval":
            return self.mCounts
        if key == "min":
            return self.mMin
        if key == "max":
            return self.mMax
        if key == "mean":
            return self.mMean
        if key == "median":
            return self.mMedian
        if key == "stddev":
            return self.mSampleStd
        if key == "sum":
            return self.mSum
        if key == "q1":
            return self.mQ1
        if key == "q3":
            return self.mQ3

        raise KeyError(key)

    def __str__(self):
        """return string representation of data."""

        if self.mMode == "int":
            format_vals = "%i"
            format_median = "%.1f"
        else:
            format_vals = self.mFormat
            format_median = self.mFormat

        return "\t".join(("%i" % self.mCounts,
                          format_vals % self.mMin,
                          format_vals % self.mMax,
                          self.mFormat % self.mMean,
                          format_median % self.mMedian,
                          self.mFormat % self.mSampleStd,
                          format_vals % self.mSum,
                          format_vals % self.mQ1,
                          format_vals % self.mQ3,
                          ))


class Summary(Result):

    """a collection of distributional parameters. Available properties
    are:

    mean, median, min, max, samplestd, sum, counts
    """

    fields = ("nval", "min", "max", "mean",
              "median", "stddev", "sum", "q1", "q3")

    def __init__(self, values=None,
                 format="%6.4f", mode="float",
                 allow_empty=True):

        Result.__init__(self)
        self._format = format
        self._mode = mode

        # note that this determintes the order of the fields at output
        self.counts, self.min, self.max, self.mean, self.median, self.samplestd, self.sum, self.q1, self.q3 = \
            (0, 0, 0, 0, 0, 0, 0, 0, 0)

        if values is not None:

            values = [x for x in values if x is not None]

            if len(values) == 0:
                if allow_empty:
                    return
                else:
                    raise ValueError("no data for statistics")

            # convert
            self._nerrors = 0
            if type(values[0]) not in (int, float):
                n = []
                for x in values:
                    try:
                        n.append(float(x))
                    except ValueError:
                        self._nerrors += 1
            else:
                n = values

            # use a non-sort algorithm?
            n.sort()
            if len(n):
                self.q1 = n[len(n) // 4]
                self.q3 = n[len(n) * 3 // 4]
            else:
                self.q1 = self.q3 = 0

            self.counts = len(n)
            self.min = min(n)
            self.max = max(n)
            self.mean = scipy.mean(n)
            self.median = scipy.median(n)
            self.samplestd = scipy.std(n)
            self.sum = reduce(lambda x, y: x + y, n)

    def getHeaders(self):
        """returns header of column separated values."""
        return self.fields

    def getHeader(self):
        """returns header of column separated values."""
        return "\t".join(self.getHeaders())

    def __str__(self):
        """return string representation of data."""

        if self._mode == "int":
            format_vals = "%i"
            format_median = "%.1f"
        else:
            format_vals = self._format
            format_median = self._format

        return "\t".join(("%i" % self.counts,
                          format_vals % self.min,
                          format_vals % self.max,
                          self._format % self.mean,
                          format_median % self.median,
                          self._format % self.samplestd,
                          format_vals % self.sum,
                          format_vals % self.q1,
                          format_vals % self.q3,
                          ))


def adjustPValuesR(pvalues, method):
    '''adjust P-Values for multiple testing using
    the p.adjust() method in R.

    Possible values of method are:

    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    '''
    return R.p_adjust(pvalues, method)


def smoothPValues(pvalues,
                  vlambda=numpy.arange(0, 0.95, 0.05),
                  smooth_df=3,
                  smooth_log_pi0=False):

    if min(pvalues) < 0 or max(pvalues) > 1:
        raise ValueError("p-values out of range")

    if len(vlambda) > 1 and len(vlambda) < 4:
        raise ValueError(
            "if length of vlambda greater than 1, you need "
            "at least 4 values.")

    if len(vlambda) > 1 and (min(vlambda) < 0 or max(vlambda) >= 1):
        raise ValueError("vlambda must be within [0, 1).")

    m = len(pvalues)

    pi0 = numpy.zeros(len(vlambda), numpy.float)

    for i in range(len(vlambda)):
        pi0[i] = numpy.mean([x >= vlambda[i]
                             for x in pvalues]) / (1.0 - vlambda[i])

    R.assign("pi0", pi0)
    R.assign("vlambda", vlambda)
    print("pi0=", pi0)

    if smooth_log_pi0:
        pi0 = math.log(pi0)

    R.assign("smooth_df", smooth_df)

    spi0 = R("""spi0 <- smooth.spline(vlambda,pi0, df = smooth_df)""")
    pi0 = R("""pi0 <- predict( spi0, x = max(vlambda) )$y""")

    print(spi0)
    if smooth_log_pi0:
        pi0 = math.exp(pi0)

    return pi0


def getPi0(pvalues,
           vlambda=numpy.arange(0, 0.95, 0.05),
           pi0_method="smoother",
           smooth_df=3,
           smooth_log_pi0=False):
    '''used within nubiscan.'''

    if min(pvalues) < 0 or max(pvalues) > 1:
        raise ValueError("p-values out of range")

    if len(vlambda) > 1 and len(vlambda) < 4:
        raise ValueError(
            "if length of vlambda greater than 1, you "
            "need at least 4 values.")

    if len(vlambda) > 1 and (min(vlambda) < 0 or max(vlambda) >= 1):
        raise ValueError("vlambda must be within [0, 1).")

    m = len(pvalues)

    # these next few functions are the various ways to estimate pi0
    if len(vlambda) == 1:
        vlambda = vlambda[0]
        if vlambda < 0 or vlambda >= 1:
            raise ValueError("vlambda must be within [0, 1).")

        pi0 = numpy.mean([x >= vlambda for x in pvalues]) / (1.0 - vlambda)
        pi0 = min(pi0, 1.0)
        R.assign("pi0", pi0)

    else:

        pi0 = numpy.zeros(len(vlambda), numpy.float)

        for i in range(len(vlambda)):
            pi0[i] = numpy.mean([x >= vlambda[i]
                                 for x in pvalues]) / (1.0 - vlambda[i])

        R.assign("pi0", pi0)
        R.assign("vlambda", vlambda)

        if pi0_method == "smoother":
            if smooth_log_pi0:
                pi0 = math.log(pi0)

            R.assign("smooth_df", smooth_df)

            spi0 = R("""spi0 <- smooth.spline(vlambda,pi0, df = smooth_df)""")
            pi0 = R("""pi0 <- predict( spi0, x = max(vlambda) )$y""")

            if smooth_log_pi0:
                pi0 = math.exp(pi0)

        elif pi0_method == "bootstrap":

            minpi0 = min(pi0)

            mse = numpy.zeros(len(vlambda), numpy.float)
            pi0_boot = numpy.zeros(len(vlambda), numpy.float)

            R.assign("pvalues", pvalues)
            pi0 = R("""
            m <- length(pvalues)
            minpi0 <- min(pi0)
            mse <- rep(0,length(vlambda))
            pi0_boot <- rep(0,length(vlambda))
            for(i in 1:100)
            {
                pvalues_boot <- sample(pvalues,size=m,replace=TRUE)
                for(i in 1:length(vlambda))
                {
                    pi0_boot[i] <- mean(pvalues_boot>vlambda[i])/(1-vlambda[i])
                }
                mse <- mse + (pi0_boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])""")
        else:
            raise ValueError(
                "'pi0_method' must be one of 'smoother' or 'bootstrap'.")

        pi0 = min(pi0, 1.0)

    if pi0 <= 0:
        raise ValueError(
            "The estimated pi0 <= 0. Check that you have valid p-values "
            "or use another vlambda method.")

    return pi0


class FDRResult:

    def __init__(self):
        pass

    def plot(self, hardcopy=None):

        if hardcopy:
            R.png(hardcopy, width=1024, height=768, type="cairo")

        R.require('qvalue')

        # build a qobj
        R.assign("pval", self.mPValues)
        R.assign("pi0", self.mPi0)
        R.assign("qval", self.mQValues)
        R.assign("lambda", self.mLambda)
        R("""qobj <-list( pi0=pi0, qvalues=qval, pvalues=pval, lambda=lambda)""")
        R(""" class(qobj) <- "qvalue" """)

        R("""qplot(qobj)""")

        if hardcopy:
            R.dev_off()


def doFDR(pvalues,
          vlambda=None,
          pi0_method="smoother",
          fdr_level=None,
          robust=False,
          smooth_df=3,
          smooth_log_pi0=False,
          plot=False):
    """modeled after code taken from
http://genomics.princeton.edu/storeylab/qvalue/linux.html.

    I did not like the error handling so I translated most to python.

    Compute FDR after method by Storey et al. (2002).

    """

    # set to default of qvalue method
    if vlambda is None:
        vlambda = numpy.arange(0, 0.95, 0.05)

    if min(pvalues) < 0 or max(pvalues) > 1:
        raise ValueError("p-values out of range")

    if type(vlambda) == float:
        vlambda = (vlambda, )

    if len(vlambda) > 1 and len(vlambda) < 4:
        raise ValueError(
            "if length of vlambda greater than 1, "
            "you need at least 4 values.")

    if len(vlambda) > 1 and (min(vlambda) < 0 or max(vlambda) >= 1):
        raise ValueError("vlambda must be within [0, 1).")

    m = len(pvalues)

    # these next few functions are the various ways to estimate pi0
    if len(vlambda) == 1:
        vlambda = vlambda[0]
        if vlambda < 0 or vlambda >= 1:
            raise ValueError("vlambda must be within [0, 1).")

        pi0 = numpy.mean([x >= vlambda for x in pvalues]) / (1.0 - vlambda)
        pi0 = min(pi0, 1.0)
        R.assign("pi0", pi0)
    else:
        pi0 = numpy.zeros(len(vlambda), numpy.float)

        for i in range(len(vlambda)):
            pi0[i] = numpy.mean([x >= vlambda[i]
                                 for x in pvalues]) / (1.0 - vlambda[i])

        R.assign("pi0", pi0)
        R.assign("vlambda", vlambda)

        if pi0_method == "smoother":
            if smooth_log_pi0:
                pi0 = math.log(pi0)

            R.assign("smooth_df", smooth_df)
            spi0 = R("""spi0 <- smooth.spline(vlambda,pi0, df = smooth_df)""")
            if plot:
                import matplotlib.pyplot as plt
                plt.figure()
                plt.plot(vlambda, pi0)
                x2 = numpy.arange(0, 1, 0.001)
                R.assign("x2", x2)
                y2 = R("""y2 <- predict( spi0, x = x2 )$y""")
                plt.plot(x2, y2)
                plt.show()

            pi0 = R("""pi0 <- predict( spi0, x = max(vlambda) )$y""")[0]

            if smooth_log_pi0:
                pi0 = math.exp(pi0)

        elif pi0_method == "bootstrap":

            minpi0 = min(pi0)

            mse = numpy.zeros(len(vlambda), numpy.float)
            pi0_boot = numpy.zeros(len(vlambda), numpy.float)

            R.assign("pvalues", pvalues)
            pi0 = R("""
            m <- length(pvalues)
            minpi0 <- min(pi0)
            mse <- rep(0,length(vlambda))
            pi0_boot <- rep(0,length(vlambda))
            for(i in 1:100)
            {
                pvalues_boot <- sample(pvalues,size=m,replace=TRUE)
                for(i in 1:length(vlambda))
                {
                    pi0_boot[i] <- mean(pvalues_boot>vlambda[i])/(1-vlambda[i])
                }
                mse <- mse + (pi0_boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])""")[0]
        else:
            raise ValueError(
                "'pi0_method' must be one of 'smoother' or 'bootstrap'.")

        pi0 = min(pi0, 1.0)
        R.assign("pi0", pi0)

    if pi0 <= 0:
        raise ValueError(
            "The estimated pi0 (%f) <= 0. Check that you have valid p-values "
            "or use another vlambda method." % pi0)

    if fdr_level is not None and (fdr_level <= 0 or fdr_level > 1):
        raise ValueError("'fdr_level' must be within (0, 1].")

    # The estimated q-values calculated here
    # u = numpy.argsort( p )

    # change by Alan
    # ranking function which returns number of observations less than or equal
    ro.globalenv['pvalues'] = ro.FloatVector(pvalues)
    R.assign("robust", robust)
    qvalues = R("""u <- order(pvalues)
    qvalues.rank <- function(x)
{
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)

      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
}

v <- qvalues.rank(pvalues)
m <- length(pvalues)
qvalues <- pi0 * m * pvalues / v

if(robust)
{
        qvalues <- pi0*m*pvalues/(v*(1-(1-pvalues)^m))
}
qvalues[u[m]] <- min(qvalues[u[m]],1)

rqvalues <- qvalues
for(i in (m-1):1)
{
   qvalues[u[i]] <- min(qvalues[u[i]],qvalues[u[i+1]],1)
}


qvalues
""")

    result = FDRResult()
    result.mQValues = qvalues

    if fdr_level is not None:
        result.mPassed = [x <= fdr_level for x in result.mQValues]
    else:
        result.mPassed = [False for x in result.mQValues]

    result.mPValues = pvalues
    result.mPi0 = pi0
    result.mLambda = vlambda

    return result


def doFDRPython(pvalues,
                vlambda=None,
                pi0_method="smoother",
                fdr_level=None,
                robust=False,
                smooth_df=3,
                smooth_log_pi0=False,
                pi0=None,
                plot=False):
    """modeled after code taken from
    http://genomics.princeton.edu/storeylab/qvalue/linux.html.

    I did not like the error handling so I translated most to python.

    Compute FDR after method by Storey et al. (2002).

    """

    if min(pvalues) < 0 or max(pvalues) > 1:
        raise ValueError("p-values out of range")

    # set to default of qvalue method
    if vlambda is None:
        vlambda = numpy.arange(0, 0.95, 0.05)

    m = len(pvalues)
    pvalues = numpy.array(pvalues, dtype=numpy.float)

    if pi0 is None:
        if type(vlambda) == float:
            vlambda = (vlambda,)

        if len(vlambda) > 1 and len(vlambda) < 4:
            raise ValueError(
                "if length of vlambda greater than 1, you need at least 4 values.")

        if len(vlambda) > 1 and (min(vlambda) < 0 or max(vlambda) >= 1):
            raise ValueError("vlambda must be within [0, 1).")

        # estimate pi0
        if len(vlambda) == 1:
            vlambda = vlambda[0]
            if vlambda < 0 or vlambda >= 1:
                raise ValueError("vlambda must be within [0, 1).")

            pi0 = numpy.mean([x >= vlambda for x in pvalues]) / (1.0 - vlambda)
            pi0 = min(pi0, 1.0)
        else:

            pi0 = numpy.zeros(len(vlambda), numpy.float)

            for i in range(len(vlambda)):
                pi0[i] = numpy.mean([x >= vlambda[i]
                                     for x in pvalues]) / (1.0 - vlambda[i])

            if pi0_method == "smoother":

                if smooth_log_pi0:
                    pi0 = math.log(pi0)

                tck = scipy.interpolate.splrep(vlambda,
                                               pi0,
                                               k=smooth_df,
                                               s=10000)

                if plot:
                    import matplotlib.pyplot as plt
                    plt.figure()
                    plt.plot(vlambda, pi0)
                    x2 = numpy.arange(0, 1, 0.001)
                    y2 = scipy.interpolate.splev(x2, tck)
                    plt.plot(x2, y2)
                    plt.show()

                pi0 = scipy.interpolate.splev(max(vlambda), tck)
                if smooth_log_pi0:
                    pi0 = math.exp(pi0)

            elif pi0_method == "bootstrap":

                minpi0 = min(pi0)

                mse = numpy.zeros(len(vlambda), numpy.float)
                pi0_boot = numpy.zeros(len(vlambda), numpy.float)

                for i in range(100):
                    # sample pvalues
                    idx_boot = numpy.random.random_integers(0, m - 1, m)
                    pvalues_boot = pvalues[idx_boot]

                    for x in range(len(vlambda)):
                        # compute number of pvalues larger than lambda[x]
                        pi0_boot[x] = numpy.mean(
                            pvalues_boot > vlambda[x]) / (1.0 - vlambda[x])
                    mse += (pi0_boot - minpi0) ** 2
                pi0 = min(pi0[mse == min(mse)])
            else:
                raise ValueError(
                    "'pi0_method' must be one of 'smoother' or 'bootstrap'.")

            pi0 = min(pi0, 1.0)

    if pi0 <= 0:
        raise ValueError(
            "The estimated pi0 <= 0. Check that you have valid p-values "
            "or use another vlambda method.")

    if fdr_level is not None and (fdr_level <= 0 or fdr_level > 1):
        raise ValueError("'fdr_level' must be within (0, 1].")

    # compute qvalues
    idx = numpy.argsort(pvalues)
    # monotonically decreasing bins, so that bins[i-1] > x >=  bins[i]
    bins = numpy.unique(pvalues)[::-1]

    # v[i] = number of observations less than or equal to pvalue[i]
    # could this be done more elegantly?
    val2bin = len(bins) - numpy.digitize(pvalues, bins)
    v = numpy.zeros(m, dtype=numpy.int)
    lastbin = None
    for x in range(m - 1, -1, -1):
        bin = val2bin[idx[x]]
        if bin != lastbin:
            c = x
        v[idx[x]] = c + 1
        lastbin = bin

    qvalues = pvalues * pi0 * m / v
    if robust:
        qvalues /= (1.0 - (1.0 - pvalues) ** m)

    # bound qvalues by 1 and make them monotonic
    qvalues[idx[m - 1]] = min(qvalues[idx[m - 1]], 1.0)
    for i in range(m - 2, -1, -1):
        qvalues[idx[i]] = min(min(qvalues[idx[i]], qvalues[idx[i + 1]]), 1.0)

    result = FDRResult()
    result.mQValues = qvalues

    if fdr_level is not None:
        result.mPassed = [x <= fdr_level for x in result.mQValues]
    else:
        result.mPassed = [False for x in result.mQValues]

    result.mPValues = pvalues
    result.mPi0 = pi0
    result.mLambda = vlambda

    result.xvalues = qvalues

    return result


class CorrelationTest:

    '''coefficient is r, not r squared'''

    def __init__(self,
                 r_result=None,
                 s_result=None,
                 method=None):
        self.mPValue = None
        self.mMethod = None

        if r_result:
            self.mCoefficient = r_result['estimate']['cor']
            self.mPValue = float(r_result['p.value'])
            self.mNObservations = r_result['parameter']['df']
            self.mMethod = r_result['method']
            self.mAlternative = r_result['alternative']
        elif s_result:
            self.mCoefficient = s_result[0]
            self.mPValue = s_result[1]
            self.mNObservations = 0
            self.mAlternative = "two-sided"
        else:
            self.mCoefficient = 0
            self.mPValue = 1
            self.mSignificance = "na"
            self.mNObservations = 0
            self.mAlternative = "na"
            self.mMethod = "na"

        if method:
            self.mMethod = method

        if self.mPValue is not None:
            self.mSignificance = getSignificance(self.mPValue)

    def __str__(self):
        return "\t".join((
            "%6.4f" % self.mCoefficient,
            "%e" % self.mPValue,
            self.mSignificance,
            "%i" % self.mNObservations,
            self.mMethod,
            self.mAlternative))

    @classmethod
    def getHeaders(cls):
        return ("coeff", "pvalue", "significance", "observations",
                "method", "alternative")


def filterMasked(xvals, yvals, missing=("na", "Nan", None, ""),
                 dtype=numpy.float):
    """convert xvals and yvals to numpy array skipping pairs with
    one or more missing values."""
    xmask = [i in missing for i in xvals]
    ymask = [i in missing for i in yvals]
    return (numpy.array([xvals[i] for i in range(len(xvals))
                         if not xmask[i]], dtype=dtype),
            numpy.array([yvals[i] for i in range(len(yvals))
                         if not ymask[i]], dtype=dtype))


def doCorrelationTest(xvals, yvals):
    """compute correlation between x and y.

    Raises a value-error if there are not enough observations.
    """

    if len(xvals) <= 1 or len(yvals) <= 1:
        raise ValueError("can not compute correlation with no data")
    if len(xvals) != len(yvals):
        raise ValueError("data vectors have unequal length")

    x, y = filterMasked(xvals, yvals)

    result = CorrelationTest(s_result=scipy.stats.pearsonr(x, y),
                             method="pearson")
    result.mNObservations = len(x)

    return result


def getPooledVariance(data):
    """return pooled variance from a
    list of tuples (sample_size, variance)."""
    t, var = 0, 0

    for n, s in data:
        t += n
        var += (n - 1) * s

    assert t > len(data), "sample size smaller than samples combined"

    return var / float(t - len(data))


def computeROC(values):
    '''return a roc curve for *values*. Values
    is a sorted list of (value, bool) pairs.

    Deprecated - use getPerformance instead

    returns a list of (FPR,TPR) tuples.
    '''
    roc = []

    npositives = len([x for x in values if x[1]])
    if npositives == 0:
        raise ValueError("no positives among values")

    ntotal = len(values)

    last_value, last_fpr = None, None
    tp, fp = 0, 0
    tn, fn = ntotal - npositives, npositives

    for value, is_positive in values:
        if is_positive:
            tp += 1
            fn -= 1
        else:
            fp += 1
            tn -= 1

        if last_value != value:

            try:
                tpr = float(tp) / (tp + fn)
            except ZeroDivisionError:
                tpr = 0

            try:
                fpr = float(fp) / (fp + tn)
            except ZeroDivisionError:
                fpr = 0

            if last_fpr != fpr:
                roc.append((fpr, tpr))
                last_fpr = fpr

        last_values = value

    return roc


class TTest:

    def __init__(self):
        pass


class WelchTTest:

    def __init__(self):
        pass

PairedTTest = collections.namedtuple("PairedTTest", "statistic pvalue")


def doPairedTTest(vals1, vals2):
    '''perform paired t-test.

    vals1 and vals2 need to contain the same number of elements.
    '''
    return PairedTTest._make(scipy.stats.ttest_rel(vals1, vals2))


def doWelchsTTest(n1, mean1, std1,
                  n2, mean2, std2,
                  alpha=0.05):
    '''Welch''s approximate t-test for the difference of two means of
    heteroscedasctic populations.

    This functions does a two-tailed test.

    see PMID: 12016052

    :Parameters:
        n1 : int
            number of variates in sample 1
        n2 : int
            number of variates in sample 2
        mean1 : float
            mean of sample 1
        mean2 : float
            mean of sample 2
        std1 : float
            standard deviation of sample 1
        std2 : float
            standard deviation of sample 2

    returns a WelchTTest
    '''
    if std1 == 0 and std2 == 0:
        raise ValueError('standard deviations are 0.')

    # convert standard deviation to sample variance
    svar1 = std1 ** 2 * n1 / float(n1 - 1)
    svar2 = std2 ** 2 * n2 / float(n2 - 1)

    # compute df and test statistic
    df = ((svar1 / n1 + svar2 / n2) ** 2) / \
        (((svar1 / n1) ** 2) / (n1 - 1) + ((svar2 / n2) ** 2) / (n2 - 1))
    denom = numpy.sqrt(svar1 / n1 + svar2 / n2)
    z = abs(mean1 - mean2) / denom

    # do the test
    pvalue = 2 * scipy.stats.t.sf(z, df)
    result = WelchTTest()
    result.mPValue = pvalue
    result.mDegreesFreedom = df
    result.mZ = z
    result.mMean1 = mean1
    result.mMean2 = mean2
    result.mSampleVariance1 = svar1
    result.mSampleVariance2 = svar2
    result.mDifference = mean1 - mean2
    result.mZLower = scipy.stats.t.ppf(alpha, df)
    result.mZUpper = scipy.stats.t.ppf(1.0 - alpha, df)
    result.mDifferenceLower = result.mZLower * denom
    result.mDifferenceUpper = result.mZUpper * denom

    return result


def getAreaUnderCurve(xvalues, yvalues):
    '''compute area under curve from a set of discrete x,y coordinates
    using trapezoids.

    This is only as accurate as the density of points.
    '''

    assert len(xvalues) == len(yvalues)
    last_x, last_y = xvalues[0], yvalues[0]
    auc = 0
    for x, y in zip(xvalues, yvalues)[1:]:
        dx = x - last_x
        assert not dx <= 0, "x not increasing: %f >= %f" % (last_x, x)
        dy = abs(last_y - y)
        my = min(last_y, y)
        # rectangle plus triangle
        auc += dx * my + dx * dy / 2
        last_x, last_y = x, y

    return auc


def getSensitivityRecall(values):
    '''return sensitivity/selectivity.

    Values is a sorted list of (value, bool) pairs.

    Deprecated - use getPerformance instead
    '''

    npositives = 0.0
    npredicted = 0.0
    l = None
    result = []
    total = float(len(values))
    for value, is_positive in values:
        npredicted += 1.0
        if is_positive > 0:
            npositives += 1.0
        if value != l:
            result.append((value, npositives / npredicted, npredicted / total))
        l = value
    if l:
        result.append((l, npositives / npredicted, npredicted / total))

    return result

ROCResult = collections.namedtuple("ROCResult",
                                   "value pred tp fp tn fn tpr fpr tnr fnr rtpr rfnr")


def getPerformance(values,
                   skip_redundant=True,
                   false_negatives=False,
                   bin_by_value=True,
                   monotonous=False,
                   multiple=False,
                   increasing=True,
                   total_positives=None,
                   total_false_negatives=None,
                   ):
    '''compute performance estimates for a list of ``(score, flag)``
    tuples in *values*.

    Values is a sorted list of (value, bool) pairs.

    If the option *false-negative* is set, the input is +/- or 1/0 for
    a true positive or false negative, respectively.

    TP: true positives
    FP: false positives
    TPR: true positive rate  = true_positives /  predicted
    P: predicted
    FPR: false positive rate = false positives  / predicted
    value: value

    '''

    true_positives = 0
    predicted = 0

    last_value = None

    binned_values = []

    for value, flag in values:
        if not bin_by_value:
            if last_value != value:
                binned_values.append((true_positives, predicted, value))
        else:
            if last_value is not None and last_value != value:
                binned_values.append((true_positives, predicted, last_value))

        predicted += 1

        if flag:
            true_positives += 1

        last_value = value

    binned_values.append((true_positives, predicted, last_value))
    binned_values.append((true_positives, predicted, value))

    if true_positives == 0:
        raise ValueError("# no true positives!")

    if total_positives is None:
        if total_false_negatives:
            positives = float(predicted)
        else:
            positives = float(true_positives)
    else:
        positives = float(total_positives)

    last_positives = None
    last_tpr = None
    result = []

    for true_positives, predicted, value in binned_values:

        if (predicted == 0):
            predicted = 1

        if total_false_negatives:
            false_negatives = predicted - true_positives
            false_positives = 0
            true_negatives = 0
        else:
            true_negatives = 0
            false_negatives = positives - true_positives
            false_positives = predicted - true_positives

        tpr = float(true_positives) / predicted
        fpr = float(false_positives) / (true_positives + false_negatives)
        fnr = float(false_negatives) / positives
        tnr = 0

        # relative rates
        rfpr = float(false_positives) / predicted
        rfnr = float(false_negatives) / predicted

        if monotonous and last_tpr and last_tpr < tpr:
            continue

        if skip_redundant and true_positives == last_positives:
            continue

        if (predicted > 0):
            result.append(ROCResult._make(
                (value,
                 predicted,
                 true_positives,
                 false_positives,
                 true_negatives,
                 false_negatives,
                 tpr, fpr, tnr, fnr,
                 rfpr, rfnr)))

        last_positives = true_positives
        last_tpr = tpr

    return result


def doMannWhitneyUTest(xvals, yvals):
    '''apply the Mann-Whitney U test to test for the difference of medians.'''

    r_result = R.wilcox_test(xvals, yvals, paired=False)

    result = Result().fromR(
        (("pvalue", 'p.value'),
         ('alternative', None),
         ('method', None)),
        r_result)

    return result


def adjustPValues(pvalues, method='fdr', n=None):
    '''returns an array of adjusted pvalues

    Reimplementation of p.adjust in the R package.

    p: numeric vector of p-values (possibly with 'NA's).  Any other
    R is coerced by 'as.numeric'.

    method: correction method. Valid values are:

    n: number of comparisons, must be at least 'length(p)'; only set
    this (to non-default) when you know what you are doing

    For more information, see the documentation of the
    p.adjust method in R.
    '''

    if n is None:
        n = len(pvalues)

    if method == "fdr":
        method = "BH"

    # optional, remove NA values
    p = numpy.array(pvalues, dtype=numpy.float)
    lp = len(p)

    assert n <= lp

    if n <= 1:
        return p
    if n == 2 and method == "hommel":
        method = "hochberg"
    if method == "bonferroni":
        p0 = n * p
    elif method == "holm":
        i = numpy.arange(lp)
        o = numpy.argsort(p)
        ro = numpy.argsort(o)
        m = numpy.maximum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "hommel":
        raise NotImplementedError("hommel method not implemented")
    # if (n > lp) p <- c(p, rep.int(1, n - lp))
    #    i = numpy.arange(n)
    #    o = numpy.argsort(p)
    #    p = p[o]
    #    ro = numpy.argsort(o)
    #
    #     q <- pa <- rep.int(min(n * p/i), n)
    #     for (j in (n - 1):2) {
    #         ij <- seq_len(n - j + 1)
    #         i2 <- (n - j + 2):n
    #         q1 <- min(j * p[i2]/(2:j))
    #         q[ij] <- pmin(j * p[ij], q1)
    #         q[i2] <- q[n - j + 1]
    #         pa <- pmax(pa, q)
    #     }
    #     pmax(pa, p)[if (lp < n) ro[1:lp] else ro]
    elif method == "hochberg":
        i = numpy.arange(0, lp)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        m = numpy.minimum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "BH":
        i = numpy.arange(1, lp + 1)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        m = numpy.minimum.accumulate(float(n) / i * p[o])
        p0 = m[ro]
    elif method == "BY":
        i = numpy.arange(1, lp + 1)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        q = numpy.sum(1.0 / numpy.arange(1, n + 1))
        m = numpy.minimum.accumulate(q * float(n) / i * p[o])
        p0 = m[ro]
    elif method == "none":
        p0 = p

    return numpy.minimum(p0, numpy.ones(len(p0)))


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.

    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only
        smoothing)

    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).

    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """

    try:
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order + 1))
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = numpy.mat([[k**i for i in order_range]
                   for k in range(-half_window, half_window + 1)])
    m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * math.factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - numpy.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + numpy.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))
    return numpy.convolve(m[::-1], y, mode='valid')
