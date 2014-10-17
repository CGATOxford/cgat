import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import scipy.stats
import numpy.ma
import Stats
import Histogram

from CGATReport.Tracker import *
from cpgReport import *

##########################################################################


class zinbaSummary(cpgTracker):
    mPattern = "_zinba_intervals$"

    def __call__(self, track, slice=None):

        data = self.getAll(
            '''SELECT count(*) as intervals FROM %(track)s_zinba_intervals''' % locals())
        return data

##########################################################################


class OverlapZinbaMacs(cpgTracker):

    """Count of intervals overlapping MACS intervals for each dataset. """

    mPattern = "_macs_zinba_intervals_shared$"

    def __call__(self, track, slice=None):
        data1 = self.getValue(
            "SELECT count(*) as intervals FROM %(track)s_macs_zinba_intervals_shared" % locals())
        data2 = self.getValue(
            "SELECT count(*) as intervals FROM %(track)s_macs_intervals" % locals())
        data3 = self.getValue(
            "SELECT count(*) as intervals FROM %(track)s_zinba_intervals" % locals())

        return dict(Overlap=data1, MACS=data2, ZINBA=data3)

##########################################################################


class zinbaIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_zinba_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT pstop-pstart as length FROM %(track)s_zinba_intervals" % locals())
        return data

##########################################################################


class zinbaMaxVal(cpgTracker):

    """return median coverage for all intervals. """

    mPattern = "_zinba_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT maxval FROM %(track)s_zinba_intervals" % locals())
        return data
