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


class SicerSummary(cpgTracker):
    pattern = "(sicer_summary)"

    def getTracks(self, subset=None):
        return self.getValues("SELECT track FROM sicer_summary ORDER BY track")

    def __call__(self, track, slice=None):

        fields = ("window_mean", "window_min", "score_threshold", "total_islands", "chip_library_size",
                  "control_library_size", "chip_island_reads", "control_island_reads", "significant_islands")
        f = ",".join(fields)
        data = self.getFirstRow(
            '''SELECT %(f)s FROM sicer_summary WHERE track="%(track)s"''' % locals())
        result = odict(zip(fields, data))
        return result

##########################################################################


class OverlapMacs(cpgTracker):

    """Count of intervals overlapping MACS intervals for each dataset. """

    mPattern = "_macs_sicer_intervals_shared$"
    #

    def __call__(self, track, slice=None):
        data1 = self.get(
            "SELECT count(*) as intervals FROM %(track)s_macs_sicer_intervals_shared" % locals())
        data2 = self.get(
            "SELECT count(*) as intervals FROM %(track)s_macs_intervals" % locals())
        data3 = self.get(
            "SELECT count(*) as intervals FROM %(track)s_sicer_intervals" % locals())

        ol = str(data1[0]).replace("(", "").replace(")", "").replace(",", "")
        macs = str(data2[0]).replace("(", "").replace(")", "").replace(",", "")
        sicer = str(data3[0]).replace(
            "(", "").replace(")", "").replace(",", "")
        return dict(Overlap=ol, MACS=macs, SICER=sicer)
# return odict( zip( ("Overlap", "MACS", "SICER" ), (data1 + data2 +
# data3)) )

##########################################################################


class SicerIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_sicer_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT stop-start as length FROM %(track)s_sicer_intervals" % locals())
        return data

##########################################################################


class SicerFoldChange(cpgTracker):

    """return fold changes for all intervals. """

    mPattern = "_sicer_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT fold FROM %(track)s_sicer_intervals" % locals())
        return data
