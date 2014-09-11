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


class IntervalsSummary(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice=None):
        data = self.getRow(
            "SELECT COUNT(*) as Intervals, round(AVG(length),0) as Mean_length, round(AVG(nprobes),0) as Mean_reads FROM %(track)s_macs_intervals" % locals())
        return data

##########################################################################


class IntervalsSummaryFiltered(cpgTracker):

    """Summary stats of intervals after filtering by fold change and merging nearby intervals. """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data = self.getRow(
            "SELECT COUNT(*) as Intervals, round(AVG(length),0) as Mean_length, round(AVG(nprobes),0) as Mean_reads FROM %(track)s_macs_merged_intervals" % locals())
        return data

##########################################################################


class IntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT length FROM %(track)s_macs_merged_intervals" % locals())
        return data

##########################################################################


class IntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT peakval FROM %(track)s_macs_merged_intervals" % locals())
        return data

##########################################################################


class IntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT avgval FROM %(track)s_macs_merged_intervals" % locals())
        return data

##########################################################################


class IntervalFoldChange(cpgTracker):

    """return fold changes for all intervals. """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT fold FROM %(track)s_macs_merged_intervals" % locals())
        return data

##########################################################################
##########################################################################
##########################################################################


class PeakLocation(cpgTracker):
    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_macs_merged_intervals" % locals())
        data2 = self.getValues(
            "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_macs_merged_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################


class PeakDistance(cpgTracker):
    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT PeakCenter - start FROM %(track)s_macs_merged_intervals" % locals())
        data2 = self.getValues(
            "SELECT end - PeakCenter FROM %(track)s_macs_merged_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################
##########################################################################
##########################################################################


class CpGDensity(cpgTracker):
    mPattern = "_composition$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT pCpG FROM %(track)s_composition" % locals())
        return data

##########################################################################


class CpGObsExp1(cpgTracker):
    pattern = "(?<!replicated)_composition$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT CpG_ObsExp1 FROM %(track)s_composition" % locals())
        return data

##########################################################################


class CpGObsExp2(cpgTracker):
    mPattern = "_composition$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT CpG_ObsExp2 FROM %(track)s_composition" % locals())
        return data

##########################################################################


class CpGNumber(cpgTracker):
    mPattern = "_composition$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT nCpG FROM %(track)s_composition" % locals())
        return data

##########################################################################


class GCContent(cpgTracker):
    mPattern = "_composition$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT pGC FROM %(track)s_composition" % locals())
        return data
