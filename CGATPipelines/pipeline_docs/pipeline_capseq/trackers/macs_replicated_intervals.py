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

from SphinxReport.Tracker import *
from cpgReport import *

##########################################################################


class replicatedIntervalSummary(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data = self.getRow(
            "SELECT COUNT(*) as Intervals, round(AVG(length),0) as Mean_length, round(AVG(nprobes),0) as Mean_reads FROM %(track)s_replicated_intervals" % locals())
        return data

##########################################################################


class replicatedIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT length FROM %(track)s_replicated_intervals" % locals())
        return data

##########################################################################


class replicatedIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT peakval FROM %(track)s_replicated_intervals" % locals())
        return data

##########################################################################


class replicatedIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT avgval FROM %(track)s_replicated_intervals" % locals())
        return data

##########################################################################


class replicatedIntervalFoldChange(cpgTracker):

    """return fold changes for all intervals. """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT fold FROM %(track)s_replicated_intervals" % locals())
        return data

##########################################################################
##########################################################################
##########################################################################


class replicatedIntervalPeakLocation(cpgTracker):
    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_replicated_intervals" % locals())
        data2 = self.getValues(
            "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_replicated_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################


class replicatedIntervalPeakDistance(cpgTracker):
    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT PeakCenter - start FROM %(track)s_replicated_intervals" % locals())
        data2 = self.getValues(
            "SELECT end - PeakCenter FROM %(track)s_replicated_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################
##########################################################################
##########################################################################


class replicatedIntervalCpGDensity(cpgTracker):
    pattern = "(.*)_replicated_composition"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT pCpG FROM %(track)s_replicated_composition" % locals())
        data2 = self.getValues(
            "SELECT pCpG FROM %(track)s_replicated_composition_control" % locals())
        data3 = self.getValues(
            "SELECT pCpG FROM %(track)s_replicated_composition_flanking5" % locals())
        data4 = self.getValues(
            "SELECT pCpG FROM %(track)s_replicated_composition_flanking3" % locals())
        return odict(zip(("CAPseq composition", "Control composition", "5` Flank Composition", "3` Flank Composition"), (data1, data2, data3, data4)))

##########################################################################


class replicatedIntervalCpGObsExp1(cpgTracker):
    pattern = "(.*)_replicated_composition"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT CpG_ObsExp1 FROM %(track)s_replicated_composition" % locals())
        data2 = self.getValues(
            "SELECT CpG_ObsExp1 FROM %(track)s_replicated_composition_control" % locals())
        data3 = self.getValues(
            "SELECT CpG_ObsExp1 FROM %(track)s_replicated_composition_flanking5" % locals())
        data4 = self.getValues(
            "SELECT CpG_ObsExp1 FROM %(track)s_replicated_composition_flanking3" % locals())
        return odict(zip(("CAPseq composition", "Control composition", "5` Flank Composition", "3` Flank Composition"), (data1, data2, data3, data4)))

##########################################################################


class replicatedIntervalCpGObsExp2(cpgTracker):
    pattern = "(.*)_replicated_composition"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT CpG_ObsExp FROM %(track)s_replicated_composition" % locals())
        data2 = self.getValues(
            "SELECT CpG_ObsExp FROM %(track)s_replicated_composition_control" % locals())
        data3 = self.getValues(
            "SELECT CpG_ObsExp FROM %(track)s_replicated_composition_flanking5" % locals())
        data4 = self.getValues(
            "SELECT CpG_ObsExp FROM %(track)s_replicated_composition_flanking3" % locals())
        return odict(zip(("CAPseq composition", "Control composition", "5` Flank Composition", "3` Flank Composition"), (data1, data2, data3, data4)))

##########################################################################


class replicatedIntervalGCContent(cpgTracker):
    pattern = "(.*)_replicated_composition"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT pGC FROM %(track)s_replicated_composition" % locals())
        data2 = self.getValues(
            "SELECT pGC FROM %(track)s_replicated_composition_control" % locals())
        data3 = self.getValues(
            "SELECT pGC FROM %(track)s_replicated_composition_flanking5" % locals())
        data4 = self.getValues(
            "SELECT pGC FROM %(track)s_replicated_composition_flanking3" % locals())
        return odict(zip(("CAPseq composition", "Control composition", "5` Flank Composition", "3` Flank Composition"), (data1, data2, data3, data4)))
