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


class MacsSummary(cpgTracker):
    pattern = "(macs_summary)"

    def getTracks(self, subset=None):
        return self.getValues("SELECT track FROM macs_summary ORDER BY track")

    def __call__(self, track, slice=None):

        resultsdir = os.path.join(EXPORTDIR, "MACS")

        fields = ("ncandidates_positive", "ncandidates_negative",
                  "called_positive", "called_negative", "min_tags",
                  "paired_peaks", "shift")

        f = ",".join(fields)
        data = self.getFirstRow(
            '''SELECT %(f)s FROM macs_summary WHERE track="%(track)s"''' % locals())
        result = odict(zip(fields, data))

        if os.path.exists(resultsdir):
            print resultsdir
            result[
                "link"] = "`pdf <%(resultsdir)s/%(track)s.macs_model.pdf>`_" % locals()
        return result

##########################################################################


class MacsSoloSummary(cpgTracker):
    pattern = "(macs_solo_summary)"

    def getTracks(self, subset=None):
        return self.getValues("SELECT track FROM macs_solo_summary ORDER BY track")

    def __call__(self, track, slice=None):

        resultsdir = os.path.join(EXPORTDIR, "MACS")

        fields = ("ncandidates_positive",
                  "called_positive", "min_tags",
                  "paired_peaks", "shift")

        f = ",".join(fields)
        data = self.getFirstRow(
            '''SELECT %(f)s FROM macs_solo_summary WHERE track="%(track)s"''' % locals())
        result = odict(zip(fields, data))

        if os.path.exists(resultsdir):
            print resultsdir
            result[
                "link"] = "`pdf <%(resultsdir)s/%(track)s.macs_model.pdf>`_" % locals()
        return result

##########################################################################


class MacsIntervalsSummary(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*), round(AVG(length),0), round(AVG(nprobes),0)  FROM %(track)s_macs_intervals" % locals())
        return odict(zip(("intervals_count", "mean_interval_length", "mean_reads_per_interval"), data))

##########################################################################


class FoldChangeThreshold(cpgTracker):

    """Count of intervals exceeding fold change threshold for each dataset. """

    mPattern = "_foldchange$"
    #

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT threshold, intervals FROM %(track)s_foldchange" % locals())
        return odict(zip(("Threshold", "Intervals"), zip(*data)))

##########################################################################


class BackgroundSummary(cpgTracker):

    """Summary stats of reads mapping inside/outside binding intervals. """

    mPattern = "_background$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow( '''select in_peaks, out_peaks, (in_peaks+out_peaks) as total, 
                                    round(((in_peaks+0.00)/(out_peaks+in_peaks+0.00))*100,1) as ratio,
                                    round(((out_peaks+0.00)/(out_peaks+in_peaks+0.00))*100,1) as ratio2 
                                    from %(track)s_background;''' % locals() )
        return odict(zip(("Reads Overlapping intervals", "Reads Outwith Intervals", "Total Reads", "Percent in Intervals", "Percent Background"), data))

##########################################################################


class MacsDiagnostics(cpgTracker):

    """Closest distance of transcript models to gene models in the reference set."""

    pattern = "(.*)_macsdiag"

    def __call__(self, track, slice=None):

        data = self.get(
            "SELECT fc,npeaks,p20,p30,p40,p50,p60,p70,p80,p90 FROM %(track)s_macsdiag" % locals())

        result = odict()
        for fc, npeaks, p20, p30, p40, p50, p60, p70, p80, p90 in data:
            result[fc] = odict()
            result[fc]["npeaks"] = npeaks
            result[fc]["proportion of reads"] = range(20, 100, 10)
            result[fc]["proportion of peaks"] = map(
                float, (p20, p30, p40, p50, p60, p70, p80, p90))

        return result
