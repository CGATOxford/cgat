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


class ExternalIntervalLists(cpgTracker):

    """Summary stats of external interval lists used for comparison. """

    mPattern = "external_interval_sets$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT bed, intervals FROM external_interval_sets" % locals())
        return odict(zip(("Dataset", "Intervals"), zip(*data)))

##########################################################################


class IntervalOverlapFoldChange(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_foldchange_shared$"

    def __call__(self, track, slice=None):
        data1 = self.get(
            "SELECT track as comparison, threshold, intervals FROM %(track)s_foldchange_shared order by track, threshold asc" % locals())
        data2 = self.get(
            "SELECT '%(track)s' as comparison, threshold, intervals FROM %(track)s_foldchange order by threshold asc" % locals())

        result = odict()

        for d in data1:
            result[d[0]] = {'threshold': (), 'intervals': ()}
        for d in data2:
            result[d[0]] = {'threshold': (), 'intervals': ()}

        for d in data1:
            result[d[0]]['threshold'] += (d[1],)
            result[d[0]]['intervals'] += (d[2],)

        for d in data2:
            result[track]['threshold'] += (d[1],)
            result[track]['intervals'] += (d[2],)
        return result

##########################################################################


class UniqueIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_merged_unique_intervals" % locals())
        return odict(zip(("Unique intervals", "mean_interval_length"), data))

##########################################################################


class SharedIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_merged_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_merged_shared_intervals" % locals())
        return odict(zip(("Shared intervals", "mean_interval_length"), data))

##########################################################################


class OverlapIntervals(cpgTracker):

    """Count of overlapping intervals for each dataset. """

    mPattern = "_merged_overlap$"
    #

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT track, overlap FROM %(track)s_merged_overlap" % locals())
        return odict(zip(("Track", "Intervals"), zip(*data)))

##########################################################################


class OverlapCpG(cpgTracker):

    """Count of intervals overlapping CpG islands for each dataset. """

    mPattern = "_merged_cgi$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_merged_cgi c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_macs_merged_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), zip(*data)))

##########################################################################


class OverlapChipseq(cpgTracker):

    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_merged_chipseq$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track,b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_merged_chipseq c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_macs_merged_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), zip(*data)))

##########################################################################


class OverlapCAPseq(cpgTracker):

    """Count of intervals overlapping CAPseq intervals for each dataset. """

    mPattern = "_merged_capseq$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_merged_capseq c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_macs_merged_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), zip(*data)))

##########################################################################


class OverlapChromatinMarks(cpgTracker):

    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_merged_chromatin$"
    #

    def __call__(self, track, slice=None):
        data = self.get( """SELECT c.track, b.a_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(b.a_intervals+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_merged_chromatin c, external_interval_sets e, (select count(*) as a_intervals from %(track)s_macs_merged_intervals) b 
                            WHERE c.track=e.bed""" % locals() )

        return odict(zip(("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B"), zip(*data)))

##########################################################################


class gatResults(cpgTracker):

    """Summary stats of GAT analysis. """

    mPattern = "gat_results$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM external_dataset_gat_results ")
        return odict(zip(("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value"), zip(*data)))
