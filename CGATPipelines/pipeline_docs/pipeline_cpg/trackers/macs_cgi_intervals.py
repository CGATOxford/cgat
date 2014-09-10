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


class cgiIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_cgi_cap_bed$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_cgi_cap_bed" % locals())
        return odict(zip(("CGI intervals", "mean_interval_length"), data))

##########################################################################


class cgiIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT (stop-start) FROM %(track)s_cgi_cap_bed" % locals())
        return {"length": data}

##########################################################################


class cgiIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_cgi_cap_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_cgi_cap_bed u, %(track_base)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"peakval": data}

##########################################################################


class cgiIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getValues( '''SELECT i.avgval FROM %(track)s_cgi_bed u, %(track_base)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"avgval": data}

##########################################################################


class cgiIntervalFoldChange(cpgTracker):

    """Distribution of fold change """

    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getValues( '''SELECT fold FROM %(track)s_cgi_bed u, %(track_base)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"Fold Change": data}

##########################################################################


class cgiIntervalTSS(cpgTracker):

    """Distribution of distance to closest TSS """

    mPattern = "_cgi_cap_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getValues( '''SELECT closest_dist FROM %(track)s_cgi_cap_bed u, 
                                  %(track_base)s_macs_intervals i, %(track_base)s_tss t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND t.gene_id=i.interval_id''' % locals() )
        return {"distance": data}

##########################################################################


class cgiIntervalCpGDensity(cpgTracker):
    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getAll( '''SELECT pCpG FROM %(track)s_cgi_bed u, 
                               %(track_base)s_macs_intervals i,%(track_base)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class cgiIntervalCpGObsExp1(cpgTracker):
    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getAll( '''SELECT CpG_ObsExp1 FROM %(track)s_cgi_bed u, 
                               %(track_base)s_macs_intervals i,%(track_base)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class cgiIntervalCpGObsExp2(cpgTracker):
    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        track_base = track.replace("_non", "").replace("_pred", "")
        data = self.getAll( '''SELECT CpG_ObsExp2 FROM %(track)s_cgi_bed u, 
                               %(track_base)s_macs_intervals i,%(track_base)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class cgiIntervalCpGNumber(cpgTracker):
    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_cgi_bed u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class cgiIntervalGCContent(cpgTracker):
    mPattern = "_cgi_bed$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_cgi_bed u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data
