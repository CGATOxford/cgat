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


class cgiIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    pattern = "(.*_replicated.*_predicted_cgi.*)"
    as_table = True

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s" % locals())
        return odict(zip(("CGI intervals", "mean_interval_length"), data))

##########################################################################


class cgiIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT (stop-start) FROM %(track)s_replicated_predicted_cgi_and_cap" % locals())
        data2 = self.getValues(
            "SELECT (stop-start) FROM %(track)s_replicated_cap_not_predicted_cgi" % locals())
        data3 = self.getValues(
            "SELECT (stop-start) FROM %(track)s_replicated_predicted_cgi_not_cap" % locals())
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}


##########################################################################
class cgiIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getValues( '''SELECT peakval FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        data2 = self.getValues( '''SELECT peakval FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        data3 = self.getValues(
            '''SELECT peakval FROM %(track)s_replicated_predicted_cgi_not_cap''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}

##########################################################################


class cgiIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getValues( '''SELECT avgval FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        data2 = self.getValues( '''SELECT avgval FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        data3 = self.getValues(
            '''SELECT avgval FROM %(track)s_replicated_predicted_cgi_not_cap''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}

##########################################################################


class cgiIntervalFoldChange(cpgTracker):

    """Distribution of fold change """
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getValues( '''SELECT fold FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        data2 = self.getValues( '''SELECT fold FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                                  %(track)s_replicated_intervals i
                                  WHERE u.interval_id=i.interval_id''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2}

##########################################################################


class cgiIntervalTSS(cpgTracker):

    """Distribution of distance to closest TSS """
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getValues( '''SELECT closest_dist FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                                  %(track)s_replicated_intervals i, %(track)s_replicated_tss t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND t.gene_id=i.interval_id''' % locals() )
        data2 = self.getValues( '''SELECT closest_dist FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                                  %(track)s_replicated_intervals i, %(track)s_replicated_tss t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND t.gene_id=i.interval_id''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2}

##########################################################################


class cgiIntervalCpGDensity(cpgTracker):
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getAll( '''SELECT pCpG FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data2 = self.getAll( '''SELECT pCpG FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data3 = self.getAll( '''SELECT pCpG FROM %(track)s_replicated_predicted_cgi_not_cap u, cgi_comp c
                               WHERE c.gene_id=u.interval_id''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}

##########################################################################


class cgiIntervalCpGObsExp(cpgTracker):
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getAll( '''SELECT CpG_ObsExp FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data2 = self.getAll( '''SELECT CpG_ObsExp FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data3 = self.getAll( '''SELECT CpG_ObsExp FROM %(track)s_replicated_predicted_cgi_not_cap u, cgi_comp c
                               WHERE c.gene_id=u.interval_id''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}

##########################################################################


class cgiIntervalGCContent(cpgTracker):
    pattern = "(.*)_replicated_predicted_cgi_and_cap"

    def __call__(self, track, slice=None):
        data1 = self.getAll( '''SELECT pGC FROM %(track)s_replicated_predicted_cgi_and_cap u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data2 = self.getAll( '''SELECT pGC FROM %(track)s_replicated_cap_not_predicted_cgi u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        data3 = self.getAll( '''SELECT pGC FROM %(track)s_replicated_predicted_cgi_not_cap u, cgi_comp c
                               WHERE c.gene_id=u.interval_id''' % locals() )
        return {"Predicted CGI & CAPseq": data1, "CAPseq not Predicted CGI": data2, "Predicted CGI not CAPseq": data3}
