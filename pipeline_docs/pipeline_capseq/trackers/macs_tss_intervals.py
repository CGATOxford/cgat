import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class tssIntervalSummary( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getFirstRow( "SELECT COUNT(t.gene_id) as tss FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getFirstRow( "SELECT COUNT(t.gene_id) as tss FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        return odict( zip( ("TSS intervals", "Non TSS intervals" ), data1 + data2) )

##################################################################################
class tssIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT i.end-i.start as length FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT i.end-i.start as length FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT i.peakval FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT i.peakval FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT i.avgval FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT i.avgval FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT i.fold FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT i.fold FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        #return odict( [("TSS interval length", data1),("Non-TSS interval length",data2)] )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }


##################################################################################
class tssIntervalCpGDensity( cpgTracker ):

    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT c.pCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT c.pCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalCpGObsExp1( cpgTracker ):
    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT c.CpG_ObsExp1 FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT c.CpG_ObsExp1 FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalCpGObsExp2( cpgTracker ):
    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT c.CpG_ObsExp2 FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT c.CpG_ObsExp2 FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalCpGNumber( cpgTracker ):
    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT c.nCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT c.nCpG FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }

##################################################################################
class tssIntervalGCContent( cpgTracker ):
    mPattern = "_replicated_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT c.pGC FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT c.pGC FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss t, %(track)s_replicated_composition c WHERE i.interval_id=t.gene_id AND c.gene_id=i.interval_id AND t.closest_dist >1000" % locals() )
        return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }


