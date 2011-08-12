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

    mPattern = "_tss$"

    def __call__(self, track, slice = None):
        data1 = self.getFirstRow( "SELECT COUNT(t.gene_id) as tss FROM %(track)s_macs_intervals i, %(track)s_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getFirstRow( "SELECT COUNT(t.gene_id) as tss FROM %(track)s_macs_intervals i, %(track)s_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        return odict( zip( ("TSS intervals", "Non TSS intervals" ), data1 + data2) )

##################################################################################
class tssIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_tss$"

    #tracks= ("<=1000",">1000")

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT i.end-i.start as length FROM %(track)s_macs_intervals i, %(track)s_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist <=1000" % locals() )
        data2 = self.getValues( "SELECT i.end-i.start as length FROM %(track)s_macs_intervals i, %(track)s_tss t WHERE i.interval_id=t.gene_id AND t.closest_dist >1000" % locals() )
        #return { "TSS interval length" : data1, "Non-TSS interval length" : data2 }
        return (("TSS interval length", data1), ("Non-TSS interval length", data2))

##################################################################################
class tssIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_shared_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "peakval" : data }

##################################################################################
class tssIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_shared_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "avgval" : data }

##################################################################################
class tssIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT fold FROM %(track)s_shared_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "Fold Change" : data }

##################################################################################
class tssIntervalCpGDensity( cpgTracker ):
    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pCpG FROM %(track)s_shared_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class tssIntervalCpGObsExp1( cpgTracker ):
    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp1 FROM %(track)s_shared_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class tssIntervalCpGObsExp2( cpgTracker ):
    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp2 FROM %(track)s_shared_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class tssIntervalCpGNumber( cpgTracker ):
    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_shared_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class tssIntervalGCContent( cpgTracker ):
    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_shared_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data


