import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *


##################################################################################
class replicatedUniqueIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_unique_intervals" % locals() )
        return odict( zip( ("Unique intervals", "mean_interval_length" ), data) )

##################################################################################
class replicatedUniqueIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT (stop-start) FROM %(track)s_unique_intervals" % locals() )
        return { "length" : data }

##################################################################################
class replicatedUniqueIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_unique_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "peakval" : data }

##################################################################################
class replicatedUniqueIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_unique_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "avgval" : data }

##################################################################################
class replicatedUniqueIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT fold FROM %(track)s_unique_intervals u, %(track)s_macs_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return odict( [("Fold Change", data)] )

##################################################################################
class replicatedUniqueIntervalTSS( cpgTracker ):
    """Distribution of distance to closest TSS """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT closest_dist FROM %(track)s_unique_intervals u, 
                                  %(track)s_macs_intervals i, %(track)s_tss t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start
                                  AND t.gene_id=i.interval_id''' % locals() )
        return { "distance" : data }

##################################################################################
class replicatedUniqueIntervalCpGDensity( cpgTracker ):
    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pCpG FROM %(track)s_unique_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedUniqueIntervalCpGObsExp1( cpgTracker ):
    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp1 FROM %(track)s_unique_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedUniqueIntervalCpGObsExp2( cpgTracker ):
    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp2 FROM %(track)s_unique_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedUniqueIntervalCpGNumber( cpgTracker ):
    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_unique_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedUniqueIntervalGCContent( cpgTracker ):
    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_unique_intervals u, 
                               %(track)s_macs_intervals i,%(track)s_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data


