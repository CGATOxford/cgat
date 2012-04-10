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
class replicatedSharedIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_shared_intervals" % locals() )
        return odict( zip( ("Shared intervals", "mean_interval_length" ), data) )

##################################################################################
class replicatedOverlapIntervals( cpgTracker ):
    """Count of overlapping intervals for each dataset. """

    mPattern = "_overlap$"
    #  
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_overlap" % locals() )
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class replicatedOverlapCpG( cpgTracker ):
    """Count of intervals overlapping CpG islands for each dataset. """

    mPattern = "_cgi$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( """SELECT c.track, count(i.start) as A_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(count(i.interval_id)+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_cgi c, external_interval_sets e, %(track)s_replicated_intervals i 
                            WHERE c.track=e.bed""" % locals() )
        
        return odict( zip( ("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B" ), zip(*data)) )

##################################################################################
class replicatedOverlapChipseq( cpgTracker ):
    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_chipseq$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( """SELECT c.track, count(i.start) as A_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(count(i.interval_id)+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_chipseq c, external_interval_sets e, %(track)s_replicated_intervals i 
                            WHERE c.track=e.bed""" % locals() )
        
        return odict( zip( ("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B" ), zip(*data)) )

##################################################################################
class replicatedOverlapCAPseq( cpgTracker ):
    """Count of intervals overlapping CAPseq intervals for each dataset. """

    mPattern = "_capseq$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( """SELECT c.track, count(i.start) as A_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(count(i.interval_id)+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_capseq c, external_interval_sets e, %(track)s_replicated_intervals i 
                            WHERE c.track=e.bed""" % locals() )
        
        return odict( zip( ("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B" ), zip(*data)) )

##################################################################################
class replicatedOverlapChromatinMarks( cpgTracker ):
    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_chromatin$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( """SELECT c.track, count(i.start) as A_intervals, e.intervals, c.overlap, 
                            round(((c.overlap+0.0)/(count(i.interval_id)+0.0))*100,2) as perca, 
                            round(((c.overlap+0.0)/(e.intervals+0.0))*100,2) as percb 
                            FROM %(track)s_chromatin c, external_interval_sets e, %(track)s_replicated_intervals i 
                            WHERE c.track=e.bed""" % locals() )
        
        return odict( zip( ("Track", "A Intervals", "B Intervals", "Overlap", "%A", "%B" ), zip(*data)) )

##################################################################################
class replicatedGatResults( cpgTracker ):
    """Summary stats of GAT analysis. """

    mPattern = "gat_results$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM external_dataset_gat_results " )
        return odict( zip( ("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value" ), zip(*data)) )


