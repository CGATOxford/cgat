import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class UniqueIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_unique_intervals" % locals() )
        return odict( zip( ("Unique intervals", "mean_interval_length" ), data) )

##################################################################################
class SharedIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_shared_intervals" % locals() )
        return odict( zip( ("Shared intervals", "mean_interval_length" ), data) )

##################################################################################
class OverlapIntervals( cpgTracker ):
    """Count of overlapping intervals for each dataset. """

    mPattern = "_overlap$"
    #  
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_overlap" % locals() )
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class OverlapCpG( cpgTracker ):
    """Count of intervals overlapping CpG islands for each dataset. """

    mPattern = "_cgi$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_cgi" % locals() )
        
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class OverlapChipseq( cpgTracker ):
    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_chipseq$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_chipseq" % locals() )
        
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class OverlapCAPseq( cpgTracker ):
    """Count of intervals overlapping CAPseq intervals for each dataset. """

    mPattern = "_capseq$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_capseq" % locals() )
        
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class OverlapChromatinMarks( cpgTracker ):
    """Count of intervals overlapping ChIPseq intervals for each dataset. """

    mPattern = "_chromatin$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_chromatin" % locals() )
        
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )

##################################################################################
class gatResults( cpgTracker ):
    """Summary stats of GAT analysis. """

    mPattern = "gat_results$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM gat_results " )
        return odict( zip( ("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value" ), zip(*data)) )


