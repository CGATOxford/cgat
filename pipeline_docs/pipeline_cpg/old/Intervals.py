import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class IntervalsSummary( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*), round(AVG(length),0), round(AVG(nprobes),0)  FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("intervals_count", "mean_interval_length", "mean_reads_per_interval" ), data) )

##################################################################################
class IntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT length FROM %(track)s_macs_intervals" % locals() )
        return { "length" : data }

##################################################################################
class IntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT peakval FROM %(track)s_macs_intervals" % locals() )
        return { "peakval" : data }

##################################################################################
class IntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT avgval FROM %(track)s_macs_intervals" % locals() )
        return { "avgval" : data }

##################################################################################
class FoldChange( cpgTracker ):
    """return fold changes for all intervals. """

    mPattern = "_macs_intervals$"
    
    def __call__(self, track, slice = None ):
        data = self.getValues( "SELECT fold FROM %(track)s_macs_intervals" % locals() )
        return { 'fold' : data }

##################################################################################
##################################################################################
##################################################################################
class PeakLocation( cpgTracker ):
    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_macs_intervals" % locals() )
        data2 = self.getValues( "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_macs_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
class PeakDistance( cpgTracker ):
    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT PeakCenter - start FROM %(track)s_macs_intervals" % locals() )
        data2 = self.getValues( "SELECT end - PeakCenter FROM %(track)s_macs_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
##################################################################################
##################################################################################
class IntervalLengthVsAverageValue( cpgTracker ):
    """Length vs average value. """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, avgval FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("length", "avgval"), zip(*data) ) )

##################################################################################
class IntervalLengthVsPeakValue( cpgTracker ):
    """Length vs peak value """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, peakval FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("length", "peakval"), zip(*data) ) )

##################################################################################
class IntervalLengthVsFoldChange( cpgTracker ):
    """Length vs fold change"""

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, fold FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("length", "foldchange"), zip(*data) ) )

##################################################################################
class IntervalAvgValVsPeakVal( cpgTracker ):
    """average value vs peak value """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT avgval, peakval FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("avgval", "peakval"), zip(*data) ) )

##################################################################################
class IntervalAvgValVsFoldChange( cpgTracker ):
    """average value vs fold change """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT avgval, fold FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("avgval", "foldchange"), zip(*data) ) )

##################################################################################
class IntervalPeakValVsFoldChange( cpgTracker ):
    """Peak value vs fold change """

    mPattern = "_macs_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT peakval, fold FROM %(track)s_macs_intervals" % locals() )
        return odict( zip( ("peakval", "foldchange"), zip(*data) ) )


