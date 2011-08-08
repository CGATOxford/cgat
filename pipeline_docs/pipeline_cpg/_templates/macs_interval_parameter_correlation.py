import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

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


