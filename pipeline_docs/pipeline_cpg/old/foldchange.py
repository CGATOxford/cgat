import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *


##################################################################################
class FoldChangeThreshold( cpgTracker ):
    """Count of intervals exceeding fold change threshold for each dataset. """

    mPattern = "_foldchange$"
    #  
    def __call__(self, track, slice = None):
        data = self.get( "SELECT threshold, intervals FROM %(track)s_foldchange" % locals() )
        return odict( zip( ("Threshold", "Intervals" ), zip(*data)) )


##################################################################################
class OverlapFoldChangeThreshold( cpgTracker ):
    """Count of overlapping intervals  in different tracks for each fold change for each dataset. """

    mPattern = "_foldchange_shared$"
    # 
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, threshold, intervals FROM %(track)s_foldchange_shared order by track, threshold" % locals() )
        
        return odict( zip( ("Track", "Threshold", "Intervals" ), zip(*data)) )



