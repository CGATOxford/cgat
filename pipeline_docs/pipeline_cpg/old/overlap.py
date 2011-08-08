import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *


##################################################################################
class OverlapIntervals( cpgTracker ):
    """Count of overlapping intervals for each dataset. """

    mPattern = "_overlap$"
    #  
    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, overlap FROM %(track)s_overlap" % locals() )
        return odict( zip( ("Track", "Intervals" ), zip(*data)) )


##################################################################################
class OverlapCpG2( cpgTracker ):
    """Count of intervals overlapping CpG islands for each dataset. """

    mPattern = "_overlap$"
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

