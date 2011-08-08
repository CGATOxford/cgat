import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class BackgroundSummary( cpgTracker ):
    """Summary stats of reads mapping inside/outside binding intervals. """

    mPattern = "_background$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( '''select in_peaks, out_peaks, (in_peaks+out_peaks) as total, 
                                    round(((in_peaks+0.00)/(out_peaks+in_peaks+0.00))*100,1) as ratio,
                                    round(((out_peaks+0.00)/(out_peaks+in_peaks+0.00))*100,1) as ratio2 
                                    from %(track)s_background;''' % locals() )
        return odict( zip( ("Reads Overlapping intervals", "Reads Outwith Intervals", "Total Reads", "Percent in Intervals", "Percent Background" ), data) )


