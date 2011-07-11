import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class gatResults( cpgTracker ):
    """Summary stats of GAT analysis. """

    mPattern = "gat_results$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM gat_results " )
        return odict( zip( ("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value" ), zip(*data)) )

