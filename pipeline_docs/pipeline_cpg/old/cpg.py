import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *


##################################################################################
class CpGDensity( cpgTracker ):
    mPattern = "_composition$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT pCpG FROM %(track)s_composition" % locals() )
        return { "CpG_density" : data1 }

##################################################################################
class CpGObsExp1( cpgTracker ):
    mPattern = "_composition$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT CpG_ObsExp1 FROM %(track)s_composition" % locals() )
        return { "CpG_ObsExp" : data1 }

##################################################################################
class CpGObsExp2( cpgTracker ):
    mPattern = "_composition$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT CpG_ObsExp2 FROM %(track)s_composition" % locals() )
        return { "CpG_ObsExp" : data1 }

##################################################################################
class CpGNumber( cpgTracker ):
    mPattern = "_composition$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT nCpG FROM %(track)s_composition" % locals() )
        return { "CpG_number" : data1 }

##################################################################################
class GCContent( cpgTracker ):
    mPattern = "_composition$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT pGC FROM %(track)s_composition" % locals() )
        return { "GC_content" : data1 }
