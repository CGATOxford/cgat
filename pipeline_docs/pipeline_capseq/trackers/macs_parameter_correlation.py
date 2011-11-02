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

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, avgval FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("length", "avgval"), zip(*data) ) )

##################################################################################
class IntervalLengthVsPeakValue( cpgTracker ):
    """Length vs peak value """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, peakval FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("length", "peakval"), zip(*data) ) )

##################################################################################
class IntervalLengthVsFoldChange( cpgTracker ):
    """Length vs fold change"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, fold FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("length", "foldchange"), zip(*data) ) )

##################################################################################
class IntervalLengthVsCpG( cpgTracker ):
    """Length vs CpG Observed/expected"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.length, c.CpG_ObsExp2 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("length", "CpG_ObsExp"), zip(*data) ) )

##################################################################################
class IntervalLengthVsGC( cpgTracker ):
    """Length vs GC content"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.length, c.pGC 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("length", "GC_content"), zip(*data) ) )

##################################################################################
class TSSDistVsLength(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT i.length, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i 
                       WHERE i.interval_id = d.gene_id """ % locals()
        data = self.getAll( statement ) 
        return data

##################################################################################
##################################################################################
class IntervalAvgValVsPeakVal( cpgTracker ):
    """average value vs peak value """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT avgval, peakval FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("avgval", "peakval"), zip(*data) ) )

##################################################################################
class IntervalAvgValVsFoldChange( cpgTracker ):
    """average value vs fold change """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT avgval, fold FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("avgval", "foldchange"), zip(*data) ) )

##################################################################################
class IntervalAvgValVsCpG( cpgTracker ):
    """average value vs CpG Observed/expected"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.avgval, c.CpG_ObsExp2 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("avgval", "CpG_ObsExp"), zip(*data) ) )

##################################################################################
class IntervalAvgValVsGC( cpgTracker ):
    """average value vs GC content"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.avgval, c.pGC 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("avgval", "GC_content"), zip(*data) ) )

##################################################################################
class TSSDistVsAvgVal(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT i.avgval, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i 
                       WHERE i.interval_id = d.gene_id """ % locals()
        data = self.getAll( statement ) 
        return data

##################################################################################
##################################################################################
class IntervalPeakValVsFoldChange( cpgTracker ):
    """Peak value vs fold change """

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT peakval, fold FROM %(track)s_macs_merged_intervals" % locals() )
        return odict( zip( ("peakval", "foldchange"), zip(*data) ) )

##################################################################################
class IntervalPeakValVsCpG( cpgTracker ):
    """Peak value vs CpG Observed/expected"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.peakval, c.CpG_ObsExp2 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("peakval", "CpG_ObsExp"), zip(*data) ) )

##################################################################################
class IntervalPeakValVsGC( cpgTracker ):
    """Peak value vs GC content"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.peakval, c.pGC 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("peakval", "GC_content"), zip(*data) ) )

##################################################################################
class TSSDistVsPeakVal(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT i.peakval, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i 
                       WHERE i.interval_id = d.gene_id """ % locals()
        data = self.getAll( statement ) 
        return data

##################################################################################
##################################################################################
class IntervalFoldChangeVsCpG( cpgTracker ):
    """Fold change vs CpG Observed/expected"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.fold, c.CpG_ObsExp2 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("Fold_Change", "CpG_ObsExp"), zip(*data) ) )

##################################################################################
class IntervalFoldChangeVsGC( cpgTracker ):
    """Fold change vs GC content"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT i.fold, c.pGC 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("Fold_Change", "GC_content"), zip(*data) ) )

##################################################################################
class TSSDistVsFoldChange(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT i.fold, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i 
                       WHERE i.interval_id = d.gene_id """ % locals()
        data = self.getAll( statement ) 
        return data

##################################################################################
##################################################################################
class IntervalCpGVsGC( cpgTracker ):
    """CpG vs GC content"""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( """SELECT c.CpG_ObsExp2 , c.pGC 
                            FROM %(track)s_macs_merged_intervals i, %(track)s_merged_composition c
                            WHERE i.interval_id=c.gene_id""" % locals() )
        return odict( zip( ("CpG_ObsExp", "GC_content"), zip(*data) ) )

##################################################################################
class TSSDistVsCpG(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT c.CpG_ObsExp2, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i, %(track)s_merged_composition c
                       WHERE i.interval_id = d.gene_id AND i.interval_id = c.gene_id""" % locals()
        data = self.getAll( statement ) 
        return data

##################################################################################
class TSSDistVsGC(cpgTracker):
    """for each interval, return peakval and the distance to the closest TSS."""

    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice = None ):

        statement = """SELECT c.pGC, d.closest_dist FROM %(track)s_merged_tss AS d, %(track)s_macs_merged_intervals AS i, %(track)s_merged_composition c
                       WHERE i.interval_id = d.gene_id AND i.interval_id = c.gene_id""" % locals()
        data = self.getAll( statement ) 
        return data


