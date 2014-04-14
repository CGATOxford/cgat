import os, sys, re, types, itertools, glob
import pandas
import pandas.rpy.common

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

################################################################################
################################################################################
## parameterization

EXPORTDIR=P['IDR_exportdir']
DATADIR  =P['IDR_datadir']
DATABASE =P['IDR_backend']

################################################################################
class ProjectTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    
################################################################################
# Retrieve tables summarizing the number of peaks called for each condition
################################################################################
class getPeakSummary_IndvidualReplicates( SingleTableTrackerRows ):
    """
    Returns a summary of number of peaks called with tissue field as tracks
    """
    table = "npeak_summary_individual_replicates"
    fields = [ "Sample_id", ]


class getPeakSummary_PseudoReplicates( SingleTableTrackerRows ):
    """
    Returns a summary of number of peaks called with tissue field as tracks
    """
    table = "npeak_summary_pseudoreplicates"
    fields = [ "Sample_id", ]


class getPeakSummary_PooledPseudoReplicates( SingleTableTrackerRows ):
    """
    Returns a summary of number of peaks called with tissue field as tracks
    """
    table = "npeak_summary_pooled_pseudoreplicates"
    fields = [ "Sample_id", ]

################################################################################
# Retrieve tables summarizing the number of peaks above specified IDR threshold 
################################################################################

class getIndividualNPeaks( SingleTableTrackerRows ):
    """
    Returns the table nPeaks with tissue field as tracks
    """
    table = "individual_replicates_nPeaks"
    fields = [ "idr_comp", ]


class getPseudoreplicateNPeaks( SingleTableTrackerRows ):
    """
    Returns the table nPeaks with tissue field as tracks
    """
    table = "pseudoreplicates_nPeaks"
    fields = [ "idr_comp", ]


class getPooledPseudoreplicateNPeaks( SingleTableTrackerRows ):
    """
    Returns the table nPeaks with tissue field as tracks
    """
    table = "pooled_pseudoreplicates_nPeaks"
    fields = [ "idr_comp", ]


################################################################################
# Retrieve tables summarizing the number of peaks at different IDR thresholds 
################################################################################
class getnPeaksAboveIDR( TrackerSQL ):
    """
    """
    def getTracks( self, subset = None ):
        return ( "pro_K4", 
                 "pre_K4" , 
                 "immature_K4", 
                 "mature_K4", 
                 "follicular_K4", 
                 "marginal_K4", 
                 "germinal_K4", 
                 "b1a_K4",
                 "pro_K4me1", 
                 "pre_K4me1" , 
                 "immature_K4me1", 
                 "mature_K4me1", 
                 "follicular_K4me1", 
                 "marginal_K4me1", 
                 "germinal_K4me1", 
                 "b1a_K4me1" ) 


class getIndividualnPeaksAboveIDR( getnPeaksAboveIDR ):
    """
    """

    def __call__( self, track ):
        statement = ( "SELECT * FROM  %(track)s_npeaks_aboveIDR"  )
        return self.getAll( statement )


class getPseudoreplicatenPeaksAboveIDR( getnPeaksAboveIDR ):
    """
    """

    def __call__( self, track ):
        statement = ( "SELECT * FROM  %(track)s_pseudoreplicates_npeaks_aboveIDR"  )
        return self.getAll( statement )


class getPooledPseudoreplicatenPeaksAboveIDR( getnPeaksAboveIDR ):
    """
    """

    def __call__( self, track ):
        statement = ( "SELECT * FROM  %(track)s_pooled_pseudoreplicates_npeaks_aboveIDR"  )
        return self.getAll( statement )

################################################################################
# load the plots output by batch-consistency-plot.r
################################################################################

class loadBatchConsistencyPlots( TrackerImages ):
    pass
