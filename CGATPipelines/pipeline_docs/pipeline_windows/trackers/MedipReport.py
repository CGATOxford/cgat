import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script

from SphinxReport.Utils import PARAMS as P

EXPORTDIR = P.get('windows_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('windows_datadir', P.get('datadir', '.'))
DATABASE = P.get('windows_backend', P.get('sql_backend', 'sqlite:///./csvdb'))

###########################################################################
###########################################################################
###########################################################################
class ProjectTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

class PicardDuplicatesMetrics( ProjectTracker, SingleTableTrackerRows ):
    table = "picard_duplicates_duplicate_metrics"

class PicardDuplicatesHistogram( ProjectTracker, SingleTableTrackerHistogram ):
    table = "picard_duplicates_duplicate_histogram"
    column = "duplicates"

class TagCountsSummary( ProjectTracker, SingleTableTrackerRows ):
    table = "counts_stats"
    fields = ("metric", )

class WindowsSummary( ProjectTracker, SingleTableTrackerRows ):
    table = "windows_stats"
    fields = ("data", )

class WindowsSizes( ProjectTracker, SingleTableTrackerHistogram ):
    table = "windows_hist"
    column = "residues"

class TagCountsCorrelations( ProjectTracker ):
    pattern = "(.*)_correlation" 
    
    def __call__(self, track ):
        return self.getDict( "SELECT * FROM %(track)s_correlation" )
    
class TagCountsDesignSummary( ProjectTracker):
    pattern = "design(.*)_stats" 

    def __call__(self,track):
        return self.getAll( "SELECT * FROM design%(track)s_stats" )
