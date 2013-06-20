import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['annotations_exportdir']
DATADIR=P['annotations_datadir']
DATABASE=P['annotations_backend']

###########################################################################
class AnnotationTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

class GenomicFeatureCoverage( TrackerSQL ):
    '''return coverage per feature.'''
    
    tablename = "annotation_summary"

    def getTracks(self):
        return self.getValues( "SELECT DISTINCT feature FROM %(tablename)s" )
    
    def __call__(self, track):
        return self.getValue( """SELECT SUM( total_percent_coverage) FROM %(tablename)s 
                               WHERE feature = '%(track)s'""" )

class ChromosomeFeatureCoverage( TrackerSQL ):
    '''return coverage per feature and chromosome.'''
    
    tablename = "annotation_summary"

    def getTracks(self):
        return self.getValues( "SELECT DISTINCT contig FROM %(tablename)s" )
    
    def __call__(self, track):
        return odict( self.get( """SELECT feature, SUM( percent_coverage) AS coverage FROM %(tablename)s 
                               WHERE contig = '%(track)s' GROUP BY feature """ ) )
