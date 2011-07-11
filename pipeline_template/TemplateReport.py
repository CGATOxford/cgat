import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['@template@_exportdir']
DATADIR  =P['@template@_datadir']
DATABASE =P['@template@_backend']

###########################################################################
class ProjectTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    
class WordFrequencies( ProjectTracker ):
    pattern = "(.*)_counts" 
    def __call__(self, track ):
        return self.getValues( "SELECT freq FROM %(track)s_counts" )
