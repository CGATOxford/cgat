import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script

from SphinxReport.Utils import PARAMS as P
EXPORTDIR=P['windows_exportdir']
DATADIR=P['windows_datadir']
DATABASE=P['windows_backend']

###########################################################################
###########################################################################
###########################################################################
class ProjectTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )


