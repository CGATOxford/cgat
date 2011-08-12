import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

# get from config file
UCSC_DATABASE="hg19"
EXPORTDIR="export"

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script

from SphinxReport.Utils import PARAMS as P
EXPORTDIR=P['medip_exportdir']
DATADIR=P['medip_datadir']
DATABASE=P['medip_backend']

###########################################################################

class MedipTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )


