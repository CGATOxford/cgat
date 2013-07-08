import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['intervals_exportdir']
DATADIR=P['intervals_datadir']
DATABASE=P['intervals_backend']

###########################################################################
class IntervalTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )


