import os
import sys
import re
import types
import itertools
import glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

import CGAT.IOTools as IOTools
import CGAT.Stats as Stats

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
# parameterization

EXPORTDIR = P['intervals_exportdir']
DATADIR = P['intervals_datadir']
DATABASE = P['intervals_backend']

###########################################################################


class IntervalTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
