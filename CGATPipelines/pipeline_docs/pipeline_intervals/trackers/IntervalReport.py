import os
import sys
import re
import types
import itertools
import glob

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from collections import OrderedDict as odict

import CGAT.IOTools as IOTools
import CGAT.Stats as Stats

from CGATReport.ResultBlock import ResultBlock, ResultBlocks
from CGATReport import Utils

###################################################################
###################################################################
# parameterization
EXPORTDIR = P.get('intervals_exportdir',
                  P.get('exportdir', 'export'))
DATADIR = P.get('intervals_datadir',
                P.get('datadir', '.'))
DATABASE = P.get('intervals_backend',
                 P.get('sql_backend', 'sqlite:///./csvdb'))


class IntervalTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
