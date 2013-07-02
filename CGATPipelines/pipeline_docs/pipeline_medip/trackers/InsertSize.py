import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from MedipReport import *

class PicardInsertSizeStats( MedipTracker, SingleTableTrackerRows ):
    table = "picard_stats_insert_size_metrics"


    
