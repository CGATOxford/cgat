import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from medipReport import *

class PicardInsertSizeStats( MedipTracker, SingleTableTrackerRows ):
    table = "picard_isize_stats"


    
