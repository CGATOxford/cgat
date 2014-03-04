import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *

class PicardInsertSizeStats( ExomeTracker, SingleTableTrackerRows ):
    table = "picard_isize_stats"


    
