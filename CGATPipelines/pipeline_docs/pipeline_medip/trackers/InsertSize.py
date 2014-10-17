import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict
from MedipReport import *


class PicardInsertSizeStats(MedipTracker, SingleTableTrackerRows):
    table = "picard_stats_insert_size_metrics"
