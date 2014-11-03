import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class PicardInsertSizeStats(ExomeTracker, SingleTableTrackerRows):
    table = "picard_isize_stats"
