import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class PicardInsertSizeStats(ExomeTracker, SingleTableTrackerRows):
    table = "picard_isize_stats"
