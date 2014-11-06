import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from rrbsReport import *


class seqStart(RrbsTracker, SingleTableTrackerRows):
    table = "read_start_summary"

class coverage(RrbsTracker, SingleTableTrackerRows):
    table = "coverage"

class readsRemaining(RrbsTracker, SingleTableTrackerRows):
    table = "reads_remaining_by_threshold"

class readsRemaining(RrbsTracker, SingleTableTrackerRows):
    table = "coverage_overlap"

