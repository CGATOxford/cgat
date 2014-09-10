import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from IntervalReport import *


class ContextSummary(IntervalTracker, SingleTableTrackerRows):
    table = "context_stats"
