import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from IntervalReport import *

class Summary( IntervalTracker, SingleTableTrackerRows ):
    table = "context_stats"
