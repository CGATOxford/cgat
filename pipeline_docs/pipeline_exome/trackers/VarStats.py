import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *

class VariantSummary( ExomeTracker, SingleTableTrackerRows ):
    table = "vcf_stats"



