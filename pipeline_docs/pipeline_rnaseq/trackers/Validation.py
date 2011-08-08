import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from RnaseqReport import *

class ExonValidationSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "exon_validation"
