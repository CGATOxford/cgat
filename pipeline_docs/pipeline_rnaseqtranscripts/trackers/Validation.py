import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from RnaseqTranscriptsReport import *

class ExonValidationSummary( RnaseqTranscriptsTracker, SingleTableTrackerRows ):
    table = "exon_validation"
