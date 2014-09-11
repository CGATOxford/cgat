import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from RnaseqTranscriptsReport import *


class ExonValidationSummary(RnaseqTranscriptsTracker, SingleTableTrackerRows):
    table = "exon_validation"
