import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from RnaseqReport import *


class ExonValidationSummary(RnaseqTracker, SingleTableTrackerRows):
    table = "exon_validation"
