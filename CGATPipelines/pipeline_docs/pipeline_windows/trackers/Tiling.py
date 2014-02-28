import os
import sys
import re
import types
import itertools
import math
import numpy

from MedipReport import *


class TileSizeHistogram(SingleTableTrackerHistogram, ProjectTracker):
    table = "tileinfo_hist"
    column = "residues"


class TileSizeStats(SingleTableTrackerRows, ProjectTracker):
    table = "tileinfo_stats"
