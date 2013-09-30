import os, sys, re, types, itertools, math, numpy

from MedipReport import *

class TileSizeHistogram( SingleTableTrackerHistogram, ProjectTracker ):
    table = "tileinfo_hist"
    column = "residues"

class TileSizeStats( SingleTableTrackerRows, ProjectTracker ):
    table = "tileinfo_stats"
