import os, sys, re, types, itertools, math, numpy

from MedipReport import *

class TileSizeHistogram( SingleTableTrackerHistogram, MedipTracker ):
    table = "tileinfo_hist"
    column = "residues"

class TileSizeStats( SingleTableTrackerRows, MedipTracker ):
    table = "tileinfo_stats"
