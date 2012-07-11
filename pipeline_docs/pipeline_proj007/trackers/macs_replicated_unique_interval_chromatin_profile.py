import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import cpgReport

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
##################################################################################
##################################################################################
class h3k4me3ProfileTracker(TrackerImages):
    """Chromatin profile per gene """

