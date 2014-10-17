import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class CoverageSummary(ExomeTracker):

    def __call__(self, track, slice=None):
        statement = '''SELECT track, MEAN_TARGET_COVERAGE, PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X, PCT_TARGET_BASES_30X FROM coverage_stats;'''
        #print (statement)
        return self.getAll(statement)


class CoveragePlot(ExomeTracker):

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM coverage_stats")
        return tuple([x[0] for x in d])

    def __call__(self, track, slice=None):
        statement = '''SELECT MEAN_TARGET_COVERAGE FROM coverage_stats where track='%(track)s';'''
        return self.getAll(statement)
