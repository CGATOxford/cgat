import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *


class CoverageSummary(ExomeTracker):

    def __call__(self, track, slice=None):
        statement = '''SELECT TRACK, round(AVG(cov_mean),2) as mean_cov, round(AVG(cov_median),2) as median_cov, round(AVG(cov_sd),2) as sd_cov, 
                        count(feature) as features, sum(covered) as covered
                        FROM coverage_stats, (select track as trk, feature as feat, case when cov_mean > 10 then 1 else 0 end as covered from coverage_stats) as covered
                        where covered.trk=coverage_stats.track and covered.feat=coverage_stats.feature
                        GROUP BY TRACK;'''
        #print (statement)
        return self.getAll(statement)


class CoveragePlot(ExomeTracker):

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM coverage_stats")
        return tuple([x[0] for x in d])

    def __call__(self, track, slice=None):
        statement = '''SELECT cov_mean FROM coverage_stats where track='%(track)s';'''
        return self.getAll(statement)
