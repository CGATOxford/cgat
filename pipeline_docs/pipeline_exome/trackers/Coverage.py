import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *

class CoverageSummary( ExomeTracker ):
    def __call__(self, track, slice = None ):
        statement = '''SELECT TRACK, round(AVG(cov_mean),2) as mean_cov, round(AVG(cov_median),2) as median_cov, round(AVG(cov_sd),2) as sd_cov, 
                        count(feature) as features, sum(covered) as covered
                        FROM coverage_stats, (select track as trk, feature as feat, case when cov_mean > 10 then 1 else 0 end as covered from coverage_stats) as covered
                        where covered.trk=coverage_stats.track and covered.feat=coverage_stats.feature
                        GROUP BY TRACK;'''
        #print (statement)
        return self.getAll( statement )

class CoveragePlot( ExomeTracker ):

    def __call__(self, track, slice = None ):
        statement = '''SELECT TRACK, cov_mean FROM coverage_stats;'''
        return self.getAll( statement )

class CoveragePlot2( ExomeTracker, SingleTableHistogram ):
    table = "coverage_stats"
    columns = "cov_mean"
    group_by = "track"
