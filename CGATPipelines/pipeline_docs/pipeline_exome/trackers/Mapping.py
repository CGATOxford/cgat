import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class MappingSummary(ExomeTracker, SingleTableTrackerRows):
    table = "bam_stats"


class PicardAlign(ExomeTracker):

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.track, TOTAL_READS, MEAN_READ_LENGTH, PCT_PF_READS_ALIGNED, PCT_READS_ALIGNED_IN_PAIRS, STRAND_BALANCE, MEDIAN_INSERT_SIZE, MEDIAN_ABSOLUTE_DEVIATION, PERCENT_DUPLICATION FROM picard_stats_alignment_summary_metrics pas, picard_stats_insert_size_metrics pis, picard_duplicate_stats_duplicate_metrics pds WHERE pas.track=pds.track AND pas.track=pis.track AND pas.CATEGORY='PAIR';'''
        #print (statement)
        return self.getAll(statement)


class PicardAlignPlot(ExomeTracker):

    @property
    def tracks(self):
        d = self.get(
            "SELECT DISTINCT track FROM picard_stats_alignment_summary_metrics")
        d = list([x[0] for x in d])
        d = [x.replace('.', '_') for x in d]
        d = [x.replace('-', '_') for x in d]
        return tuple(d)

    def __call__(self, track, slice=None):
        statement = '''SELECT insert_size, %(track)s FROM picard_stats_insert_size_histogram;'''
        return self.getAll(statement)


# class AlignmentSummary( SingleTableTrackerRows ):
#    table = "alignment_stats"

# class MappingFlagsMismatches( SingleTableTrackerHistogram ):
#    table = "bam_stats_nm"
#    column = "nm"

# class MappingFlagsHits( SingleTableTrackerHistogram ):
#    table = "bam_stats_nh"
#    column = "nh"

# class AlignmentQualityByCycle( SingleTableTrackerHistogram ):
#    table = "alignment_stats_quality_by_cycle_metrics"
#    column = "cycle"

# class AlignmentQualityDistribution( SingleTableTrackerHistogram ):
#    table = "alignment_stats_quality_distribution_metrics"
#    column = "quality"
