import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


# class MappingSummary(ExomeTracker, SingleTableTrackerRows):
#    table = "bam_stats"


class PicardAlignmentSummaryMetrics(ExomeTracker, SingleTableTrackerRows):
    table = "picard_stats_alignment_summary_metrics"

#    def __call__(self, track, slice=None):
#        statement = '''SELECT pas.track, TOTAL_READS, MEAN_READ_LENGTH, PCT_PF_READS_ALIGNED, PCT_READS_ALIGNED_IN_PAIRS, STRAND_BALANCE, MEDIAN_INSERT_SIZE, MEDIAN_ABSOLUTE_DEVIATION, PERCENT_DUPLICATION FROM picard_stats_alignment_summary_metrics pas, picard_stats_insert_size_metrics pis, picard_duplicate_stats_duplicate_metrics pds WHERE pas.track=pds.track AND pas.track=pis.track AND pas.CATEGORY='PAIR';'''
# print (statement)
#        return self.getAll(statement)


class PicardInsertSizeMetrics(ExomeTracker, SingleTableTrackerRows):
    table = "picard_stats_insert_size_metrics"


class PicardDuplicatesMetrics(ExomeTracker, SingleTableTrackerRows):
    table = "picard_duplicate_stats_duplicate_metrics"


class PicardInsertSizeHistogram(ExomeTracker, SingleTableTrackerHistogram):
    table = "picard_stats_insert_size_histogram"
    column = "insert_size"


class PicardDuplicatesHistogram(ExomeTracker, SingleTableTrackerHistogram):
    table = "picard_duplicate_stats_duplicate_histogram"
    column = "duplicates"


class PicardQualityByCycleHistogram(ExomeTracker, SingleTableTrackerHistogram):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"


class PicardQualityDistributionHistogram(ExomeTracker,
                                         SingleTableTrackerHistogram):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"


class PicardCoverageStats(ExomeTracker, SingleTableTrackerRows):
    table = "coverage_stats"

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
