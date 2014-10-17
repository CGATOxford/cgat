import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict
from MedipReport import *


class MappingSummary(MedipTracker, SingleTableTrackerRows):
    table = "bam_stats"


class BamSummary(MedipTracker, SingleTableTrackerRows):
    table = "bam_stats"


class MappingFlagsMismatches(MedipTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nm"
    column = "nm"


class MappingFlagsHits(MedipTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nh"
    column = "nh"


class PicardAlignmentSummaryMetrics(MedipTracker, SingleTableTrackerRows):
    table = "picard_stats_alignment_summary_metrics"


class PicardInsertSizeMetrics(MedipTracker, SingleTableTrackerRows):
    table = "picard_stats_insert_size_metrics"


class PicardDuplicatesMetrics(MedipTracker, SingleTableTrackerRows):
    table = "picard_duplicates_duplicate_metrics"


class PicardInsertSizeHistogram(MedipTracker, SingleTableTrackerHistogram):
    table = "picard_stats_insert_size_histogram"
    column = "insert_size"


class PicardDuplicatesHistogram(MedipTracker, SingleTableTrackerHistogram):
    table = "picard_duplicates_duplicate_histogram"
    column = "duplicates"


class PicardQualityByCycleHistogram(MedipTracker, SingleTableTrackerHistogram):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"


class PicardQualityDistributionHistogram(MedipTracker, SingleTableTrackerHistogram):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"


class MappingContext(MedipTracker, SingleTableTrackerRows):
    table = "context_stats"


class PicardAlign(MedipTracker):

    mPattern = "picard_stats_alignment_summary_metrics"

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.TRACK, total_reads/2 as Total_read_pairs, round(pct_reads_aligned_in_pairs*pct_pf_reads_aligned*100,2) as Aligned_pairs, 
                       round(strand_balance*100,2) as Strand_balance, round(pds.percent_duplication*100,2) as Duplicates, round(pct_adapter,2) as Adapter
                       FROM picard_stats_alignment_summary_metrics pas, picard_duplicates_duplicate_metrics pds
                       where pas.track=pds.track and pas.category='PAIR';'''

        return self.getAll(statement)


class PicardAlignPlot(MedipTracker):

    mPattern = "picard_stats_alignment_summary_metrics"

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.TRACK, pas.total_reads/2 as Total_read_pairs, pas.reads_aligned_in_pairs/2 as Aligned_pairs, 
                       ROUND((pas.reads_aligned_in_pairs/2)*(1-strand_balance),0) as reverse, pds.read_pair_duplicates as duplicate_pairs
                       FROM picard_stats_alignment_summary_metrics pas, picard_duplicates_duplicate_metrics pds
                       WHERE pas.track=pds.track and pas.category='PAIR';'''
        data = self.get(statement)
        result = odict()

        # Define tracks as first column
        for d in data:
            result[d[0]] = odict()

        # Define slices as other columns
        for d in data:
            for s, v in zip(("total", "mapped", "reverse", "duplicates"), d[1:]):
                result[d[0]][s] = v
        # print result
        return result


class CpGCoverage(MedipTracker, SingleTableTrackerHistogram):
    table = "cpg_coverage"
    column = "bin"


class CpGCoverageHistogram(MedipTracker, SingleTableTrackerHistogram):
    table = "cpg_coverage"
    column = "bin"

    def __call__(self):
        pass
