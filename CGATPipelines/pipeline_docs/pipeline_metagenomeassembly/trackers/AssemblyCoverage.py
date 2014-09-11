from CGATReport.Tracker import *
import sqlite3
import collections
import numpy as np


class AlignmentSummary(TrackerSQL):

    '''
    class to collect the alignment statistics from
    picard
    '''

    def __call__(self, track, slice=None):
        return self.getAll("""SELECT * FROM picard_stats_alignment_summary_metrics""")


class ReadsAligned(TrackerSQL):

    '''
    percent reads aligned
    '''

    def __call__(self, track, slice=None):

        result = {}
        for data in self.execute("""SELECT track, PCT_PF_READS_ALIGNED FROM picard_stats_alignment_summary_metrics"""):
            result[data[0]] = data[1]
        return result


class ReadsAlignedInPairs(TrackerSQL):

    '''
    percent reads aligned in pairs
    '''

    def __call__(self, track, slice=None):
        result = {}
        for data in self.execute("""SELECT track, PCT_READS_ALIGNED_IN_PAIRS FROM picard_stats_alignment_summary_metrics"""):
            result[data[0]] = data[1]
        return result


class MismatchRate(TrackerSQL):

    '''
    percent reads aligned in pairs
    '''

    def __call__(self, track, slice=None):

        result = {}
        for data in self.execute("""SELECT track, PF_MISMATCH_RATE FROM picard_stats_alignment_summary_metrics"""):
            result[data[0]] = data[1]
        return result


class InsertSizeSummary(TrackerSQL):

    '''
    class to collect the insert size statistics from
    picard
    '''

    def __call__(self, track, slice=None):
        return self.getAll("""SELECT * FROM picard_stats_insert_size_metrics""")


class CoverageSd(TrackerSQL):

    '''
    class to collect data on the standard deviation of base coverage
    across contigs
    '''
    pattern = "(.*)_coverage_stats"

    def __call__(self, track, slice=None):

        result = []
        for data in self.execute("""SELECT cov_sd FROM %(track)s_coverage_stats WHERE cov_sd > 0""" % locals()).fetchall():
            result.append(data[0])
        return np.log2(result)


class CoverageMean(TrackerSQL):

    '''
    class to collect data on the standard deviation of base coverage
    across contigs
    '''
    pattern = "(.*)_coverage_stats"

    def __call__(self, track, slice=None):

        result = []
        for data in self.execute("""SELECT cov_mean FROM %(track)s_coverage_stats WHERE cov_mean > 0""" % locals()).fetchall():
            result.append(data[0])
        return np.log2(result)
