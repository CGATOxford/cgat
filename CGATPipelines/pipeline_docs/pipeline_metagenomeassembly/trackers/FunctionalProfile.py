from CGATReport.Tracker import *
import sqlite3
import collections
import numpy as np


class AlignmentCounts(TrackerSQL):

    def __call__(self, track, slice=None):
        '''return counts of unique sequence
        alignments to protein database
        '''
        result = collections.defaultdict(list)
        for data in self.execute("""SELECT track, proportion FROM rpsblast_alignment_stats""").fetchall():
            result["track"].append(data[0])
            result["proportion"].append(data[1])
        return result


class CogCounts(TrackerSQL):

    pattern = "(.*)_cog_counts"

    def __call__(self, track, slice=None):
        '''
        return the distribution of reads across
        COG categories
        '''
        result = {}
        for data in self.execute("""SELECT funtion, proportion FROM %s_cog_counts""" % track).fetchall():
            result[data[0]] = data[1]
        return result
