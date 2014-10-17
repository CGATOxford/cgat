import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict
from titrationReport import *


class MappingSummary(benchmarkTracker, SingleTableTrackerRows):
    table = "bam_stats"


class MappingFlagsMismatches(benchmarkTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nm"
    column = "nm"


class MappingFlagsHits(benchmarkTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nh"
    column = "nh"


class PicardAlign(benchmarkTracker):

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.TRACK, total_reads/2 as Total_read_pairs, round(pct_reads_aligned_in_pairs*100,2) || '\%' as Aligned_pairs, 
                       round(strand_balance*100,2) as Strand_balance, round(pds.percent_duplication*100,2) as Duplicates  
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track and pas.category='PAIR';'''
        #print (statement)
        return self.getAll(statement)


class PicardAlignPlot(benchmarkTracker):

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM picard_align_stats")
        return tuple([x[0] for x in d])

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.total_reads/2 as Total_read_pairs, pas.reads_aligned_in_pairs/2 as mapped_pairs, pds.read_pair_duplicates as duplicate_pairs
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track and pas.category='PAIR' and pas.track='%(track)s';'''
        return self.getAll(statement)
