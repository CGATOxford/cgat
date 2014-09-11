import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict
from cpgReport import *

##########################################################################


class MappingSummary(cpgTracker, SingleTableTrackerRows):
    table = "bam_stats"

##########################################################################


class MappingFlagsMismatches(cpgTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nm"
    column = "nm"

##########################################################################


class PicardAlignPerc(TrackerSQL):

    mPattern = "picard_align_stats"

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.TRACK as sample, total_reads as Total, 
                       round(pct_pf_reads_aligned*100,2) as Mapped, 
                       round(pds.percent_duplication*100,2) as Duplicates,
                       round((1-strand_balance)*100,2) as Reverse, 
                       round(pct_adapter*100,2) as Adapter,
                       round(((pas.pf_reads_aligned - pds.unpaired_read_duplicates+0.0) / (total_reads+0.0))*100,2) as Unique_reads
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        return self.getAll(statement)

##########################################################################


class PicardAlignCount(TrackerSQL):

    mPattern = "picard_align_stats"

    def __call__(self, track, slice=None):
        statement = '''SELECT pas.track as sample, pas.total_reads as Total, 
                       pas.pf_reads_aligned as Mapped, 
                       pds.unpaired_read_duplicates as Duplicates,
                       round(pas.pf_reads_aligned*(1-strand_balance),0) as Reverse, 
                       pas.pf_reads_aligned-pds.unpaired_read_duplicates as Unique_reads
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        data = self.getAll(statement)
        return data

##########################################################################


class PicardAlignPlot(TrackerSQL):

    mPattern = "picard_align_stats"

    def getTracks(self):
        return self.getValues("select distinct track from picard_align_stats")

    def __call__(self, track, slice=None):

        statement = '''SELECT pas.total_reads as Total, 
                       pas.pf_reads_aligned as Mapped, 
                       pds.unpaired_read_duplicates as Duplicates,
                       round(pas.pf_reads_aligned*(1-strand_balance),0) as Reverse, 
                       pas.pf_reads_aligned-pds.unpaired_read_duplicates as Unique_reads
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track and pas.track="%(track)s";'''
        data = self.getRow(statement)

        return data
