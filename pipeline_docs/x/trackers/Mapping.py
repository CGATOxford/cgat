import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from cpgReport import *

class MappingSummary( cpgTracker, SingleTableTrackerRows ):
    table = "bam_stats"

class MappingFlagsMismatches( cpgTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nm"
    column = "nm"

class PicardAlign( TrackerSQL ):

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.TRACK, total_reads, round(pct_pf_reads_aligned*100,2) || '\%' as reads_mapped, round(pds.percent_duplication*100,2) || '\%' as Duplicates,
                       round(strand_balance*100,2) as Strand_balance, round(pct_adapter*100,2) || '\%' as Adapter
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        #print (statement)
        return self.getAll( statement )

class PicardAlignPlot( TrackerSQL ):

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.total_reads, pas.pf_reads_aligned as Mapped, pds.unpaired_read_duplicates as Duplicates
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        return self.getAll( statement )

