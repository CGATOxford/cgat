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

    mPattern = "picard_align_stats"

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.TRACK, total_reads as total, round(pct_pf_reads_aligned*100,2) as mapped, round(pds.percent_duplication*100,2) as Duplicates,
                       round((1-strand_balance)*100,2) as reverse, round(pct_adapter*100,2) as Adapter
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        #print (statement)
        return self.getAll( statement )

class PicardAlignPlot( TrackerSQL ):

    mPattern = "picard_align_stats"

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.track, pas.total_reads, pas.pf_reads_aligned as Mapped, 
                       round(pas.pf_reads_aligned*(1-strand_balance),0) as reverse, pds.unpaired_read_duplicates as Duplicates
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track;'''
        data = self.get( statement )
        result = odict()

        # Define tracks as first column
        for d in data:             
            result[d[0]] = odict()

        # Define slices as other columns
        for d in data:
            for s, v in zip( ("total", "mapped", "reverse", "duplicates"), d[1:]):
                result[d[0]][s] = v
        print result
        return result

