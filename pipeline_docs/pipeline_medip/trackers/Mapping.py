import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from medipReport import *

class MappingSummary( MedipTracker, SingleTableTrackerRows ):
    table = "bam_stats"

class MappingFlagsMismatches( MedipTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nm"
    column = "nm"

class PicardAlign( MedipTracker ):

    mPattern = "picard_align_stats"

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.TRACK, total_reads/2 as Total_read_pairs, round(pct_reads_aligned_in_pairs*pct_pf_reads_aligned*100,2) as Aligned_pairs, 
                       round(strand_balance*100,2) as Strand_balance, round(pds.percent_duplication*100,2) as Duplicates, round(pct_adapter,2) as Adapter
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track and pas.category='PAIR';'''
        #print (statement)
        return self.getAll( statement )

class PicardAlignPlot( MedipTracker ):

    mPattern = "picard_align_stats"

    def __call__(self, track, slice = None ):
        statement = '''SELECT pas.TRACK, pas.total_reads/2 as Total_read_pairs, pas.reads_aligned_in_pairs/2 as Aligned_pairs, 
                       round((pas.reads_aligned_in_pairs/2)*(1-strand_balance),0) as reverse, pds.read_pair_duplicates as duplicate_pairs
                       FROM picard_align_stats pas, picard_duplicate_stats pds
                       where pas.track=pds.track and pas.category='PAIR';'''
        data = self.get( statement )
        result = odict()

        # Define tracks as first column
        for d in data:             
            result[d[0]] = odict()

        # Define slices as other columns
        for d in data:
            for s, v in zip( ("total", "mapped", "reverse", "duplicates"), d[1:]):
                result[d[0]][s] = v
        #print result
        return result

