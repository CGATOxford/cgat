import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from MappingReport import *

class MappingSummary( MappingTracker, SingleTableTrackerRows ):
    table = "view_mapping"

class TophatSummary( MappingTracker, SingleTableTrackerRows ):
    table = "tophat_stats"

class StarSummary( MappingTracker, SingleTableTrackerRows ):
    table = "star_stats"

class BamSummary( MappingTracker, SingleTableTrackerRows ):
    table = "bam_stats"

class PicardSummary( MappingTracker, SingleTableTrackerRows ):
    table = "picard_stats_alignment_summary_metrics"

class PicardDuplicationSummary( MappingTracker, SingleTableTrackerRows ):
    table = "picard_duplication_stats_duplication_metrics"

class PicardAlignmentSummaryMetrics( MappingTracker, SingleTableTrackerRows ):
    table = "picard_stats_alignment_summary_metrics"

class PicardInsertSizeMetrics( MappingTracker, SingleTableTrackerRows ):
    table = "picard_stats_insert_size_metrics"

class PicardDuplicatesMetrics( MappingTracker, SingleTableTrackerRows ):
    table = "picard_duplicates_duplicate_metrics"

class PicardInsertSizeHistogram( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_stats_insert_size_histogram"
    column = "insert_size"

class PicardDuplicatesHistogram( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_duplicates_duplicate_histogram"
    column = "duplicates"

class PicardQualityByCycleHistogram( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"

class PicardQualityDistributionHistogram( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"

class DuplicationMetricsTable( MappingTracker, SingleTableTrackerHistogram ):

    table = "picard_duplication_stats_duplication_histogram"

    def __call__(self, track = None, slice = None ):
        cols = self.getColumns(self.table)
        fields = ", ".join(cols)
        data = self.getAll( "SELECT %(fields)s FROM %(table)s ORDER BY coverage_multiple")
        return data


class MappingFlagsMismatches( MappingTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nm"
    column = "nm"

class MappingFlagsHits( MappingTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nh"
    column = "nh"

class AlignmentQualityByCycle( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_stats_quality_by_cycle_histogram"
    column = "cycle"

class DuplicationMetrics( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_duplication_stats_duplication_histogram"
    column = "coverage_multiple"

class AlignmentQualityDistribution( MappingTracker, SingleTableTrackerHistogram ):
    table = "picard_stats_quality_distribution_histogram"
    column = "quality"

class MappingContext( MappingTracker, SingleTableTrackerRows ):
    table = "context_stats"

class FilteringSummary( MappingTracker, SingleTableTrackerRows ):
    table = "mapping_stats"

class BigwigSummary( MappingTracker, SingleTableTrackerRows ):
    table = "bigwig_stats"

##############################################################
##############################################################
##############################################################
class BamReport( MappingTracker ):
    tracks = [ "all" ]

    slices = ("genome", "accepted", "mismapped" )

    def __call__(self, track, slice = None ):
        edir = EXPORTDIR

        toc_text = []
        link_text = []
        
        filenames = sorted( [x.asFile() for x in TRACKS ] )
        
        for fn in filenames:
            link = "%(slice)s_%(fn)s" % locals() 
            toc_text.append( "* %(link)s_" % locals() )
            link_text.append( ".. _%(link)s: %(edir)s/bamstats/%(fn)s.%(slice)s.html" % locals() )
            
        toc_text = "\n".join(toc_text)
        link_text =  "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict( (("text", rst_text),) )
    

##############################################################
##############################################################
##############################################################
class FastQCReport( MappingTracker ):
    tracks = [ "all" ]

    def __call__(self, track, slice = None ):
        edir = EXPORTDIR

        toc_text = []
        link_text = []
        
        filenames = sorted( [x.asFile() for x in TRACKS ] )
        
        for fn in filenames:
            toc_text.append( "* %(fn)s_" % locals()) 
            link_text.append( ".. _%(fn)s: %(edir)s/fastqc/%(fn)s.genome_fastqc/fastqc_report.html" % locals() )
            
        toc_text = "\n".join(toc_text)
        link_text =  "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict( (("text", rst_text),) )
