import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from RnaseqReport import *

class MappingSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "view_mapping"

class TophatSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "tophat_stats"

class BamSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "bam_stats"

class AlignmentSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "alignment_stats"

class MappingContext( RnaseqTracker, SingleTableTrackerRows ):
    table = "context_stats"

class FilteringSummary( RnaseqTracker, SingleTableTrackerRows ):
    table = "mapping_stats"

class MappingFlagsMismatches( RnaseqTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nm"
    column = "nm"

class MappingFlagsHits( RnaseqTracker, SingleTableTrackerHistogram ):
    table = "bam_stats_nh"
    column = "nh"

class AlignmentQualityByCycle( RnaseqTracker, SingleTableTrackerHistogram ):
    table = "alignment_stats_quality_by_cycle_metrics"
    column = "cycle"

class AlignmentQualityDistribution( RnaseqTracker, SingleTableTrackerHistogram ):
    table = "alignment_stats_quality_distribution_metrics"
    column = "quality"
    
