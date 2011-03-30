import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

from RnaseqReport import *

class TophatSummary( SingleTableTrackerRows ):
    table = "tophat_stats"
class MappingSummary( SingleTableTrackerRows ):
    table = "bam_stats"
class AlignmentSummary( SingleTableTrackerRows ):
    table = "alignment_stats"
class MappingContext( SingleTableTrackerRows ):
    table = "context_stats"
class FilteringSummary( SingleTableTrackerRows ):
    table = "mapping_stats"

class MappingFlagsMismatches( SingleTableTrackerHistogram ):
    table = "bam_stats_nm"
    column = "nm"

class MappingFlagsHits( SingleTableTrackerHistogram ):
    table = "bam_stats_nh"
    column = "nh"

class AlignmentQualityByCycle( SingleTableTrackerHistogram ):
    table = "alignment_stats_quality_by_cycle_metrics"
    column = "cycle"

class AlignmentQualityDistribution( SingleTableTrackerHistogram ):
    table = "alignment_stats_quality_distribution_metrics"
    column = "quality"
    
