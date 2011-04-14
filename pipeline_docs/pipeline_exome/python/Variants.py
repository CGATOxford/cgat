import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *

class VariantSummary( SingleTableTrackerRows ):
    table = "vcf_stats"
#class AlignmentSummary( SingleTableTrackerRows ):
#    table = "alignment_stats"

#class MappingFlagsMismatches( SingleTableTrackerHistogram ):
#    table = "bam_stats_nm"
#    column = "nm"

#class MappingFlagsHits( SingleTableTrackerHistogram ):
#    table = "bam_stats_nh"
#    column = "nh"

#class AlignmentQualityByCycle( SingleTableTrackerHistogram ):
#    table = "alignment_stats_quality_by_cycle_metrics"
#    column = "cycle"

#class AlignmentQualityDistribution( SingleTableTrackerHistogram ):
#    table = "alignment_stats_quality_distribution_metrics"
#    column = "quality"
    
