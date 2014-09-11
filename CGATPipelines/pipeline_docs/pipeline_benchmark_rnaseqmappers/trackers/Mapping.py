import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from BenchmarkReport import *


class MappingSummary(BenchmarkTracker, SingleTableTrackerRows):
    table = "view_mapping"


class TophatSummary(BenchmarkTracker, SingleTableTrackerRows):
    table = "tophat_stats"


class BamSummary(BenchmarkTracker, SingleTableTrackerRows):
    table = "bam_stats"


class AlignmentSummary(BenchmarkTracker, SingleTableTrackerRows):
    table = "alignment_stats"


class MappingContext(BenchmarkTracker, SingleTableTrackerRows):
    table = "context_stats"


class FilteringSummary(BenchmarkTracker, SingleTableTrackerRows):
    table = "mapping_stats"


class MappingFlagsMismatches(BenchmarkTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nm"
    column = "nm"


class MappingFlagsHits(BenchmarkTracker, SingleTableTrackerHistogram):
    table = "bam_stats_nh"
    column = "nh"


class AlignmentQualityByCycle(BenchmarkTracker, SingleTableTrackerHistogram):
    table = "alignment_stats_quality_by_cycle_metrics"
    column = "cycle"


class AlignmentQualityDistribution(BenchmarkTracker, SingleTableTrackerHistogram):
    table = "alignment_stats_quality_distribution_metrics"
    column = "quality"

##############################################################
##############################################################
##############################################################


class BamReport(BenchmarkTracker):
    tracks = ["all"]

    slices = ("genome", "accepted", "mismapped")

    def __call__(self, track, slice=None):
        edir = EXPORTDIR

        toc_text = []
        link_text = []

        filenames = sorted([x.asFile() for x in TRACKS])

        for fn in filenames:
            link = "%(slice)s_%(fn)s" % locals()
            toc_text.append("* %(link)s_" % locals())
            link_text.append(
                ".. _%(link)s: %(edir)s/bamstats/%(fn)s.%(slice)s.html" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))


##############################################################
##############################################################
##############################################################
class FastQCReport(BenchmarkTracker):
    tracks = ["all"]

    def __call__(self, track, slice=None):
        edir = EXPORTDIR

        toc_text = []
        link_text = []

        filenames = sorted([x.asFile() for x in TRACKS])

        for fn in filenames:
            toc_text.append("* %(fn)s_" % locals())
            link_text.append(
                ".. _%(fn)s: %(edir)s/fastqc/%(fn)s_fastqc/fastqc_report.html" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))
