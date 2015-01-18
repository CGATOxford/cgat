from MedipReport import *


class CpGCoverageByReads(ProjectTracker, SingleTableTrackerHistogram):
    '''CpG coverage in regions in the genome that are covered
    by tags.'''

    table = "cpg_coverage_by_reads"
    column = "bin"


class CpGContext(ProjectTracker, SingleTableTrackerRows):
    table = "cpg_context"
    fields = ('category',)


class CpGDistributionInTags(ProjectTracker, SingleTableTrackerHistogram):
    '''distribution of percent CpG in regions covered by reads.'''
    table = "pcpg_in_coveredregions"
    column = "bin"

