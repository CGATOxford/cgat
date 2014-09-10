from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from collections import OrderedDict as odict

###################################################################
# parameterization

EXPORTDIR = P.get(
    'annotations_exportdir',
    P.get('exportdir', 'export'))
DATADIR = P.get(
    'annotations_datadir',
    P.get('datadir', '.'))
DATABASE = P.get(
    'annotations_backend',
    P.get('sql_backend', 'sqlite:///./csvdb'))


class AnnotationTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class BedSummaryIntervals(TrackerSQL):
    pattern = '(.*)_bedsummary'

    def __call__(self, track):
        return self.getValue(
            "SELECT sum(nintervals) FROM %(track)s_bedsummary")


class BedSummaryIntervalsPerContig(TrackerSQL):
    pattern = '(.*)_bedsummary'

    def __call__(self, track):
        return self.getAll(
            "SELECT track as contig, nintervals FROM %(track)s_bedsummary")


class BedSummaryBases(TrackerSQL):
    pattern = '(.*)_bedsummary'

    def __call__(self, track):
        return self.getValue(
            "SELECT SUM(nbases) FROM %(track)s_bedsummary")


class BedSummaryBasesPerContig(TrackerSQL):
    pattern = '(.*)_bedsummary'

    def __call__(self, track):
        return self.getAll(
            "SELECT track as contig, nbases FROM %(track)s_bedsummary")


class GTFSummaryPerSource(TrackerSQL):
    pattern = '(.*)_(gtf|gff)summary'

    def __call__(self, track, slice):
        if not self.hasTable("%s_%ssummary" % (track, slice)):
            return None

        return self.getAll(
            """SELECT source,
            SUM(bases) as bases,
            SUM(intervals) as intervals,
            SUM(total_percent_coverage) as percent_coverage
            FROM %(track)s_%(slice)ssummary
            GROUP BY source""")


class GTFStats(SingleTableTrackerRows):
    table = 'gtf_stats'


class GFFStats(SingleTableTrackerRows):
    table = 'gff_stats'


class InfoTable(TrackerSQLCheckTables):
    pattern = "(.*_info|.*_gtf)$"
    slices = ["gene_id", "transcript_id", "gene_name"]
