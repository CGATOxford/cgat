from MedipReport import ProjectTracker
from SphinxReport.Tracker import SingleTableTrackerRows
from SphinxReport.Tracker import SingleTableTrackerHistogram


class TileSizeHistogram(SingleTableTrackerHistogram, ProjectTracker):
    table = "tileinfo_hist"
    column = "residues"


class TileSizeStats(SingleTableTrackerRows, ProjectTracker):
    table = "tileinfo_stats"
