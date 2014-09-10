from MedipReport import ProjectTracker
from CGATReport.Tracker import SingleTableTrackerRows
from CGATReport.Tracker import SingleTableTrackerHistogram


class TileSizeHistogram(SingleTableTrackerHistogram, ProjectTracker):
    table = "tileinfo_hist"
    column = "residues"


class TileSizeStats(SingleTableTrackerRows, ProjectTracker):
    table = "tileinfo_stats"
