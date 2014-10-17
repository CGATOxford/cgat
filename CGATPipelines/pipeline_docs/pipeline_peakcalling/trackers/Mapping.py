from PeakcallingReport import *
from CGATReport.Tracker import *


class BamSummary(CallingTracker, SingleTableTrackerRows):
    table = "bam_stats"
