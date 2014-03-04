from PeakcallingReport import *
from SphinxReport.Tracker import *


class BamSummary(CallingTracker, SingleTableTrackerRows):
    table = "bam_stats"
