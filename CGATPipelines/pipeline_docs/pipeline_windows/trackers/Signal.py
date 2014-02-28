from MedipReport import *


class SignalMedian(ProjectTracker, SingleTableTrackerColumns):
    table = "counts_l2foldchange_median"
    column = None


class SignalInput(ProjectTracker, SingleTableTrackerColumns):
    table = "counts_l2foldchange_input"
    column = None
