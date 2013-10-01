from MedipReport import *

class CpGCoverage( ProjectTracker, SingleTableTrackerHistogram):
    table = "cpg_coverage"
    column = "bin"
    

class CpGContext( ProjectTracker, SingleTableTrackerRows ):
    table = "cpg_context"
    fields = ('category',)

