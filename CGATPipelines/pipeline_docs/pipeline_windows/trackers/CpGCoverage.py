from MedipReport import *

class CpGCoverage( ProjectTracker, SingleTableTrackerHistogram):
    table = "cpg_coverage"
    column = "bin"
    

class CpGContext( ProjectTracker, SingleTableTrackerRows ):
    table = "cpg_context"
    fields = ('category',)

class CpGDistribution( ProjectTracker ):
    '''distribution of percent CpG in regions covered by reads.'''

    pattern = "(.*)_covered_composition"
    
    def __call__(self,track):
        return self.getValues( "SELECT pCpG FROM %(track)s_covered_composition" )

