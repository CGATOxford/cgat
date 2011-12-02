from IntervalReport import *

import IOTools

class BindingPatterns(IntervalTracker):
    '''output summary counts of binding patterns.
    
    The empty pattern is excluded.
    '''
    pattern = "(.*)_binding"
    
    def __call__(self, track ):
        return self.getAll( """SELECT pattern, COUNT(*) as counts 
                             FROM %(track)s_binding  
                             WHERE CAST(pattern AS INT) != 0
                             GROUP BY pattern""" )

class BindingSummary(IntervalTracker):
    '''output summary counts of binding patterns.'''
    pattern = "(.*)_binding"
    
    def __call__(self, track ):
        
        data = odict()

        # data["no binding"]  = self.getValue( """
        #                       SELECT COUNT(*) FROM %(track)s_binding 
        #                       WHERE CAST(pattern AS INT) == 0""" )

        for section in ( "flank5", "utr5", "cds", "intron", "utr3", "flank3" ):
            data[section] = self.getValue( """
                              SELECT COUNT(*) FROM %(track)s_binding 
                              WHERE %(section)s_overlap > 0""" )
        return data
