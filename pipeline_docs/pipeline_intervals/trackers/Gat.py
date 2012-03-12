import os, sys, re, types, itertools, math

from IntervalReport import *

class GatTracker( IntervalTracker ):
    pass

class GatGenomicContextTable( GatTracker ):
    pattern = "gat_context_(.*)$"

    def __call__(self, track):
        return self.getAll( "SELECT * FROM gat_context_%(track)s" )

class GatGenomicAnnotationTable( GatTracker ):
    pattern = "gat_annotation_(.*)$"

    def __call__(self, track):
        return self.getAll( "SELECT * FROM gat_annotation_%(track)s" )


##################################################################################
##################################################################################
##################################################################################
## GAT results
##################################################################################
class GatResults( IntervalTracker, SingleTableTrackerRows ):
    '''All gat results.'''
    fields = ('track', 'annotation')
    extra_columns = { "colour" : "CASE WHEN qvalue < 0.05 THEN 'red' ELSE 'blue' END" }
    sort = 'l2fold'

class GatLogFold( IntervalTracker ):
    pattern = "gat_(.*)"

    fdr = 2.0

    def __call__(self, track ):
        return self.getDict( """SELECT annotation, fold, 
                                       CASE WHEN qvalue < %(fdr)f THEN 'red' ELSE 'blue' END AS colour
                               FROM gat_%(track)s ORDER BY fold""")
