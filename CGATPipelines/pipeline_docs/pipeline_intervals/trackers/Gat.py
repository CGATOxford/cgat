from IntervalReport import *


# class GatGenomicContextTable( GatTracker ):
#     pattern = "gat_context_(.*)$"

#     def __call__(self, track):
#         return self.getAll( "SELECT * FROM gat_context_%(track)s" )

# class GatGenomicAnnotationTable( GatTracker ):
#     pattern = "gat_annotations_(.*)$"

#     def __call__(self, track):
#         return self.getAll( "SELECT * FROM gat_annotations_%(track)s" )

##########################################################################
##########################################################################
##########################################################################
# GAT results
##########################################################################
# class GatResults( IntervalTracker, SingleTableTrackerRows ):
#     '''All gat results.'''
#     fields = ('track', 'annotation')
#     extra_columns = { "colour" : "CASE WHEN qvalue < 0.05 THEN 'red' ELSE 'blue' END" }
#     sort = 'l2fold'


class GatFold(IntervalTracker, SingleTableTrackerEdgeList):

    '''fold change matrix.'''
    row = "track"
    column = "annotation"
    value = "fold"
    where = "pvalue < 0.05"


class GatLogFold(IntervalTracker):

    '''logfold - colour is signficance'''
    fdr = 0.05
    as_tables = True

    def __call__(self, track):
        return self.getDict("""SELECT annotation, l2fold,
        CASE WHEN qvalue < %(fdr)f THEN 'red' ELSE 'gray' END AS colour
        FROM %(track)s ORDER BY fold""")


class GatResults(IntervalTracker):
    as_tables = True

    def __call__(self, track):
        return self.getAll("SELECT * FROM %(track)s")


class GatSignificantResults(IntervalTracker):
    as_tables = True
    fdr = 0.05

    def __call__(self, track):
        return self.getAll("SELECT * FROM %(track)s WHERE qvalue < %(fdr)f")


class GatTableAnnotations:
    pattern = "gat_annotations_(.*)"


class GatTableContext:
    pattern = "gat_context_(.*)"


class GatTableFunctions:
    pattern = "gat_functions_(.*)"

_gat_analysis = {"Results": GatResults,
                 "SignificantResults": GatSignificantResults,
                 "Fold": GatLogFold,
                 "LogFold": GatLogFold}

_gat_sets = {"Annotations": GatTableAnnotations,
             "Context": GatTableContext,
             "Functions": GatTableFunctions,
             }

for a, aa in _gat_analysis.items():
    for b, bb in _gat_sets.items():
        n = "Gat%s%s" % (a, b)
        globals()[n] = type(n, (bb, aa), {})
