import math

from GeneSetsReport import *


class GeneSetsSummaryForeground(GeneSetsTracker):
    '''Size of genesets used in the analyses.'''
    table = "genelists_foreground"

    def getTracks(self):
        return [x for x in self.getColumns(self.table) if x != "id"]

    def __call__(self, track):
        return self.getValue("SELECT SUM(%(track)s) FROM %(table)s")


class GeneSetsSummaryBackground(GeneSetsSummaryForeground):
    '''Size of background genesets used in the analyses.'''
    table = "genelists_background"


class GOAnalysisComparison(GeneSetsTracker):

    limit = 100
    tracks = ('biol_process', 'mol_function')
    # only select with fdr > 1% (-10logqvalue >= 20)

    # 5% threshold
    fdr_threshold = 0.05

    def __call__(self, track, slice):
        columns = ",".join(["fold.%s AS %s_fold" % (x, x)
                           for x in self.columns])

        # remove the GO identifier
        # convert to phred like fdr scores
        fdr = - 10.0 * math.log(self.fdr_threshold, 10)

        return self.getAll('''SELECT substr( fold.category,12) AS description, 
        %(columns)s
        FROM %(table_prefix)s_%(track)s_l2fold AS fold,
        %(table_prefix)s_%(track)s_l10qvalue AS qvalue
        WHERE fold.category = qvalue.category
        AND qvalue.%(slice)s >= %(fdr)f
        ORDER BY fold.%(slice)s DESC
        LIMIT %(limit)i
        ''')


class GOAnalysisComparisonPromotor(GOAnalysisComparison):

    '''gene lists built by promotor definition.'''

    columns = ('irf5_rela_upstream_all',
               #'irf5_rela_overlap_upstream_all',
               #'irf5_rela_nooverlap_upstream_all',
               'rela_noirf5_upstream_all',
               'irf5_norela_upstream_all')

    slices = ('irf5_rela_upstream_all',
              'irf5_rela_overlap_upstream_all',
              'irf5_rela_nooverlap_upstream_all',
              'rela_noirf5_upstream_all',
              'irf5_norela_upstream_all')

    table_prefix = "genelists_goall_alldesc"


