from GeneSetsReport import *

import CGAT.IOTools as IOTools


class HypergeometricResults(GeneSetsTracker):
    pattern = "(\S+)_results"

    def __call__(self, track):
        return self.getAll(
            """SELECT
            track as geneset,
            category, goid, ratio, pvalue, fdr, description
            FROM %(track)s_results
            ORDER BY fdr
            """)


class FoldMatrix(GeneSetsTracker):

    tracks = ("pathways", "go")
    slices = ("biol_process", "cell_location")

    def __call__(self, track, slice):
        fold_matrix, rownames, colnames = IOTools.readMatrix(
            IOTools.openFile(
                "hypergeometric.dir/%s.dir/alldesc.%s.fold" %
                (track, slice)))

        return odict((('matrix', fold_matrix),
                      ('rows', rownames),
                      ('columns', colnames)))


class GOComparison(GeneSetsTracker):
    pattern = "pairs_(.*)"

    qvalue = 0.05

    def __call__(self, track):
        return self.getAll(
            "SELECT * FROM pairs_%(track)s WHERE qvalue < %(qvalue)f")


