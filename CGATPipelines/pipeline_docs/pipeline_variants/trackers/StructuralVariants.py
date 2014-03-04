import os
import sys
import re
import types
from VariantsReport import *

SV_SLICES = ["all", "other", "del", "ins"]


class StructuralVariantsSummary(VariantsTracker):
    mPattern = "_sv_stats$"

    def __call__(self, track, slice=None):

        columns = "event", "name", "counts", "nval", "mean", "median", "sum"
        cols = ",".join(columns)
        data = self.get( '''SELECT %(cols)s
                               FROM %(track)s_sv_stats''' % locals() )
        return odict(zip(columns, zip(*data)))


class StructuralVariantsSummaryGenes(VariantsTracker):
    mPattern = "_sv_genes$"

    def __call__(self, track, slice=None):

        columns = ("naffected_ins",
                   "naffected_del",
                   "naffected_other",
                   "naffected_all",
                   "ndeleted",
                   "is_ins",
                   "is_del",
                   "is_other",
                   "is_all",
                   "is_deleted")

        cols = ",".join(["SUM(%s) AS %s" % (x, x) for x in columns])
        data = self.get( '''SELECT %(cols)s
                               FROM %(track)s_sv_genes''' % locals() )
        return odict(zip(columns, zip(*data)))


class StructuralVariantsDeletedGenesExons(VariantsTracker):

    '''return (max) number of exons for deleted genes.'''
    mPattern = "_sv_genes$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT max(nval) 
                                   FROM annotations.transcript_info as i, 
                                   %(track)s_sv_genes as s, 
                                   transcript_stats as a 
                                  WHERE i.gene_id = s.gene_id and 
                                        s.is_deleted and 
                                        a.transcript_id = i.transcript_id 
                                  GROUP by s.gene_id''' % locals() )
        return odict((("exons", data),))


class StructuralVariantsDeletedGenes(VariantsTracker):

    '''return genes that have been deleted through structural variants.

    Exclude genes that 
    * have only a single exon
    * have the mgi type "Pseudogene"

    '''
    mPattern = "_sv_genes$"

    def __call__(self, track, slice=None):
        statement = '''
        SELECT DISTINCT i.gene_id, i.gene_name, b.description 
        FROM %(track)s_sv_genes as g, 
        annotations.transcript_info AS i, 
        transcript_stats AS s, 
        gene_status as b 
        LEFT JOIN mgi_markers as m on m.symbol = i.gene_name 
        WHERE is_deleted AND 
        s.nval > 1 AND 
        i.transcript_id = s.transcript_id AND 
        i.gene_id = g.gene_id AND 
        b.gene_id = i.gene_id AND 
        (m.type != 'Pseudogene' OR m.type IS NULL)''' % locals()

        data = self.get(statement)

        return odict(zip(("gene_id",
                          "gene_name",
                          "description"),
                         zip(*data)))


class StructuralVariantsTracker(VariantsTracker):

    def getSlices(self, subset=None):
        return SV_SLICES


class StructuralVariantsSummaryTranscripts(StrainTracker):

    def __call__(self, track, slice=None):

        if track == "all":
            return

        result = odict()
        for slice in SV_SLICES:
            statement = """SELECT DISTINCT COUNT(*) FROM %(track)s_sv_overlap_%(slice)s
                           WHERE nover > 0
                        """ % locals()
            result[slice] = self.getValue(statement)

        return result


class StructuralVariantsTranscripts(StrainTracker, StructuralVariantsTracker):

    def __call__(self, track, slice=None):

        if track == "all":
            return
        if slice == "all":
            return

        statement = """SELECT gene_id, gene_name, i.transcript_id, nover, nover1, nover2, pover1, pover2
                           FROM %(track)s_sv_overlap_%(slice)s AS a,
                                annotations.transcript_info AS i
                           WHERE nover > 0 AND i.transcript_id = a.transcript_id
                        """ % locals()

        data = self.get(statement)

        return odict(zip(("transcript_id", "nover", "nover1", "nover2", "pover1", "pover2"),
                         zip(*data)))
