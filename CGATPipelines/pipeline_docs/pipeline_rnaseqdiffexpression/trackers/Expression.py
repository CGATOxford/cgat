
from CGATReport.Tracker import *
from collections import OrderedDict as odict

from RnaseqDiffExpressionReport import *

import numpy


class TrackerGeneCounts(ProjectTracker):

    '''tracker with gene counts.'''
    pattern = "(\S+)_gene_counts$"


class TrackerRealEstate(TrackerGeneCounts):
    '''return proportion of reads mapped to genes stratified by
    expression level.'''

    slices = ("counted_all",)

    def __call__(self, track, slice):
        statement = '''
        SELECT count(gene_id) AS ngenes,
               %(slice)s * count(gene_id) AS reads_times_genes
        FROM %(track)s_gene_counts
        GROUP BY %(slice)s
        ORDER BY %(slice)s DESC'''

        data = self.getAll(statement)
        r = numpy.array(data["reads_times_genes"], dtype=numpy.float)
        data["reads_times_genes"] = r.cumsum() / r.sum() * 100.0

        r = numpy.array(data["ngenes"], dtype=numpy.float)
        data["ngenes"] = r.cumsum()

        return data


class TrackerGeneFPKM(ProjectTracker):
    min_fpkm = 1
    pattern = "(.*)_genefpkm$"
    as_tables = True


class TrackerTranscriptFPKM(ProjectTracker):
    min_fpkm = 1
    pattern = "(.*)_fpkm$"
    as_tables = True


class ExpressionFPKM(TrackerGeneFPKM):

    def __call__(self, track, slice=None):
        statement = '''SELECT fpkm FROM %(track)s WHERE fpkm > %(min_fpkm)f'''
        return odict((("fpkm", self.getValues(statement)),))


class ExpressionGenesLog2FPKM(TrackerGeneFPKM):
    '''log2 transformed FPKM values with a +1 pseudocount'''
    def __call__(self, track, slice=None):
        statement = '''SELECT fpkm FROM %(track)s'''
        data = numpy.array(self.getValues(statement),
                           dtype=numpy.float)
        data = data + 1.0
        data = numpy.log2(data)
        return odict((("fpkm", data),))


class ExpressionTranscriptsLog2FPKM(TrackerTranscriptFPKM):
    '''log2 transformed FPKM values with a +1 pseudocount'''
    def __call__(self, track, slice=None):
        statement = '''SELECT fpkm FROM %(track)s'''
        data = numpy.array(self.getValues(statement),
                           dtype=numpy.float)
        data = data + 1.0
        data = numpy.log2(data)
        return odict((("fpkm", data),))


class ExpressionNormalizedFPKM(TrackerGeneFPKM):

    def __call__(self, track, slice=None):
        max_fpkm = float(self.getValue( '''SELECT max(fpkm) 
        FROM %(track)s'''))
        statement = '''SELECT CAST( fpkm AS FLOAT) / %(max_fpkm)f FROM %(track)s WHERE fpkm > %(min_fpkm)f'''
        return odict((("percent of max(fpkm)", self.getValues(statement)),))


class ExpressionFPKMConfidence(TrackerGeneFPKM):

    def __call__(self, track, slice=None):
        # divide by two to get relative error
        statement = '''SELECT (fpkm_conf_hi - fpkm_conf_lo ) / fpkm / 2
                       FROM %(track)s WHERE fpkm > %(min_fpkm)f'''
        return odict((("relative_error", self.getValues(statement)),))


class ExpressionFPKMConfidenceCorrelation(TrackerGeneFPKM):

    def __call__(self, track, slice=None):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns(table):
            return None
        statement = '''SELECT %(slice)s_fpkm AS fpkm,
        (%(slice)s_conf_hi - %(slice)s_conf_lo ) AS confidence
        FROM %(table)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getAll(statement)
        return data


class ExpressionHighestExpressedGenes(TrackerGeneFPKM):

    '''output ten highest expressed genes.'''
    limit = 10

    def __call__(self, track, slice=None):
        statement = '''SELECT gene_id, gene_short_name, locus, fpkm
                              FROM %(track)s
                              ORDER BY fpkm DESC LIMIT %(limit)i'''
        data = self.getAll(statement)

        data['gene_id'] = [linkToEnsembl(x) for x in data["gene_id"]]
        data["locus"] = [linkToUCSC(*splitLocus(x)) for x in data["locus"]]
        return data


class ExpressionHighestExpressedGenesDetailed(TrackerGeneFPKM):

    '''output ten highest expressed genes.'''
    limit = 10

    def __call__(self, track, slice=None):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        geneset = track[:track.index("_")]
        if c not in self.getColumns(table):
            return None
        statement = '''SELECT tracking_id, gene_short_name, locus, source,
        class, sense,
        %(slice)s_fpkm AS fpkm,
        mappability.median AS median_mappability
        FROM %(table)s,
        %(geneset)s_class AS class,
        %(geneset)s_mappability AS mappability
        WHERE tracking_id = class.transcript_id AND
        tracking_id = mappability.transcript_id
        ORDER BY %(slice)s_fpkm DESC LIMIT %(limit)i'''

        data = self.getAll(statement)
        data['tracking_id'] = [linkToEnsembl(x) for x in data["tracking_id"]]
        data["locus"] = [linkToUCSC(*splitLocus(x)) for x in data["locus"]]

        return data
