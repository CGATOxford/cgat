from RnaseqDiffExpressionReport import *


class TopDifferentiallyExpressedGenes(ProjectTracker):

    '''output differentially expressed genes.'''
    limit = 10

    pattern = '(.*)_gene_diff'
    sort = ''

    def __call__(self, track, slice=None):
        statement = '''SELECT a.gene_name,
                              a.gene_id,
                              a.gene_biotype,
                              t.l2fold,
                              t.treatment_mean,
                              t.control_mean,
                              t.pvalue,
                              t.qvalue,
                              s.contig, s.start, s.end
                              FROM %(track)s_gene_diff as t,
                                   annotations.gene_info as a,
                                   annotations.gene_stats as s
                              WHERE a.gene_id = t.test_id AND
                                    s.gene_id = t.test_id AND
                                    t.significant
                              ORDER BY %(sort)s
                              LIMIT %(limit)i'''
        data = self.getAll(statement)
        if len(data) == 0:
            return data

        data['gene_id'] = [linkToEnsembl(x) for x in data["gene_id"]]
        data["locus"] = [linkToUCSC(*x) for x in zip(
            data["contig"],
            data["start"],
            data["end"])]
        return data


class TopUpRegulatedGenes(TopDifferentiallyExpressedGenes):
    sort = 't.l2fold DESC'


class TopDownRegulatedGenes(TopDifferentiallyExpressedGenes):
    sort = 't.l2fold Asc'


class AllDifferentiallyExpressedGenes(ProjectTracker):

    '''output differentially expressed genes.'''
    limit = 1000

    pattern = '(.*)_gene_diff'

    def __call__(self, track, slice=None):
        statement = '''SELECT a.gene_name, 
                              a.gene_id,
                              a.gene_biotype, 
                              t.l2fold,
                              t.treatment_mean,
                              t.control_mean,
                              t.pvalue,
                              t.qvalue,
                              s.contig, s.start, s.end
                              FROM %(track)s_gene_diff as t, 
                                   annotations.gene_info as a,
                                   annotations.gene_stats as s
                              WHERE a.gene_id = t.test_id AND
                                    s.gene_id = t.test_id AND
                                    t.significant
                              ORDER BY t.l2fold DESC LIMIT %(limit)i'''
        data = self.getAll(statement)
        if len(data) == 0:
            return data

        data['gene_id'] = [linkToEnsembl(x) for x in data["gene_id"]]
        data["locus"] = [linkToUCSC(*x) for x in zip(
            data["contig"],
            data["start"],
            data["end"])]
        return data
