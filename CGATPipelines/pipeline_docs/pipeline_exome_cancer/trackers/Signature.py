from SphinxReport.Tracker import *
from exomeReport import *


#class signature(ExomeTracker):
#
#    pattern = "(.*)_mutect_snp_annotated_tsv$"
#
#    def __call__(self, track, slice=None):
#
#        statement = '''
#        SELECT A.REF, A.ALT, COUNT(*)
#        FROM %(track)s_mutect_snp_annotated_tsv AS A
#        JOIN %(track)s_call_stats_out AS B
#        ON A.CHROM = B.contig AND A.POS = B.position
#        WHERE A.FILTER!='REJECT' AND B.t_alt_count > 4
#        GROUP BY A.REF, A.ALT;
#        ''' % locals()
#
#        return self.getAll(statement)


class SnpSummary(ExomeTracker):

    def __call__(self, track, slice=None):
        table = "mutational_signature"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        return self.getAll(statement)


class SnpSummaryTable(ExomeTracker):

    def __call__(self, track, slice=None):
        table = "mutational_signature_table"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        return self.getAll(statement)
