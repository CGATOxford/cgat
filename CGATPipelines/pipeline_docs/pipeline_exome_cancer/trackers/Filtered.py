import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class snp(ExomeTracker):

    pattern = "(.*)_mutect_snp_annotated_tsv$"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.CHROM AS Chr, A.POS AS Pos,
        A.SNPEFF_GENE_NAME AS Gene,
        A.SNPEFF_EXON_ID AS Exon,
        A.REF, A.ALT,
        A.SNPEFF_IMPACT AS Impact, A.SNPEFF_GENE_BIOTYPE AS Biotype,
        SNPEFF_AMINO_ACID_CHANGE AS AA_change,
        SNPEFF_CODON_CHANGE AS Codon_change,
        C.function, C.type as NCG, C.cancer_type,
        C.duplicability,
        B.n_ref_count AS Normal_Ref, B.n_alt_count AS Normal_Alt,
        B.t_ref_count AS Tumor_Ref, B.t_alt_count AS Tumor_Alt
        FROM %(track)s_mutect_snp_annotated_tsv AS A
        JOIN %(track)s_call_stats_out AS B
        JOIN cancergenes as C
        ON A.CHROM = B.contig AND A.POS = B.position
        AND A.SNPEFF_GENE_NAME = C.symbol
        WHERE A.FILTER!="REJECT" AND B.t_alt_count > 4;
        ''' % locals()

        return self.getAll(statement)


class indel(ExomeTracker):

    pattern = "(.*)_indels_annotated_tsv$"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT A.CHROM AS Chr, A.POS AS Pos,
        A.SNPEFF_GENE_NAME AS Gene,
        A.SNPEFF_EXON_ID AS Exon,
        A.REF, A.ALT,
        A.SNPEFF_IMPACT AS Impact, A.SNPEFF_GENE_BIOTYPE AS Biotype,
        A.SNPEFF_AMINO_ACID_CHANGE AS AA_change,
        A.SNPEFF_CODON_CHANGE AS Codon_change,
        B.function, B.type as NCG, B.cancer_type,
        B.duplicability,
        A.NORMAL_DP AS Normal_depth,
        A.TUMOR_DP AS Tumor_depth,
        A.NORMAL_TAR as Normal_Ref, A.NORMAL_TIR as Normal_Alt,
        A.TUMOR_TAR as Tumor_Ref, A.TUMOR_TIR as Tumor_Alt
        FROM %(track)s_indels_annotated_tsv AS A
        JOIN cancergenes as B
        ON A.SNPEFF_GENE_NAME = B.symbol
        WHERE A.QSI_NT > 20 AND A.IHP < 12
        AND A.RC < 12 AND A.IC < 12;
        ''' % locals()

        return self.getAll(statement)


class filterSummary(ExomeTracker):

    pattern = "(.*)_mutect_filtering_summary$"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT * FROM %(track)s_mutect_filtering_summary
        ;
        ''' % locals()

        return self.getAll(statement)
