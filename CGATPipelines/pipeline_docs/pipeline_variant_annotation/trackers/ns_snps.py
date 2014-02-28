import os
import sys
import re
import types
from VariantsReport import *

from Effects import TrackerEffects

#####################################################
#####################################################
#####################################################


class TranscriptListSNP(TrackerEffects):

    '''output a genelist of genes with nonsynonymous SNPs'''

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "transcript_id", "chr", "snp_position",
                   "reference", "ref_aa", "Variant_type", "genotype", "variant_aa")

        #field = self.getPrefix(slice) + "stop_min"
        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.contig,
            e.snp_position,
            e.reference,
            e.orig_codons as ref_aa,
            e.variant_type,
            e.variant_bases as genotype, 
            e.variant_codons as variant_aa
        FROM
            %(track)s_effects_cds AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        AND e.code='N'
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))
