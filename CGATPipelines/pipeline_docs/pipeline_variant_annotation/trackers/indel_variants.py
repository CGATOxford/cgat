import os
import sys
import re
import types
from VariantsReport import *

#####################################################
#####################################################
#####################################################


class TrackerEffects(VariantsTracker):

    # minimum number of truncated codons
    min_truncated = 5

    mPattern = "_effects$"

    def getPrefix(self, slice):
        if slice == None or slice == "all":
            prefix = ""
        else:
            prefix = "%s_" % slice
        return prefix

    def getSlices(self, subset=None):
        if subset == None:
            return []
        elif "separate" in subset:
            return ("all", "splice", "cds")
        return subset


#####################################################
#####################################################
#####################################################
class TranscriptListFrameshift(TrackerEffects):

    '''output a genelist of genes with frameshifts'''

    mPattern = "_effects_cds$"

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "contig",
                   "snp_position", "exon_id", "nexons", "reference", "Variant_type", "indel", "effect")

        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            e.contig,
            e.snp_position,
            e.exon_id,
            e.nexons,
            e.reference,
            e.variant_type,
            e.variant_bases as indel,
            'Frameshift' as effect
        FROM
            %(track)s_effects_cds AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        AND e.code='F'
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))

#####################################################
#####################################################
#####################################################


class TranscriptListDeletions(TrackerEffects):

    '''output a genelist of genes with in frame deletions'''

    mPattern = "_effects_cds$"

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "contig",
                   "snp_position", "exon_id", "nexons", "reference", "Variant_type", "indel", "effect")

        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            e.contig,
            e.snp_position,
            e.exon_id,
            e.nexons,
            e.reference,
            e.variant_type,
            e.variant_bases as indel,
            'in-frame deletion' as effect
        FROM
            %(track)s_effects_cds AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        AND e.code='D'
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))


#####################################################
#####################################################
#####################################################
class TranscriptListInsertions(TrackerEffects):

    '''output a genelist of genes with in frame deletions'''

    mPattern = "_effects_cds$"

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "contig",
                   "snp_position", "exon_id", "nexons", "reference", "Variant_type", "indel", "effect")

        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            e.contig,
            e.snp_position,
            e.exon_id,
            e.nexons,
            e.reference,
            e.variant_type,
            e.variant_bases as indel,
            'in-frame insertion' as effect
        FROM
            %(track)s_effects_cds AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        AND e.code='I'
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))
