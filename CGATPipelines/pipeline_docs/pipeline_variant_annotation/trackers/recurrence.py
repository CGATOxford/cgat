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
        if slice is None or slice == "all":
            prefix = ""
        else:
            prefix = "%s_" % slice
        return prefix

    def getSlices(self, subset=None):
        if subset is None:
            return []
        elif "separate" in subset:
            return ("all", "splice", "cds")
        return subset


#####################################################
#####################################################
#####################################################
class RecurrentVariants(TrackerEffects):

    '''output a list of variants which occur in >1 sample'''

    mPattern = "^annotations$"

    def __call__(self, track, slice=None):

        headers = ("chromosome", "position", "genotype",
                   "reference", "Variant_type", "code", "samples")

        statement = '''
        SELECT * FROM (SELECT
            a.chromosome,
            a.position,
            a.genotype,
            a.reference_base,
            a.variant_type,
            a.code,
            count(a.track) as samples
        FROM annotations a
        GROUP BY a.chromosome, a.position
        ORDER BY samples desc, a.chromosome, a.position)
        WHERE samples >1
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))

#####################################################
#####################################################
#####################################################


class RecurrentEffects(TrackerEffects):

    '''output a list of genes which are disrupted in >1 sample'''

    mPattern = "^effects_genes$"

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "samples")

        statement = '''
        SELECT * FROM (SELECT
            i.gene_id,
            i.gene_name,
            count( distinct e.track) as samples
        FROM
            effects_genes AS e,
            annotations.transcript_info AS i
        WHERE i.gene_id = e.gene_id
        AND min_nalleles >0
        GROUP BY i.gene_id
        ORDER BY i.gene_name)
        WHERE samples >1
        ''' % self.members(locals())

        return odict(zip(headers, zip(*self.get(statement))))
