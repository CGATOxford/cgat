import os
import sys
import re
import types
from VariantsReport import *

coding_codes = "abceABCE"

##########################################################################
##########################################################################
##########################################################################
# Exon positions
##########################################################################


class ExonEssentiality(VariantsTracker):

    """types of snps.
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        codes = "','".join(coding_codes)
        data = self.get(
            "SELECT exons_nused, exons_ntranscripts FROM %(track)s_annotations WHERE code IN ('%(codes)s')" % locals())
        return odict((("used", "ntranscript"), zip(*data)))


##########################################################################
##########################################################################
##########################################################################
# Exon positions
##########################################################################
class EffectCounts(VariantsTracker):

    """counts of effects. Substitutions are counted as:
    synonymous
    non-synonymous
    truncating
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        data = []
        data.append(("synonymous",
                     self.getValue(
                         "SELECT COUNT(*) FROM %(track)s_annotations WHERE variant_type = 'S' AND reference_aa = consensus_aa" % locals())
                     ))
        data.append(("nonsynonymous",
                     self.getValue(
                         "SELECT COUNT(*) FROM %(track)s_annotations WHERE variant_type = 'S' AND reference_aa != consensus_aa" % locals())
                     ))

        return odict((("used", "ntranscript"), zip(*data)))
