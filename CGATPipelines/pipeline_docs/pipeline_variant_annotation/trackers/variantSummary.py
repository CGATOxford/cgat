import os
import sys
import re
import types
import collections
from VariantsReport import *

SNP_CODES = ("E", "O")

##########################################################################
##########################################################################
##########################################################################
# Variant types
##########################################################################


class VariantTypeCounts(VariantsTracker):

    """return counts for each variant type. The counts are symmetrized, i.e., W,D = D,W"""

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        data = self.get(
            "SELECT variant_type, COUNT(*) FROM %(track)s_annotations GROUP BY variant_type" % locals())
        combined = collections.defaultdict(int)
        for key, counts in data:
            try:
                if "," in key:
                    key = ",".join(sorted(key.split(",")))
            except TypeError:
                pass
            combined[key] += counts
        return combined

##########################################################################
##########################################################################
##########################################################################
# Genomic locations of variants
##########################################################################


class VariantLocation(VariantsTracker):

    """Genomic location of sequence variants"""

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):
        data = self.get(
            """SELECT UPPER(code), COUNT(*) FROM %(track)s_annotations GROUP BY code""" % locals() )
        return odict(data)
