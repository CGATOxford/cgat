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


class CDSCountsTranscripts(VariantsTracker):

    '''number of transcripts with certain features.'''
    mPattern = "_effects_translation$"

    def __call__(self, track, slice=None):
        columns = ("ncorrected_frames",
                   "ndeleted_bases",
                   "ndeleted_codons",
                   "nincomplete_codons",
                   "ninserted_bases",
                   "ninserted_codons",
                   "nnonsynonymous_codons",
                   "noframe_codons",
                   "nstop_codons",
                   "nstops",
                   "nsynonymous_codons",
                   "nunaffected_codons",
                   "nwrong_frames")

        result = odict()
        for column in columns:
            result[column] = self.getValue('''SELECT COUNT(*) FROM %(track)s_effects_translation
                                    WHERE %(column)s > 0''' % locals() )
        return result

#####################################################
#####################################################
#####################################################


class CDSCountsVariants(VariantsTracker):

    '''number of transcripts with certain features.'''
    mPattern = "_effects_translation$"

    def __call__(self, track, slice=None):
        columns = ("ncorrected_frames",
                   "ndeleted_bases",
                   "ndeleted_codons",
                   "nincomplete_codons",
                   "ninserted_bases",
                   "ninserted_codons",
                   "nnonsynonymous_codons",
                   "noframe_codons",
                   "nstop_codons",
                   "nstops",
                   "nsynonymous_codons",
                   "nunaffected_codons",
                   "nwrong_frames")

        select = ",".join(["SUM(%s)" % x for x in columns])
        statement = '''SELECT %(select)s FROM %(track)s_effects_translation''' % locals(
        )

        return odict(zip(columns, self.getFirstRow(statement)))

#####################################################
#####################################################
#####################################################


class TrackerVariants(VariantsTracker):

    '''type of cds variants.'''

    mPattern = "_effects_cds$"
    column = None
    aggregate = None

    def process(self, data):
        return data

    def __call__(self, track, slice=None):
        if not self.column:
            raise NotImplementedError

        if self.aggregate == None or self.aggregate == "transcript":
            statement = 'SELECT %(column)s, COUNT(*) FROM %(track)s_effects_cds GROUP BY %(column)s'
        elif self.aggregate == "position":
            statement = 'SELECT %(column)s, COUNT(DISTINCT snp_position) FROM %(track)s_effects_cds GROUP BY %(column)s'

        data = odict(self.get(statement %
                              dict(self.members(), **locals())))

        return self.process(data)


class VariantsCDSEffectCodes(TrackerVariants):
    column = "code"


class VariantsCDSEffectCodesPerStrain(VariantsCDSEffectCodes):
    aggregate = "position"


class VariantsCDSEffectCodesPerPosition(VariantsTracker):

    '''count effects per SNP.'''
    mPattern = "^effects_cds$"
    mAsTable = True

    def __call__(self, track, slice=None):

        statement = """SELECT 
                  COUNT(*) AS 'all',
                  SUM(CASE WHEN X > 0 AND N = 0 AND S = 0 THEN 1 ELSE 0 END) AS 'X', 
                  SUM(CASE WHEN N > 0 AND X = 0 AND S = 0 THEN 1 ELSE 0 END) AS 'N',
                  SUM(CASE WHEN S > 0 AND X = 0 AND N = 0 THEN 1 ELSE 0 END) AS 'S' 
                  FROM %(track)s"""

        result = self.getFirstRow(statement % locals())

        r = odict(zip(("all", "X", "N", "S"), result))
        r["ambiguous"] = result[0] - sum(result[1:])
        del r["all"]
        return r


class VariantsCDSVariantCodes(TrackerVariants):

    '''type of cds variants according variant code.'''
    column = "variant_code"

    def process(self, data):
        '''quote `+` and `-`.'''
        result = odict()
        for key, counts in data.iteritems():
            result["``%s``" % key] = counts
        return result


class VariantsCDSVariantCodesPerPosition(VariantsCDSVariantCodes):
    aggregate = "position"


class VariantsCDSVariantTypes(TrackerVariants):

    '''type of cds variants according variant code.

    The counts are symmetrized, i.e., W,D = D,W.

    Note that the counts are not normalized by transcripts.
    If there is a variant affecting multiple transcripts,
    it will be counted multiple times.
    '''
    column = "variant_type"

    def process(self, data):
        '''symmetrize counts.'''
        result = collections.defaultdict(int)
        for key, counts in data.iteritems():
            try:
                if "," in key:
                    key = ",".join(sorted(key.split(",")))
            except TypeError:
                pass
            result[key] += counts
        return result
