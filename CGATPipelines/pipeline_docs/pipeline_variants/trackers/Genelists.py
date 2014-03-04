import os
import sys
import re
import types
import collections
import math
from VariantsReport import *

#####################################################
#####################################################
#####################################################


class TrackerGenelists(VariantsTracker):

    '''default tracker for allele analysis.'''
    pattern = "(.*)_gla$"
    qvalue_threshold = 0.05


class GenelistPValues(TrackerGenelists):

    def __call__(self, track, slice=None):

        statement = '''SELECT pvalue FROM %(track)s_gla''' % locals()
        return odict((('pvalues', self.getValues(statement)), ))


class GenelistSummary(TrackerGenelists):

    def getSlices(self, subset):
        return StrainTracker().getTracks(subset)

    def __call__(self, track, slice=None):

        statement = '''SELECT genelist, COUNT(*) 
                              FROM %(track)s_gla 
                              WHERE fdr < %(qvalue_threshold)f
                              AND genelist LIKE '%(slice)s%%'
                              GROUP BY genelist ''' % self.members(locals())
        data = self.get(statement)
        return odict(data)


class GenelistSignificantCategories(TrackerGenelists):

    # def getSlices( self, subset ):
    #     return StrainTracker().getTracks( subset )

    def __call__(self, track, slice=None):

        statement = '''SELECT genelist, go_description, ratio, pvalue
                              FROM %(track)s_gla 
                              WHERE fdr < %(qvalue_threshold)f
                              ''' % self.members(locals())
        data = self.get(statement)
        r = collections.defaultdict(dict)
        for genelist, desc, ratio, pvalue in data:
            if ratio > 0:
                ratio = math.log10(ratio)
            pvalue = -math.log10(pvalue)

            r[genelist][desc] = {'pvalue': pvalue,
                                 'fold': ratio}

        return r
