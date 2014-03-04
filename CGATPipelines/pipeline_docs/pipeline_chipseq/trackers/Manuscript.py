import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma

from ChipseqReport import *


class ReproducibilityBetweenSamples(DefaultTracker):

    def __call__(self, track, slice=None):

        set1 = track + "R1"
        set2 = track + "R2"

        statement = '''
        SELECT MAX(pexons_ovl1, pexons_ovl2) FROM overlap
        WHERE set1 IN ('%(set1)s', '%(set2)s') AND set2 IN ('%(set1)s', '%(set2)s') ''' % locals()

        return odict((("overlap", self.getValue(statement)),))


class AnnotatorTracksComparisonBetweenMotifs(DefaultTracker):

    '''get fold enrichment for AnnotatorTracks comparing
    intervals with and without motifs.'''

    prefix = "annotator_tracks"

    def getSlices(self, subset=None):
        return self.getValues("SELECT DISTINCT category FROM %s" % self.prefix)

    def __call__(self, track, slice=None):

        statement = '''
        SELECT fold, pvalue, qvalue FROM %(tablename)s WHERE track = '%(track)s' AND category = '%(slice)s' 
        '''

        tablename = self.prefix + "_with_motif"
        with_motif, with_pvalue, with_qvalue = self.getFirstRow(
            statement % locals())
        tablename = self.prefix + "_without_motif"
        without_motif, without_pvalue, without_qvalue = self.getFirstRow(
            statement % locals())

        try:
            ratio = "%5.2f" % (100.0 * (with_motif / without_motif))
        except ZeroDivisionError:
            ratio = "na"

        return odict((("with motif - fold", with_motif),
                      ("without motif", without_motif),
                      ("with/without motif [%]", ratio),
                      ("with motif - pvalue", with_pvalue),
                      ("with motif - qvalue", with_qvalue),
                      ("without motif - pvalue", without_pvalue),
                      ("without motif - qvalue", without_qvalue)))


class AnnotatorRegionsOfInterestComparisonBetweenMotifs(AnnotatorTracksComparisonBetweenMotifs):

    '''get fold enrichment for AnnotatorRegionsOfInterest comparing
    intervals with and without motifs.'''

    prefix = "annotator_roi"
