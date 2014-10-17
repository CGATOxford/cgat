"""specialized trackers for the songbird project.
"""

import math
import sys
import os
import collections

from CGATReport.Tracker import *
import CGATReport.Stats


class IntronicEnrichment(TrackerSQL):

    """test for intronic enrichment of transcripts using a chi-squared test.

    Collect the proportion of intronic transcripts among unknown transcriptios
    in an expression set and the set unchanged.
    """

    def getSlices(self, subset=None):
        return []

    def getTracks(self):
        return ["FastDown", "SlowDown", "FastUp", "SlowUp"]

    def __call__(self, track, slice=None):

        master = trackers_derived_slices[track]
        ref = "unchanged"

        statement_in = """SELECT COUNT(*) FROM %(master)s_annotation AS m, %(track)s AS t WHERE t.gene_id = m.gene_id AND m.is_unknown AND m.is_intronic"""
        statement_out = """SELECT COUNT(*) FROM %(master)s_annotation AS m, %(track)s AS t WHERE t.gene_id = m.gene_id AND m.is_unknown AND NOT m.is_intronic"""

        matrix = numpy.matrix(
            ((self.getValue(statement_in % {"master": master, "track": track}),
                self.getValue(statement_out % {"master": master, "track": track})),
             (self.getValue(statement_in % {"master": master, "track": ref}),
              self.getValue(statement_out % {"master": master, "track": ref}))))

        result = CGATReport.Stats.doChiSquaredTest(matrix)

        return odict((("P", "%e" % result.mProbability),
                      ("df", "%i" % result.mDegreesFreedom),
                      ("chi2", "%5.2f" % result.mChiSquaredValue),
                      ("passed", result.mPassed),
                      ("phi", "%5.2f" % result.mPhi)))
