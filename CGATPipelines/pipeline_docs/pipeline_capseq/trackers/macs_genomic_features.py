import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import cpgReport

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##########################################################################


class mergedIntervalEnsemblTranscriptOverlap(cpgReport.cpgTracker):

    """return overlap of interval with genomic features """

    mPattern = "_merged_ensembl_transcript_overlap$"

    def __call__(self, track, slice=None):
        data = self.getValues( """SELECT count(distinct gene_id) as intervals FROM (
                               SELECT gene_id,
                               CASE WHEN  tss_transcript_extended_pover1 > 0  THEN 'TSS'
                               WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                               WHEN genes_pover1 > 0 THEN 'Gene'
                               WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                               ELSE 'Intergenic'
                               END AS feature_class
                               FROM %(track)s_merged_ensembl_transcript_overlap)
                               group by feature_class
                               order by feature_class asc""" % locals() )

        result = odict(
            zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data))
        return result

##########################################################################


class mergedIntervalEnsemblGeneOverlap(cpgReport.cpgTracker):

    """return overlap of interval with genomic features """

    mPattern = "_merged_ensembl_gene_overlap$"

    def __call__(self, track, slice=None):
        data = self.getValues( """SELECT count(distinct gene_id) as intervals FROM (
                               SELECT gene_id,
                               CASE WHEN  tss_gene_extended_pover1 > 0  THEN 'TSS'
                               WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                               WHEN genes_pover1 > 0 THEN 'Gene'
                               WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                               ELSE 'Intergenic'
                               END AS feature_class
                               FROM %(track)s_merged_ensembl_gene_overlap)
                               group by feature_class
                               order by feature_class asc""" % locals() )

        result = odict(
            zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data))
        return result

##########################################################################


class RepeatOverlap(cpgReport.cpgTracker):

    """Overlap with repeats."""
    mPattern = "_merged_repeats$"

    def __call__(self, track, slice=None):
        statement = '''SELECT SUM(CASE WHEN nover>0 THEN 1 ELSE 0 END) as with, SUM(CASE WHEN nover=0 THEN 1 ELSE 0 END) AS without
                       FROM %(track)s_merged_repeats '''
        return odict(zip(("with", "without"), self.getFirstRow(statement)))
