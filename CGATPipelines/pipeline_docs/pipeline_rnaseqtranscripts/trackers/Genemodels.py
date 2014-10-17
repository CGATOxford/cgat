import os
import sys
import re
import types
import itertools

from CGATReport.Tracker import *
from collections import OrderedDict as odict

from RnaseqTranscriptsReport import *


class TrackerGenemodels(RnaseqTranscriptsTracker):
    mPattern = "_gene_expression"

    def getTracks(self, subset=None):
        return self.getValues("SELECT DISTINCT track FROM %s" % self.name)


class GeneModelsBenchmark(TrackerGenemodels):
    name = "agg_agg_agg_cuffcompare_benchmark"

    def getSlices(self, subset=None):
        '''slice by contig'''
        return ('all',)  # self.getValues( "SELECT DISTINCT contig FROM transcripts_compare" )

    def __call__(self, track, slice=None):
        return self.getRow( """
               SELECT baselevel_sp, baselevel_sn,
                      exonlevel_sp, exonlevel_sn,
                      transcriptlevel_sp, transcriptlevel_sn,
                      intronlevel_sp, intronlevel_sn,
                      locuslevel_sp, locuslevel_sn,
                      100.0 * missedexons_counts / missedexons_total AS missed_exons,
                      100.0 * missedloci_counts / missedloci_total AS missed_loci,
                      100.0 * novelexons_counts / novelexons_total AS wrong_exons,
                      100.0 * novelloci_counts / novelloci_total AS wrong_loci
                      FROM %(name)s WHERE track = '%(track)s' AND contig = '%(slice)s'
               """ % self.members(locals()))


class GeneModelsCodes(RnaseqTranscriptsTracker):
    pattern = "(.*)_cuffcompare_tracking"
    as_tables = True

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice=None):
        return self.getValue( """SELECT COUNT(*) AS ntransfrags FROM %(track)s WHERE code = '%(slice)s'""" )


class GeneModelsSharedLoci(RnaseqTranscriptsTracker):

    '''number of times a locus appears in experiments.'''
    mPattern = "_cuffcompare_loci"
    mAsTables = True

    def __call__(self, track, slice=None):
        return self.getDict("SELECT nexperiments, count(*) AS nloci FROM %(track)s group by nexperiments")


class GeneModelsSharedTransfrags(RnaseqTranscriptsTracker):

    '''number of times a transfrag appears in experiments.'''
    mPattern = "_cuffcompare_tracking"
    mAsTables = True

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice=None):
        data = self.getDict(
            "SELECT nexperiments, count(*) AS ntransfrags FROM %(track)s WHERE code = '%(slice)s' group by nexperiments")
        return data


class ExpressionByClass(RnaseqTranscriptsTracker):

    '''number of times a transfrag appears in experiments.'''
    mPattern = "_cuffcompare_tracking"
    mAsTables = False

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice=None):
        vals = self.getValues( """SELECT avg(FPKM)
                                      FROM %(track)s_cuffcompare_tracking AS t,
                                           %(track)s_cuffcompare_transcripts AS a
                                      WHERE code = '%(slice)s' AND 
                                      a.transfrag_id = t.transfrag_id
                                      GROUP BY a.transfrag_id""" % locals() )
        return odict((("fpkm", vals), ))


class TransfragCorrelation(RnaseqTranscriptsTracker):

    '''return correlation table 
    '''
    mPattern = "_reproducibility"

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice=None):
        data = self.getAll( """SELECT track1, track2, 
                                      coeff, pvalue, significance,
                                      pairs, both_null, null1, null2, 
                                      method, alternative
                               FROM %(track)s_reproducibility 
                               WHERE code = '%(slice)s'""" )
        return data


class TransfragReproducibility2(RnaseqTranscriptsTracker):

    '''return proportion of transfrags present in a pair of replicates.
    '''

    pattern = "(.*)_reproducibility"

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice):
        data = self.getAll( """SELECT track1, track2, 
                                      ROUND( CAST( not_null AS FLOAT) / (pairs-both_null),2) AS pcalled,
                                      ROUND( coeff, 2) as correlation
                               FROM %(track)s_reproducibility 
                               WHERE code = '%(slice)s'""" )
        return data


class TransfragReproducibility(RnaseqTranscriptsTracker):

    '''return proportion of transfrags present in a pair of replicates.
    '''

    mPattern = "_reproducibility"

    def getSlices(self, subset=None):
        return tuple("=cjeiopruxs.*")

    def __call__(self, track, slice=None):
        data = self.getDict( """SELECT track1 || '_x_' || track2 AS tracks, 
                                      ROUND( CAST( not_null AS FLOAT) / (pairs-both_null),2) AS pcalled,
                                      ROUND( coeff, 2) as correlation
                               FROM %(track)s_reproducibility 
                               WHERE code = '%(slice)s'""" )
        return data


class GenesetSummary(RnaseqTranscriptsTracker, SingleTableTrackerRows):

    '''summary properties of genesets.'''
    table = "geneset_stats"
    column = "track"


class GenesetMappability(RnaseqTranscriptsTracker):

    '''return average mappability for all transcripts.'''
    mPattern = "_mappability"

    def __call__(self, track, slice=None):
        return odict( (('mean', self.getValues( '''SELECT mean FROM %(track)s_mappability''' )), ))


class TranscriptClassCounts(RnaseqTranscriptsTracker):

    '''return number of transcripts within each class.'''
    pattern = "(.*)_class"

    def getSlices(self, subset=None):
        '''slice by contig'''
        return self.getValues("SELECT DISTINCT source FROM agg_agg_agg_class") + ["-"]

    def __call__(self, track, slice=None):
        if slice == "-":
            stmt = "is NULL"
        else:
            stmt = "= '%(slice)s'" % locals()

        return self.getDict( '''SELECT class || '-' || CASE WHEN sense = 's' THEN 'sense' WHEN sense = 'a' THEN 'antisense' ELSE 'anysense' END, 
                                       COUNT(*) AS ntranscripts
                                       FROM %(track)s_class 
                                       WHERE source %(stmt)s
                                       GROUP BY class, sense''' )


class TranscriptClassCountsSummaryBySource(RnaseqTranscriptsTracker):

    '''return number of transcripts within each class.'''
    pattern = "(.*)_class$"

    def __call__(self, track, slice=None):
        return self.getDict( '''SELECT source, COUNT(*) AS ntranscripts FROM %(track)s_class GROUP BY source''')


class TranscriptClassCountsSummaryByClass(RnaseqTranscriptsTracker):

    '''return number of transcripts within each class.'''
    pattern = "(.*)_class$"

    def __call__(self, track, slice=None):
        return self.getDict( '''SELECT class, COUNT(*) AS ntranscripts FROM %(track)s_class GROUP BY class''')


class BuildSummary(RnaseqTranscriptsTracker):

    '''summary of gene set construction.'''
    pattern = "(.*)_build_summary"

    def __call__(self, track, slice=None):
        return self.getAll( '''SELECT category, transcripts FROM %(track)s_build_summary''')
