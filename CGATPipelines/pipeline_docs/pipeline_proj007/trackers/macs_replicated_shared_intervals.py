import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import scipy.stats
import numpy.ma
import Stats
import Histogram

from CGATReport.Tracker import *
from cpgReport import *

##########################################################################


class replicatedSharedIntervals(cpgTracker):

    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_replicated_shared_intervals" % locals())
        return odict(zip(("Shared intervals", "mean_interval_length"), data))

##########################################################################


class replicatedsharedIntervalLengths(cpgTracker):

    """Distribution of interval length. """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT (stop-start) FROM %(track)s_replicated_shared_intervals" % locals())
        return {"length": data}

##########################################################################


class replicatedSharedIntervalPeakValues(cpgTracker):

    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"peakval": data}

##########################################################################


class replicatedSharedIntervalAverageValues(cpgTracker):

    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return {"avgval": data}

##########################################################################


class replicatedSharedIntervalFoldChange(cpgTracker):

    """Distribution of fold change """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues( '''SELECT i.fold FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return odict([("Fold Change", data)])

##########################################################################


class replicatedSharedIntervalTSS(cpgTracker):

    """Distribution of distance to closest TSS """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( '''SELECT closest_dist FROM %(track)s_replicated_shared_intervals u, 
                                  %(track)s_replicated_intervals i, %(track)s_replicated_%(ANNOTATIONS_NAME)s_transcript_tss_distance t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND t.gene_id=i.interval_id''' % locals() )
        return {"distance": data}

##########################################################################


class replicatedSharedIntervalCpGDensity(cpgTracker):
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT pCpG FROM %(track)s_replicated_shared_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedSharedIntervalCpGObsExp(cpgTracker):
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT CpG_ObsExp FROM %(track)s_replicated_shared_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedSharedIntervalCpGNumber(cpgTracker):
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_replicated_shared_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################


class replicatedSharedIntervalGCContent(cpgTracker):
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_replicated_shared_intervals u, 
                               %(track)s_replicated_intervals i,%(track)s_replicated_capseq_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##########################################################################
##########################################################################
##########################################################################


class replicatedSharedIntervalLengthVsAverageValue(cpgTracker):

    """Length vs average value. """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT length, avgval FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("length", "avgval"), zip(*data)))

##########################################################################


class replicatedSharedIntervalLengthVsPeakValue(cpgTracker):

    """Length vs peak value """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT length, peakval FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("length", "peakval"), zip(*data)))

##########################################################################


class replicatedSharedIntervalLengthVsFoldChange(cpgTracker):

    """Length vs fold change"""

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT length, fold FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("length", "foldchange"), zip(*data)))

##########################################################################


class replicatedSharedIntervalAvgValVsPeakVal(cpgTracker):

    """average value vs peak value """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT avgval, peakval FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("avgval", "peakval"), zip(*data)))

##########################################################################


class replicatedSharedIntervalAvgValVsFoldChange(cpgTracker):

    """average value vs fold change """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT avgval, fold FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("avgval", "foldchange"), zip(*data)))

##########################################################################


class replicatedSharedIntervalPeakValVsFoldChange(cpgTracker):

    """Peak value vs fold change """

    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        data = self.get( '''SELECT peakval, fold FROM %(track)s_replicated_shared_intervals u, %(track)s_replicated_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict(zip(("peakval", "foldchange"), zip(*data)))

##########################################################################


class replicatedSharedIntervalTranscriptOverlap(featureOverlap):

    """return overlap of interval with  protein-coding transcripts """
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN  tss_transcript_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap o, %(track)s_replicated_shared_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )

        return odict(zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data))

##########################################################################


class replicatedSharedIntervalGeneOverlap(featureOverlap):

    """return overlap of interval with  protein-coding genes """
    mPattern = "_replicated_shared_intervals$"

    def __call__(self, track, slice=None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN tss_gene_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap o, %(track)s_replicated_shared_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )

        return odict(zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data))
