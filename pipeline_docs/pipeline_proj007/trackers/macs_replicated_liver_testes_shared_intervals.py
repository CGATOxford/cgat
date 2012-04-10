import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class replicatedSharedIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(end-start),0) as length FROM liver_testes_shared_intervals" % locals() )
        return odict( zip( ("Shared intervals", "mean_interval_length" ), data) )

##################################################################################
class replicatedsharedIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT (end-start) FROM liver_testes_shared_intervals" % locals() )
        return { "length" : data }

##################################################################################
class replicatedSharedIntervalTSS( cpgTracker ):
    """Distribution of distance to closest TSS """

    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( '''SELECT closest_dist FROM liver_testes_shared_intervals u, 
                                  liver_testes_merged_intervals i, liver_testes_merged_%(ANNOTATIONS_NAME)s_transcript_tss_distance t
                                  WHERE u.interval_id=i.interval_id 
                                  AND t.gene_id=i.interval_id''' % locals() )
        return { "distance" : data }

##################################################################################
class replicatedSharedIntervalCpGDensity( cpgTracker ):
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pCpG FROM liver_testes_shared_intervals u, 
                               liver_testes_merged_intervals i, liver_testes_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedSharedIntervalCpGObsExp( cpgTracker ):
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp FROM liver_testes_shared_intervals u, 
                               liver_testes_merged_intervals i,liver_testes_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedSharedIntervalCpGNumber( cpgTracker ):
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT nCpG FROM liver_testes_shared_intervals u, 
                               liver_testes_merged_intervals i, liver_testes_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedSharedIntervalGCContent( cpgTracker ):
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pGC FROM liver_testes_shared_intervals u, 
                               liver_testes_merged_intervals i, liver_testes_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class replicatedSharedIntervalTranscriptOverlap(featureOverlap):
    """return overlap of interval with  protein-coding transcripts """
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN  tss_transcript_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM liver_testes_merged_%(ANNOTATIONS_NAME)s_overlap o, liver_testes_shared_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )
               
        return odict(zip(("Downstream","Gene","Intergenic","TSS","Upstream"),data))

##################################################################################
class replicatedSharedIntervalGeneOverlap(featureOverlap):
    """return overlap of interval with  protein-coding genes """
    mPattern = "liver_testes_shared_intervals$"

    def __call__(self, track, slice = None):
        ANNOTATIONS_NAME = P['annotations_name']
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN tss_gene_extended_pover1 > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM liver_testes_merged_%(ANNOTATIONS_NAME)s_overlap o, liver_testes_shared_intervals u
                                   WHERE u.interval_id=o.gene_id)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )
               
        return odict(zip(("Downstream","Gene","Intergenic","TSS","Upstream"),data))
        
