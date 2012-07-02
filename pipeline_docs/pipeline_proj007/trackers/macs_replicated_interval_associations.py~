import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
from cpgReport import *

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from macs_annotations import *

##################################################################################
class replicatedAssociationsHierarchy(cpgTracker):
    """ """

    pattern = "(.*)_replicated_intervals$"
    
    def __call__(self, track, slice = None ):

        ANNOTATIONS_NAME = P['annotations_name']
        statement =  """SELECT count(distinct interval_id) as NMIs, feature_class
                        FROM (
                        SELECT interval_id,
                        CASE WHEN coding_tss > 0  THEN 'Protein-coding gene TSS'
                        WHEN lincrna_tss > 0 THEN 'lncRNA gene TSS'
                        WHEN short_rna_tss > 0 THEN 'short RNA TSS'
                        WHEN pseudogene_tss > 0 THEN 'Pseudogene TSS'
                        WHEN processed_transcript_tss > 0 THEN 'Processed transcript TSS'
                        WHEN enhancer >0 THEN 'Enhancer (H3K4Me1)'
                        WHEN rnaseq >0 THEN 'Novel RNAseq transcript TSS'
                        ELSE 'Intergenic'
                        END AS feature_class FROM (
                        SELECT i.interval_id, a.coding_tss, b.lincrna_tss, c.short_rna_tss, d.pseudogene_tss, e.processed_transcript_tss, f.enhancer, g.rnaseq FROM %(track)s_replicated_intervals i 
                        left join 
                        (SELECT distinct gene_id, 1 as coding_tss 
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap
                        where (genes_nover>0 OR downstream_flank_nover>0 OR upstream_flank_nover>0) ) a 
                        on a.gene_id=i.interval_id
                        left join 
                        (SELECT distinct gene_id, 1 as lincrna_tss 
                        FROM %(track)s_replicated_lncrna_tss_distance
                        where closest_dist < 1000) b 
                        on i.interval_id=b.gene_id
                        left join 
                        (SELECT distinct n.gene_id as gene_id, 1 as short_rna_tss
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n,
                        annotations.transcript_info t, %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m
                        where n.gene_id=m.interval_id
                        AND m.gene_id=t.gene_id
                        AND t.gene_biotype IN ("miRNA","snRNA","snoRNA","rRNA")
                        AND n.closest_dist < 1000) c
                        on i.interval_id=c.gene_id
                        left join 
                        (SELECT distinct n.gene_id as gene_id, 1 as pseudogene_tss
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n,
                        annotations.transcript_info t, %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m
                        where n.gene_id=m.interval_id
                        AND m.gene_id=t.gene_id
                        AND t.gene_biotype="pseudogene"
                        AND n.closest_dist < 1000) d
                        on i.interval_id=d.gene_id
                        left join 
                        (SELECT distinct n.gene_id as gene_id, 1 as processed_transcript_tss
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n,
                        annotations.transcript_info t, %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m
                        where n.gene_id=m.interval_id
                        AND m.gene_id=t.gene_id
                        AND t.gene_biotype="processed_transcript"
                        AND n.closest_dist < 1000) e
                        on i.interval_id=e.gene_id
                        left join 
                        (SELECT distinct interval_id, 1 as enhancer 
                        FROM %(track)s_replicated_h3k4me1_intervals) f
                        on i.interval_id=f.interval_id
                        left join 
                        (SELECT distinct gene_id, 1 as rnaseq 
                        FROM %(track)s_replicated_rnaseq_tss_distance
                        WHERE closest_dist < 1000) g
                        on i.interval_id=g.gene_id))
                        group by feature_class
                        order by feature_class asc;""" % locals()
        data = self.getAll(statement)
        return data

##################################################################################
class replicatedAssociationsHierarchy3(cpgTracker):
    """ """

    pattern = "(.*)_replicated_intervals$"
    
    def __call__(self, track, slice = None ):

        ANNOTATIONS_NAME = P['annotations_name']
        statement =  """SELECT count(distinct interval_id) as NMIs, feature_class
                        FROM (
                        SELECT interval_id,
                        CASE WHEN coding_tss > 0  THEN 'Protein-coding gene TSS'
                        WHEN lincrna_tss > 0 THEN 'lncRNA gene TSS'
                        WHEN short_rna_tss > 0 THEN 'short RNA TSS'
                        WHEN enhancer >0 THEN 'Enhancer (H3K4Me1)'
                        WHEN rnaseq >0 THEN 'Novel RNAseq transcript TSS'
                        ELSE 'Intergenic'
                        END AS feature_class FROM (
                        SELECT i.interval_id, a.coding_tss, b.lincrna_tss, c.short_rna_tss, f.enhancer, g.rnaseq FROM %(track)s_replicated_intervals i 
                        left join 
                        (SELECT distinct gene_id, 1 as coding_tss 
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap
                        where (genes_nover>0 OR downstream_flank_nover>0 OR upstream_flank_nover>0) ) a 
                        on a.gene_id=i.interval_id
                        left join 
                        (SELECT distinct gene_id, 1 as lincrna_tss 
                        FROM %(track)s_replicated_lncrna_tss_distance
                        where closest_dist < 1000) b 
                        on i.interval_id=b.gene_id
                        left join 
                        (SELECT distinct n.gene_id as gene_id, 1 as short_rna_tss
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance n,
                        annotations.transcript_info t, %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_noncoding_mapping m
                        where n.gene_id=m.interval_id
                        AND m.gene_id=t.gene_id
                        AND t.gene_biotype IN ("miRNA","snRNA","snoRNA","rRNA")
                        AND n.closest_dist < 1000) c
                        on i.interval_id=c.gene_id
                        left join 
                        (SELECT distinct interval_id, 1 as enhancer 
                        FROM %(track)s_replicated_h3k4me1_intervals) f
                        on i.interval_id=f.interval_id
                        left join 
                        (SELECT distinct gene_id, 1 as rnaseq 
                        FROM %(track)s_replicated_rnaseq_tss_distance
                        WHERE closest_dist < 1000) g
                        on i.interval_id=g.gene_id))
                        group by feature_class
                        order by feature_class asc;""" % locals()
        data = self.getAll(statement)
        return data
        
##################################################################################
class replicatedAssociationsHierarchy2(cpgTracker):
    """ """

    pattern = "(.*)_replicated_intervals$"
    
    def __call__(self, track, slice = None ):

        ANNOTATIONS_NAME = P['annotations_name']
        statement =  """SELECT count(distinct interval_id) as NMIs, feature_class
                        FROM (
                        SELECT interval_id,
                        CASE WHEN coding_tss > 0  THEN 'Protein-coding gene TSS'
                        WHEN noncoding_tss > 0 THEN 'Non-coding gene TSS'
                        ELSE 'Intergenic'
                        END AS feature_class FROM (
                        SELECT i.interval_id, a.coding_tss, b.noncoding_tss 
                        FROM %(track)s_replicated_intervals i 
                        left join 
                        (SELECT distinct gene_id, 1 as coding_tss 
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap
                        where (genes_nover>0 OR downstream_flank_nover>0 OR upstream_flank_nover>0) ) a 
                        on a.gene_id=i.interval_id
                        left join 
                        (SELECT distinct gene_id, 1 as noncoding_tss 
                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance
                        where closest_dist < 1000) b 
                        on i.interval_id=b.gene_id))
                        group by feature_class
                        order by feature_class asc;""" % locals()
        data = self.getAll(statement)
        return data
        

##################################################################################
class replicatedAssociations(cpgTracker):
    """ """

    pattern = "(.*)_replicated_intervals$"
    
    def __call__(self, track, slice = None ):
        ANNOTATIONS_NAME = P['annotations_name']
        try: 
            data1 = self.getValue( """SELECT count(distinct gene_id) as intervals
                                        FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_overlap
                                        where (genes_nover>0 OR downstream_flank_nover>0 OR upstream_flank_nover>0)""" % locals() )
        except:
            data1 = "0"
        try: 
            data2 = self.getValue( """SELECT count(distinct gene_id) as intervals
                                       FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_noncoding_tss_distance
                                       where closest_dist < 1000""" % locals() )
        except:
            data2 = "0"
        try:
            data3 = self.getValue( """SELECT distinct count(distinct interval_id) as intervals, "enhancer" as feature_class 
                                       FROM %(track)s_replicated_h3k4me1_intervals""" % locals() )
        except:
            data3 = "0"
        try:
            data4 = self.getValue( """SELECT count(distinct gene_id) as intervals
                                        FROM %(track)s_replicated_rnaseq_tss_distance
                                        where closest_dist < 1000""" % locals() )
        except:
            data4 = "0"
        try: 
            data5 = self.getValue( """SELECT count(distinct gene_id) as intervals
                                        FROM %(track)s_replicated_lncrna_tss_distance
                                        where closest_dist < 1000""" % locals() )
        except:
            data5 = "0"
        return odict( zip(("Protein-coding TSS","Non-coding TSS","H3K4Me1 Enhancer", "RNAseq transcript", "lincRNA TSS"), (data1, data2, data3, data4, data5)) )
        
        
