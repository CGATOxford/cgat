import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
class cgiAnnotations(cpgTracker):
    """Breakdown of overlap of predicted CGIs with genomic regions """

    mPattern = "cgi_annotations"

    def __call__(self, track, slice = None ):

        data = self.getFirstRow( """SELECT 
                                sum(is_cds) AS cds, 
                                sum(is_utr) AS utr, 
                                sum(is_upstream) AS upstream, 
                                sum(is_downstream) AS downstream,
                                sum(is_intronic) AS intronic, 
                                sum(is_intergenic) AS intergenic, 
                                sum(is_flank) AS flank, 
                                sum(is_ambiguous) AS ambiguous 
                                FROM cgi_annotations""" )
        mColumns = [ "cds", 
                 "utr", 
                 "upstream", 
                 "downstream", 
                 "intronic", 
                 "intergenic", 
                 "flank",
                 "ambiguous" ]

        return odict( zip(mColumns, data) )

##################################################################################
class cgiGenomicFeatures(cpgTracker):
    """return overlap of interval with genomic features """

    mPattern = "cgi_genomic_features$"
    
    def __call__(self, track, slice = None ):
        data = self.getValues( """SELECT count(distinct gene_id) as intervals FROM (
                               SELECT gene_id,
                               CASE WHEN  tss_extended_pover1 > 0  THEN 'TSS'
                               WHEN genes_pover1 > 0 THEN 'Gene'
                               WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                               WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                               ELSE 'Intergenic'
                               END AS feature_class
                               FROM cgi_genomic_features)
                               group by feature_class
                               order by feature_class asc""" % locals() )
               
        result = odict(zip(("Downstream","Gene","Intergenic","TSS","Upstream"),data))
        return result

##################################################################################
class cgitssoverlap(cpgTracker):
    """overlap of predicted CGIs with TSS """

    mPattern = "tss_cgi_venn"

    def __call__(self, track, slice = None ):
        data = self.getAll("SELECT track, intervals from tss_cgi_venn")
        return data

##################################################################################
class CGI_CpGObsExp2( cpgTracker ):
    mPattern = "_comp$"

    def __call__(self, track, slice = None):
        data = self.getAll( "SELECT CpG_ObsExp FROM %(track)s_comp" % locals() )
        return data

##################################################################################
class CGI_GCContent( cpgTracker ):
    mPattern = "_comp$"

    def __call__(self, track, slice = None):
        data = self.getAll( "SELECT pGC FROM %(track)s_comp" % locals() )
        return data



