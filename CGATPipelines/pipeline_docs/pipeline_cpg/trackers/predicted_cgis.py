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
        data = self.getAll( "SELECT CpG_ObsExp2 FROM %(track)s_comp" % locals() )
        return data

##################################################################################
class CGI_GCContent( cpgTracker ):
    mPattern = "_comp$"

    def __call__(self, track, slice = None):
        data = self.getAll( "SELECT pGC FROM %(track)s_comp" % locals() )
        return data



