import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import ChipseqReport

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class AnnotationSlicer:
    """Returns the default slices.

    The option *subset* returns a group of slices. The available
    groups are
       * default: only default slices ("all", "known", "unknown")
       * complete: all slices 
       * derived: only derived slices %s, see conf.py
       * "track": only the derived slices for a "track", see conf.py
       * mixture: combination of default slices and derived slices
       * set.track: combine one of the default slices with all slices of track "track"
       """
    def getSlices( self, subset = "default" ):

        if type(subset) in ( types.ListType, types.TupleType ):
            if len(subset) > 1:
                return subset
            else:
                subset = subset[0]
        
        all_slices = trackers_default_slices + trackers_derived_slices.keys()
        if subset == None or subset == "default":
            return trackers_default_slices
        elif subset == "complete":
            return all_slices
        elif subset == "derived":
            return trackers_derived_slices.keys()
        elif subset == "mixture":
            slices = []
            for x in trackers_default_slices:
                for y in trackers_derived_slices.keys():
                    slices.append("%s.%s" % (x,y) )
            return slices
        elif "." in subset:
            slices = []
            set1,track = subset.split(".")
            for y in [x[0] for x in trackers_derived_slices.items() if x[1] == track]:
                slices.append("%s.%s" % (set1,y) )
            return slices
        elif subset in all_slices:
            # a single slice
            return [subset,]
        else:
            # slices by track
            return [x[0] for x in trackers_derived_slices.items() if x[1] == subset ]

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class AnnotationsAssociated(ChipseqReport.DefaultTracker, AnnotationSlicer):
    """simple join between a data table and table defining slices.

    :attr:`mTable`
       table to join with
    :attr:`mColums`
       columns to output
    """
    mPattern = "_annotations$"
    mTable = None
    mColumns = None
    mWhere = "1"
    mSelectAll = "SELECT %(columns)s FROM %(track)s_%(table)s AS t WHERE %(where)s"
    mSelectSubset = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_annotation AS a WHERE a.gene_id = t.gene_id AND a.is_%(slice)s AND %(where)s"
    mSelectSlice = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_%(slice)s AS s WHERE s.gene_id = t.gene_id AND %(where)s"
    mSelectMixture = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(subset)s AS s, %(track)s_annotation AS a WHERE a.gene_id = t.gene_id AND a.is_%(slice)s AND s.gene_id = t.gene_id AND %(where)s"

    def getStatement( self, track, slice = None ):
        columns = self.mColumns
        table = self.mTable
        where = self.mWhere
        if not table or not columns: raise NotImplementedError
        if slice and "." in slice:
            slice, subset = slice.split(".")
            if subset in trackers_derived_slices and track != trackers_derived_slices[subset]:
                return None
            else:
                return self.mSelectMixture % locals()
        elif slice in trackers_derived_slices:
            if track == trackers_derived_slices[slice]:
                return self.mSelectSlice % locals()
            else:
                return None
        elif slice == "all" or slice == None:
            return self.mSelectAll % locals()
        else:
            return self.mSelectSubset % locals()


##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class Annotations(ChipseqReport.DefaultTracker, AnnotationSlicer):
    """Base class for trackers getting info from the annotations tables.

    Derived Trackers should define the two attributes :attr:`mSelect` and 
    :attr:`mColumns`.
    """
    pattern = "(.*)_annotations$"
    mTable = "annotations"
    mSelect = None
    mColumns = None
    mWhere = "1"

    def __call__(self, track, slice = None ):

        where = self.mWhere
        select = self.mSelect
        table = self.mTable

        if slice in trackers_derived_slices:
            if track == trackers_derived_slices[slice]:
                data = self.getFirstRow( """%(select)s
                        FROM %(track)s_%(table)s AS a, %(track)s_%(slice)s AS s WHERE %(where)s AND a.gene_id = s.gene_id""" % locals() )
            else:
                return []
        elif slice == "all" or slice == None:
            data = self.getFirstRow( "%(select)s FROM %(track)s_%(table)s WHERE %(where)s" % locals() )
        else:
            data = self.getFirstRow( "%(select)s FROM %(track)s_%(table)s WHERE %(where)s AND is_%slices" % locals() )
      
        return odict( zip(self.mColumns, data) )

class AllAnnotations(Annotations):
    """Annotations of all transcript models."""

    mColumns = [ "cds", 
                 "utr", 
                 "upstream", 
                 "downstream", 
                 "intronic", 
                 "intergenic", 
                 "flank",
                 "ambiguous", 
                 "unclassified" ]

    mSelect = """SELECT 
			sum(is_cds) AS cds, 
                        sum(is_utr) AS utr, 
                        sum(is_upstream) AS upstream, 
			sum(is_downstream) AS downstream,
                        sum(is_intronic) AS intronic, 
                        sum(is_intergenic) AS intergenic, 
                        sum(is_flank) AS flank, 
                        sum(is_ambiguous) AS ambiguous, 
			count(*) - (sum(is_cds+is_utr+is_upstream+is_downstream+is_intronic+is_intergenic+is_flank+is_ambiguous)) as unclassified"""


class AnnotationsBases(Annotations):
    """Annotations as bases."""
    mColumns = [ "total", "CDS", "UTRPromotor", "intronic", "intergenic" ]
    mSelect = """SELECT 
                 sum( exons_sum) AS total,
		 sum( nover_CDS ) AS cds,
                 sum( nover_UTR + nover_UTR3 + nover_UTR5 + nover_flank + nover_5flank + nover_3flank) AS utr, 
                 sum( nover_intronic) AS intronic,
                 sum( nover_intergenic) AS intergenic
                 """
