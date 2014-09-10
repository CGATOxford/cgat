import os
import sys
import re
import types
import itertools
import PeakcallingReport

from CGATReport.Tracker import *
from collections import OrderedDict as odict

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class AnnotationSlicer:

    """Returns the default slices.
    """

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class AnnotationsAssociated(PeakcallingReport.DefaultTracker, AnnotationSlicer):

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

    def getStatement(self, track, slice=None):
        columns = self.mColumns
        table = self.mTable
        where = self.mWhere
        if not table or not columns:
            raise NotImplementedError
        if slice and "." in slice:
            slice, subset = slice.split(".")
            return self.mSelectMixture % locals()
        elif slice == "all" or slice is None:
            return self.mSelectAll % locals()
        else:
            return self.mSelectSubset % locals()


##########################################################################
##########################################################################
##########################################################################
##
##########################################################################
class Annotations(PeakcallingReport.DefaultTracker, AnnotationSlicer):

    """Base class for trackers getting info from the annotations tables.

    Derived Trackers should define the two attributes :attr:`mSelect` and 
    :attr:`mColumns`.
    """
    pattern = "(.*)_annotations$"
    table = "annotations"
    select = None
    columns = None
    where = "1"

    def __call__(self, track, slice=None):

        where = self.where
        select = self.select
        table = self.table

        if slice == "all" or slice is None:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_%(table)s WHERE %(where)s" % locals())
        else:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_%(table)s WHERE %(where)s AND is_%slices" % locals())

        return odict(zip(self.columns, data))


class AllAnnotations(Annotations):

    """Annotations of all transcript models."""

    columns = ["cds",
               "utr",
               "upstream",
               "downstream",
               "intronic",
               "intergenic",
               "flank",
               "ambiguous",
               "unclassified"]

    select = """SELECT 
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
    columns = ["total", "CDS", "UTRPromotor", "intronic", "intergenic"]
    select = """SELECT 
                 sum( exons_sum) AS total,
		 sum( nover_CDS ) AS cds,
                 sum( nover_UTR + nover_UTR3 + nover_UTR5 + nover_flank + nover_5flank + nover_3flank) AS utr, 
                 sum( nover_intronic) AS intronic,
                 sum( nover_intergenic) AS intergenic
                 """
