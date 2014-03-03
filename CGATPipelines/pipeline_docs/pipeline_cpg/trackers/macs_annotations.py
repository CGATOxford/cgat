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


class Annotations(cpgReport.cpgTracker):

    """Base class for trackers getting info from the annotations tables.
    Derived Trackers should define the two attributes :attr:`mSelect` and :attr:`mColumns`. """

    pattern = "(.*)_annotations$"
    mTable = "annotations"
    mSelect = None
    mColumns = None
    mWhere = "1"

    def __call__(self, track, slice=None):

        where = self.mWhere
        select = self.mSelect
        table = self.mTable

        if slice == "all" or slice is None:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_%(table)s WHERE %(where)s" % locals())
        else:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_%(table)s WHERE %(where)s AND is_%slices" % locals())

        return odict(zip(self.mColumns, data))

##########################################################################


class AllAnnotations(Annotations):

    """Annotations of all transcript models."""

    mColumns = ["cds",
                "utr",
                "upstream",
                "downstream",
                "intronic",
                "intergenic",
                "flank",
                "ambiguous"]

    mSelect = """SELECT 
                  sum(is_cds) AS cds, 
                  sum(is_utr) AS utr, 
                  sum(is_upstream) AS upstream, 
                  sum(is_downstream) AS downstream,
                  sum(is_intronic) AS intronic, 
                  sum(is_intergenic) AS intergenic, 
                  sum(is_flank) AS flank, 
                  sum(is_ambiguous) AS ambiguous"""

##########################################################################


class AnnotationsBases(Annotations):

    """Annotations as bases."""
    mColumns = ["total", "CDS", "UTRPromotor", "intronic", "intergenic"]
    mSelect = """SELECT 
                 sum( exons_sum) AS total,
		             sum( nover_CDS ) AS cds,
                 sum( nover_UTR + nover_UTR3 + nover_UTR5 + nover_flank + nover_5flank + nover_3flank) AS utr, 
                 sum( nover_intronic) AS intronic,
                 sum( nover_intergenic) AS intergenic """

##########################################################################


class AnnotationsAssociated(cpgReport.cpgTracker):

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


class RepeatOverlap(AnnotationsAssociated):

    """Overlap with repeats."""
    mPattern = "_repeats$"
    mColumns = "SUM(CASE WHEN nover>0 THEN 1 ELSE 0 END) as with, SUM(CASE WHEN nover=0 THEN 1 ELSE 0 END) AS without"
    mTable = "repeats"

    def __call__(self, track, slice=None):
        statement = self.getStatement(track, slice)
        if not statement:
            return []
        return odict(zip(("with", "without"), self.getFirstRow(statement)))

##########################################################################
##########################################################################
##########################################################################


class TSSOverlap(cpgReport.cpgTracker):

    '''number of TSS that an interval overlaps.'''
    mPattern = "_tss$"
    mAnnotations = "annotations"
    mTable = "tss"
    mColumn = "d.is_overlap"
    mWhere = "d.is_overlap < 5 "

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        hist, bins = numpy.histogram(
            data, bins=numpy.arange(0, max(data) + 1, 1))
        return odict(zip(map(str, bins[:-1]), hist))

##########################################################################


class TSSClosest(cpgReport.cpgTracker):

    """for each interval, return the distance to the closest TSS."""

    mXLabel = "distance / bases"
    mPattern = "_tss$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "annotations"
    mTable = "tss"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.get(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.get( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        return data

##########################################################################


class TSSClosestUpstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mColumn = "d.dist5"
    mWhere = "d.dist5 > 0"

##########################################################################


class TSSClosestDownstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mColumn = "d.dist3"
    mWhere = "d.dist3 > 0"

##########################################################################


class TSSProfile(cpgReport.cpgTracker):

    """Get profile around TSS"""
    mPattern = "_tss$"

    def __call__(self, track, slice=None):
        statement1 = """SELECT (closest_dist*-1) as d from %(track)s_tss where closest_dist=dist5 """
        statement2 = """SELECT closest_dist as d from %(track)s_tss where closest_dist=dist3 """

        data1 = self.getValues(statement1)
        data2 = self.getValues(statement2)

        return {"Genomic_distance": data1 + data2}

##########################################################################


class TTSProfile(cpgReport.cpgTracker):

    """Get profile around TTS"""
    mPattern = "_tts$"

    def __call__(self, track, slice=None):
        statement1 = """SELECT (closest_dist*-1) as d from %(track)s_tts where closest_dist=dist5 """
        statement2 = """SELECT closest_dist as d from %(track)s_tts where closest_dist=dist3 """

        data1 = self.getValues(statement1)
        data2 = self.getValues(statement2)

        return {"Genomic_distance": data1 + data2}
