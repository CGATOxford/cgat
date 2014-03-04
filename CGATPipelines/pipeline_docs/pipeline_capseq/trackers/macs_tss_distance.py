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
##########################################################################
##########################################################################


class transcriptTSSOverlap(cpgReport.cpgTracker):

    '''number of transcript TSSs that an interval overlaps.'''
    mPattern = "_merged_tss$"
    mAnnotations = "merged_annotations"
    mTable = "merged_tss"
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
    mPattern = "_merged_tss$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "merged_annotations"
    mTable = "merged_tss"

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


class transcriptTSSClosestUpstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist5 is null THEN 1000000 ELSE dist5 END as dist5 "

##########################################################################


class transcriptTSSClosestDownstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist3 is null THEN 1000000 ELSE dist3 END as dist3 "

##########################################################################
##########################################################################
##########################################################################


class geneTSSOverlap(cpgReport.cpgTracker):

    '''number of gene TSSs that an interval overlaps.'''
    mPattern = "_merged_gene_tss$"
    mAnnotations = "merged_annotations"
    mTable = "merged_gene_tss"
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


class geneTSSClosest(TSSClosest):

    """for each interval, return the distance to the closest TSS."""
    mPattern = "_merged_gene_tss$"
    mColumn = "d.closest_dist"
    mTable = "merged_gene_tss"

##########################################################################


class geneTSSClosestUpstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mPattern = "_merged_gene_tss$"
    mTable = "merged_gene_tss"
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist5 is null THEN 1000000 ELSE dist5 END as dist5 "

##########################################################################


class geneTSSClosestDownstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mPattern = "_merged_gene_tss$"
    mTable = "merged_gene_tss"
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist3 is null THEN 1000000 ELSE dist3 END as dist3 "
