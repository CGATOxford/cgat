import os
import sys
import re
import types
import itertools
import Annotations

from SphinxReport.Tracker import *
from IntervalReport import *

##########################################################################
##########################################################################
##########################################################################
# Looking at distance
##########################################################################


class TSSClosest(Annotations.AnnotationSlicer, IntervalTracker):

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
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        return odict(((self.mColumn, data),))


class TSSDistanceVersusPeakVal(Annotations.AnnotationSlicer, IntervalTracker):

    """for each interval, return peakval and the distance to the closest TSS."""
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
            print """SELECT %(column)s, i.peakval FROM %(track)s_%(table)s AS d, %(track)s_intervals AS i 
                                      WHERE i.interval_id = d.gene_id AND %(where)s""" % locals()
            data = self.get( """SELECT %(column)s, i.peakval FROM %(track)s_%(table)s AS d, %(track)s_intervals AS i 
                                      WHERE i.interval_id = d.gene_id AND %(where)s""" % locals() )
        else:
            data = self.get( """SELECT %(column)s, i.peakval FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a, %(track)s_intervals AS i 
                                      WHERE i.interval_id = d.gene_id AND d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        return odict(zip(("distance", "peakval"), zip(*data)))


class TSSClosestUpstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mColumn = "d.upstream_dist"
    mWhere = "d.upstream_dist > 0"


class TSSClosestDownstream(TSSClosest):

    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mColumn = "d.downstream_dist"
    mWhere = "d.downstream_dist > 0"


class TSSOverlap(Annotations.AnnotationSlicer, IntervalTracker):

    '''number of intervals that overlap 0 or more than one transcription start sites.'''
    mPattern = "_tss$"
    mAnnotations = "annotations"
    mTable = "tss"
    mColumn = "d.is_overlap"
    mWhere = 1

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d 
                                      WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, 
                                      %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id 
                                      AND a.is_%(slice)s AND %(where)s""" % locals() )

        if len(data) == 0:
            return None

        return odict((("no overlap", len([x for x in data if x == 0])),
                      ("overlap", len([x for x in data if x > 0]))))


class TSSDistances(IntervalTracker):
    pattern = "(.*)_annotations$"

    def getSlices(self, subset=None):
        return ["5", "3"]

    def getResultsDir(self, track):
        return os.path.join(EXPORTDIR, "annotator_distance", "%s-tss.gtf.annodist" % track)

    def __call__(self, track, slice=None):

        resultsdir = self.getResultsDir(track)
        if not os.path.exists(resultsdir):
            return []

        data = []
        headers = None
        for line in open(os.path.join(resultsdir, "proximity.table")):
            d = line[:-1].split("\t")
            label = d[0]

            if line.startswith("label"):
                if not headers:
                    headers = d
            elif line.startswith(slice):
                data.append(
                    ("plots", "`plot <%(resultsdir)s/%(slice)s.png>`_" % locals()))
                data.extend([(x, y) for x, y in zip(headers, d)[1:]])

        return odict(data)


class IntergenicDistances(TSSDistances):

    def getResultsDir(self, track):
        return os.path.join(EXPORTDIR, "annotator_distance", "%s-intergenic.annodist" % track)


class IntronicDistances(TSSDistances):

    def getResultsDir(self, track):
        return os.path.join(EXPORTDIR, "annotator_distance", "%s-intronic.annodist" % track)


class GenicDistances(TSSDistances):

    def getResultsDir(self, track):
        return os.path.join(EXPORTDIR, "annotator_distance", "%s-genic.annodist" % track)
