import math
import sys
import os
import collections

from CGATReport.Tracker import *

##########################################################################
##########################################################################
##########################################################################
# Display results from annotator_distance.py
##########################################################################


class AnnotatorDistance(TrackerSQL):
    mPattern = "_annotation"
    mTableName = "annotatordistance"
    mSelect = "SELECT label, CAST (observed AS FLOAT) / expected AS fold, observed, expected, CIlower, CIupper, pvalue, qvalue FROM %s WHERE qvalue < 0.05"
    mColumns = ("label", "fold", "observed", "expected",
                "cilower", "ciupper", "pvalue", "qvalue")
    mOrder = "fold"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track FROM %s" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT slice || ':' || subset || ':' || counter || ':' || workspace FROM %s" % self.mTableName)

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")

        select = self.mSelect % self.mTableName
        order = self.mOrder
        if slice == "all" or slice is None:
            data = list(self.execute(
                """%s AND track = '%s' ORDER BY %s""" % (select, track, order)).fetchone() )
        else:
            slice, subset, counter, workspace = slice.split(":")
            data = self.getAll( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s' 
                        AND subset='%(subset)s' AND counter = '%(counter)s' AND workspace = '%(workspace)s' ORDER BY %(order)s""" % locals())

        return odict(zip(self.mColumns, zip(*data)))


class AnnotatorDistanceVolcano(TrackerSQL):

    """return log(fold change) versus log(P-Value)."""
    mPattern = None
    mTableName = "annotatordistance"

    mXLabel = "fold change"
    mYLabel = "log P-Value"

    def getTracks(self):
        return ["intronic", "intergenic"]

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):

        table = self.mTableName
        data = self.getAll(
            "SELECT CAST( observed AS FLOAT)/ expected, max( pvalue, 0.00001)  FROM %(table)s WHERE workspace = '%(track)s'" % locals())

        data = [(x[0], math.log10(x[1])) for x in data]
        return odict(zip(("fold", "log(pvalue)"), zip(*data)))


class AnnotatorDistanceSummary(AnnotatorDistance):
    mSelect = "SELECT COUNT(*), SUM( CASE code WHEN '+' THEN 1 ELSE 0 END) FROM %s WHERE passed AND code != '?'"
    mColumns = ("total", "over", "under")

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")

        select = self.mSelect % self.mTableName

        if slice == "all" or slice is None:
            data = self.getFirstRow(
                """%s AND track = '%s'""" % (select, track, order))
        else:
            slice, subset, counter, workspace = slice.split(":")
            data = self.getFirstRow( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s' 
                        AND workspace = '%(workspace)s' 
                        AND subset='%(subset)s' AND counter = '%(counter)s'""" % locals())

        data.append(data[0] - data[1])
        return odict(zip(self.mColumns, data))


class AnnotatorDistanceEnrichment(AnnotatorDistance):
    mSelect = "SELECT label, 100.0 * (ratio - 1), pover FROM %s WHERE passed AND code = '+' "
    mColumns = ("category", "enrichment/%", "pvalue")
    mOrder = "ratio DESC"


class AnnotatorDistanceDepletion(AnnotatorDistance):
    mSelect = "SELECT label, 100.0 * (1 - ratio), punder FROM %s WHERE passed AND code = '-' "
    mColumns = ("category", "depletion/%", "pvalue")
    mOrder = "ratio"


class AnnotatorDistancePower(TrackerSQL):
    mPattern = "_annotation"
    mTableName = "annotatordistance"
    mSelect = "SELECT 100.0 * ( 1.0 - (cast(cilower as float) / expected) ) FROM %s WHERE 1"
    mXLabel = "Power / %"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || counter || ':' || workspace FROM %s ORDER BY track,slice,subset,counter" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not self.mSelect:
            raise NotImplementedError(
                "invalid use of base class - only use derived classes")

        select = self.mSelect % self.mTableName
        xtrack, xslice, subset, counter, workspace = track.split(":")

        stmt = """%(select)s AND track = '%(xtrack)s' 
                        AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' 
                        AND counter = '%(counter)s'
                        AND workspace = '%(workspace)s' 
                        ORDER BY label""" % locals()

        data = self.getValues(stmt)

        return odict((("power", data),))


class AnnotatorDistanceMatrix(TrackerSQL):

    """Display annotator distance results in a matrix.

    All tracks are output, but they are masked if they do
    not pass the FDR test (qvalue < 0.5).

    Slicing is possible using a keyword:value syntax. Keyword
    can be one of 'counter', 'slice', 'subset' or 'track'. For example,
    the slices counter:ensembl,slice:unknown will only display
    results where counter is ensembl and slice is unknown.
    """
    mTableName = "annotatordistance"
    mSelect = "SELECT label, CASE qvalue < 0.05 WHEN 1 THEN (100.0 * (1- (cast(observed as float) / expected)) ) ELSE 0 END FROM %s WHERE 1 "

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || counter || ':' || workspace FROM %s ORDER BY slice,track,subset,counter" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")
        select = self.mSelect % self.mTableName
        xtrack, xslice, subset, counter, workspace = track.split(":")

        if slice:
            o = collections.defaultdict(list)
            pairs = slice.split(",")
            for pair in pairs:
                key, value = pair.split(":")
                o[key] = value
            if "counter" in o and counter not in o["counter"]:
                return []
            if "subset" in o and subset not in o["subset"]:
                return []
            if "slice" in o and xslice not in o["slice"]:
                return []
            if "track" in o and xtrack not in o["track"]:
                return []
            if "workspace" in o and workspace not in o["workspace"]:
                return []

        stmt = """%(select)s 
                        AND track = '%(xtrack)s' 
                        AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' 
                        AND counter='%(counter)s'
                        AND workspace='%(workspace)s'
                        ORDER BY label""" % locals()

        data = self.getAll(stmt)

        return odict(data)
