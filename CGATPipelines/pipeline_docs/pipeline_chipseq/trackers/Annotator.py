import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
from ChipseqReport import *
import Annotations

from SphinxReport.Tracker import *

##########################################################################
##########################################################################
##########################################################################
# generic implementations of Annotator results
##########################################################################


class AnnotatorSlicer(object):

    def getSlices(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def addFilter(self, statement, track, slice):

        statement += " AND track = '%(track)s' " % locals()
        order = self.order

        if slice == "all" or slice == None:
            statement += " ORDER BY %(order)s""" % locals()
        elif slice:
            o = collections.defaultdict(list)
            if ":" in slice:
                pairs = slice.split(",")
                for pair in pairs:
                    key, value = pair.split(":")
                    o[key] = value

                if "workspace" in o and workspace not in o["workspace"]:
                    return None
                if "subset" in o and subset not in o["subset"]:
                    return None
                if "slice" in o and xslice not in o["slice"]:
                    return None
                if "track" in o and xtrack not in o["track"]:
                    return None

                statement += """ AND slice = '%(slice)s' 
                        AND subset='%(subset)s' AND workspace = '%(workspace)s' ORDER BY %(order)s""" % locals()

            else:
                categories = slice.split("-")
                statement += "AND category IN ('%s')" % ("','".join(categories))

        return statement


class Annotator(AnnotatorSlicer, TrackerSQL):

    '''
    observed
       observed percentage of bases in :term:`segments` that overlap annotations
    expected
       expected percentage of bases in :term:`segments` that overlap annotations
    '''

    mPattern = "_annotations$"
    tablename = None
    columns = ("over", "category", "fold", "pvalue",
               "qvalue", "observed", "expected")
    order = "pvalue"
    fdr = 0.10

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")
        all_tracks = self.getValues(
            "SELECT DISTINCT track FROM %s" % self.tablename)
        return selectTracks(all_tracks, subset)

    def __call__(self, track, slice=None):

        if not self.tablename:
            raise NotImplementedError("table not specified")

        tablename = self.tablename
        fdr = self.fdr
        order = self.order

        statement = self.select % locals()
        statement = self.addFilter(statement, track, slice)

        data = self.get(statement)
        return odict(zip(self.columns, zip(*data)))


class AnnotatorTable(Annotator):
    select = """SELECT CASE over WHEN '+' THEN 'over' ELSE 'under' END, 
                        category, fold, pvalue, qvalue, observed, expected 
                 FROM %(tablename)s WHERE qvalue < %(fdr)f """


class AnnotatorFullTable(Annotator):
    select = """SELECT CASE over WHEN '+' THEN 'over' ELSE 'under' END, 
                        category, fold, pvalue, qvalue, observed, expected 
                 FROM %(tablename)s WHERE 1 """


class AnnotatorSummary(Annotator):
    select = """SELECT COUNT(*), SUM( CASE over WHEN '+' THEN 1 ELSE 0 END) 
                       FROM %(tablename)s WHERE qvalue < %(fdr)f """
    columns = ("total", "over", "under")

    def __call__(self, track, slice=None):

        if not self.tablename:
            raise NotImplementedError("table not specified")

        tablename = self.tablename
        select = self.select
        fdr = self.fdr

        statement = self.select % locals()
        statement = self.addFilter(statement, track, slice)

        if not statement:
            return []

        data = self.getFirstRow(statement)
        data.append(data[0] - data[1])
        return odict(zip(self.columns, data))


class AnnotatorEnrichment(Annotator):
    select = """SELECT category, 100.0 * (fold - 1), pvalue
                 FROM %(tablename)s WHERE qvalue <= %(fdr)f AND over = '+' """
    columns = ("category", "enrichment/%", "pvalue")
    order = "fold DESC"


class AnnotatorDepletion(Annotator):
    select = """SELECT category, 100.0 * (1 - fold), pvalue
              FROM %(tablename)s WHERE qvalue <= %(fdr)f AND over = '-' """
    columns = ("category", "depletion/%", "pvalue")
    order = "fold"


class AnnotatorPower(TrackerSQL):
    mPattern = "_annotation"
    tablename = None

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || workspace FROM %s ORDER BY track,slice,subset,workspace" % self.tablename)

    def getSlices(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def __call__(self, track, slice=None):

        if not self.tablename:
            raise NotImplementedError("table not specified")
        if not self.select:
            raise NotImplementedError(
                "invalid use of base class - only use derived classes")

        select = self.select
        tablename = self.tablename
        fdr = self.fdr

        xtrack, xslice, subset, workspace = track.split(":")

        stmt = self.select % locals() +\
            """ AND track = '%(xtrack)s' AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' AND workspace = '%(workspace)s'
                        ORDER BY category""" % locals()

        data = self.getValues(stmt)

        return odict((("data", data),))


class AnnotatorPowerEnrichment(AnnotatorPower):

    """display power of go analyses.

    Show are the fold enrichment (in %) that could be detected
    given the setup of the test (size of sample and size of workspace).
    """
    select = "SELECT 100.0 * ( (ci95high / expected) - 1) FROM %(tablename)s WHERE 1"
    mXLabel = "Power : fold enrichment / %"


class AnnotatorPowerDepletion(AnnotatorPower):

    """display power of go analyses.

    Show are the fold depletion (in %) that could be detected
    given the setup of the test (size of sample and size of workspace).
    """
    select = "SELECT 100.0 * ( 1 - (ci95low / expected)) FROM %(tablename)s WHERE 1"
    mXLabel = "Power : fold depletion / %"


class AnnotatorMatrix(TrackerSQL):

    '''return enrichment/depletion for each category and track.

    Enrichment is returned as percent (100 = 2-fold enrichment) and
    depletion is returned as negative percent of the depletion (-50 = 2-fold
    depletion). Depletion is bounded by -100, while enrichment is unbounded.

    Slicing is possible using a keyword:value syntax. Keyword
    can be one of 'workspace', 'slice', 'subset' or 'track'. For example,
    the slices workspace:interegenic,slice:unknown will only display
    results where workspace is intergenic and slice is unknown.

    Only results passing the FDR rate test are shown. 
    '''

    tablename = None
    select = '''SELECT category, CASE over WHEN '+' THEN (100.0*(fold-1)) ELSE (-100.0 * (1-fold)) END 
                 FROM %(tablename)s WHERE qvalue <= %(fdr)f'''
    columns = ("enrichment", "pvalue")
    fdr = 0.10

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getSlices(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def getTracks(self, subset=None):
        if not self.tablename:
            raise NotImplementedError("table not specified")

        if subset:
            tracks = DefaultTracker().getTracks(subset)
            return self.getValues( """SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || workspace 
                                    FROM %s WHERE track IN ('%s') ORDER BY track,slice,subset,workspace""" %
                                   (self.tablename, "','".join(tracks)))
        else:
            return self.getValues( """SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || workspace 
                                      FROM %s ORDER BY track,slice,subset,workspace""" % self.tablename )

    def __call__(self, track, slice=None):

        if not self.tablename:
            raise NotImplementedError("table not specified")

        select = self.select
        tablename = self.tablename
        fdr = self.fdr

        xtrack, xslice, subset, workspace = track.split(":")

        category_sql = "AND category NOT like '%%:unknown' "

        if slice:
            o = collections.defaultdict(list)
            if ":" in slice:
                pairs = slice.split(",")
                for pair in pairs:
                    key, value = pair.split(":")
                    o[key] = value

                if "workspace" in o and workspace not in o["workspace"]:
                    return []
                if "subset" in o and subset not in o["subset"]:
                    return []
                if "slice" in o and xslice not in o["slice"]:
                    return []
                if "track" in o and xtrack not in o["track"]:
                    return []
            else:
                categories = slice.split("-")
                category_sql = "AND category IN ('%s')" % (
                    "','".join(categories))

        # remove the unknown tracks
        stmt = self.select % locals() +\
            """ AND track = '%(xtrack)s' AND slice = '%(xslice)s' 
        AND subset='%(subset)s' AND workspace = '%(workspace)s'
        %(category_sql)s
        ORDER BY category""" % locals()

        data = self.get(stmt)
        return odict(data)


class AnnotatorFullMatrix(AnnotatorMatrix):

    '''return enrichment/depletion for each category and track.

    Enrichment is returned as percent (100 = 2-fold enrichment) and
    depletion is returned as negative percent of the depletion (-50 = 2-fold
    depletion). Depletion is bounded by -100, while enrichment is unbounded.
    A value of -100 generally means `no-overlap`.

    Slicing is possible using a keyword:value syntax. Keyword
    can be one of 'workspace', 'slice', 'subset' or 'track'. For example,
    the slices workspace:interegenic,slice:unknown will only display
    results where workspace is intergenic and slice is unknown.

    All results are returned (no matter of significance).
    '''

    select = """SELECT category, CASE over WHEN '+' THEN (100.0*(fold-1)) ELSE (-100.0 * (1-fold)) END 
                 FROM %(tablename)s WHERE 1"""

# =================================================================
# mixin classes to pair annotator results with a table.
# =================================================================


class AnnotatorPromotorsGO(object):
    tablename = "annotator_promotors_go"


class AnnotatorPromotorsGOSlim(object):
    tablename = "annotator_promotors_goslim"


class AnnotatorTerritoriesGO(object):
    tablename = "annotator_territories_go"


class AnnotatorTerritoriesGOSlim(object):
    tablename = "annotator_territories_goslim"


class AnnotatorArchitecture(object):
    tablename = "annotator_architecture"


class AnnotatorArchitectureWithMotif(object):
    tablename = "annotator_architecture_with_motif"


class AnnotatorArchitectureWithoutMotif(object):
    tablename = "annotator_architecture_without_motif"


class AnnotatorTracks(object):
    tablename = "annotator_tracks"


class AnnotatorTracksWithMotif(object):
    tablename = "annotator_tracks_with_motif"


class AnnotatorTracksWithoutMotif(object):
    tablename = "annotator_tracks_without_motif"


class AnnotatorRegionsOfInterest(object):
    tablename = "annotator_roi"


class AnnotatorRegionsOfInterestWithMotif(object):
    tablename = "annotator_roi_with_motif"


class AnnotatorRegionsOfInterestWithoutMotif(object):
    tablename = "annotator_roi_without_motif"


class AnnotatorRegionsOfInterestTop(object):
    tablename = "annotator_roi_top"


class AnnotatorRegionsOfInterestBottom(object):
    tablename = "annotator_roi_bottom"


class AnnotatorRegionsOfInterestGO(object):
    tablename = "annotator_roi_go"


class AnnotatorRegionsOfInterestGOSlim(object):
    tablename = "annotator_roi_goslim"


class AnnotatorRegionsOfInterestOverlapGO(object):
    tablename = "annotator_roi_overlap_go"


class AnnotatorRegionsOfInterestOverlapGOSlim(object):
    tablename = "annotator_roi_overlap_goslim"


class AnnotatorGWAS(object):
    tablename = "annotator_gwas"


class AnnotatorSelection(object):
    tablename = "annotator_selection"

# =================================================================
# specific implementations of Annotator results
# =================================================================
_annotator_analysis = {
    "Enrichment": AnnotatorEnrichment,
    "Depletion": AnnotatorDepletion,
    "Summary": AnnotatorSummary,
    "Table": AnnotatorTable,
    "FullTable": AnnotatorFullTable,
    "Matrix": AnnotatorMatrix,
    "FullMatrix": AnnotatorFullMatrix,
    "PowerEnrichment": AnnotatorPowerEnrichment,
    "PowerDepletion": AnnotatorPowerDepletion,
}

_annotator_territories = {
    "PromotorsGO": AnnotatorPromotorsGO,
    "PromotorsGOSlim": AnnotatorPromotorsGOSlim,
    "TerritoriesGO": AnnotatorTerritoriesGO,
    "TerritoriesGOSlim": AnnotatorTerritoriesGOSlim,
    "GWAS": AnnotatorGWAS,
    "Selection": AnnotatorSelection,
    "Architecture": AnnotatorArchitecture,
    "ArchitectureWithMotif": AnnotatorArchitectureWithMotif,
    "ArchitectureWithoutMotif": AnnotatorArchitectureWithoutMotif,
    "Tracks": AnnotatorTracks,
    "TracksWithMotif": AnnotatorTracksWithMotif,
    "TracksWithoutMotif": AnnotatorTracksWithoutMotif,
    "RegionsOfInterest": AnnotatorRegionsOfInterest,
    "RegionsOfInterestWithMotif": AnnotatorRegionsOfInterestWithMotif,
    "RegionsOfInterestWithoutMotif": AnnotatorRegionsOfInterestWithoutMotif,
    "RegionsOfInterestTop": AnnotatorRegionsOfInterestTop,
    "RegionsOfInterestBottom": AnnotatorRegionsOfInterestBottom,
    "RegionsOfInterestGO": AnnotatorRegionsOfInterestGO,
    "RegionsOfInterestGOSlim": AnnotatorRegionsOfInterestGOSlim,
    "RegionsOfInterestOverlapGO": AnnotatorRegionsOfInterestOverlapGO,
    "RegionsOfInterestOverlapGOSlim": AnnotatorRegionsOfInterestOverlapGOSlim,
}

# the order of the base classes is important
# also: make sure that these are new-style classes
for a, aa in _annotator_analysis.items():
    for b, bb in _annotator_territories.items():
        n = "Annotator%s%s" % (a, b)
        globals()[n] = type(n, (bb, aa), {})
