import os
import sys
import re
import types
import itertools
import numpy

import CGAT.Stats as Stats
import CGAT.IndexedGenome as IndexedGenome

from CGATReport.Tracker import *
from PeakcallingReport import *

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class OverlapsBase(DefaultTracker):

    """Overlap between sets.

    This tracker returns the overlap between a track and 
    all other tracks. Only one attribute is returned 
    given by :attr:`mColumn`. As the table is not symmetrized,
    the mColumn attribute should be specified without the 
    suffix, i.e. ('nbases_unique' instead of 'nbases_unique1').

    Possible values for mColumn are combinations of 'A_B'.

    A is one of ('nbases', 'nexons', 'ngenes', 'pbases', 'pexons', 'pgenes') 
    where the prefix ''n'' or ''p'' denote the counts or the percent, respectively.

    B is one of ('total', 'uniq', 'ovl')
    """
    tablename = None
    column = None
    pattern = "(.*)_intervals$"

    def getSlices(self, subset=None):
        if subset is not None:
            return subset
        else:
            return []

    def __call__(self, track, slice=None):

        tablename = self.tablename
        if self.column is None:
            raise NotImplementedError("column not set in derived class.")
        column = self.column

        if slice is None:
            result = odict(self.get("SELECT set2, %(column)s1 FROM %(tablename)s WHERE set1 = '%(track)s'" % locals()) +
                           self.get("SELECT set1, %(column)s2 FROM %(tablename)s WHERE set2 = '%(track)s'" % locals()))
        elif "-" in slice:
            slices = "','".join(slice.split("-"))
            result = odict(self.get("SELECT set2, %(column)s1 FROM %(tablename)s WHERE set1 = '%(track)s' AND set2 IN ('%(slices)s')" % locals()) +
                           self.get("SELECT set1, %(column)s2 FROM %(tablename)s WHERE set2 = '%(track)s' AND set1 IN ('%(slices)s')" % locals()))
        else:
            result = odict(self.get("SELECT set2, %(column)s1 FROM %(tablename)s WHERE set1 = '%(track)s' AND set2 = '%(slice)s'" % locals()) +
                           self.get("SELECT set1, %(column)s2 FROM %(tablename)s WHERE set2 = '%(track)s' AND set1 = '%(slice)s'" % locals()))

        return result


class ExonsCounts(OverlapsBase):
    column = "nexons_ovl"


class ExonsPercent(OverlapsBase):
    column = "pexons_ovl"

    def __call__(self, track, slice=None):
        # add the diagonal element of 100%
        x = OverlapsBase.__call__(self, track, slice)
        x[track] = 100.0
        return x


class BasesCounts(OverlapsBase):
    column = "nbases_ovl"


class BasesPercent(OverlapsBase):
    column = "pbases_ovl"

    def __call__(self, track, slice=None):
        x = OverlapsBase.__call__(self, track, slice)
        x[track] = 100.0
        return x


class BasesNormalized(OverlapsBase):
    column = "nbases_ovl"

    def __call__(self, track, slice=None):
        tablename = self.tablename
        if self.column is None:
            raise NotImplementedError("column not set in derived class.")
        column = self.column
        return odict(self.get("SELECT set2, 100.0 * nbases_ovl1 / (nbases_ovl1 + nbases_unique1 + nbases_unique2) FROM %(tablename)s WHERE set1 = '%(track)s'" % locals()) +
                     self.get("SELECT set1, 100.0 * nbases_ovl2 / (nbases_ovl2 + nbases_unique1 + nbases_unique2) FROM %(tablename)s WHERE set2 = '%(track)s'" % locals()) +
                     [(track, 100.0)])

# =================================================================
# mixin classes for table
# =================================================================


class Overlaps(object):
    tablename = "overlap"


class OverlapsUCSC(object):
    tablename = "ucsc_overlap"

# =================================================================
# specific implementations of Annotator results
# =================================================================
_overlap_analysis = {
    "ExonsCounts": ExonsCounts,
    "ExonsPercent": ExonsPercent,
    "BasesCounts": BasesCounts,
    "BasesPercent": BasesPercent,
    "BasesNormalized": BasesNormalized,
}

_overlap_tables = {
    "Overlaps": Overlaps,
    "UCSCOverlaps": OverlapsUCSC,
}

# the order of the base classes is important
# also: make sure that these are new-style classes
for a, aa in _overlap_analysis.items():
    for b, bb in _overlap_tables.items():
        n = "%s%s" % (b, a)
        globals()[n] = type(n, (bb, aa), {})

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class OverlapVersusPeakval(DefaultTracker):

    """Overlap between experiments
    """
    tablename = "reproducibility"
    pattern = "(.*)_reproducibility$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT pexons_union, pexons_ovl FROM %(track)s_%(tablename)s" % locals())
        return odict(zip(("recall", "reproducibility"), zip(*data)))


class OverlapROC(DefaultTracker):

    """Overlap between experiments. 

    This tracker computes ROC curves examining various 
    interval variables to see if they improve reproducibility.

    True positives are those intervals which reproducible, i.e.,
    appear in all biological replicates.
    """
    pattern = "(.*)_reproducibility$"
    mFields = ("peakval", "avgval", "length")
    mXLabel = "FPR"
    mYLabel = "TPR"

    def __call__(self, track, slice=None):

        result = odict()

        merged = None
        rocs = []

        for field in self.mFields:
            data = []
            for replicate in EXPERIMENTS.getTracks(track):
                statement = "SELECT contig, start, end,%(field)s FROM %(replicate)s_intervals" % locals(
                )
                data.append(self.get(statement))

            idx = []
            for x in range(len(data)):
                i = IndexedGenome.IndexedGenome()
                for contig, start, end, peakval in data[x]:
                    i.add(contig, start, end, peakval)
                idx.append(i)

            def _iter(all):
                all.sort()
                last_contig, first_start, last_end, last_value = all[0]
                for contig, start, end, value in all[1:]:
                    if contig != last_contig or last_end < start:
                        yield (last_contig, first_start, last_end)
                        last_contig, first_start, last_end = contig, start, end
                    else:
                        last_end = max(last_end, end)
                yield (last_contig, first_start, last_end)

            if not merged:
                all = [x for x in itertools.chain(*data)]
                merged = list(_iter(all))

            roc_data = []
            for contig, start, end in merged:
                intervals = []
                for i in idx:
                    try:
                        intervals.append(list(i.get(contig, start, end)))
                    except KeyError:
                        continue

                if len(intervals) == 0:
                    continue

                is_repro = len([x for x in intervals if x != []]) == len(data)
                value = max([x[2] for x in itertools.chain(*intervals)])

                # fpr, tpr
                roc_data.append((value, is_repro))

            roc_data.sort()
            roc_data.reverse()

            roc = zip(*Stats.computeROC(roc_data))
            result[field] = odict((("FPR", roc[0]), (field, roc[1])))

        return result


class OverlapMatrix(DefaultTracker):
    pattern = "(.*)_reproducibility$"
    field = "pexons_ovl"

    def __call__(self, track):

        data = self.get(
            "SELECT set1, set2, %(field)s1, %(field)s2 FROM %(track)s_reproducibility")

        rows = sorted(
            list(set([x[0] for x in data]).union(set([x[1] for x in data]))))

        map_row2index = dict([(x[1], x[0]) for x in enumerate(rows)])

        matrix = numpy.zeros((len(rows), len(rows)))
        for row, col, value1, value2 in data:
            matrix[map_row2index[row]][map_row2index[col]] = value1
            matrix[map_row2index[col]][map_row2index[row]] = value2

        return odict((('matrix', matrix),
                      ('rows', rows),
                      ('columns', rows)))
