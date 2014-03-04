import os
import sys
import re
import types
import itertools
from IntervalReport import *

import Annotations

##########################################################################
##########################################################################
##########################################################################
# Make histograms of the interval GC contents
##########################################################################


class pGC(Annotations.AnnotationSlicer, IntervalTracker):

    """for each interval, return the GC percentage"""

    mXLabel = "distance / bases"
    mPattern = "_intervals$"
    mColumn = "pGC"
    mTable = "nuc"
    mAnnotations = "annotations"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table = self.mTable
        column = self.mColumn
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_%(table)s AS d """ % locals() )

        return odict(((self.mColumn, data),))

##########################################################################
##########################################################################
##########################################################################
# Make histograms of the interval GC contents
##########################################################################


class TSS(Annotations.AnnotationSlicer, IntervalTracker):

    """for each interval, return the GC percentage"""

    mXLabel = "distance / bases"
    mPattern = "_intervals$"
    mColumn = "pGC"
    mTable1 = "tss_withoverlap_nuc"
    mTable2 = "tss_woutoverlap_nuc"
    mAnnotations = "annotations"

    def getNormValues(self, track, table):
        table1 = self.mTable1
        table2 = self.mTable2
        column = self.mColumn
        rawdata = self.getValues(
            """SELECT %(column)s FROM %(track)s_%(table)s AS d """ % locals() )
        data = [float(x) for x in rawdata if x != "None"]
        return(data)

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table1 = self.mTable1
        table2 = self.mTable2
        column = self.mColumn
        if not slice or slice == "all":
            d1 = self.getNormValues(track, table1)
            d2 = self.getNormValues(track, table2)

        return odict((("with overlap", d1), ("without overlap", d2)))


##########################################################################
##########################################################################
##########################################################################
# Make scatter plots of the interval GC vs interval average counts!
##########################################################################
class GCvsCounts(Annotations.AnnotationSlicer, IntervalTracker):

    """for each interval, return the GC percentage and the number of counts"""

    mXLabel = "distance / bases"
    mPattern = "_intervals$"
    mColumn1 = "pGC"
    mColumn2 = "avgval"
    mTable1 = "nuc"
    mTable2 = "intervals"
    mAnnotations = "annotations"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table1 = self.mTable1
        table2 = self.mTable2
        column1 = self.mColumn1
        column2 = self.mColumn2
        if not slice or slice == "all":

            data = self.getAll(
                """SELECT %(column1)s,%(column2)s FROM %(track)s_%(table1)s AS t1, %(track)s_%(table2)s AS t2 WHERE t1.contig==t2.contig and t1.start==t2.start and t1.end==t2.end """ % locals() )

        return data  # odict( ((self.mColumn, data),) )

# select avgval,pGC from macs_nuc as n, macs_intervals as i where
# n.end==i.end and n.start=i.start and n.contig==i.contig limit 1;
