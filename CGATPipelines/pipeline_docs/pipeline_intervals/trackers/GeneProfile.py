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
# Make line plots of the gene read profiles .....
##########################################################################


class MultiProfilePlot(Annotations.AnnotationSlicer, IntervalTracker):

    """for each interval, return the GC percentage and the number of counts"""

    mXLabel = "distance / bases"
    mPattern = "_intervals$"
    mColumn1 = "upstream"
    mColumn2 = "exons"
    mColumn3 = "downstream"
    mTable1 = "withoverlap_geneprofile_counts"
    mTable2 = "woutoverlap_geneprofile_counts"
    mAnnotations = "annotations"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations
        table1 = self.mTable1
        table2 = self.mTable2
        column1 = self.mColumn1
        column2 = self.mColumn2
        column3 = self.mColumn3
        if not slice or slice == "all":
            withdata = self.getAll(
                """SELECT %(column1)s,%(column2)s,%(column3)s FROM %(track)s_%(table1)s""" % locals() )
            woutdata = self.getAll(
                """SELECT %(column1)s,%(column2)s,%(column3)s FROM %(track)s_%(table2)s""" % locals() )
            withd = withdata[column1] + withdata[column2] + withdata[column3]
            witho = woutdata[column1] + woutdata[column2] + woutdata[column3]
            wd = odict(
                (("position", range(0, len(withd))), ("density", withd)))
            wo = odict(
                (("position", range(0, len(witho))), ("density", witho)))

            data = odict((("with_overlap", wd), ("without_overlap", wo)))

        return data
