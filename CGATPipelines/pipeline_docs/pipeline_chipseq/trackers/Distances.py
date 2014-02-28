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

##########################################################################
##########################################################################
##########################################################################
# Looking at distance
##########################################################################


class ClosestDistance(Annotations.AnnotationSlicer, DefaultTracker):

    """Closest distance of transcript models to gene models in the reference set."""
    mXLabel = "distance / bases"
    mPattern = "_distances$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "annotations"

    def __call__(self, track, slice=None):

        annotations = self.mAnnotations

        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_distances AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_distances AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )
        return odict(data)


class ClosestDistanceTo5(ClosestDistance):

    """Closest distance of transcript models at 5' end of a gene."""
    mColumn = "d.amin5"
    mWhere = "d.amin5 > 0"


class ClosestDistanceTo3(ClosestDistance):

    """Closest distance of transcript models at 3' end of a gene."""
    mColumn = "d.amin3"
    mWhere = "d.amin3 > 0"
