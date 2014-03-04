import os
import sys
import re
import types
import itertools
import math

from IntervalReport import *

##########################################################################
##########################################################################
##########################################################################
# peak shape
##########################################################################


class PeakShapeTracker(Tracker):

    '''return peakshape data.

    Only 1000 rows are returned.
    '''

    tracks = [os.path.basename(x)[:-len(".peakshape.tsv.gz")]
              for x in glob.glob(os.path.join(DATADIR, "peakshapes", "*.peakshape.tsv.gz"))]
    slices = ["peak_height", "peak_width"]

    scale = 1000

    def __call__(self, track, slice=None):
        fn = os.path.join(
            DATADIR, "peakshapes/%(track)s.peakshape.tsv.gz.matrix_%(slice)s.gz" % locals())
        if not os.path.exists(fn):
            return

        matrix, rownames, colnames = IOTools.readMatrix(IOTools.openFile(fn))
        nrows = len(rownames)
        if nrows == 0:
            return

        if nrows > self.scale:
            take = numpy.array(
                numpy.floor(numpy.arange(0, nrows, float(nrows + 1) / self.scale)), dtype=int)
            rownames = [rownames[x] for x in take]
            matrix = matrix[take]

        return odict((('matrix', matrix),
                      ('rows', rownames),
                      ('columns', colnames)))


class PeakShapeSummary(Tracker):

    '''summary information about peak shapes.'''
    pattern = "(.*)_peakshape"
