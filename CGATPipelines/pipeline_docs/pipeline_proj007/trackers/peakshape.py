import os
import sys
import re
import types
import itertools
import math
import numpy
import IOTools

from CGATReport.Tracker import *
from cpgReport import *

##########################################################################
##########################################################################
##########################################################################
# peak shape


class PeakShapeTracker(Tracker):

    '''return peakshape data.  '''
    pattern = ".reads.peakshape.gz"
    scale = 2000

    def getTracks(self):
        pattern = self.pattern
        return [os.path.basename(x)[:-len(pattern)] for x in glob.glob(os.path.join(DATADIR, "liver_vs_testes", "*" + pattern))]

    slices = ["peak_height", "peak_width", "interval_score", "interval_width"]

    def __call__(self, track, slice=None):
        pattern = self.pattern
        fn = os.path.join(
            DATADIR, "liver_vs_testes/%(track)s%(pattern)s.matrix_%(slice)s.gz" % locals())
        if not os.path.exists(fn):
            return

        x = IOTools.openFile(fn)
        matrix, rownames, colnames = IOTools.readMatrix(x)

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

##########################################################################


class PeakShapeTrackerCentre(PeakShapeTracker):

    '''return peakshape data.  '''
    pattern = ".centre.peakshape.gz"
    scale = 2000

##########################################################################


class PeakShapeTrackerCentreNoScale(PeakShapeTracker):

    '''return peakshape data.  '''
    pattern = ".centre.peakshape.gz"
    scale = 1000000

##########################################################################


class PeakShapeSummary(Tracker):

    '''summary information about peak shapes.'''
    pattern = "(.*)_peakshape"
