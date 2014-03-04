import os
import sys
import re
import types
import itertools
import math
import numpy

from MedipReport import *


class MedipsTracker(MedipTracker):
    pass

##############################################################
##############################################################
##############################################################


class MedipsPlots(MedipsTracker):

    tracks = [x.asFile() for x in TRACKS]

    slices = ("calibration",
              "cpg_coverage",
              "saturation")

    def __call__(self, track, slice=None):

        # note there are spaces behind the %(image)s directive to accomodate
        # for path substitution
        block = '''
.. figure:: %(image)s                                     
   :height: 300 
'''
        filename = os.path.join("%s.dir" % track,
                                "%s.prep.medips_%s.png" % (track, slice))

        image = os.path.abspath(filename)

        b = block % locals()

        return odict((("rst", b),))
