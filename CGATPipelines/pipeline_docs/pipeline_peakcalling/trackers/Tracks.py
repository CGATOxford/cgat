'''Trackers for analysing chromatin tracks.'''

import os
import sys
import re
import types
import itertools
import glob
from PeakcallingReport import *
from CGATReport.Tracker import *


class TrackerTracks(TrackerSQL):

    mPattern = "_tracks$"

    def getTracks(self, subset=None):
        return [x for x in TrackerSQL.getTracks(self, subset) if "annotator" not in x]

    def getSlices(self, subset):
        return self.getValues("SELECT DISTINCT set2 FROM ucsc_overlap")


class TracksPeakval(TrackerTracks):

    '''return peakval for intervals overlapping and non-overlapping
    a track.'''

    def __call__(self, track, slice):

        result = odict()
        result["overlapping"] = self.getValues( '''
        SELECT i.peakval FROM 
        %(track)s_intervals AS i, 
        %(track)s_tracks AS t
        WHERE t.gene_id = i.interval_id 
        AND t.%(slice)s_nover > 0''' % locals())

        result["non-overlapping"] = self.getValues( '''
        SELECT i.peakval FROM 
        %(track)s_intervals AS i, 
        %(track)s_tracks AS t
        WHERE t.gene_id = i.interval_id 
        AND t.%(slice)s_nover = 0''' % locals())

        return result
