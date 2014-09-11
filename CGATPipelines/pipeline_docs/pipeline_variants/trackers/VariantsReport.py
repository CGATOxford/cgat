import os
import sys
import re
import types
import collections
import numpy
import numpy.ma
import Stats

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from CGATReport.odict import OrderedDict as odict

SNP_CODES = ("E", "O")

###################################################################
###################################################################
# parameterization

EXPORTDIR = P['variants_exportdir']
DATADIR = P['variants_datadir']
DATABASE = P['variants_backend']

ANNOTATIONS_DB = P['annotations_database']

###########################################################################


class VariantsTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

    def connect(self):
        '''connect to database.'''

        if not self.db:
            # call base class connect
            TrackerSQL.connect(self)

            # there might a better way of doing this, see for example:
            # http://www.mail-archive.com/sqlalchemy@googlegroups.com/msg03606.html

            statement = "ATTACH DATABASE '%s' as annotations" % (
                ANNOTATIONS_DB)
            self.db.execute(statement)


class StrainTracker(VariantsTracker):

    '''tracker returning results per strain.'''
    pattern = "(.*)_effects$"

    def getTracks(self, subset=None):

        tracks = TrackerSQL.getTracks(self)

        if subset is None:
            return ["all"] + tracks
        else:
            return tracks
