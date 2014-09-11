from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from collections import OrderedDict as odict

###################################################################
# parameterization

EXPORTDIR = P.get(
    'genesets_exportdir',
    P.get('exportdir', 'export'))
DATADIR = P.get(
    'genesets_datadir',
    P.get('datadir', '.'))
DATABASE = P.get(
    'genesets_backend',
    P.get('sql_backend', 'sqlite:///./csvdb'))


class GeneSetsTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
