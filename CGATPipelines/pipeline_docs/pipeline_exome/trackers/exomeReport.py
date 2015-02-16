import re
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

# get from config file
UCSC_DATABASE = "hg19"
EXPORTDIR = "export"

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script
EXPORTDIR = P.get('exome_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('exome_datadir', P.get('datadir', '.'))
DATABASE = P.get('exome_backend', P.get('sql_backend', 'sqlite:///./csvdb'))

TRACKS = ['WTCHG_10997_01', 'WTCHG_10997_02']


###########################################################################
def splitLocus(locus):
    if ".." in locus:
        contig, start, end = re.match("(\S+):(\d+)\.\.(\d+)", locus).groups()
    elif "-" in locus:
        contig, start, end = re.match("(\S+):(\d+)\-(\d+)", locus).groups()

    return contig, int(start), int(end)


def linkToUCSC(contig, start, end):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link

###########################################################################


class ExomeTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class SingleTableHistogram(TrackerSQL):
    columns = None
    table = None
    group_by = None

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT %(group_by)s, %(columns)s FROM %(table)s")
        return data
