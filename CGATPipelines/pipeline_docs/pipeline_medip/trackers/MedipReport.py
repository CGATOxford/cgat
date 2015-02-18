import itertools
import glob
import Pipeline
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import PipelineTracks

# get from config file
UCSC_DATABASE = "hg19"
EXPORTDIR = "export"

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script


EXPORTDIR = P.get('medip_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('medip_datadir', P.get('datadir', '.'))
DATABASE = P.get('medip_backend', P.get('sql_backend', 'sqlite:///./csvdb'))

###################################################################
# cf. pipeline_medip.py
# This should be automatically gleaned from pipeline_chipseq.py
###################################################################
PARAMS_PIPELINE = Pipeline.peekParameters(".",
                                          "pipeline_medip.py")

Sample = PipelineTracks.Sample3

suffixes = ["export.txt.gz",
            "sra",
            "fastq.gz",
            "fastq.1.gz",
            "csfasta.gz"]

TRACKS = sum(itertools.chain([PipelineTracks.Tracks(Sample).loadFromDirectory(
    [x for x in glob.glob("%s/*.%s" % (DATADIR, s)) if "input" not in x],
    "%s/(\S+).%s" % (DATADIR, s)) for s in suffixes]),
    PipelineTracks.Tracks(Sample))

Sample.setDefault("asTable")

ALL = PipelineTracks.Aggregate(TRACKS)
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))


class MedipTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
