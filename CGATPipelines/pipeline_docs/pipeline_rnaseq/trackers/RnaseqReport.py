import re
import glob

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import PipelineTracks

# get from config file
UCSC_DATABASE = "hg19"
ENSEMBL_DATABASE = "Homo_sapiens"
RX_ENSEMBL_GENE = re.compile("ENSG")
RX_ENSEMBL_TRANSCRIPT = re.compile("ENST")

REFERENCE = "refcoding"

###################################################################
###################################################################
# parameterization

EXPORTDIR = P['rnaseq_exportdir']
DATADIR = P['rnaseq_datadir']
DATABASE = P['rnaseq_backend']

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################


TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
    glob.glob("%s/*.sra" % DATADIR), "(\S+).sra") +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % DATADIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % DATADIR), "(\S+).fastq.1.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        glob.glob("*.csfasta.gz"), "(\S+).csfasta.gz")

ALL = PipelineTracks.Aggregate(TRACKS)
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

###########################################################################
# tracks for the gene sets


class GenesetTrack(PipelineTracks.Sample):
    attributes = ("geneset",)

GENESET_TRACKS = PipelineTracks.Tracks(GenesetTrack).loadFromDirectory(
    glob.glob("%s/*.cuffdiff" % DATADIR),
    "%s/(\S+).cuffdiff" % DATADIR)

CUFFDIFF_LEVELS = ("gene", "isoform", "cds", "tss")

###########################################################################
# shorthand
MAP_TRACKS = {
    'default': EXPERIMENTS,
    'experiments': EXPERIMENTS,
    'conditions': CONDITIONS,
    'tissues': TISSUES,
    'merged': ALL,
    'geneset-summary': GENESET_TRACKS}

###########################################################################


def selectTracks(subset):
    '''select tracks from *all_tracks* according to *subset*.
    '''
    if subset is None or subset == "default":
        return MAP_TRACKS["default"]
    elif subset in MAP_TRACKS:
        return MAP_TRACKS[subset]

    return subset

###########################################################################


def splitLocus(locus):
    if ".." in locus:
        contig, start, end = re.match("(\S+):(\d+)\.\.(\d+)", locus).groups()
    elif "-" in locus:
        contig, start, end = re.match("(\S+):(\d+)\-(\d+)", locus).groups()

    return contig, int(start), int(end)

###########################################################################


def linkToUCSC(contig, start, end):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i-%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link

###########################################################################


def linkToEnsembl(id):
    ensembl_database = ENSEMBL_DATABASE
    if RX_ENSEMBL_GENE.match(id):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_database)s/Gene/Summary?g=%(id)s>`_" \
            % locals()
    elif RX_ENSEMBL_TRANSCRIPT.match(id):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_database)s/Transcript/Summary?t=%(id)s>`_" \
            % locals()
    else:
        link = id
    return link

###########################################################################


class RnaseqTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)
