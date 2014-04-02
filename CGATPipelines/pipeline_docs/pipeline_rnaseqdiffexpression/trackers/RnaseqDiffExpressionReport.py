import os
import sys
import re
import types
import itertools
import glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

# get from config file
UCSC_DATABASE = P["genome"]
ENSEMBL_DATABASE = P["ensembl_database"]
RX_ENSEMBL_GENE = re.compile(P["ensembl_gene_prefix"])
RX_ENSEMBL_TRANSCRIPT = re.compile(P["ensembl_transcript_prefix"])

REFERENCE = "refcoding"

###################################################################
###################################################################
# parameterization

EXPORTDIR = P.get(
    'rnaseqdiffexpression_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('rnaseqdiffexpression_datadir', P.get('datadir', '.'))
DATABASE = P.get('rnaseqdiffexpression_backend', P.get(
    'sql_backend', 'sqlite:///./csvdb'))

DATABASE_ANNOTATIONS = P['annotations_database']

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.bam" % DATADIR), "(\S+).bam")

ALL = PipelineTracks.Aggregate(TRACKS)
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

GENESETS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"), "(\S+).gtf.gz")

DESIGNS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("design*.tsv"), "(\S+).tsv")

METHODS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("*_stats.tsv"), "(\S+)_stats.tsv")

###########################################################################
CUFFDIFF_LEVELS = ("gene", "isoform", "cds", "tss")

###########################################################################
# shorthand
MAP_TRACKS = {
    'default': EXPERIMENTS,
    'experiments': EXPERIMENTS,
    'conditions': CONDITIONS,
    'tissues': TISSUES,
    'merged': ALL,
    'geneset-summary': GENESETS}

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


class ProjectTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self,
                            *args,
                            backend=DATABASE,
                            attach=[(DATABASE_ANNOTATIONS, 'annotations')],
                            **kwargs)
