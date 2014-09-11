import os
import sys
import re
import types
import itertools
import glob
import PipelineTracks

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict
from CGATReport.Utils import PARAMS as P

EXPORTDIR = P['cpg_exportdir']
DATADIR = P['cpg_datadir']
DATABASE = P['cpg_backend']
UCSC_GENOME = P['ucsc_genome']
ANNOTATIONS_DB = P['annotations_db']

##########################################################################
Sample = PipelineTracks.Sample3
Sample.setDefault("asTable")
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    [x for x in glob.glob("*.export.txt.gz") if P["tracks_control"] not in x],
    "(\S+).export.txt.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.sra") if P["tracks_control"] not in x],
        "(\S+).sra" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.fastq.gz") if P["tracks_control"] not in x],
        "(\S+).fastq.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.fastq.1.gz") if P["tracks_control"] not in x],
        "(\S+).fastq.1.gz" ) +\
    PipelineTracks.Tracks(PipelineTracks.Sample3).loadFromDirectory(
        [x for x in glob.glob("*.csfasta.gz") if P["track_control"] not in x],
        "(\S+).csfasta.gz")

Sample.setDefault("asTable")

for X in TRACKS:
    print "TRACK=", X, "\n"

ALL = PipelineTracks.Aggregate(TRACKS)
EXPERIMENTS = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))
CONDITIONS = PipelineTracks.Aggregate(TRACKS, labels=("condition", ))
TISSUES = PipelineTracks.Aggregate(TRACKS, labels=("tissue", ))

MAP_TRACKS = {
    'master': map(str, list(EXPERIMENTS) + list(CONDITIONS)),
    'replicates': map(str, list(TRACKS)),
    'default': map(str, list(EXPERIMENTS)),
    'experiments': map(str, list(EXPERIMENTS)),
    'conditions': map(str, list(CONDITIONS)),
    'tissues': map(str, list(TISSUES)),
    'merged': map(str, list(EXPERIMENTS)),
}

for x, y in MAP_TRACKS.iteritems():
    print "MAP_TRACK=", x, "--", y, "\n"


##########################################################################
class cpgTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

    def getTracks(self, subset=None):
        if subset:
            for key, tracks in MAP_TRACKS.iteritems():
                if key in subset:
                    return tracks

        return TrackerSQL.getTracks(self)

    def connect(self):
        '''connect to database.'''

        if not self.db:
            # call base class connect
            TrackerSQL.connect(self)

            statement = "ATTACH DATABASE '%s' as annotations" % (
                ANNOTATIONS_DB)
            self.db.execute(statement)
