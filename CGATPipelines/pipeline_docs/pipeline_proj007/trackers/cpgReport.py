import os
import sys
import re
import types
import itertools
import glob
import PipelineTracks

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from SphinxReport.Utils import PARAMS as P

EXPORTDIR = P['cpg_exportdir']
DATADIR = P['cpg_datadir']
DATABASE = P['cpg_backend']
UCSC_GENOME = P['ucsc_genome']
ANNOTATIONS_DB = P['annotations_db']
ANNOTATIONS_NAME = P['annotations_name']

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

        # issuing the ATTACH DATABASE into the sqlalchemy ORM (self.db.execute( ... ))
        # does not work. The database is attached, but tables are not accessible in later
        # SELECT statements.
        if not self.db:
            def _create():
                import sqlite3
                conn = sqlite3.connect(re.sub("sqlite:///", "", DATABASE))
                statement = "ATTACH DATABASE '%s' AS annotations; " % (
                    ANNOTATIONS_DB)
                conn.execute(statement)
                return conn

            self.connect(creator=_create)

    def getTracks(self, subset=None):
        if subset:
            for key, tracks in MAP_TRACKS.iteritems():
                if key in subset:
                    return tracks

        return TrackerSQL.getTracks(self)

##########################################################################


class featureOverlap(cpgTracker):

    """return overlap of interval with genomic features """

    mPattern = "_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_transcript_extended_pover1"

    def __call__(self, track, slice=None):
        table = self.mTable
        where = self.mWhere
        data = self.getValues( """ SELECT count(distinct gene_id) as intervals FROM (
                                   SELECT gene_id,
                                   CASE WHEN  %(where)s > 0  THEN 'TSS'
                                   WHEN genes_pover1 > 0 THEN 'Gene'
                                   WHEN upstream_flank_pover1 >0 THEN 'Upstream'
                                   WHEN downstream_flank_pover1 >0 THEN 'Downstream'
                                   ELSE 'Intergenic'
                                   END AS feature_class
                                   FROM %(track)s%(table)s)
                                   group by feature_class
                                   order by feature_class asc""" % locals() )

        result = odict(
            zip(("Downstream", "Gene", "Intergenic", "TSS", "Upstream"), data))
        return result
