import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

# get from config file
UCSC_DATABASE="hg19"
EXPORTDIR="export"

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script

from SphinxReport.Utils import PARAMS as P
EXPORTDIR=P['medip_exportdir']
DATADIR=P['medip_datadir']
DATABASE=P['medip_backend']

###################################################################
# cf. pipeline_medip.py
# This should be automatically gleaned from pipeline_chipseq.py
###################################################################
import Pipeline
PARAMS_PIPELINE = Pipeline.peekParameters( ".",
                                           "pipeline_medip.py" )

import PipelineTracks

Sample = PipelineTracks.Sample3

suffixes = ["export.txt.gz",
            "sra",
            "fastq.gz",
            "fastq.1.gz",
            "csfasta.gz" ]

TRACKS = sum( itertools.chain( [ PipelineTracks.Tracks( Sample ).loadFromDirectory( 
                [ x for x in glob.glob( "%s/*.%s" % (DATADIR, s) ) if "input" not in x],
                "%s/(\S+).%s" % (DATADIR, s) ) for s in suffixes ] ), 
              PipelineTracks.Tracks( Sample ) )

Sample.setDefault( "asTable" )

ALL = PipelineTracks.Aggregate( TRACKS )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###########################################################################
###########################################################################
###########################################################################
class MedipTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )


