import os, sys, re, types, itertools, glob
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import sqlalchemy

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

###################################################################
###################################################################
## parameterization

EXPORTDIR=P.get('calling_exportdir', P.get( 'exportdir', 'export'))
DATADIR=P.get('calling_datadir', P.get( 'datadir', '.'))
DATABASE=P.get('calling_backend', P.get( 'sql_backend', 'sqlite:///./csvdb'))

###################################################################
# cf. pipeline_chipseq.py
# This should be automatically gleaned from pipeline_chipseq.py
###################################################################
import CGAT.Pipeline as Pipeline
PARAMS_PIPELINE = Pipeline.peekParameters( ".",
                                           "pipeline_chipseq.py" )

import CGATPipelines.PipelineTracks as PipelineTracks

Sample = PipelineTracks.Sample3

suffixes = ["export.txt.gz",
            "sra",
            "fastq.gz",
            "fastq.1.gz",
            "csfasta.gz" ]

TRACKS = sum( itertools.chain( [ PipelineTracks.Tracks( Sample ).loadFromDirectory( 
        [ x for x in glob.glob( "%s/*.%s" % (DATADIR, s) ) if "input" not in x ],
        "%s/(\S+).%s" % (DATADIR, s) ) for s in suffixes ] ), 
              PipelineTracks.Tracks( Sample ) )

Sample.setDefault( "asTable" )

ALL = PipelineTracks.Aggregate( TRACKS )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

############################################################################
# The folllowing need to be parameterized in a config file
# TISSUES=["GM00855", "GM00861" ]
# CONDITIONS=["D3", "unstim" ]
# REPLICATES=["R1", "R2" ]
TAG_UNSTIM = PARAMS_PIPELINE["tracks_unstimulated"]
UCSC_GENOME = PARAMS_PIPELINE["genome"]

if "motifs_plot" in P and P["motifs_plot"]:
    MOTIFS = [x.strip() for x in P["motifs_plot"].split(",")]
else:
    MOTIFS = None

###########################################################################
## shorthand
# use list to convert trackers to strings
MAP_TRACKS = {
    'master' : map( str, list(EXPERIMENTS) + list( CONDITIONS )),
    'replicates' : map( str, list(TRACKS) ),
    'default' : map(str, list(EXPERIMENTS)),
    'experiments' : map( str, list(EXPERIMENTS)),
    'conditions' : map( str, list(CONDITIONS)),
    'tissues' : map(str, list(TISSUES)),
    'merged' : map(str, list(EXPERIMENTS)), 
    }

# MAP_TRACKS = { "default":
#                    [ "%s_%s" % x for x in itertools.product( TISSUES, CONDITIONS ) ] +\
#                    [ "agg_D3", "agg_unstim" ],
#                "master":
#                    [ "%s_%s" % x for x in itertools.product( TISSUES, CONDITIONS ) ] +\
#                    [ "agg_D3", "agg_unstim"  ],
#                "replicates" :
#                    [ "%s_%s_%s" % x for x in itertools.product(  TISSUES, CONDITIONS, REPLICATES ) ],
#                "merged" :
#                    [ "%s_%s" % x for x in itertools.product(  TISSUES, CONDITIONS ) ],
#                }

###########################################################################

def selectTracks( all_tracks, subset ):
    '''select tracks from *all_tracks* according to *subset*.
    '''
    if subset == None:
        return MAP_TRACKS["default"]
    elif "all" in subset: 
        return sorted(all_tracks)

    for key, tracks in MAP_TRACKS.iteritems():
        if key in subset: return tracks
        
    # user specified tracks
    tracks = subset

    all_tracks = set(all_tracks)
    tracks = [ x for x in tracks if x in all_tracks ]

    return tracks

def getReplicates( track ):
    '''return replicates for a track.'''
    

def linkToUCSC( contig, start, end ):
    '''build URL for UCSC.'''

    link = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=%(contig)s:%(start)i..%(end)i>`_" \
          % locals()
    return link

###########################################################################
###########################################################################
###########################################################################
## Trackers
###########################################################################
class CallingTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

class DefaultTracker( CallingTracker ):
    '''Define convenience tracks for plots'''
    # def __init__(self, *args, **kwargs ):
    #     ChipseqTracker.__init__(self, *args, **kwargs)


