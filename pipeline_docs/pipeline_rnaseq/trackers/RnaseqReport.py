import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

# get from config file
UCSC_DATABASE="hg19"

REFERENCE="refcoding"

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['rnaseq_exportdir']
DATADIR=P['rnaseq_datadir']
DATABASE=P['rnaseq_backend']

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.sra" % DATADIR), "%s/(\S+).sra" % DATADIR) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.gz" % DATADIR), "%s/(\S+).fastq.gz" % DATADIR ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.1.gz" % DATADIR), "%s/(\S+).fastq.1.gz" % DATADIR ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" )

ALL = PipelineTracks.Aggregate( TRACKS )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###########################################################################
## tracks for the gene sets
class GenesetTrack( PipelineTracks.Sample ):
    attributes = ("geneset",)

GENESET_TRACKS = PipelineTracks.Tracks( GenesetTrack ).loadFromDirectory( 
    glob.glob( "%s/*.cuffdiff" % DATADIR ), 
    "%s/(\S+).cuffdiff" % DATADIR )

CUFFDIFF_LEVELS= ("gene", "isoform", "cds", "tss")

###########################################################################
## shorthand
MAP_TRACKS = {
    'default' : EXPERIMENTS,
    'experiments' : EXPERIMENTS,
    'conditions' : CONDITIONS,
    'tissues' : TISSUES,
    'merged' : ALL,
    'geneset-summary': GENESET_TRACKS }

###########################################################################
def selectTracks( subset ):
    '''select tracks from *all_tracks* according to *subset*.
    '''
    if subset == None or subset == "default":
        return MAP_TRACKS["default"]
    elif subset in MAP_TRACKS:
        return MAP_TRACKS[subset]

    return subset

###########################################################################
def splitLocus( locus ):
    if ".." in locus:
        contig, start, end = re.match("(\S+):(\d+)\.\.(\d+)", locus ).groups()
    elif "-" in locus:
        contig, start, end = re.match("(\S+):(\d+)\-(\d+)", locus ).groups()
        
    return contig, int(start), int(end)

def linkToUCSC( contig, start, end ):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link

###########################################################################
class RnaseqTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    
