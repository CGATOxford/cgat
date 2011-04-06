import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

# get from config file
UCSC_DATABASE="hg19"
EXPORTDIR="export"
REFERENCE="refcoding"

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("conf.py"): 
    execfile("conf.py")

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.sra" % datadir), "%s/(\S+).sra" % datadir) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.gz" % datadir), "%s/(\S+).fastq.gz" % datadir ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.1.gz" % datadir), "%s/(\S+).fastq.1.gz" % datadir ) +\
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
    glob.glob( "%s/*.cuffdiff" % datadir ), 
    "%s/(\S+).cuffdiff" % datadir )

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
class DefaultTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''

    def getTracks( self, subset = None):
        return selectTracks( self.tracks, subset )

class SingleTableTrackerRows( TrackerSQL ):
    '''Tracker to interrogate a single table.

    Rows are unique in *field*.
    '''
    exclude_columns = ()
    table = None
    fields = ("track",)

    @property
    def tracks( self ):
        d = self.get( "SELECT DISTINCT %s FROM %s" % (",".join(self.fields), self.table ))
        if len(self.fields) == 1:
            return tuple( [x[0] for x in d ] )
        else:
            return tuple( [tuple(x) for x in d ] )

    @property
    def slices( self ):
        columns = self.getColumns( self.table )
        return [ x for x in columns if x not in self.exclude_columns and x not in self.fields ]

    def __call__(self, track, slice = None ):
        if len(self.fields) == 1: track = (track,)
        wheres = " AND ".join([ "%s = '%s'" % (x,y) for x,y in zip( self.fields, track ) ] )
        return self.getValue( "SELECT %(slice)s FROM %(table)s WHERE %(wheres)s" )

class SingleTableTrackerColumns( TrackerSQL ):
    exclude_columns = ("track,")
    table = None
    column = None

    @property
    def tracks(self):
        columns = self.getColumns( self.table )
        return [ x for x in columns if x not in self.exclude_columns and x != self.column ]

    def __call__(self, track, slice = None ):
        data = self.getAll( "SELECT %(column)s, %(track)s FROM %(table)s" )
        return data

class SingleTableTrackerHistogram( TrackerSQL ):
    exclude_columns = ("track,")
    table = None
    column = None

    @property
    def tracks(self):
        columns = self.getColumns( self.table )
        return [ x for x in columns if x not in self.exclude_columns and x != self.column ]

    def __call__(self, track, slice = None ):
        data = self.getAll( "SELECT %(column)s, %(track)s FROM %(table)s" )
        return data

