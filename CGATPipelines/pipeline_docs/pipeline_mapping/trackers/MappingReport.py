import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

# get from config file
UCSC_DATABASE="hg19"
ENSEMBL_DATABASE="Homo_sapiens"
RX_ENSEMBL_GENE = re.compile("ENSG")
RX_ENSEMBL_TRANSCRIPT = re.compile("ENST")

REFERENCE="refcoding"

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['mapping_exportdir']
DATADIR=P['mapping_datadir']
DATABASE=P['mapping_backend']

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "%s/*.sra" % DATADIR), "%s/(\S+).sra" % DATADIR) +\
    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.gz" % DATADIR), "%s/(\S+).fastq.gz" % DATADIR ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "%s/*.fastq.1.gz" % DATADIR), "%s/(\S+).fastq.1.gz" % DATADIR ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" )


###########################################################################
## tracks for the gene sets
class GenesetTrack( PipelineTracks.Sample ):
    attributes = ("geneset",)

GENESET_TRACKS = PipelineTracks.Tracks( GenesetTrack ).loadFromDirectory( 
    glob.glob( "%s/*.cuffdiff" % DATADIR ), 
    "%s/(\S+).cuffdiff" % DATADIR )

CUFFDIFF_LEVELS= ("gene", "isoform", "cds", "tss")

###########################################################################
def splitLocus( locus ):
    if ".." in locus:
        contig, start, end = re.match("(\S+):(\d+)\.\.(\d+)", locus ).groups()
    elif "-" in locus:
        contig, start, end = re.match("(\S+):(\d+)\-(\d+)", locus ).groups()
        
    return contig, int(start), int(end)

###########################################################################
def linkToUCSC( contig, start, end ):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i-%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link

###########################################################################
def linkToEnsembl( id ):
    ensembl_database = ENSEMBL_DATABASE
    if RX_ENSEMBL_GENE.match( id ):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_database)s/Gene/Summary?g=%(id)s>`_" \
                % locals()
    elif RX_ENSEMBL_TRANSCRIPT.match( id ):
        link = "`%(id)s <http://www.ensembl.org/%(ensembl_database)s/Transcript/Summary?t=%(id)s>`_" \
                % locals()
    else:
        link = id
    return link

###########################################################################
class MappingTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    
