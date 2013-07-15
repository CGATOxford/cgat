import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['readqc_exportdir']
DATADIR=P['readqc_datadir']
DATABASE=P['readqc_backend']

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
class ReadqcTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

##############################################################
##############################################################
##############################################################
class TrackerFastQC( ReadqcTracker ):
    tracks = [ "all" ]

    def __call__(self, track, slice = None ):

        edir = EXPORTDIR

        toc_text = []
        link_text = []
        
        tracks = sorted( [x.asFile() for x in TRACKS ] )
        print tracks
        print "hello"
        for track in tracks:
            
              for x, fn in enumerate( glob.glob( os.path.join( EXPORTDIR, "fastqc", "%s*_fastqc" % track ) )):
                y = x + 1
                toc_text.append( "* %(track)s-%(y)i_" % locals()) 
                link_text.append( ".. _%(track)s-%(y)i: %(fn)s/fastqc_report.html" % locals() )
        
        print toc_text
        toc_text = "\n".join(toc_text)
        print toc_text
        print link_text
        link_text =  "\n".join(link_text)
        print link_text
        rst_text = "\n%(toc_text)s\n\n%(link_text)s\n" % locals()

        print rst_text
        return odict( (("text", rst_text),) )

##########################################################
##########################################################
##########################################################
class FilteringSummary( SingleTableTrackerRows ):
    table = "filtering_summary"

##############################################################
##############################################################
##############################################################
class FastQCDetails( ReadqcTracker ):
    tracks = [ "all" ]
    slices = ("duplication_levels",
              "kmer_profiles",
              "per_base_gc_content",
              "per_base_n_content",
              "per_base_quality",
              "per_base_sequence_content",
              "per_sequence_gc_content",
              "per_sequence_quality",
              "sequence_length_distribution" )

    def __call__(self, track, slice = None ):

        # note there are spaces behind the %(image)s directive to accomodate
        # for path substitution
        block = '''
.. figure:: %(image)s                                     
   :height: 300 
'''

        blocks = ResultBlocks()
        tracks = sorted( [x.asFile() for x in TRACKS ] )
        
        for track in tracks:
            
            files = glob.glob( os.path.join( EXPORTDIR, "fastqc", "%s*_fastqc" % track ) )
            for x, fn in enumerate( sorted(files) ):
                y = x + 1
                
                image = os.path.abspath(os.path.join( fn, "Images", "%s.png" % slice ))
                if not os.path.exists( image ): continue

                blocks.append( ResultBlock( text = block % locals(),
                                            title = os.path.basename(fn) ) )

        return odict( (("rst", "\n".join( Utils.layoutBlocks( blocks, layout = "columns-2"))),))


class FastqcSummary( ReadqcTracker ):
    pattern = "(.*)_Basic_Statistics"
    slices = ("File type", "Filename", "Encoding", "Total Sequences", "Sequence Length", "%GC" )
    def __call__(self, track, slice ):
        return self.getAll( "SELECT * FROM %(track)s_Basic_Statistics WHERE measure = '%(slice)s'" )

class ProcessingDetails( ReadqcTracker ):
    '''return summary of the read processing steps.'''
    pattern = "(.*)_processed$"
    def __call__(self, track ):
        return self.getAll( """SELECT pair,input,output,pair, 100.0 * output / input as percent 
                              FROM %(track)s_processed""" )
    
class ProcessingSummary( ReadqcTracker, SingleTableTrackerRows ):
    table = "processing_summary"
