import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['readqc_exportdir']
DATADIR=P['readqc_datadir']
DATABASE=P['readqc_backend']
PE=P.get('readqc_pe',False)


###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import PipelineTracks

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
        
        filenames = sorted( [x.asFile() for x in TRACKS ] )
        
        for fn in filenames:
            if PE=="True":
                fn1 = fn + ".1"
                fn2 = fn + ".2"
                toc_text.append( "* %(fn1)s_" % locals()) 
                toc_text.append( "* %(fn2)s_" % locals()) 
                link_text.append( ".. _%(fn1)s: %(edir)s/fastqc/%(fn1)s_fastqc/fastqc_report.html" % locals() )
                link_text.append( ".. _%(fn2)s: %(edir)s/fastqc/%(fn2)s_fastqc/fastqc_report.html" % locals() )
            else:
                toc_text.append( "* %(fn)s_" % locals()) 
                link_text.append( ".. _%(fn)s: %(edir)s/fastqc/%(fn)s_fastqc/fastqc_report.html" % locals() )
            
        toc_text = "\n".join(toc_text)
        link_text =  "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

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

        filenames = sorted( [x.asFile() for x in TRACKS ] )

        blocks = ResultBlocks()

        # note there are spaces behind the %(image)s directive to accomodate
        # for path substitution
        block = '''
.. figure:: %(image)s                                     
   :height: 300 
'''

        def _add( fn ):
            image = os.path.abspath(os.path.join( EXPORTDIR, "fastqc", "%s_fastqc" % fn, "Images", "%s.png" % slice ))
            if not os.path.exists( image ): return

            blocks.append( ResultBlock( text = block % locals(),
                                        title = fn ) )
            
        for fn in filenames:

            if PE=="True":
                _add( fn + ".1" )
                _add( fn + ".2" )
            else:
                _add( fn )

        return odict( (("rst", "\n".join( Utils.layoutBlocks( blocks, layout = "columns-2"))),))

