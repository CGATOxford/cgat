import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['readqc_exportdir']
DATADIR=P['readqc_datadir']
DATABASE=P['readqc_backend']
PE=P['readqc_pe']

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
    
