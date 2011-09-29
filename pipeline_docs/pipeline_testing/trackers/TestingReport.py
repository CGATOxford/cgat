import os, sys, re, types, itertools, glob

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['testing_exportdir']
DATADIR=P['testing_datadir']
DATABASE=P['testing_backend']

###########################################################################
class TestingTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

##############################################################
##############################################################
##############################################################
class PipelineStatus( Status ):

    tracks = [ os.path.splitext(x, ".dir")[0] for x in glob.glob( "*.dir" ) ]

    def testCompletion( self, track ):
        '''check if pipeline ran to completion.
        '''
        
        lines = open( track ).readlines()
        
        x = re.search( "# job finished in (\d+) seconds", lines[-1] )
        if not x: 
            return 'FAIL', 0
        else:
            return 'PASS', x.groups()[0]

class X:
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
    


