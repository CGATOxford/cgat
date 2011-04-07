import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

from RnaseqReport import *

##############################################################
##############################################################
##############################################################
class TrackerDESeqFit( Tracker ):
    method = "deseq"
    tracks = [ x.asFile() for x in EXPERIMENTS ]
    slices = [ x.asFile() for x in GENESET_TRACKS ]

    def __call__(self, track, slice = None ):
        edir = exportdir
        level = "gene"
        geneset = slice
        method = self.method

        rst_text = '''
%(level)s %(track)s %(geneset)s
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: %(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_fit_%(track)s.png

.. figure:: %(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_residuals_%(track)s.png

''' % locals()

        return odict( (("text", rst_text),) )

##############################################################
##############################################################
##############################################################
class TrackerDESummaryDESeq( SingleTableTrackerRows ):
    table = "deseq_stats"
    fields = ("level", "geneset", "track1", "track2" )

class TrackerDESummaryCuffdiff( SingleTableTrackerRows ):
    table = "cuffdiff_stats"
    fields = ("level", "geneset", "track1", "track2" )

##############################################################
##############################################################
##############################################################
class TrackerDESummaryPlots( Tracker ):

    tracks = [ x.asFile() for x in GENESET_TRACKS ]

    def __call__(self, track, slice = None ):
        edir = exportdir
        geneset = track
        method = self.method

        rst_text = []

        for level in self.levels:
            for fn in (
                "%(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_heatmap.png",
                "%(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_scvplot.png",
                "%(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_pvalue_vs_length.png" ):
                f = fn % locals()
                if not os.path.exists(f): continue
                rst_text.append( ".. figure:: %(f)s" % locals() )
                
        if rst_text:
            rst_text = '''
%(geneset)s
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

''' % locals() + "\n\n".join( rst_text )
        else:
            rst_text = ""

        return odict( (("text", rst_text),) )

class TrackerDESummaryPlotsDESeq(TrackerDESummaryPlots ):
    method = "deseq"
    levels = ("gene",)

class TrackerDESummaryPlotsCuffdiff(TrackerDESummaryPlots ):
    method = "cuffdiff"
    levels = CUFFDIFF_LEVELS

################################################################
################################################################
################################################################
class TrackerDifferentialExpression( Tracker ):

    tracks = [ x.asFile() for x in GENESET_TRACKS ]

    def __call__(self, track, slice = None ):

        edir, method = exportdir, self.method
        rst_text = []
        
        geneset = track
        for level in self.levels:
            for x,y in itertools.combinations( EXPERIMENTS, 2 ):
                filename = "%(edir)s/%(method)s/%(geneset)s_%(method)s_%(level)s_%(x)s_vs_%(y)s_significance.png" % locals()
                if not os.path.exists( filename ): continue

                rst_text.append('''
%(geneset)s %(level)s %(x)s vs %(y)s 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: %(filename)s

''' % locals())
                            
        return odict( (("text", "\n".join( rst_text)),) )

class TrackerDEPairwiseDESeq( TrackerDifferentialExpression ):
    method = "deseq"
    levels = ("gene",)
    
class TrackerDEPairwiseCuffdiff( TrackerDifferentialExpression ):
    method = "cuffdiff"
    levels = CUFFDIFF_LEVELS

                       
