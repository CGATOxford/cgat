import os, sys, re, types, itertools, math
import numpy
import IOTools

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
##################################################################################
##################################################################################
## peak shape
class PeakShapeTracker( Tracker ):
    '''return peakshape data.  '''

    def getTracks(self):
        return [ os.path.basename( x )[:-len(".peakshape.gz")] for x in glob.glob( os.path.join( DATADIR , "liver_vs_testes", "*.liver.testes.merge.peakshape.gz" )) ]

    slices = ["peak_height", "peak_width", "interval_score", "interval_width" ]
    
    scale = 2000

    def __call__(self, track, slice = None):

        fn = os.path.join( DATADIR, "liver_vs_testes/%(track)s.peakshape.gz.matrix_%(slice)s.gz" % locals() )
        if not os.path.exists( fn ): 
            return
        
        x = IOTools.openFile( fn )
        matrix, rownames, colnames = IOTools.readMatrix( x )
        
        nrows = len(rownames)
        if nrows == 0: return
        if nrows > self.scale:
            take = numpy.array( numpy.floor( numpy.arange( 0, nrows, float(nrows + 1) / self.scale ) ), dtype = int )
            rownames = [ rownames[x] for x in take ]
            matrix = matrix[ take ]
            
        return odict( (('matrix', matrix),
                       ('rows', rownames),
                       ('columns', colnames)) )

##################################################################################        
class PeakShapeSummary( Tracker ):
    '''summary information about peak shapes.'''
    pattern = "(.*)_peakshape"

