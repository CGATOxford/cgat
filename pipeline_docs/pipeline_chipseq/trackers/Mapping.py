import os, sys, re, types, itertools, glob
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from ChipseqReport import *

from SphinxReport.Tracker import *

class TrackerReadStats( Tracker ):

    def getFilename( self, track ):
        return os.path.join( DATADIR, track) + ".readstats"

    def getTracks( self, subset = None ):
        tracks = [ os.path.basename(x)[:-len(".readstats")] for x in glob.glob( os.path.join(DATADIR,"*.readstats") ) ]
        return tracks

class MappingDuplicates( TrackerReadStats ):

    def __call__(self, track, slice = None ):

        data = [ x[:-1].split("\t") for x in open(self.getFilename(track) + ".histogram") ]
        data = [ map(int, x) for x in data[1:]]
        return odict( zip( ("duplicates", "counts"), zip(*data)))


class MappingSummary( TrackerReadStats ):

    def __call__(self, track, slice = None ):

        data = [ x[:-1].split("\t") for x in open(self.getFilename(track)) ]
        
        result = odict()

        for x,y in data[1:]:
            y = re.sub( "\(.*\)", "", y)
            y = re.sub( "^+\s\d+\s", "", y)
            result[y] = int(x)

        return result
