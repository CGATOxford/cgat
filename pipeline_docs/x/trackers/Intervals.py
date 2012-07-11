import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class IntervalsSummary( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*), AVG(length), SUM(nprobes)  FROM %(track)s_intervals" % locals() )
        return odict( zip( ("nintervals", "<length>", "nprobes" ), data) )

##################################################################################
class IntervalLengths( cpgTracker ):
    """Distribution of interval sizes. """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT length FROM %(track)s_intervals" % locals() )
        return { "length" : data }

##################################################################################
class IntervalPeakValues( cpgTracker ):
    """Distribution of peak values (the number of reads at peak).
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT peakval FROM %(track)s_intervals" % locals() )
        return { "peakval" : data }

##################################################################################
class IntervalAverageValues( cpgTracker ):
    """Distribution of average values (the average number of reads within the interval) """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT avgval FROM %(track)s_intervals" % locals() )
        return { "avgval" : data }

##################################################################################
class IntervalLengthVsAverageValue( cpgTracker ):
    """Length vs average value. """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, avgval FROM %(track)s_intervals" % locals() )
        return odict( zip( ("length", "avgval"), zip(*data) ) )

##################################################################################
class IntervalLengthVsPeakValue( cpgTracker ):
    """Length vs peak value """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, peakval FROM %(track)s_intervals" % locals() )
        return odict( zip( ("length", "peakval"), zip(*data) ) )

##################################################################################
class PeakLocation( cpgTracker ):
    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals() )
        data2 = self.getValues( "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
class PeakDistance( cpgTracker ):
    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT PeakCenter - start FROM %(track)s_intervals" % locals() )
        data2 = self.getValues( "SELECT end - PeakCenter FROM %(track)s_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
class IntervalList( cpgTracker ):
    '''list of intervals.'''

    nresults = 10
    mColumnsFixed = ("pos", "length" )
    mColumnsVariable= ( "peakval", "avgval" )
    mPattern = "_intervals$"

    def getSQLStatement( self, track, slice = None ):

        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, i.peakval, i.avgval
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC''' % locals()

        if self.nresults:
            statement += " LIMIT %i" % self.nresults

        return statement

    def __call__(self, track, slice = None ):

        statement = self.getSQLStatement( track, slice )
        data = self.get( statement  )
        ucsc_genome = UCSC_GENOME
        n = odict()
        for d in data:
            id, contig, start, end = d[:4]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_genome)s&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[str(id)] = odict( zip(self.mColumnsFixed + self.mColumnsVariable, (pos, end-start,) + d[4:]))
            
        return n

##################################################################################
class IntervalListFull( cpgTracker ):
    '''list of all intervals.

    Table for export.
    '''

    nresults = None
    mPattern = "_intervals$"

    def __call__(self, track, slice = None ):

        statement = '''
          SELECT 
               i.contig, i.start, i.end, i.peakval, i.avgval
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC''' % locals()

        data = self.get( statement )
        return odict( zip( 
                ("contig", "start", "end", "peakval", "avgval"),
                zip(*data ) ))

##################################################################################
class IntervalListPeakval( IntervalList ):
    '''list of intervals.'''

    def getSQLStatement( self, track, slice = None ):

        nresults = self.nresults
        
        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, i.peakval, i.avgval, i.length
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC
          LIMIT %(nresults)s''' % locals()

        return statement

##################################################################################
class IntervalListFoldChange( IntervalList ):
    '''list of intervals.'''

    mColumnsVariable= ( "fg_mean", "bg_mean", "fold_mean", "fg_max", "bg_max", "fold_max" )
    mMinFold = 1.5
    mPattern = "_readcounts$"

    def getSQLStatement( self, track, slice = None ):
        nresults = self.nresults
        minfold = self.mMinFold

        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, 
               fg_mean, bg_mean,
               fg_mean / bg_mean AS fold_mean,
               fg_max, bg_max,
               cast(fg_max as float)/ bg_max AS fold_max
          FROM 
               %(track)s_intervals AS i,
               %(track)s_readcounts AS c
          WHERE
               i.interval_id = c.gene_id AND 
               cast(fg_max as float)/ bg_max > %(minfold)f
          ORDER BY fold_max DESC
          LIMIT %(nresults)s''' % locals()

        return statement

##################################################################################
class FoldChange( cpgTracker ):
    """return fold changes for all intervals.
    """

    mPattern = "_readcounts$"
    
    def __call__(self, track, slice = None ):
        data = self.getValues( "SELECT CAST(fg_max AS float)/ bg_max FROM %(track)s_readcounts" % locals() )
        return { 'fold' : data }

##################################################################################
class FoldChangeCounts( cpgTracker ):
    """Correlation between all sets.
    """

    pattern = "(.*)_readcounts$"
    mMinFoldChange = 2.0

    def __call__(self, track, slice = None ):
        data = []
        upfold = self.mMinFoldChange
        downfold = 1.0 / upfold
        data.append( ("> %5.2f fold" % upfold, self.getValue( "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_max AS float)/ bg_max > %(upfold)f " % locals() )) )
        data.append( ("unchanged", self.getValue( "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_max AS float)/ bg_max between %(downfold)f and %(upfold)f" % locals() )) )
        data.append( ("< %5.2f fold" % downfold, self.getValue( "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_max AS float)/ bg_max < %(downfold)f " % locals() )) )

        return odict(data)



