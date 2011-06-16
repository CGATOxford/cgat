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
        data = self.getFirstRow( "SELECT COUNT(*), round(AVG(length),0), round(AVG(nprobes),0)  FROM %(track)s_intervals" % locals() )
        return odict( zip( ("intervals_count", "mean_interval_length", "mean_reads_per_interval" ), data) )

##################################################################################
class IntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT length FROM %(track)s_intervals" % locals() )
        return { "length" : data }

##################################################################################
class IntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT peakval FROM %(track)s_intervals" % locals() )
        return { "peakval" : data }

##################################################################################
class IntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

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

    nresults = 20
    mColumnsFixed = ("pos", "length" )
    mColumnsVariable= ( "peakval", "avgval", "fold" )
    mPattern = "-?_intervals$"

    def getSQLStatement( self, track, slice = None ):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_intervals AS i
                       ORDER BY i.avgval DESC''' % locals()

        if self.nresults:
            statement += " LIMIT %i" % self.nresults

        return statement

    def __call__(self, track, slice = None ):

        statement = self.getSQLStatement( track, slice )
        data = self.get( statement  )
        ucsc_genome = UCSC_GENOME
        n = odict()
        for d in data:
            id, contig, start, end, length = d[:5]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_genome)s&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[str(id)] = odict( zip(self.mColumnsFixed + self.mColumnsVariable, (pos, length,) + d[5:]))
            
        return n

##################################################################################
class IntervalListFull( cpgTracker ):
    '''list of all intervals. Table for export. '''

    nresults = None
    mPattern = "_intervals$"

    def __call__(self, track, slice = None ):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_intervals AS i
                       ORDER BY i.peakval DESC''' % locals()

        data = self.get( statement )
        return odict( zip( ("contig", "start", "end", "peakval", "avgval"),  zip(*data ) ))

##################################################################################
class IntervalListPeakval( IntervalList ):
    '''list of intervals.'''

    def getSQLStatement( self, track, slice = None ):
        nresults = self.nresults
        
        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_intervals AS i
                       ORDER BY i.peakval DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##################################################################################
class IntervalListFoldChange( IntervalList ):
    '''list of intervals.'''

    def getSQLStatement( self, track, slice = None ):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_intervals AS i
                       ORDER BY fold DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##################################################################################
class IntervalListCDS( IntervalList ):
    '''list of intervals overlapping CDS.'''
    nresults = 100
    mColumnsVariable= ( "peakval", "avgval", "fold", "nover_CDS", "pover1_CDS", "pover2_CDS", "closest_id", "gene_id", "gene_name" )

    def getSQLStatement( self, track, slice = None ):
        nresults = self.nresults

        statement = '''SELECT distinct i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold, 
                       a.nover_CDS, a.pover1_CDS, a.pover2_CDS, 
                       tss.closest_id, tr.gene_id, tr.gene_name
                       FROM %(track)s_intervals AS i, %(track)s_annotations AS a,
                       %(track)s_tss AS tss, annotations.transcript_info AS tr
                       WHERE i.interval_id=a.gene_id 
                       AND i.interval_id=tss.gene_id
                       AND tr.transcript_id=tss.closest_id
                       AND a.nover_CDS>0
                       ORDER BY a.pover2_CDS DESC, a.pover1_CDS DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##################################################################################
class FoldChange( cpgTracker ):
    """return fold changes for all intervals. """

    mPattern = "_intervals$"
    
    def __call__(self, track, slice = None ):
        data = self.getValues( "SELECT fold FROM %(track)s_intervals" % locals() )
        return { 'fold' : data }

##################################################################################
class FoldChangeCounts( cpgTracker ):
    """Correlation between all sets. """

    pattern = "(.*)_intervals$"
    mMinFoldChange = 2.0

    def __call__(self, track, slice = None ):
        data = []
        upfold = self.mMinFoldChange
        downfold = 1.0 / upfold
        data.append( ("> %5.2f fold" % upfold, self.getValue( "SELECT COUNT(*) FROM %(track)s_intervals WHERE fold > %(upfold)f " % locals() )) )
        data.append( ("unchanged", self.getValue( "SELECT COUNT(*) FROM %(track)s_intervals WHERE fold between %(downfold)f and %(upfold)f" % locals() )) )
        data.append( ("< %5.2f fold" % downfold, self.getValue( "SELECT COUNT(*) FROM %(track)s_intervals WHERE fold < %(downfold)f " % locals() )) )

        return odict(data)

