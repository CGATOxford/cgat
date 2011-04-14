import os, sys, re, types, itertools, glob
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import sqlalchemy

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['chipseq_exportdir']
DATADIR=P['chipseq_datadir']
DATABASE=P['chipseq_backend']

###################################################################
# cf. pipeline_chipseq.py
# This should be automatically gleaned from pipeline_chipseq.py
###################################################################
import PipelineTracks

Sample = PipelineTracks.Sample3

TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x for x in glob.glob( "%s/*_export.txt.gz" % DATADIR) if "input" not in x ],
      "%s/(\S+)_export.txt.gz" % DATADIR )

Sample.setDefault( "asTable" )

ALL = PipelineTracks.Aggregate( TRACKS )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

############################################################################
# The folllowing need to be parameterized in a config file
# TISSUES=["GM00855", "GM00861" ]
# CONDITIONS=["D3", "unstim" ]
# REPLICATES=["R1", "R2" ]
TAG_UNSTIM = "unstim"
UCSC_GENOME = "hg19"

MOTIFS = ( "rxrvdr",
           "GM00855-D3",
           "GM00861-D3" )

###########################################################################
## shorthand
# use list to convert trackers to strings
MAP_TRACKS = {
    'master' : map( str, list(EXPERIMENTS) + list( CONDITIONS )),
    'replicates' : map( str, list(TRACKS) ),
    'default' : map(str, list(EXPERIMENTS)),
    'experiments' : map( str, list(EXPERIMENTS)),
    'conditions' : map( str, list(CONDITIONS)),
    'tissues' : map(str, list(TISSUES)),
    'merged' : map(str, list(EXPERIMENTS)), 
    }

# MAP_TRACKS = { "default":
#                    [ "%s_%s" % x for x in itertools.product( TISSUES, CONDITIONS ) ] +\
#                    [ "agg_D3", "agg_unstim" ],
#                "master":
#                    [ "%s_%s" % x for x in itertools.product( TISSUES, CONDITIONS ) ] +\
#                    [ "agg_D3", "agg_unstim"  ],
#                "replicates" :
#                    [ "%s_%s_%s" % x for x in itertools.product(  TISSUES, CONDITIONS, REPLICATES ) ],
#                "merged" :
#                    [ "%s_%s" % x for x in itertools.product(  TISSUES, CONDITIONS ) ],
#                }

###########################################################################

def selectTracks( all_tracks, subset ):
    '''select tracks from *all_tracks* according to *subset*.
    '''
    if subset == None:
        return MAP_TRACKS["default"]
    elif "all" in subset: 
        return sorted(all_tracks)

    for key, tracks in MAP_TRACKS.iteritems():
        if key in subset: return tracks
        
    # user specified tracks
    tracks = subset

    all_tracks = set(all_tracks)
    tracks = [ x for x in tracks if x in all_tracks ]

    return tracks

def getReplicates( track ):
    '''return replicates for a track.'''
    

def linkToUCSC( contig, start, end ):
    '''build URL for UCSC.'''

    link = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=%(contig)s:%(start)i..%(end)i>`_" \
          % locals()
    return link

###########################################################################
###########################################################################
###########################################################################
## Trackers
###########################################################################
class ChipseqTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )

    def getTracks( self, subset = None):
        if subset:
            for key, tracks in MAP_TRACKS.iteritems():
                if key in subset: return tracks

        return TrackerSQL.getTracks( self )

class DefaultTracker( ChipseqTracker ):
    '''Define convenience tracks for plots'''
    # def __init__(self, *args, **kwargs ):
    #     ChipseqTracker.__init__(self, *args, **kwargs)

class FoldChangeTracker( TrackerSQL ):
    '''the fold change tracker ignores unstimulated tracks.'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    
    def getTracks( self, subset = None ):
        tracks = TrackerSQL.getTracks( self , subset = subset )
        return [ x for x in tracks if TAG_UNSTIM not in x ]

##################################################################################
##################################################################################
##################################################################################
## Annotation of bases with SNPs
##################################################################################
class IntervalsSummary( DefaultTracker ):
    """Summary stats of intervals called by the peak finder.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*), AVG(length), SUM(nprobes)  FROM %(track)s_intervals" % locals() )
        return odict( zip( ("nintervals", "<length>", "nprobes" ), data) )

##################################################################################
##################################################################################
##################################################################################
## Distribution of interval lengths
##################################################################################
class IntervalLengths( DefaultTracker ):
    """Distribution of interval sizes.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT length FROM %(track)s_intervals" % locals() )
        return { "length" : data }

##################################################################################
##################################################################################
##################################################################################
## Distribution of peak values
##################################################################################
class IntervalPeakValues( DefaultTracker ):
    """Distribution of peak values (the number of reads at peak).
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT peakval FROM %(track)s_intervals" % locals() )
        return { "peakval" : data }

##################################################################################
##################################################################################
##################################################################################
## Distribution of peak values
##################################################################################
class EnrichmentOverUnstim( DefaultTracker ):
    """For each peakval, present the fold enrichment of a track
    compared to the Unstim set.

    Useful for diagnosing cutoffs.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):

        cellline, condition, replicate = splitTrack( track )        
        if condition == "Unstim": return []

        if replicate == None: replicate = ""
        other_track = "run" +  cellline + "Unstim" + replicate
        
        data_fg = self.getValues( "SELECT peakval FROM %(track)s_intervals ORDER BY peakval DESC" % locals() )
        data_bg = self.getValues( "SELECT peakval FROM %(other_track)s_intervals ORDER BY peakval DESC" % locals() )

        mi = min( data_bg + data_fg )
        ma = max( data_bg + data_fg )

        n_fg = len(data_fg)
        n_bg = len(data_bg)

        
        hist_fg, bin_edges = numpy.histogram( data_fg, bins = numpy.arange( mi, ma + 1.0, 1.0 ) )
        hist_bg, bin_edges = numpy.histogram( data_bg, bins = numpy.arange( mi, ma + 1.0, 1.0 ) )

        hist_fg = hist_fg[::-1].cumsum()[::-1]
        hist_bg = hist_bg[::-1].cumsum()[::-1]
        
        fold = []
        for fg, bg in zip( hist_fg, hist_bg ):
            fold.append( float(bg) / fg )

        return odict( ( ( "peakval" , bin_edges[:-1]),
                        ( "fold" , fold ) ) )

##################################################################################
##################################################################################
##################################################################################
## Distribution of average values
##################################################################################
class IntervalAverageValues( DefaultTracker ):
    """Distribution of average values (the average number of reads within the interval)
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT avgval FROM %(track)s_intervals" % locals() )
        return { "avgval" : data }

##################################################################################
##################################################################################
##################################################################################
## Distribution of average values
##################################################################################
class IntervalLengthVsAverageValue( DefaultTracker ):
    """Length vs average value.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, avgval FROM %(track)s_intervals" % locals() )
        return odict( zip( ("length", "avgval"), zip(*data) ) )

##################################################################################
##################################################################################
##################################################################################
## Distribution of average values
##################################################################################
class IntervalLengthVsPeakValue( DefaultTracker ):
    """Length vs peak value
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( "SELECT length, peakval FROM %(track)s_intervals" % locals() )
        return odict( zip( ("length", "peakval"), zip(*data) ) )

##################################################################################
##################################################################################
##################################################################################
## peak location within interval
##################################################################################
class PeakLocation( DefaultTracker ):
    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals() )
        data2 = self.getValues( "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
##################################################################################
##################################################################################
## distance of peak to end of interval
##################################################################################
class PeakDistance( DefaultTracker ):
    mPattern = "_intervals$"

    def __call__(self, track, slice = None):
        data1 = self.getValues( "SELECT PeakCenter - start FROM %(track)s_intervals" % locals() )
        data2 = self.getValues( "SELECT end - PeakCenter FROM %(track)s_intervals" % locals() )
        return { "distance" : data1 + data2 }

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class IntervalList( DefaultTracker ):
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
##################################################################################
##################################################################################
## 
##################################################################################
class IntervalListFull( DefaultTracker ):
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
##################################################################################
##################################################################################
## 
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
##################################################################################
##################################################################################
## 
##################################################################################
class IntervalListFoldChange( FoldChangeTracker, IntervalList ):
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
##################################################################################
##################################################################################
## correlations
##################################################################################
class Correlations( ChipseqTracker ):
    """Correlation between all sets.
    """

    pattern = "(.*)_intervals$"
    mSkipColumns = ("id", "contig", "start", "end" )
    
    def __call__(self, track, slice = None ):
        table = "%s_correlation" % self.mField
        data = self.getValues( "SELECT %s FROM %s AS e ORDER BY id" % (track, table))
        return odict( ((self.mField,data),) )

class CorrelationsPeakval( Correlations ):
    mField = "peakval"

class CorrelationsAvgval( Correlations ):
    mField = "avgval"

class CorrelationsLength( Correlations ):
    mField = "length"

##################################################################################
##################################################################################
##################################################################################
## fold change
##################################################################################
class FoldChange( FoldChangeTracker ):
    """return fold changes for all intervals.
    """

    mPattern = "_readcounts$"
    
    def __call__(self, track, slice = None ):
        data = self.getValues( "SELECT CAST(fg_max AS float)/ bg_max FROM %(track)s_readcounts" % locals() )
        return { 'fold' : data }

##################################################################################
##################################################################################
##################################################################################
## fold change
##################################################################################
class FoldChangeCounts( FoldChangeTracker ):
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



