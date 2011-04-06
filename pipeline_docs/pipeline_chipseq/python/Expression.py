import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram
import ChipseqReport


from SphinxReport.Tracker import *

# for trackers_derived_sets and trackers_master
if not os.path.exists("conf.py"):
    raise IOError( "could not find conf.py" )

execfile( "conf.py" )

def splitExpressionTrack( track ):
    '''split track from expression experiments into fields cellline, condition, replicate'''
    assert track.startswith( "exp" )
        
    cellline, timepoint, stimulus = track[3:].split("_")
    return cellline, timepoint, stimulus

def mapExpression2Chipseq( track ):
    '''map an expression track to a chip-seq track.'''
    cellline, timepoint, stimulus = splitExpressionTrack( track )
    if cellline in ('GM00855', 'GM00861'):
        cellline = cellline[-3:]
    elif cellline.startswith( 'All' ):
        # use aggregate
        cellline = ""
    elif cellline.startswith( 'LCL' ):
        # use aggregate
        cellline = ""
    else:
        # others: aggregate celllines
        cellline = ""

    # ignore timepoint
    if timepoint == "00": stimulus = "Unstim"
    
    return "run%s%s" % (cellline,stimulus)

def selectTracks( all_tracks, subset ):
    '''select tracks from *all_tracks* according to *subset*.

    master
       tracks with replicates combined. This is the default.
    all
       all tracks
    replicates
       all tracks with replicates
    celllines
       tracks combined by celll lines
    conditions
       tracks combined by condition
    subtract
       tracks build after subtracting Unstim.
    '''

    if subset == None or "master" in subset:
        return sorted(all_tracks)
    elif "conditions" in subset:
        return sorted(all_tracks)
    elif "all" in subset: 
        return sorted(all_tracks)

    for key,tracks in MAP_TRACKS.iteritems():
        if key in subset: return tracks

    # user specified tracks
    tracks = subset

    all_tracks = set(all_tracks)
    tracks = [ x for x in tracks if x in all_tracks ]

    return tracks

class DefaultTracker( TrackerSQL ):
    '''Define convenience tracks for plots.
    '''
    mPattern = "exp.*_levels"
    
    def getTracks( self, subset = None):

        all_tracks = TrackerSQL.getTracks( self )
        return selectTracks( all_tracks, subset )

class ExpressionTracker( TrackerSQL ):
    '''define convenience tracks for differential expression plots.
    '''

    def getTracks( self, subset = None ):

        all_tracks = TrackerSQL.getTracks( self )
        return selectTracks( all_tracks, subset )
    
######################################################################
######################################################################
######################################################################
class Correlations( DefaultTracker ):
    """Correlation between all sets.
    """
    mSkipColumns = ("id", "contig", "start", "end" )
    
    def __call__(self, track, slice = None ):
        field = self.mField
        data = self.getValues( "SELECT %(field)s FROM %(track)s AS e ORDER BY cluster_id" % locals() )
        return odict( ((self.mField,data),))

######################################################################
######################################################################
######################################################################
class CorrelationsMean( Correlations ):
    mField = "mean"


######################################################################
######################################################################
######################################################################
class ExpressionAffymetrixSourceCounts( TrackerSQL ):
    '''returns the numbers of cluster_ids and genes in the
    expression data grouped by Affymetrix source.
    '''

    mPattern = "expression$"

    def __call__(self, track, slice=None ):
        statement = '''SELECT ee.source,
                              COUNT(distinct e.cluster_id),
                              COUNT(distinct e.transcript_id),
                              COUNT(distinct e.gene_id) 
                              FROM %(track)s AS ee
                    LEFT JOIN probeset2transcript as e ON e.cluster_id = ee.cluster_id
                    GROUP BY ee.source''' % locals()

        data = self.get( statement ) 
        
        statement = '''SELECT 'all',
                              COUNT(distinct e.cluster_id),
                              COUNT(distinct e.transcript_id),
                              COUNT(distinct e.gene_id) 
                    FROM %(track)s AS ee 
                    LEFT JOIN probeset2transcript as e ON e.cluster_id = ee.cluster_id''' % locals()

        data.extend( self.get( statement ) )

        result = odict()
        for source, ncluster_ids, ntranscripts, ngenes in data:
            if source == None: source = "na"
            result[source] = odict( ( ('cluster_ids', ncluster_ids), 
                                      ('transcripts', ntranscripts),
                                      ('ngenes', ngenes)) )

        

        return result

######################################################################
######################################################################
######################################################################
class ExpressionSourceCounts( TrackerSQL ):
    '''returns the numbers of cluster_ids and genes in the
    expression data grouped by the ENSEMBL source.
    '''

    mPattern = "probeset2transcript"

    def __call__(self, track, slice=None ):


        statement = '''SELECT e.source,
                              COUNT(distinct e.cluster_id),
                              COUNT(distinct e.transcript_id),
                              COUNT(distinct e.gene_id) 
                    FROM %(track)s AS e
                    GROUP BY e.source''' % locals()

        data = self.get( statement ) 

        statement = '''SELECT 'all',
                              COUNT(distinct e.cluster_id),
                              COUNT(distinct e.transcript_id),
                              COUNT(distinct e.gene_id) 
                    FROM %(track)s AS e''' % locals() 

        data.extend( self.get( statement ) )

        result = odict()
        for source, ncluster_ids, ntranscripts, ngenes in data:
            if source == None: source = "na"
            result[source] = odict( ( ('cluster_ids', ncluster_ids), 
                                      ('transcripts/with', ntranscripts),
                                      ('ngenes/with', ngenes)) )

        statement = '''SELECT t.source,
                              COUNT(distinct t.transcript_id),
                              COUNT(distinct t.gene_id) 
                    FROM transcripts AS t 
                    LEFT JOIN %(track)s AS e ON t.transcript_id = e.transcript_id
                    WHERE e.transcript_id IS NULL
                    GROUP BY t.source''' % locals()

        data = self.get( statement ) 

        statement = '''SELECT 'all',
                              COUNT(distinct t.transcript_id),
                              COUNT(distinct t.gene_id) 
                    FROM transcripts AS t 
                    LEFT JOIN %(track)s AS e ON t.transcript_id = e.transcript_id
                    WHERE e.transcript_id IS NULL
                    ''' % locals()

        data.extend( self.get( statement ) )

        for source, ntranscripts, ngenes in data:
            if source == None: source = "na"
            if source not in result:
                result[source] = odict( ( ('cluster_ids', 0), 
                                      ('transcripts/with', 0),
                                      ('ngenes/with', 0)) )

            result[source]['transcripts/without'] = ntranscripts
            result[source]['genes/without'] = ngenes

        return result

######################################################################
######################################################################
######################################################################
class ExpressionCounts( DefaultTracker ):
    '''returns the numbers of cluster_ids and genes in the
    expression data.
    '''

    def __call__(self, track, slice=None ):

        table = "%s_levels" % (track)
        headers = ('clusters', 'transcripts', 'genes')
        statement = '''SELECT 
                              COUNT(distinct m.cluster_id),
                              COUNT(distinct m.transcript_id),
                              COUNT(distinct m.gene_id) 
                    FROM probeset2transcript AS m, %(table)s AS t
                    WHERE m.cluster_id = t.cluster_id
                    ''' % locals()
            
        data = self.getFirstRow( statement ) 
        result = odict( zip( headers, data) )

        statement = '''SELECT 
                              COUNT(distinct cluster_id),
                              COUNT(distinct transcript_id),
                              COUNT(distinct gene_id) 
                    FROM probeset2transcript'''
            
        data = self.getFirstRow( statement ) 
        
        for x,y in zip(headers,data): result["%s_all" % x] = y

        return result

######################################################################
######################################################################
######################################################################
class ExpressionCountsPerGroup( DefaultTracker ):
    '''returns the numbers of cluster_ids and genes in the
    expression data.
    '''

    def __call__(self, track, slice=None ):

        table = "%s_levels" % (track)
        headers = ('cluster_id', 'clusters', 'transcripts', 'genes')
        statement = '''SELECT 
                              COUNT(distinct m.cluster_id), 
                              COUNT(distinct m.transcript_id),
                              COUNT(distinct m.gene_id) 
                    FROM probeset2transcript AS m, %(table)s AS t
                    WHERE m.cluster_id = t.cluster_id''' % locals()
            
        data = self.getFirstRow( statement ) 
        result = odict( zip( headers, data) )

        statement = '''SELECT 
                              COUNT(distinct cluster_id), 
                              COUNT(distinct transcript_id),
                              COUNT(distinct gene_id) 
                    FROM probeset2transcript'''
            
        data = self.getFirstRow( statement ) 
        
        for x,y in zip(headers,data): result["%s_all" % x] = y

        return result

######################################################################
######################################################################
######################################################################
class ExpressionProbesetsPerGene( DefaultTracker ):
    '''returns the number of probesets per cluster, mrna and genes.
    '''

    mPattern = "expression$"
    def getSlices(self, subset = None ):
        return ("cluster", "transcript", "gene")

    def __call__(self,track, slice=None ):

        if slice == "cluster":
            statement = '''
                    SELECT COUNT(distinct e.probeset)
                    FROM probeset2transcript AS e 
                    GROUP BY e.cluster_id''' % locals()
        elif slice == "transcript":
            statement = '''
                    SELECT COUNT(distinct e.probeset)
                    FROM probeset2transcript AS e
                    GROUP BY e.transcript_id''' % locals()
        elif slice == "gene":
            statement = '''
                    SELECT COUNT(distinct e.probeset)
                    FROM probeset2transcript AS e
                    GROUP BY e.gene_id''' % locals()

        return odict( ((slice,
                        self.getValues(statement)), ) )

######################################################################
######################################################################
######################################################################
class ExpressionReplicateCorrelation( DefaultTracker ):
    '''returns the correlation between replicates.
    '''
    
    def __call__(self, track, slice=None ):

        table = "expression_correlation"
        if track.endswith("_levels"): track=track[:-len("_levels")]
                          
        headers = ('replicate1', 'replicate2', 'coeff', 'observations', 'pvalue', 'significance', 'alternative', 'method')
        data = self.get('''SELECT replicate1, replicate2, coeff, observations, pvalue, significance, alternative, method
        FROM %(table)s WHERE track = "%(track)s" ''' % locals() )

        result = odict()
        for x,y in enumerate(data):
            result[str(x)] = odict( zip( headers, y ))
        return result

######################################################################
######################################################################
######################################################################
class ExpressionFullCorrelation( DefaultTracker ):
    '''returns the correlation between replicates.
    '''

    tablename = "expression_full_correlation"

    def getSubset( self, subset = None):

        all_tracks = self.get( "SELECT DISTINCT track1, replicate1 FROM %s ORDER BY track1,replicate1" % self.tablename)
        if not subset: 
            return ["-".join( x) for x in all_tracks ]
        else:
            subset = set(subset)
            return ["-".join( x) for x in all_tracks if x[0] in subset ]

    def getTracks( self, subset = None):
        return self.getSubset( subset )

    def __call__(self, track, slice=None ):

        table = self.tablename
        xtrack, replicate = track.split("-")

        data = self.get('''SELECT track2, replicate2, coeff
        FROM %(table)s WHERE track1 = "%(xtrack)s" AND replicate1 = "%(replicate)s" ''' % locals() )

        data.extend( self.get('''SELECT track1, replicate1, coeff
        FROM %(table)s WHERE track2 = "%(xtrack)s" AND replicate2 = "%(replicate)s" ''' % locals() ) )

        data = [ ("-".join((x[0],x[1])), x[2]) for x in data]
        data.append( (track,1.0) )

        subset = set(self.dispatcher.tracks)
        data = [ x for x in data if x[0] in subset ] 

        result = odict()
        return odict( sorted(data))

######################################################################
######################################################################
######################################################################
class ExpressionDifferencesSummaryPerGroup( ExpressionTracker ):
    '''return the number of genes called to be differentially expressed.
    The counts are returned per group.
    '''
    def getSlices( self, subset = None ):
        return self.getValues( "SELECT DISTINCT source FROM probeset2transcript" )

    def __call__(self, track, slice=None ):

        result = odict()
        foreground = re.sub("_vs_.*", "", track)

        cutoff = "a.%s %s" % (self.mCutoffField, self.mCutoffThreshold)
        result["tested"] = self.getValue("SELECT COUNT(distinct gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE e.cluster_id = a.cluster_id AND e.source = '%(slice)s'" % locals() )
        result["called"] = self.getValue("SELECT COUNT(distinct gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE %(cutoff)s AND e.cluster_id = a.cluster_id AND e.source = '%(slice)s'" % locals() )
        result["upregulated"] = self.getValue("SELECT COUNT(distinct gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE %(cutoff)s AND fold > 1 AND e.cluster_id = a.cluster_id AND e.source = '%(slice)s'" % locals() )
        result["downregulated"] = self.getValue("SELECT COUNT(distinct gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE %(cutoff)s AND fold < 1 AND e.cluster_id = a.cluster_id AND e.source = '%(slice)s'" % locals() )
        return result

######################################################################
######################################################################
######################################################################
class ExpressionDifferencesCalibrationTTest1( ExpressionTracker ):
    '''calibrate P-Value'''
    mPattern = "exp(\S+)_vs_exp(\S+)_ttest$"
    mAsTables = True

    def __call__(self, track, slice=None ):

        tracks = self.getTracks( subset = None )

        result = odict()

        exponents = []
        for exponent in xrange(0,10):
            pvalue = math.pow( 10, -exponent)
            total = self.getValue( '''SELECT COUNT(*) FROM %(track)s AS a WHERE a.pvalue < %(pvalue)s''' % locals() )
            if total == 0: break
            exponents.append( (pvalue, total, exponent) )

        for other in tracks:
            if other == track: continue
            pvalues, overlaps = [],[]
            for pvalue, total, exponent in exponents:
                val = self.getValue( '''SELECT COUNT(*) 
                           FROM %(track)s AS a, %(other)s AS b
                           WHERE a.cluster_id = b.cluster_id AND 
                                 a.pvalue < %(pvalue)s AND b.pvalue < %(pvalue)s ''' % locals())
                pvalues.append( exponent )
                overlaps.append( 100.0 * val / total )
                
            data = odict( ( ('-log(pvalue)', pvalues), ('overlap', overlaps)) )
            result[other] = data
        return result

######################################################################
######################################################################
######################################################################
class ExpressionDifferencesCalibrationTTest2( ExpressionTracker ):
    '''calibrate P-Value.
    
    The idea was to look at the overlap of the three tests D3, EST and D3EST
    to get an idea of the noise.
    '''
    mPattern = "exp(\S+)_vs_exp(\S+)_ttest$"
    mAsTables = True

    def getTracks( self, subset = None ):
        return ChipseqReport.CELLLINES + ("",)

    def __call__(self, track, slice=None ):

        tablename_ref = "exp%(track)sD3_vs_exp%(track)sUnstim_ttest" % locals()
        tablename_fg = "exp%(track)sD3Est_vs_exp%(track)sUnstim_ttest" % locals()
        tablename_bg = "exp%(track)sEst_vs_exp%(track)sUnstim_ttest" % locals()

        counts = []

        for exponent in xrange(0,10):
            pvalue = math.pow( 10, -exponent)
            total = self.getValue( '''SELECT COUNT(*) FROM %(tablename_ref)s AS a WHERE a.pvalue < %(pvalue)s''' % locals() )
            c = []
            for tablename in (tablename_fg, tablename_bg):
                val = self.getValue( '''SELECT COUNT(*) 
                           FROM %(tablename_ref)s AS a, %(tablename)s AS b
                           WHERE a.cluster_id = b.cluster_id AND 
                                 a.pvalue < %(pvalue)s AND b.pvalue < %(pvalue)s ''' % locals())
                c.append( val )
            counts.append( [exponent,total]+c )
            
        data = []
        for exponent, ref, fg, bg in counts:
            if ref > 0:
                b = 100.0 * (ref - bg) / ref
                f = 100.0 * (ref - fg) / ref 
                data.append( (exponent, f, b) )
        
        return odict( zip( ( "-10log(P)", "fg", "bg"), zip(*data) ))

######################################################################
######################################################################
######################################################################
class PairTrackerSAM( ExpressionTracker ):
    '''Parameters for SAM filtering by qvalue.
    '''
    mPattern = "exp(\S+)_vs_exp(\S+)_sam$"
    mTablePattern = "exp%(track1)s_vs_exp%(track2)s_sam"
    mAsTables = True
    mCutoffField = "called"
    # for some reason the actual qvalue output by siggenes is higher than
    # the qvalue threshold chosen. I artifically bind them by the chosen
    # threshold. Hence use a threshold that is slightly larger.

    # called is a 0/1 flag
    mCutoffThreshold = "> 0.5"
    mWithResultsDir = True

    def called( self, val ): return val > 0.5

class PairTrackerAll( ExpressionTracker ):
    '''Parameters for no filtering. It uses the ttest tables.
    '''
    mPattern = "exp(\S+)_vs_exp(\S+)_ttest$"
    mTablePattern = "exp%(track1)s_vs_exp%(track2)s_ttest"
    mAsTables = True
    mCutoffField = "pvalue"
    mCutoffThreshold = "< 2.0"
    mWithResultsDir = False

    def called( self, val ): return True

class PairTrackerTTest( ExpressionTracker ):
    '''Parameters for ttest filtering by pvalue
    '''
    mPattern = "exp(\S+)_vs_exp(\S+)_ttest$"
    mTablePattern = "exp%(track1)s_vs_exp%(track2)s_ttest"
    mAsTables = True
    mCutoffField = "pvalue"
    mCutoffThreshold = "<= 0.01"
    mWithResultsDir = False

    def called( self, val ): return val <= 0.01

class PairTrackerTTestQValue( ExpressionTracker ):
    '''Parameters for ttest filtering by qvalue.
    '''

    mPattern = "exp(\S+)_vs_exp(\S+)_ttest$"
    mTablePattern = "exp%(track1)s_vs_exp%(track2)s_ttest"
    mAsTables = True
    mCutoffField = "qvalue"
    mCutoffThreshold = "<= 0.05"
    mWithResultsDir = False

    def called( self, val ): return val <= 0.05

######################################################################
######################################################################
######################################################################
class ExpressionDifferencesSummary( ExpressionTracker ):
    '''return the number of genes that are differentially expressed.'''

    def getSlices( self, subset = None ):
        if subset == "all":
            return self.getValues( "SELECT DISTINCT source FROM probeset2transcript" )
        else:
            return ("protein_coding",)

    def __call__(self, track, slice=None ):

        result = odict()
        foreground = re.sub("_vs_.*", "", track)

        cutoff = "a.%s %s" % (self.mCutoffField, self.mCutoffThreshold)
        result["tested"] = self.getValue("SELECT COUNT(DISTINCT e.gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE a.cluster_id = e.cluster_id AND source = '%(slice)s'" % locals() )
        result["called"] = self.getValue("SELECT COUNT(DISTINCT e.gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE a.cluster_id = e.cluster_id AND %(cutoff)s AND source = '%(slice)s'" % locals() )
        result["upregulated"] = self.getValue("SELECT COUNT(DISTINCT e.gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE a.cluster_id = e.cluster_id AND %(cutoff)s AND fold > 1 AND source = '%(slice)s'" % locals() )
        result["downregulated"] = self.getValue("SELECT COUNT(DISTINCT e.gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE a.cluster_id = e.cluster_id AND %(cutoff)s AND fold < 1 AND source = '%(slice)s'" % locals() )
        if self.mWithResultsDir:
            resultdir = os.path.join( exportdir, "SAM" )
            if os.path.exists( resultdir ): 
                result["pdf"] = "`pdf <%(resultdir)s/%(foreground)s.sam.expdiffsam.pdf>`_" % locals()
        return result

######################################################################
######################################################################
######################################################################
class ExpressionDifferencesOverlap( ExpressionTracker ):
    '''overlap between sets of differently expressed genes.'''

    def __call__(self, track, slice=None ):

        tracks = self.getTracks( subset = None )

        cutoff1 = "a.%s %s" % (self.mCutoffField, self.mCutoffThreshold)
        cutoff2 = "b.%s %s" % (self.mCutoffField, self.mCutoffThreshold)
        result = odict()
        for other in tracks:
            result[other] = \
                self.getValue( '''SELECT COUNT(*) FROM %(track)s AS a, %(other)s AS b
                           WHERE a.cluster_id = b.cluster_id AND 
                                 %(cutoff1)s AND %(cutoff2)s''' % locals() )
            
        return result

######################################################################
######################################################################
######################################################################
class DifferentiallyExpressedGenesFoldChange( ExpressionTracker ):

    mField = "fold"

    def __call__(self, track, slice = None ):
        field = self.mField
        return odict( ((self.mField, self.getValues("SELECT %(field)s FROM %(track)s" % locals())),))

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class DifferentiallyExpressedGenesFoldChange95( DifferentiallyExpressedGenesFoldChange ):
    mField = "fold95upper"

######################################################################
######################################################################
######################################################################
class DifferentiallyExpressedGenesFoldChangeCategories( ExpressionTracker ):
    '''return the number of cluster_ids that are up-regulated, down-regulated
    and intermediate depending on the fold change.'''
    
    mField = "fold"
    mThreshold = 2
    
    def __call__(self, track, slice = None ):
        field = self.mField
        
        r = odict()
        t1 = 1.0 / self.mThreshold 
        t2 = self.mThreshold 

        statement = "SELECT COUNT(DISTINCT e.gene_id) FROM %(track)s AS a, probeset2transcript AS e WHERE e.cluster_id = a.cluster_id AND %%s" % locals()

        r[">%s fold down" % self.mThreshold] = self.getValue( statement % ("fold < %(t1)f" % locals()) )
        r["intermediate"] = self.getValue( statement % ("%(t1)f and %(t2)f" % locals() ))
        r[">%s fold up" % self.mThreshold] = self.getValue( statement % ("fold > %(t2)f" % locals()) )
        
        return r

######################################################################
######################################################################
######################################################################
class DifferentiallyExpressedGenesDistanceToInterval( ExpressionTracker ):
    '''compute histogram of distance to closest interval.

    Note that annotator_distance.py will provide a statistical
    significance of an enrichment.

    The slices are

    responsive
       gene shows significant differential expression
    unresponsive
       gene shows no significant differential expression
    upregulated
       gene shows significant differential expression and is upregulated 
    downregulated
       gene shows significant differential expression and is upregulated 
    all
       all genes
    
    Up/down regulation means fold change is larger/smaller than 1.

    The measure of significance depends on the test used.
    '''

    min_distance = 0
    
    def getSlices(self, subset = None ):
        return ["responsive", "unresponsive", "up", "down", "all" ]

    def getNames( self, track ):

        track1, track2 = re.search( self.mPattern, track).groups()
        master = mapExpression2Chipseq( "exp" + track1 )
        return track1, track2, master

    def __call__(self, track, slice = None ):

        track1, track2, master = self.getNames( track )
        
        threshold = self.mCutoffThreshold
        field = self.mCutoffField

        statement = '''
        SELECT min( abs(a.dists) )
        FROM %(track)s as d, probeset2transcript as e, %(master)s_assoc AS a, %(master)s_readcounts AS c
        WHERE 
            e.cluster_id = d.cluster_id AND 
            e.transcript_id = a.transcript_id AND 
            c.gene_id = a.id AND 
            %%s
            GROUP BY c.gene_id''' % locals()

        statement = '''SELECT MIN(distance) 
        FROM probeset2transcript as e,
             %(track)s as d,
             %(master)s_distance AS a
        WHERE 
            e.cluster_id = d.cluster_id AND 
            e.transcript_id = a.transcript_id AND 
            a.distance != '' AND
            a.distance >= %(min_distance)i AND 
             %%s
            GROUP BY e.gene_id''' % self.members( locals() )

        if slice == "responsive":
            s = statement % ( "%(field)s %(threshold)s" % locals() )
        elif slice == "unresponsive":
            s = statement % ( "not( %(field)s %(threshold)s )" % locals() ) 
        elif slice == "up":
            s = statement % ( "%(field)s %(threshold)s AND fold >= 1" % locals() )
        elif slice == "down":
            s = statement % ( "%(field)s %(threshold)s AND fold <= 1" % locals() )
        elif slice == "all":
            s = statement % ( "1" % locals() )

        data = self.getValues( s )
        if len(data) == 0: return []
        return odict( (("distance", data),) )

######################################################################
######################################################################
######################################################################
class DifferentiallyExpressedGenesDistanceToIntervalDistal( DifferentiallyExpressedGenesDistanceToInterval ):
    '''compute histogram of distance to closest interval.

    the distances here are filtered excluding those greater than min_distance.
    '''
    
    min_distance = 5000

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class GeneList( ExpressionTracker ):
    '''a table with candidate genes.
    
    '''
    mColumnsFixed = ( "gene_name", "pos", )
    mColumnsVariable= ()
    mMaxEntries = None

    def getNames( self, track ):

        track1, track2 = re.search( self.mPattern, track).groups()
        master = mapExpression2Chipseq( "exp" + track1 )
        return track1, track2, master

    def __call__(self, track, slice = None ):

        statement = self.getSQLStatement( track, slice )
        
        if self.mMaxEntries: statement += " LIMIT %i" % self.mMaxEntries

        data = self.get( statement  )
        
        n = odict()
        for d in data:
            gene_id, gene_name, contig, start, end = d[:5]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[gene_id] = odict( zip (self.mColumnsFixed + self.mColumnsVariable, (gene_name, pos,) + d[5:] ) )
            
        return n
    
##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class CandidateGenesMotifs( GeneList ):
    '''a table with candidate genes.
    
    List all transcripts with
       * differential expression (P-Value < :attr:`mPValue`)
       * associated intervals with Vitamin D3 motif (P-Value < attr:`mPValue`).
       * only consider intervals after subtracting Unstim

    The columns are:

    transcript_id
       the transcript id
    gene_id
       the gene id
    pos
       the genomic location of the transcript
    nprobes
       the number of probes associated with this transcript
    nintervals
       the number of intervals in the neighbourhood of the probes
    exp-pvalue
       the minimum P-Value among probes in the test for differential expression.
    motif-pvalue
       the minimum P-Value among intervals in the neighborhood containing one of the canonical motifs.
  
    '''
    mColumnsFixed = ( "pos", )
    mColumnsVariable= ( "nprobes", "nintervals", "exp-pvalue", "motif-pvalue", "description" )
    mMotif = "rxrvdr"
    mEValueThreshold = 0.1

    # flag set if pair needs to inverted

    def getSQLStatement( self, track, slice = None ):
        motif = self.mMotif
        track1, track2, master = self.getNames( track )
        
        threshold = self.mCutoffThreshold
        evalue_threshold = self.mEValueThreshold
        field = self.mCutoffField

        # cross join necessary to use all indices        
        statement = '''
          SELECT 
               e.gene_id, ii.gene_name,
               e.contig, e.start, e.end,
               COUNT(DISTINCT e.cluster_id), 
               COUNT(DISTINCT a.id),
               min(d.%(field)s) AS p1, 
               min(m.evalue) AS p2,
               ee.description
          FROM %(master)s_assoc as a CROSS JOIN
               %(master)s_mast as m,
               %(track)s AS d,
               probeset2transcript AS e,
               gene_info AS ii,
               expression as ee
          WHERE 
                i.gene_id = e.gene_id AND
                ii.gene_id = e.gene_id AND 
                a.transcript_id = e.transcript_id AND 
                e.cluster_id = d.cluster_id AND
                a.id = m.id AND 
                e.cluster_id = ee.cluster_id AND 
                m.motif = '%(motif)s'
          GROUP BY e.gene_id 
          HAVING p1 %(threshold)s and p2 <= %(evalue_threshold)f''' % locals()

        return statement
        
##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class CandidateGenesPeakval( CandidateGenesMotifs ):
    mColumnsVariable= ( "nprobes", "nintervals", "exp-pvalue", "min-peakval", "max-peakval", "description" )

    def getSQLStatement( self, track, slice = None ):

        track1, track2, master = self.getNames( track )

        threshold = self.mCutoffThreshold
        field = self.mCutoffField

        peakvals = self.getValues( "SELECT peakval FROM %(master)s_intervals" % (locals() ))
        minval = scipy.stats.scoreatpercentile(peakvals, 75)

        # cross join necessary to use all indices
        statement = '''
          SELECT 
               e.gene_id, ii.gene_name,
               e.contig, min(e.start), max(e.end),
               COUNT(DISTINCT e.cluster_id), 
               COUNT(DISTINCT i.interval_id),
               MIN(d.%(field)s) AS p1,
               MIN(i.peakval) AS minp,
               MAX(i.peakval) AS maxp,
               ee.description
          FROM %(master)s_assoc as a CROSS JOIN
               %(master)s_intervals AS i,
               probeset2transcript AS e,
               expression AS ee,
               gene_info AS ii,
               %(track)s AS d
          WHERE 
                a.transcript_id = e.transcript_id AND
                ii.gene_id = e.gene_id AND 
                e.cluster_id = d.cluster_id AND 
                e.cluster_id = ee.cluster_id AND
                i.interval_id = a.id AND
                i.peakval >= %(minval)i
          GROUP BY e.gene_id having p1 %(threshold)s''' % locals()

        return statement

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class DifferentiallyExpressedGenes( GeneList ):
    '''output a list of differentially expressed genes.'''
    mColumnsVariable= ( "nprobes", "min-fold", "max-fold", "exp-pvalue", "description" )
    mMaxEntries = None

    def getSQLStatement( self, track, slice = None ):    

        track1, track2, master = self.getNames( track )

        threshold = self.mCutoffThreshold
        field = self.mCutoffField

        peakvals = self.getValues( "SELECT peakval FROM %(master)s_intervals" % (locals() ))
        minval = scipy.stats.scoreatpercentile(peakvals, 75)

        # cross join necessary to use all indices
        statement = '''
          SELECT 
               e.gene_id, ii.gene_name, 
               e.contig, min(e.start), max(e.end),
               COUNT(DISTINCT e.cluster_id),
               MIN(d.fold),
               MAX(d.fold),
               MIN(d.%(field)s) AS p1,
               ee.description
          FROM probeset2transcript AS e,
               %(track)s AS d,
               gene_info AS ii,
               expression AS ee
          WHERE 
                e.cluster_id = d.cluster_id AND 
                e.cluster_id = ee.cluster_id AND
                ii.gene_id = e.gene_id 
          GROUP BY e.gene_id having p1 %(threshold)s
          ORDER by p1''' % locals()
        
        return statement

##################################################################################
##################################################################################
##################################################################################
## counts 
##################################################################################
class ExpressionIntervalCounts( ExpressionTracker ):
    '''return counts of genes that:

    1. do or do not overlap with an interval at the TSS
    2. are responsive or are not responsive 

    Also, those within 5k of the TSS are further stratified by
    whether or not they contain a motif.
    '''

    mMaxClosestDist = 5000

    def getSlices( self, subset = None ):
        if subset: return subset
        return ChipseqReport.MOTIFS

    def getNames( self, track ):

        track1, track2 = re.search( self.mPattern, track).groups()
        master = mapExpression2Chipseq( "exp" + track1 )
        return track1, track2, master

    def __call__(self, track, slice = None ):

        motif = slice
        track1, track2, master = self.getNames( track )

        threshold = self.mCutoffThreshold
        field = self.mCutoffField
        
        # cross join necessary to use all indices
        statement = '''
        SELECT 
               e.gene_id,
               MIN(tss.closest_dist),
               COUNT(DISTINCT e.cluster_id), 
               COUNT(DISTINCT i.interval_id),
               MIN(d.%(field)s) AS p1,
               MIN(i.peakval) AS minp,
               MAX(i.peakval) AS maxp,
               ee.description,
               MAX(mast.nmatches)
        FROM 
               probeset2transcript AS e
               LEFT JOIN %(track)s AS d ON d.cluster_id = e.cluster_id
               LEFT JOIN %(master)s_tss as tss ON tss.closest_id  = e.gene_id 
               LEFT JOIN %(master)s_intervals AS i ON i.interval_id = tss.gene_id
               LEFT JOIN expression AS ee ON ee.cluster_id = e.cluster_id
               LEFT JOIN %(master)s_mast AS mast on mast.id = i.interval_id AND mast.motif = '%(motif)s' 
        WHERE e.source = 'protein_coding' 
        GROUP BY e.gene_id''' % locals()
               
        data = self.get( statement  )

        gene_ids = []

        xlabels = ("na", "responsive", "unresponsive")
        ylabels = ("overlap with motif", "overlap without motif", "no overlap" )
        counts = numpy.zeros( (len(xlabels),len(ylabels)), numpy.int )
        
        for d in data:
            gid, closest_dist, nprobeset, nids, pvalue, minpeakval, maxpeakval, description, nmatches = d

            if pvalue == None: x = 0
            elif self.called( pvalue ): x = 1
            else: x = 2

            if closest_dist != None and closest_dist <= self.mMaxClosestDist:
                if nmatches > 0:
                    y = 0
                else: y = 1
            else: 
                y = 2

            counts[x,y] += 1

        result = odict()
        for x,xlabel in enumerate(xlabels):
            for y,ylabel in enumerate(ylabels):
                result[ "%s/%s" % (xlabel,ylabel) ] = counts[x,y]
        return result

##################################################################################
##################################################################################
##################################################################################
## counts 
##################################################################################
class ExpressionNumMotifs( ExpressionIntervalCounts ):
    '''return the number of motifs within intervals that
    are close to responsive or non-responsive genes.
    
    The tracker returns the following values for each gene:

    nintervals
       number of intervals
    nmotifs
       number of intervals with motifs
    nmatches
       number of matches to motif over all intervals

    '''

    def __call__(self, track, slice = None ):

        motif = slice
        track1, track2, master = self.getNames( track )

        threshold = self.mCutoffThreshold
        field = self.mCutoffField
        max_dist = self.mMaxClosestDist

        statement = '''
        SELECT
               (SELECT COUNT(*) 
                       FROM %(master)s_tss as tss 
                       WHERE tss.closest_id = e.gene_id AND 
                             tss.closest_dist < %(max_dist)i) 
               AS nintervals,
               (SELECT CASE WHEN COUNT(*) THEN COUNT(*) ELSE 0 END 
                       FROM %(master)s_tss as tss, 
                            %(master)s_mast as mast
                       WHERE tss.closest_id = e.gene_id AND 
                             tss.closest_dist < %(max_dist)i AND 
                             mast.id = tss.gene_id AND 
                             mast.motif = '%(motif)s' )
               AS nmotifs,
               (SELECT CASE WHEN SUM(mast.nmatches) THEN SUM(mast.nmatches) ELSE 0 END 
                       FROM %(master)s_tss as tss, 
                            %(master)s_mast as mast
                       WHERE tss.closest_id = e.gene_id AND 
                             tss.closest_dist < %(max_dist)i AND 
                             mast.id = tss.gene_id AND 
                             mast.motif = '%(motif)s' ) 
               AS nmatches
        FROM 
               probeset2transcript AS e,
               %(track)s AS d ON d.cluster_id = e.cluster_id
        WHERE e.source = 'protein_coding'
              AND %(where)s
        GROUP BY e.gene_id
        '''

        where = "d.%(field)s %(threshold)s" % locals()
        responsive_statement = statement % locals()
        where = "not( d.%(field)s %(threshold)s )" % locals()
        unresponsive_statement = statement % locals()

        def _getData( statement ):
            data = self.get( statement )
            return odict( zip( ( "nintervals", 
                                 "nmotifs", 
                                 "nmatches" ),
                               zip(*data)))
            
        return odict( ( ( "responsive", 
                          _getData(responsive_statement) ),
                        ( "non_responsive", 
                          _getData(unresponsive_statement) ) ) )
    


##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class CandidateGenesOverlap( GeneList ):
    '''return a list of candidate genes. These are genes that are
    responsive to expression stimulation AND contain an interval
    within 5kb of the transcription start site.'''

    mColumnsVariable= ( "nprobes", "nintervals", "exp-pvalue", "min-peakval", "max-peakval", "description" )
    mMaxClosestDist = 5000

    def getSQLStatement( self, track, slice = None ):
        track1, track2, master = self.getNames( track )

        threshold = self.mCutoffThreshold
        field = self.mCutoffField
        closest_dist = self.mMaxClosestDist

        # cross join necessary to use all indices
        statement = '''
          SELECT 
               e.gene_id, e.contig, min(e.start), max(e.end),
               COUNT(DISTINCT e.cluster_id), 
               COUNT(DISTINCT i.interval_id),
               MIN(d.%(field)s) AS p1,
               MIN(i.peakval) AS minp,
               MAX(i.peakval) AS maxp,
               ee.description
          FROM 
               probeset2transcript AS e CROSS JOIN
               %(track)s AS d ON d.cluster_id = e.cluster_id,
               %(master)s_tss as tss ON tss.closest_id  = e.gene_id,
               %(master)s_intervals AS i ON i.interval_id = tss.gene_id,
               expression AS ee ON ee.cluster_id = e.cluster_id
          WHERE e.source = 'protein_coding' AND
                tss.closest_dist < %(closest_dist)i
          GROUP BY e.gene_id
          HAVING p1 %(threshold)s
          ''' % locals()
        
        return statement

###########################################################################
###########################################################################
###########################################################################
## setup the trackers for each method
###########################################################################
class SummarySAM( ExpressionDifferencesSummary, PairTrackerSAM ): pass
class SummaryTTest( ExpressionDifferencesSummary, PairTrackerTTest ): pass
class SummaryTTestQValue( ExpressionDifferencesSummary, PairTrackerTTestQValue ): pass
class OverlapSAM( ExpressionDifferencesOverlap, PairTrackerSAM ): pass
class OverlapTTest( ExpressionDifferencesOverlap, PairTrackerTTest ): pass
class FoldChangeAll( DifferentiallyExpressedGenesFoldChange, PairTrackerAll ): pass
class FoldChangeCategories( DifferentiallyExpressedGenesFoldChangeCategories, PairTrackerAll ): pass
class FoldChangeSAM( DifferentiallyExpressedGenesFoldChange, PairTrackerSAM ): pass
class FoldChangeTTest( DifferentiallyExpressedGenesFoldChange, PairTrackerTTest ): pass
class FoldChange95CISAM( DifferentiallyExpressedGenesFoldChange, PairTrackerSAM ): pass
class FoldChange95CITTest( DifferentiallyExpressedGenesFoldChange, PairTrackerTTest ): pass
class DistanceToIntervalSAM( DifferentiallyExpressedGenesDistanceToInterval, PairTrackerSAM ): pass
class DistanceToIntervalDistalSAM( DifferentiallyExpressedGenesDistanceToIntervalDistal, PairTrackerSAM ): pass
class DistanceToIntervalTTest( DifferentiallyExpressedGenesDistanceToInterval, PairTrackerTTest ): pass
class CandidateGenesMotifsTTest( CandidateGenesMotifs, PairTrackerTTest ): pass
class CandidateGenesMotifsSAM( CandidateGenesMotifs, PairTrackerSAM ): pass
class CandidateGenesPeakvalTTest( CandidateGenesPeakval, PairTrackerTTest ): pass
class CandidateGenesPeakvalSAM( CandidateGenesPeakval, PairTrackerSAM ): pass
class DifferentiallyExpressedGenesSAM( DifferentiallyExpressedGenes, PairTrackerSAM ): pass
class DifferentiallyExpressedGenesTTest( DifferentiallyExpressedGenes, PairTrackerTTest ): pass
class CandidateGenesPeakvalSAM( CandidateGenesPeakval, PairTrackerSAM ): pass
class SummaryPerGroupSAM( ExpressionDifferencesSummaryPerGroup, PairTrackerSAM ): pass
class SummaryPerGroupTTest( ExpressionDifferencesSummaryPerGroup, PairTrackerTTest ): pass
class ExpressionIntervalCountsTTest ( ExpressionIntervalCounts, PairTrackerTTest ): pass
class ExpressionIntervalCountsSAM ( ExpressionIntervalCounts, PairTrackerSAM ): pass
class CandidateGenesOverlapTTest ( CandidateGenesOverlap, PairTrackerTTest ): pass
class CandidateGenesOverlapSAM ( CandidateGenesOverlap, PairTrackerSAM ): pass
class ExpressionNumMotifsSAM ( ExpressionNumMotifs, PairTrackerSAM ): pass

##################################################################################
##################################################################################
##################################################################################
## generic implementations of GO results
##################################################################################
class GO(object):
    mPattern = "exp.*_vs.*_go$"
    suffix = "go"
    
class GOSlim(object):
    mPattern = "exp.*_vs.*_goslim$"
    suffix = "goslim"
    
class GOTable( ExpressionTracker ):
    
    select = '''SELECT CASE code WHEN '+' THEN 'over' ELSE 'under' END,
    goid, description, ratio, CASE code WHEN '+' THEN pover ELSE punder END AS pvalue
    FROM %(track)s WHERE passed AND code != '?'
         %(where)s
    ORDER by %(order)s '''
    
    columns = ("over", "goid", "category", "fold", "pvalue")
    order = "pvalue"
    
    def __init__(self, *args, **kwargs ): 
        ExpressionTracker.__init__(self, *args, **kwargs)

    def getTracks( self, subset = None ):
        tracks = ExpressionTracker.getTracks( self, subset)
        if tracks[0].endswith( self.suffix ):
            return tracks
        else:
            return [ "%s_%s" % (x, self.suffix) for x in tracks ]

    def getSlices( self, subset = None ):
        if subset == None: return ("biol_process", "cell_location", "mol_function" )
        else: return subset
        
    def __call__(self, track, slice = None ):

        order = self.order
        if slice:
            where = "AND category == '%(slice)s' " % locals()
        else:
            where = ""

        statement = self.select % (locals() )

        data = self.get( statement )
        return odict( zip(self.columns, zip(*data) ) )

class GOSummary( GOTable ):

    select = '''SELECT COUNT(*),
                SUM( CASE code WHEN '+' THEN 1 ELSE 0 END),
                SUM( CASE code WHEN '-' THEN 1 ELSE 0 END)
    FROM %(track)s WHERE passed AND code != '?' %(where)s'''
    
    columns = ("total", "over", "under" )

    def __call__(self, track, slice = None ):

        order = self.order
        if slice:
            where = "AND category == '%(slice)s' " % locals()
        else:
            where = ""

        statement = self.select % (locals() )

        data = self.getFirstRow( statement )
        return odict( zip(self.columns, data) )

class GOEnrichment( GO):
    mSelect = "SELECT description, 100.0 * (ratio - 1), pover FROM %s WHERE passed AND code = '+' "
    mColumns = ("category", "enrichment/%", "pvalue" )
    mOrder = "ratio DESC"

class GODepletion( GO ):
    mSelect = "SELECT description, 100.0 * (1 - ratio), punder FROM %s WHERE passed AND code = '-' "
    mColumns = ("category", "depletion/%", "pvalue" )
    mOrder = "ratio"

class GOPower( ExpressionTracker ):
    mPattern = "_annotation"
    mTableName = None

    def __init__(self, *args, **kwargs ): 
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice = None ):

        if not self.mTableName: raise NotImplementedError("table not specified")
        if not self.mSelect: raise NotImplementedError("invalid use of base class - only use derived classes")

        select = self.mSelect % self.mTableName
        xtrack, xslice, subset, background = track.split(":")            

        stmt = """%(select)s AND track = '%(xtrack)s' AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' AND background = '%(background)s'
                        AND code != '?'
                        ORDER BY category""" % locals()

        data = self.getValues( stmt )
      
        return data

class GOPowerEnrichment( GOPower ):
    """display power of go analyses.

    Show are the fold enrichment (in %) that could be detected
    given the setup of the test (size of sample and size of background).
    """
    mSelect = "SELECT 100.0 * ( (ci95upper / mean) - 1) FROM %s WHERE 1"
    mXLabel = "Power : fold enrichment / %"

class GOPowerDepletion( GOPower ):
    """display power of go analyses.

    Show are the fold depletion (in %) that could be detected
    given the setup of the test (size of sample and size of background).
    """
    mSelect = "SELECT 100.0 * ( 1 - (ci95lower / mean)) FROM %s WHERE 1"
    mXLabel = "Power : fold depletion / %"

class GOMatrix( ExpressionTracker ):
    """Display GO results in a matrix.
    """
    select = '''SELECT description,
                CASE passed
                WHEN 1 THEN
                          CASE code
                          WHEN '+' THEN (100.0 * (ratio -1))
                          ELSE - (100.0 * (1-ratio) )
                          END
                ELSE 0
                END
                FROM %(track)s WHERE code != '?' %(where)s ORDER BY %(order)s'''

    order = "category"
    def __init__(self, *args, **kwargs ): 
        ExpressionTracker.__init__(self, *args, **kwargs)

    def getTracks( self, subset = None ):
        tracks = ExpressionTracker.getTracks( self, subset)
        if tracks[0].endswith( self.suffix ):
            return tracks
        else:
            return [ "%s_%s" % (x, self.suffix) for x in tracks ]
                 

    def getSlices( self, subset = None ):
        if subset == None:
            return ("biol_process", "cell_location", "mol_function" )
        else: return subset

    def __call__(self, track, slice = None ):

        if not self.hasTable( track ):
            return []
        
        order = self.order
        if slice:
            where = "AND category = '%(slice)s' " % locals()
        else:
            where = ""
            
        statement = self.select % locals()
        data = self.get( statement )
        return odict( data) 

##=================================================================
## specific implementations of GO results
##=================================================================
_go_analysis = { "Enrichment" : GOEnrichment, 
                 "Depletion" : GODepletion, 
                 "Summary" : GOSummary, 
                 "Matrix" : GOMatrix,
                 "PowerEnrichment" : GOPowerEnrichment, 
                 "PowerDepletion" : GOPowerDepletion}

_go_analysis = { 
    "Summary" : GOSummary,
    "Table" : GOTable,
    "Matrix" : GOMatrix, 
    }


_go_territories = { "GO" : GO, 
                    "GOSlim" : GOSlim,
                    }

# the order of the base classes is important
# also: make sure that these are new-style classes
for a, aa in _go_analysis.items():
    for b, bb in _go_territories.items():
        n = "GO%s%s" % (a,b)
        globals()[n] = type( n, (bb,aa), {})

