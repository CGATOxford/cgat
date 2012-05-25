import os, sys, re, types, itertools, math
import numpy
import numpy.ma
import Stats
import Histogram
import xml.etree.ElementTree

from IntervalReport import *

##################################################################################
##################################################################################
##################################################################################
## compute MAST curve
##################################################################################
def computeMastCurve( evalues ):
    '''compute a MAST curve.

    see http://www.nature.com/nbt/journal/v26/n12/extref/nbt.1508-S1.pdf

    returns a tuple of arrays (evalues, with_motifs, explained )
    '''
    
    if len(evalues) == 0:
        raise ValueError( "no data" )

    mi, ma = math.floor(min(evalues)), math.ceil(max(evalues))

    if mi == ma:
        raise ValueError( "not enough data" )
        
    hist, bin_edges = numpy.histogram(evalues, bins = numpy.arange( mi, ma, 1.0) )
    with_motifs = numpy.cumsum( hist )
    explained = numpy.array( with_motifs )

    for x, evalue in enumerate(bin_edges[:-1]):
        explained[x] -= evalue
        
    explained[ explained < 0 ] = 0

    return bin_edges[:-1], with_motifs, explained

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
def getFDR( samples, control, num_bins = 1000 ):
    '''return the score cutoff at a certain FDR threshold using scores and control scores.
    Note that this method assumes that a higher score is a better result.

    The FDR is defined as fdr = expected number of false positives (FP) / number of positives (P)

    Given a certain score threshold s , the following will be used as approximations:

    FP: the number of controls with a score of less than or equal to s. These
        are all assumed to be false positives. Both samples and control should 
        contain rougly equal number of entries, but FP is scaled to be equivalent to P. 

    P: the number of samples with a score of less than or equal to s. These
       are a mixture of both true and false positives.
    
    returns the score cutoff at FDR threshold.
    '''

    if len(samples) == 0 or len(control) == 0: return None, None

    bins = 100
    mi1, ma1 = min(samples), max(samples)
    mi2, ma2 = min(control), max(control)
    mi,ma = min(mi1,mi2), max(ma1,ma2)
    hist_samples, bin_edges_samples = numpy.histogram( samples, range=(mi,ma), bins = num_bins )
    hist_control, bin_edges_control = numpy.histogram( control, range=(mi,ma), bins = num_bins )
    hist_samples = hist_samples[::-1].cumsum()
    hist_control = hist_control[::-1].cumsum()
    bin_edges = bin_edges_samples[::-1]

    ## correct for unequal size in the two sets
    correction = float(len(samples)) / len(control)

    fdrs = []
    m = 0
    for s, p, fp in zip( bin_edges[:-1], hist_samples, hist_control):
        if p != 0: 
            fdr = min(1.0, correction * float(fp) / p)
            m = max( m, fdr )
        fdrs.append( m )

    return bin_edges[:-1][::-1], fdrs[::-1]

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
def getMastEvalueCutoff( evalues, control_evalues, fdr = 0.1 ):
    '''return the E-Value cutoff at a certain FDR threshold
    using the control tracks.

    returns the evalue cutoff at FDR threshold.
    '''
    bin_edges, fdrs = getFDR( [-x for x in evalues], [-x for x in control_evalues] )

    for bin, f in zip(bin_edges, fdrs):
        if f < fdr: return -bin
        
    return 0


##################################################################################
##################################################################################
##################################################################################
## Base class for mast analysis
##################################################################################
class Mast( IntervalTracker ):
    pattern = "(.*)_mast$"

    def getSlices( self ):
        return self.getValues( "SELECT DISTINCT motif FROM motif_info" )

##################################################################################
##################################################################################
##################################################################################
## 
##################################################################################
class MastFDR( Mast ):
    """return arrays of glam2scan scores
    """
    pattern = "(*.)_mast$"

    def __call__(self, track, slice = None):
        evalues = self.getValues( "SELECT -evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        control_evalues = self.getValues( "SELECT -min_evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        bin_edges, fdrs = getFDR( evalues, control_evalues )
        if bin_edges == None: return odict()
        bin_edges = [-x for x in bin_edges ]
        print len(bin_edges), len(fdrs)

        return odict( (("score", bin_edges),
                       ("fdr", fdrs)))

##################################################################################
##################################################################################
##################################################################################
## Annotation of bases with SNPs
##################################################################################
class MastSummary( Mast ):
    """return summary of mast results.

    Return for each track the number of intervals in total,
    the number of intervals submitted to mast, 

    The evalue used as a MAST curve cutoff, the number and % explained using the Mast cutoff.

    """

    mEvalueCutoff = 1
    mFDR = 0.1

    def __call__(self, track, slice = None):

        data = []
        nintervals = self.getValue( "SELECT COUNT(*) FROM %(track)s_intervals" % locals() )
        data.append( ("nintervals", nintervals ) )
        data.append( ("nmast", self.getValue( "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )) )

        evalues = self.getValues( "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        if len(evalues) <= 1: return odict()

        try:
            bin_edges, with_motifs, explained = computeMastCurve( evalues )
        except ValueError, msg:
            return odict( ( ("msg", msg),) )
        
        if len(explained) == 0:
            return odict((("msg", "no data"), ) )

        am = numpy.argmax( explained )
        evalue = bin_edges[am]
        
        intervals_and_peakvals = \
                              self.get("""
                              SELECT peakval, evalue
                              FROM %(track)s_mast AS m, %(track)s_intervals AS i
                              WHERE i.interval_id = m.id AND motif = '%(slice)s'
                              ORDER BY peakval""" % locals() )

        intervals_with_motifs = len( [x for x in intervals_and_peakvals if x[1] <= evalue ] )
        ntop = nintervals / 4
        top25_with_motifs = len( [x for x in intervals_and_peakvals[-ntop:] if x[1] <= evalue ] )
        bottom25_with_motifs = len( [x for x in intervals_and_peakvals[:ntop] if x[1] <= evalue ] )
        
        data.append( ("MC-Evalue", bin_edges[am] ) )
        data.append( ("MC-explained", explained[am]) )
        data.append( ("MC-explained / %", "%5.2f" % (100.0 * explained[am] / nintervals)) )
        data.append( ("MC-with-motif", intervals_with_motifs ))
        data.append( ("MC-with-motif / %", "%5.2f" % (100.0 * intervals_with_motifs / nintervals) ) )
        data.append( ("MC-top25-with-motif", top25_with_motifs ))
        if ntop == 0: ntop = 1
        data.append( ("MC-top25-with-motif / %", "%5.2f" % (100.0 * top25_with_motifs / ntop) ) )
        data.append( ("MC-bottom25-with-motif", bottom25_with_motifs ))
        data.append( ("MC-bottom25-with-motif / %", "%5.2f" % (100.0 * bottom25_with_motifs / ntop) ) )

        # use control intervals to compute FDR              
        control_evalues = self.getValues( "SELECT min_evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        evalue_cutoff = getMastEvalueCutoff( evalues, control_evalues, self.mFDR )
        nexplained = len([x for x in evalues if x <= evalue_cutoff])
        data.append( ("FDR-Evalue", evalue_cutoff ) )
        data.append( ("FDR-explained", nexplained) )
        data.append( ("FDR-explained / %", "%5.2f" % (100.0 * nexplained / nintervals )) )

        # use a pre-defined E-value threshold
        threshold = self.mEvalueCutoff
        data.append( ("Evalue", self.mEvalueCutoff ) )
        n = self.getValue( "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s' AND evalue <= %(threshold)f" % locals())
        data.append( ("threshold-explained", n ) )
        data.append( ("threshold-explained / %", "%5.2f" % (100.0 * n / nintervals) ) )

        # use no threshold (nmatches > 1)
        n = self.getValue( "SELECT COUNT(*) FROM %(track)s_mast WHERE motif = '%(slice)s' AND nmatches > 0" % locals())
        data.append( ("nmatches-explained", n ) )
        data.append( ("nmatches-explained / %", "%5.2f" % (100.0 * n / nintervals) ) )

        intervals_and_peakvals = \
                               self.get("""
                               SELECT peakval, nmatches
                               FROM %(track)s_mast AS m, %(track)s_intervals AS i
                               WHERE i.interval_id = m.id AND motif = '%(slice)s'
                               ORDER BY peakval""" % locals() )

        intervals_with_motifs = len( [x for x in intervals_and_peakvals if x[1] > 0 ] )
        ntop = nintervals / 4
        top25_with_motifs = len( [x for x in intervals_and_peakvals[-ntop:] if x[1] > 0 ] )
        bottom25_with_motifs = len( [x for x in intervals_and_peakvals[:ntop] if x[1] > 0] )
        
        if ntop == 0: ntop = 1
        data.append( ("nmatches-top25-with-motif", top25_with_motifs ))
        data.append( ("nmatches-top25-with-motif / %", "%5.2f" % (100.0 * top25_with_motifs / ntop) ) )
        data.append( ("nmatches-bottom25-with-motif", bottom25_with_motifs ))
        data.append( ("nmatches-bottom25-with-motif / %", "%5.2f" % (100.0 * bottom25_with_motifs / ntop) ) )

        return odict(data)

class MastMotifEvalues( Mast ):
    '''distribution of evalues.'''
    
    def __call__(self, track, slice = None ):
        r = odict()
        for x in ("evalue", "l_evalue", "r_evalue" ):
            r[x] = self.getValues( "SELECT %(x)s FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        return r

# class MastNumberOfMotifs( Mast ):
#     '''number of motifs matching within intervals.'''
    
#     def __call__(self, track, slice = None ):
#         data = self.getValues( "SELECT nmatches FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
#         return odict( (("nmatches", data), ))

# class MastAllCorrelations( Mast ):
#     '''correlating all measures.'''
    
#     def __call__(self, track, slice = None ):
#         field = "length"
#         data = self.get( """SELECT m.evalue, m.nmatches, i.length, i.peakval, i.avgval
#                                  FROM %(track)s_mast as m, %(track)s_intervals as i 
#                                  WHERE i.interval_id = m.id AND motif = '%(slice)s'
#                                  ORDER BY i.%(field)s DESC""" 
#                               % locals() )
#         return odict( zip( ("evalue", "nmatches", "length", "peakval", "avgval"),
#                            zip(*data)))

# class MastPairwiseCorrelation( Mast ):
#     '''base class for correlating two measures.'''
    
#     def __call__(self, track, slice = None ):

#         field1 = self.mField1
#         field2 = self.mField2
#         data = self.get( """SELECT %(field1)s as a, %(field2)s AS b
#                                  FROM %(track)s_mast as m, %(track)s_intervals as i 
#                                  WHERE i.interval_id = m.id AND motif = '%(slice)s'"""
#                               % locals() )

#         return odict(zip( (field1, field2), zip(*data) ))

# class MastEvalueVersusLength( MastPairwiseCorrelation ):
#     '''correlate evalue with interval length.'''
#     mField1 = "evalue"
#     mField2 = "i.length"

# class MastEvalueVersusNumberOfMatches( MastPairwiseCorrelation ):
#     '''correlate evalue with number of motifs found.'''
#     mField1 = "evalue"
#     mField2 = "nmatches"

# class MastEvalueVersusPeakVal( MastPairwiseCorrelation ):
#     '''correlate evalue with peak value.'''
#     mField2 = "evalue"
#     mField1 = "peakval"

# class MastNMatchesVersusPeakVal( MastPairwiseCorrelation ):
#     '''correlate evalue with peak value.'''
#     mField2 = "nmatches"
#     mField1 = "peakval"

# class MastPeakValPerNMatches( MastPairwiseCorrelation ):
#     '''correlate evalue with peak value.'''
#     mField1 = "nmatches"
#     mField2 = "peakval"

#     def __call__(self, track, slice = None ):

#         data = MastPairwiseCorrelation.__call__( self, track, slice )

#         n = odict()
#         for x in sorted( data["nmatches"] ):
#             n[x] = []
        
#         for nmatches, peakval in sorted(zip( data["nmatches"], data["peakval"] )):
#             n[nmatches].append( peakval )

#         return odict( n )
        
class MastMotifLocation( Mast ):
    '''plot median position of motifs versus the peak location.

    The position is centered around 1kb.
    '''
    def __call__(self, track, slice = None ):

        data = numpy.array ( self.getValues( """SELECT (i.peakcenter - (m.start + (m.end - m.start) / 2)) / 
                                          500.0
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.nmatches = 1"""
                            % locals() ), numpy.float )
        
        data[ data < -1.0 ] = -1.0
        data[ data > 1.0 ] = 1.0
        
        return odict( (("distance", data),))

class MastMotifLocationMiddle( Mast ):
    '''plot median position of motifs versus the center of the interval.'''
    def __call__(self, track, slice = None ):

        # difference between
        #   middle of interval: i.start + i.length / 2
        #   middle of motif: m.start + (m.end - m.start) / 2
        # divide by (intervalsize - motifsize) / 2
        #
        # only take single matches (multiple matches need not be centered)
        data = self.getValues( """SELECT ((i.start + i.length / 2) - (m.start + (m.end - m.start) / 2)) 
                                         / ((CAST(i.length AS FLOAT) - (m.end - m.start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.nmatches = 1"""
                               % locals() )
        return odict( (("distance", data),))

class MastControlLocationMiddle( Mast ):
    '''plot median position of controls versus the center of the interval.'''
    def __call__(self, track, slice = None ):

        data1 = self.getValues( """SELECT ( (m.r_length / 2) - (m.r_start + (m.r_end - m.r_start) / 2) ) / ((CAST( m.r_length as float) - (m.r_end - m.r_start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.r_nmatches = 1"""
                            % locals() )
        data2 = self.getValues( """SELECT ( (m.l_length / 2) - (m.l_start + (m.l_end - m.l_start) / 2) ) / ((CAST( m.l_length as float) - (m.l_end - m.l_start))/2)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'
                                 AND m.l_nmatches = 1"""
                            % locals() )
        
        return odict( (("distance", data1 + data2),))

class MastCurve( Mast ):
    """Summary stats of mast results.
    """
    pattern = "(.*)_mast$"

    def __call__(self, track, slice = None):

        evalues = self.getValues( "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )
        
        if len(evalues) == 0: return odict()
        try:
            bin_edges, with_motifs, explained = computeMastCurve( evalues )
        except ValueError, msg:
            return odict()

        data = odict()
        data["with_motifs"] = odict( (("evalue", bin_edges),("with_motifs", with_motifs)) )
        data["explained"] = odict( (("evalue", bin_edges),("explained", explained)) )
        
        return data
 
# class MastROC( Mast ):
#     '''return a ROC curve. The ROC tests various peak parameters
#     whether they are good descriptors of a motif.

#     True/false positives are identified by the presence/absence
#     of a motif. The presence is tested using the E-value of
#     a motif.
#     '''

#     mPattern = "_mast$"
#     mFields = ("peakval", "avgval", "length" )

#     def __call__(self, track, slice = None):
#         data = []

#         # obtain evalue distribution
#         evalues = self.getValues( "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )

#         if len(evalues) == 0: return odict()
        
#         try:
#             bin_edges, with_motifs, explained = computeMastCurve( evalues )
#         except ValueError, msg:
#             return odict()

#         # determine the e-value cutoff as the maximum of "explained"
#         cutoff = bin_edges[numpy.argmax( explained )]

#         # retrieve values of interest together with e-value
#         result = odict()
        
#         for field in self.mFields:

#             values = self.get( """SELECT i.%(field)s, m.evalue 
#                                  FROM %(track)s_mast as m, %(track)s_intervals as i 
#                                  WHERE i.interval_id = m.id AND motif = '%(slice)s'
#                                  ORDER BY i.%(field)s DESC""" 
#                               % locals() )

#             try:
#                 roc = Stats.computeROC( [ (x[0], x[1] <= cutoff) for x in values ])
#             except ValueError, msg:
#                 # ignore results where there are no positives among values.
#                 continue

#             result[field] = odict( zip( ("FPR", "TPR"), zip(*roc)) )
#         return result

# class MastROCNMatches( Mast ):
#     '''return a ROC curve. The ROC tests various peak parameters
#     whether they are good descriptors of a motif.

#     True/false positives are identified by the presence/absence
#     of a motif using the field nmatches
#     '''

#     mPattern = "_mast$"
#     mFields = ("peakval", "avgval", "length" )

#     def __call__(self, track, slice = None):
#         data = []

#         # retrieve values of interest together with e-value
#         result = odict()
        
#         for field in self.mFields:

#             values = self.get( """SELECT i.%(field)s, m.nmatches
#                                  FROM %(track)s_mast as m, %(track)s_intervals as i 
#                                  WHERE i.interval_id = m.id AND motif = '%(slice)s'
#                                  ORDER BY i.%(field)s DESC""" 
#                               % locals() )

#             try:
#                 roc = Stats.computeROC( [ (x[0], x[1] > 0) for x in values ])
#             except ValueError, msg:
#                 # ignore results where there are no positives among values.
#                 continue

#             result[field] = odict( zip( ("FPR", "TPR"), zip(*roc)) )
            
#         return result

# class MastAUC( MastROC ):
#     '''return AUC for a ROC curve. The ROC tests various peak parameters
#     whether they are good descriptors of a motif.

#     True/false positives are identified by the presence/absence
#     of a motif.
#     '''

#     mPattern = "_mast$"
#     mFields = ("peakval", "avgval", "length" )

#     def __call__(self, track, slice = None):

#         data = MastROC.__call__(self, track, slice )
#         for k, d in data.iteritems():
#             data[k] = Stats.getAreaUnderCurve( d['FPR'], d['TPR'] )
#         return data

# class MastEvalues( Mast ):
#     """return arrays of mast evalues.
#     """
#     mPattern = "_mast$"

#     def __call__(self, track, slice = None):

#         return odict( (("evalue", self.getValues( "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )),
#                        ("evalue - control", self.getValues( "SELECT min_evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )) ) )


class MastPeakValWithMotif( Mast ):
    '''return for each peakval the proportion of intervals
    that have a motif.

    This tracker uses the nmatches indicator as cutoff.
    '''

    def __call__(self, track, slice ):

        # obtain evalue distribution
        data = self.get( '''
        SELECT i.peakval, m.nmatches
        FROM %(track)s_intervals AS i,
             %(track)s_mast AS m
        WHERE m.id = i.interval_id \
           AND m.motif = '%(slice)s' ORDER BY i.peakval DESC''' % locals() )

        result = Stats.getSensitivityRecall( [ (int(x[0]), x[1] > 0) for x in data ] )
        
        return odict( zip( ("peakval", "proportion with motif", "recall" ), zip( *result ) ) )

# class MastPeakValWithMotifEvalue( Mast ):
#     '''return for each peakval the proportion of intervals
#     that have a motif.

#     This class uses the ROC Evalue as cutoff.
#     '''

#     def __call__(self, track, slice = None ):

        
#         # obtain evalue distribution
#         evalues = self.getValues( "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals() )

#         if len(evalues) == 0: return odict()
        
#         try:
#             bin_edges, with_motifs, explained = computeMastCurve( evalues )
#         except ValueError, msg:
#             return odict()

#         # determine the e-value cutoff as the maximum of "explained"
#         cutoff = bin_edges[numpy.argmax( explained )]
        
#         data = self.get( '''
#         SELECT i.peakval, m.evalue
#         FROM %(track)s_intervals AS i,
#              %(track)s_mast AS m
#         WHERE m.id = i.interval_id \
#            AND m.motif = '%(slice)s' ORDER BY i.peakval DESC''' % locals() )

#         result = Stats.getSensitivityRecall( [ (int(x[0]), x[1] < cutoff) for x in data ] )
        
#         return odict( zip( ("peakval", "proportion with motif", "recall" ), zip( *result ) ) )
    
class MemeInputSequenceComposition( IntervalTracker ):
    '''distribution of sequence composition in sequences
       submitted to motif searches.'''
    pattern = "(.*)_motifseq_stats"
    slices = ('nA','nAT','nC','nG','nGC','nN','nT','nUnk','pA','pAT','pC','pG','pGC','pN','pT')

    def __call__(self, track, slice ):
        return self.getValues( '''SELECT %(slice)s FROM %(track)s_motifseq_stats''' )
    
class MemeRuns( IntervalTracker ):
    
    def getTracks( self ):
        return self.getValues( "SELECT DISTINCT track FROM meme_summary" ) 
    
    def __call__(self, track, slice = None ):
        
        memedir = os.path.abspath( os.path.join( EXPORTDIR, "meme", "%s.meme" % track ) )
        if not os.path.exists( memedir ): return None

        tomtomdir = os.path.abspath( os.path.join( EXPORTDIR, "tomtom", "%s.tomtom" % track ) )

        data = []

        tree = xml.etree.ElementTree.ElementTree()
        tree.parse( os.path.join( memedir, "meme.xml" ) )
        model =  tree.find( "model" )
        data.append( ("nsequences", int(model.find( "num_sequences" ).text) ) )
        data.append( ("nbases", int(model.find( "num_positions" ).text) ) )
        
        motifs = tree.find( "motifs" )
        nmotifs = 0
        for motif in motifs.getiterator( "motif" ):
            nmotifs += 1

        data.append( ("nmotifs", nmotifs ) )
        data.append( ("link", "`meme_%s <%s/meme.html>`_" % (track, memedir) ) )
        data.append( ("tomtom", "`tomtom_%s <%s/tomtom.html>`_" % (track, tomtomdir) ) )

        return odict(data)

class MemeResults( IntervalTracker ):

    def getTracks( self ):
        return self.getValues( "SELECT DISTINCT track FROM meme_summary" ) 
    
    def __call__(self, track, slice = None ):
        
        resultsdir = os.path.abspath( os.path.join( EXPORTDIR, "meme", "%s.meme" % track ) )
        if not os.path.exists( resultsdir ): return []

        tree = xml.etree.ElementTree.ElementTree()
        tree.parse( os.path.join( resultsdir, "meme.xml" ) )
        model =  tree.find( "model" )
        # data.append( ("nsequences", int(model.find( "num_sequences" ).text) ) )
        # data.append( ("nbases", int(model.find( "num_positions" ).text) ) )
        
        motifs = tree.find( "motifs" )
        nmotif = 0
        result = odict()
        for motif in motifs.getiterator( "motif" ):
            nmotif += 1
            result[str(nmotif)] = odict( (\
                    ("width", motif.get("width" )),
                    ("evalue", motif.get("e_value" )),
                    ("information content",motif.get("ic")),
                    ("sites", motif.get("sites" )),
                    ("link", "`meme_%s_%i <%s/meme.html#summary%i>`_" % (track, nmotif, resultsdir, nmotif)),
                    ("img", ".. image:: %s/logo%i.png" % (resultsdir,nmotif)),
                    ("rev", ".. image:: %s/logo_rc%i.png" % (resultsdir,nmotif)),
                    ))

        return result

class TomTomResults( IntervalTracker ):
    '''overview of tomtom results.'''


    pattern = "(.*)_tomtom$"

    def __call__(self, track, slice = None ):

        data = self.getAll( """SELECT query_id, target_id, target_name,
                               optimal_offset,pvalue,qvalue,evalue, overlap, query_consensus,
                               target_consensus, orientation 
                               FROM %(track)s_tomtom""" % locals())
        
        resultsdir = os.path.abspath( os.path.join( EXPORTDIR, "tomtom", "%s.tomtom" % re.sub("_", "-", track ) ))
        if not os.path.exists( resultsdir ): return []
        
        # format is: match_q_3_t_2_M01904
        # q_3: query_id
        # t_2: target database
        # M01904: target_id
        data['link'] = [ "`tomtom <%s/tomtom.html#match_q_%s_t_2_%s>`_" % (resultsdir, target_id, target_name)
                         for target_id, target_name in zip( data['query_id'], data['target_id']) ]
        return data


# class AnnotationsMatrix( DefaultTracker ):

    
#     def getSlices( self, subset = None ):
#         if subset: return subset
#         return []

#     def __call__(self, track, slice = None ):
 
#         result = odict()
#         rows = ("intergenic", "intronic", "upstream", "downstream", "utr", "cds", "other" )

#         statement = self.getStatement( slice )
#         data = self.get( statement % locals() )
#         levels = sorted(list(set( [ x[7] for x in data ] )))

#         for row in rows:
#             m = odict()
#             for l in levels: m[l] = 0
#             result[row] = m

#         map_level2col = dict( [(y,x) for x,y in enumerate(levels)] )
#         for intergenic, intronic, upstream, downstream, utr, coding, ambiguous, level in data:
#             col = level 
#             for x,v in enumerate( (intergenic, intronic, upstream, downstream, utr, coding, ambiguous)):
#                 if v: 
#                     row=rows[x]
#                     break
#             else:
#                 row = rows[-1]

#             result[row][col] += 1

#         return result

# class AnnotationsMotifs( AnnotationsMatrix ):
#     '''return a matrix with intervals stratified by motif presence
#     and location of the interval.
#     '''

#     mPattern = "_mast$"
    
#     def getStatement( self, slice = None ):

#         statement = '''
#         SELECT a.is_intergenic, a.is_intronic, a.is_upstream, a.is_downstream, a.is_utr, a.is_cds, a.is_ambiguous, 
#           CASE WHEN m.nmatches > 0 THEN motif || '+' ELSE motif || '-' END
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id,
#         %(track)s_mast AS m ON m.id = i.interval_id'''

#         if slice != None:
#             statement += " AND motif = '%(slice)s'"
#         return statement

# class AnnotationsPeakVal( AnnotationsMatrix ):
#     '''return a matrix with intervals stratified by peakval
#     and location of the interval.
#     '''
#     mPattern = "_annotations$"

#     def getStatement( self, slice = None ):

#         statement = '''
#         SELECT a.is_intergenic, a.is_intronic, a.is_upstream, a.is_downstream, a.is_utr, a.is_cds, a.is_ambiguous,
#         peakval
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id'''

#         return statement

# class AnnotationsPeakValData( DefaultTracker ):
#     '''return peakval for intervals falling into various regions.'''
    
#     def getSlices( self, subset = None ):
#         if subset: return subset
#         return ("intergenic", "intronic", "upstream", "downstream", "utr", "cds", "other" )

#     def __call__(self, track, slice = None ):

#         if slice == "other": slice = "ambiguous"
 
#         statement = '''
#         SELECT peakval
#         FROM %(track)s_intervals AS i,
#         %(track)s_annotations AS a ON a.gene_id = i.interval_id AND a.is_%(slice)s ''' % locals()
            
#         return odict( (("peakval", self.getValues( statement )),) )



