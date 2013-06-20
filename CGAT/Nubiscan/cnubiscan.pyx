'''implement NUBISCAN like motif searches.'''
# import numpy
import math, re, sys
import collections, bisect
import CGAT.Stats as Stats

cimport cython

# note both import and cimport are necessary
import numpy
cimport numpy
DTYPE_INT = numpy.int
ctypedef numpy.int_t DTYPE_INT_t
DTYPE_FLOAT = numpy.float
ctypedef numpy.float_t DTYPE_FLOAT_t

@cython.boundscheck(False) # turn of bounds-checking for entire function
def computeScores( numpy.ndarray[DTYPE_INT_t, ndim=1]seq,
                   numpy.ndarray[DTYPE_FLOAT_t, ndim=2]matrix,
                   numpy.ndarray[DTYPE_FLOAT_t, ndim=1]weights ):
    '''return a score vector. Each position in seq is scored
    against the matrix.
    
    The matrix should be sorted by ACGT (rows) and motif positions
    (columns).

    The sequence is a numeric vector.
    '''

    assert seq.dtype == DTYPE_INT and matrix.dtype == DTYPE_FLOAT and weights.dtype == DTYPE_FLOAT
    cdef int lalphabet = matrix.shape[0]
    cdef int lmotif = matrix.shape[1]
    cdef int l = len(seq) - lmotif
    cdef numpy.ndarray[ DTYPE_FLOAT_t, ndim=1] result = numpy.zeros( l, dtype = DTYPE_FLOAT )
    cdef DTYPE_FLOAT_t score = 0
    cdef int x, y, i
    # nested loop implementation, could made more elegant and quicker using
    # numpy annotation? Convert to cython?
    for 0 <= x < l:
        score = 0
        i = 0
        for x <= y < x+lmotif:
            score += weights[i] * matrix[seq[y],i]
            i += 1
        result[x] = score

    return result

@cython.boundscheck(False) # turn of bounds-checking for entire function
def combineScores(
    numpy.ndarray[ DTYPE_FLOAT_t, ndim=1 ]sense_scores,
    numpy.ndarray[ DTYPE_FLOAT_t, ndim=1 ]antisense_scores,
    int lmotif,
    arrangement = "DR3" ):

    assert sense_scores.dtype == DTYPE_FLOAT and antisense_scores.dtype == DTYPE_FLOAT
    cdef int distance, x, y
    cdef int l = len(sense_scores)
        
    orientation, distance = arrangement[:2], int(arrangement[2:])
    distance += lmotif

    cdef numpy.ndarray[ DTYPE_FLOAT_t, ndim=1] sense_results = numpy.zeros( l, dtype = DTYPE_FLOAT )
    cdef numpy.ndarray[ DTYPE_FLOAT_t, ndim=1] antisense_results = numpy.zeros( l, dtype = DTYPE_FLOAT )
    
    if orientation == "DR":
        for 0 <= x < l - distance:
#             sense_results[x] = (sense_scores[x] + sense_scores[x+distance]) / 2.0
#             antisense_results[x] = (antisense_scores[x] + antisense_scores[x+distance]) / 2.0
              sense_results[x] = min(sense_scores[x],sense_scores[x+distance])
              antisense_results[x] = min(antisense_scores[x],antisense_scores[x+distance])
    elif orientation == "ER":
        for 0 <= x < l - distance:
#             sense_results[x] = (sense_scores[x] + antisense_scores[x+distance] ) / 2.0
#             antisense_results[x] = (antisense_scores[x] + sense_scores[x+distance]) / 2.0
            sense_results[x] = min(sense_scores[x], antisense_scores[x+distance])
            antisense_results[x] = min(antisense_scores[x], sense_scores[x+distance])
    else:
        raise ValueError("invalid orientation `%s`" % orientation )
            
    return sense_results, antisense_results

def mapSequence( seq ):
    '''mapSequence to numbers.'''
    chars = "ACGTN"
    result = []
    for x in re.sub("\s", "", seq).upper(): result.append( chars.find( x ) )
    assert min(result) >= 0
    return numpy.array( result, dtype=numpy.int )

def computeWeights( matrix ):
    '''compute weights from matrix.

    The formula in the NUBISCAN paper is Wi = 100 / ln(4) * ( sum_over_ACGT (pa * log(pa)) + ln(4) )
    '''

    lalphabet, lmotif = matrix.shape
    # ignore dummy row for N
    lalphabet -= 1
    g = math.log(lalphabet)
    weights = numpy.zeros( lmotif, numpy.float )
    for y in range( lmotif ):
        w = 0
        for x in range(lalphabet):
            v = matrix[x,y]
            if v != 0: i =  v * math.log( v ) 
            else: i = 0
            w += i 
        weights[y] = 100.0 / g * (w + g) 
    return weights

# first half site for rxrvdr
#      0.020000  0.100000  0.020000  0.860000
#      0.160000  0.000000  0.840000  0.000000
#      0.820000  0.120000  0.040000  0.020000
#      0.660000  0.340000  0.000000  0.000000
#      0.000000  1.000000  0.000000  0.000000
#      0.000000  0.320000  0.000000  0.680000
# second half site for rxrvdr (not used)
#     0.020000  0.000000  0.000000  0.980000
#     0.040000  0.040000  0.920000  0.000000
#     0.920000  0.020000  0.060000  0.000000
#     0.260000  0.740000  0.000000  0.000000
#     0.000000  1.000000  0.000000  0.000000
#     0.200000  0.380000  0.320000  0.100000

def getAntisenseMatrix( matrix, weights = None ):
    
    as_matrix = numpy.vstack( (matrix[3,::-1],
                               matrix[2,::-1],
                               matrix[1,::-1],
                               matrix[0,::-1],
                               matrix[4,::-1], ) )
    
    if weights != None:
        return as_matrix, weights[::-1]
    else:
        return as_matrix
    
class Matcher(object):

    def __init__(self, matrix, *args, **kwargs ):
        self.matrix = matrix
        self.weights = computeWeights( matrix )

        self.antisense_matrix, self.antisense_weights = \
                               getAntisenseMatrix( matrix,
                                                   self.weights )

        max_scores = numpy.array( self.matrix.max( 0 ).flat)
        self.max_score = numpy.dot( self.weights, max_scores)

    def run( self, sequences, arrangements):
        pass

    def parseArrangement( self, arrangement ):
        '''parse arrangment into tuple of orientation, distance.'''
        orientation, distance = arrangement[:2], int(arrangement[2:])
        return orientation, distance
    
NubiscanMatch = collections.namedtuple ("NubiscanMatch",
                                        '''id, length, start, end, strand,
                                        arrangement, 
                                        score, zscore, pvalue, qvalue,
                                        alternatives''' )
    
class Nubiscan( Matcher ):
    '''Find motif using the standard NUBISCAN procedure.

    Motifs are filtered using a z-score threshold.

    Z-scores are compute individually for each repeat arrangement.
    '''

    def run( self, sequences, arrangements, zscore_threshold = 6.5 ):

        lalphabet, lmotif = self.matrix.shape
        assert lalphabet == 5

        sense_matrix = self.matrix
        sense_weights = self.weights

        antisense_matrix, antisense_weights = self.antisense_matrix, \
                                              self.antisense_weights
                 
        max_scores = numpy.array( sense_matrix.max( 0 ).flat)
        max_score = numpy.dot( sense_weights, max_scores)

        # compute all scores for all sequence positions
        sense_halfsite_scores, antisense_halfsite_scores = [], []
        for sequence in sequences:
            s = mapSequence(sequence)
            sense_halfsite_scores.append( computeScores( s, sense_matrix, sense_weights ) / max_score )
            antisense_halfsite_scores.append( computeScores( s, antisense_matrix, antisense_weights ) / max_score )

        results = []                

        # iterate over arrangements
        for arrangement in arrangements:

            orientation, distance = arrangement[:2], int(arrangement[2:])
            # offset to end
            distance += 2 * lmotif 

            sense_combined_scores, antisense_combined_scores, all_scores = [], [], []
            for s_hs, as_hs in zip( sense_halfsite_scores, antisense_halfsite_scores ):
                s_cs, as_cs = combineScores(
                    s_hs,
                    as_hs,
                    lmotif,
                    arrangement )
                sense_combined_scores.append( s_cs )
                antisense_combined_scores.append( as_cs )

            all_scores = numpy.concatenate( sense_combined_scores + antisense_combined_scores )
            mean = numpy.mean( all_scores)
            std = numpy.std( all_scores )

            for x, scores in enumerate( sense_combined_scores ):
                for start, score in enumerate( scores ):
                    zscore = (score - mean) / std
                    if zscore >= zscore_threshold:
                        results.append( NubiscanMatch._make( (x,
                                                              0,
                                                              start,
                                                              start + distance,
                                                              "+",
                                                              arrangement,
                                                              score,
                                                              zscore,
                                                              0, 0, None
                                                               ) ) )

            for x, scores in enumerate( antisense_combined_scores ):
                for start, score in enumerate( scores ):
                    zscore = (score - mean) / std
                    if zscore >= zscore_threshold:
                        results.append( NubiscanMatch._make( (x,
                                                              start,
                                                              start + distance,
                                                              "-",
                                                              arrangement,
                                                              score,
                                                              zscore,
                                                              0, 0, None
                                                              ) ) )

        return results
        
class MatcherRandomisationSequence( Matcher ):
    '''Find motifs using the NUBISCAN procedure for a single sequence.

    Results are filtered using a FDR procedure.

    The procedure is as follows:

    The sequence is randomised x times and motif scores for
    all arrangements and locations are collected. From these scores an
    empirical P-Value distribution is derived. The number of samples
    can be small, as the number of positions is likely to be large.
    This randomisation ensures that sequences 

    Motif scores are computed for all arrangements and possible locations
    within the original sequence. Each motif score is assigned a pvalue based
    on the randomized positions.

    From the pvalues for motif scores for the original sequence, qvalues
    are computed using the method from Storey et al. (2002).

    Only values less than the q-value threshold are returned.
    '''

    def __init__(self, matrix, samples = 100, *args, **kwargs ):
        Matcher.__init__( self, matrix, *args, **kwargs )
        self.samples = samples
        self.binsize = 0.001
        
    def computeBackground( self, seq, arrangements ):
        '''randomise sequences and compute P-Values for scores.

        return an array of P-Values.
        '''
        
        lalphabet, lmotif = self.matrix.shape
        # create a copy
        r = numpy.array(seq)
        lseq = len(seq)
        # nbins = int(math.ceil(1.0 / self.binsize)) + 1
        bg_scores = []
        t = 0
        for x in range( self.samples ):
            numpy.random.shuffle(r)
            # compute sense(s)/antisense(a) halfsite (hs) scores
            s_hs, a_hs = computeScores( r, self.matrix, self.weights ) / self.max_score, \
                         computeScores( r, self.antisense_matrix, self.antisense_weights ) / self.max_score

            for arrangement in arrangements:
                s_cs, as_cs = combineScores(
                    s_hs,
                    a_hs,
                    lmotif,
                    arrangement )
                # for x in s_cs: hist[int(x/self.binsize)] += 1
                # for x in as_cs: hist[int(x/self.binsize)] += 1
                
                # bg_scores.extend( s_cs )
                # bg_scores.extend( as_cs )
                bg_scores.append( s_cs )
                bg_scores.append( as_cs )
                
        #bg_scores = numpy.array( bg_scores )
        #bg_mean = numpy.mean( bg_scores)
        #bg_std = numpy.std( bg_scores )
        #t = len(bg_scores)
        hist, bin_edges = numpy.histogram( numpy.concatenate(bg_scores),
                                           bins = numpy.arange(0,1.0+2*self.binsize,self.binsize),
                                           normed = False, new = True )
        bg_mean = 0.5
        bg_std = 0.1
        t = sum(hist)
        # t = lseq * self.samples * 2 
        pvalues = (hist[::-1].cumsum() / float(t))[::-1]
        return pvalues, t, bg_mean, bg_std

    def collectMatches( self, seq, arrangements, pvalues ):
        '''compute matching score and pvalue for each site and arrangement.
        '''

        lalphabet, lmotif = self.matrix.shape
        sense_halfsite_scores, anti_halfsite_scores = [], []
        s_hs = computeScores( seq, self.matrix, self.weights ) \
               / self.max_score 
        a_hs = computeScores( seq, self.antisense_matrix, self.antisense_weights ) \
               / self.max_score 

        results = []
        # iterate over arrangements
        append = results.append
        for arrangement in arrangements:

            s_cs, as_cs = combineScores( s_hs,
                                         a_hs,
                                         lmotif,
                                         arrangement )
            
            for start, score in enumerate( s_cs ):
                pvalue = pvalues[ int(score / self.binsize) ]
                append( (arrangement, start, "+", pvalue, score ) )
                
            for start, score in enumerate( as_cs ):
                pvalue = pvalues[ int(score / self.binsize) ]
                append( (arrangement, start, "-", pvalue, score ) )

        return results
        
    def run( self, sequence, arrangements, qvalue_threshold = None):
        
        lalphabet, lmotif = self.matrix.shape
        assert lalphabet == 5

        seq = mapSequence(sequence)

        # do the sampling
        bg_pvalues, t, bg_mean, bg_std = self.computeBackground( seq, arrangements )

        # compute all matches
        matches = self.collectMatches( seq, arrangements, bg_pvalues )

        if qvalue_threshold != None:
            # collect pvalues
            pvalues = [ x[3] for x in matches ]

            # estimate qvalues for matches
            fdr = Stats.doFDR( pvalues )
            qvalues = fdr.mQValues
        else:
            qvalues = [0] * len(matches)
            qvalue_threshold = 1

        results = []
        # filter by qvalue
        for 0 <= x < len(matches):
            match = matches[x]
            qvalue = qvalues[x]
            if qvalue > qvalue_threshold: continue
            arrangement, start, strand, pvalue, score = match
            orientation, distance = self.parseArrangement( arrangement )
            zscore = (score - bg_mean) / bg_std
            results.append( NubiscanMatch._make( ("",
                                                  start,
                                                  start + distance + 2*lmotif,
                                                  strand,
                                                  arrangement,
                                                  score,
                                                  zscore,
                                                  pvalue,
                                                  qvalue,
                                                  None
                                                  ) ) )
        
        return results

class MatcherRandomisationSequences( MatcherRandomisationSequence ):
    '''Find motifs using the NUBISCAN procedure for multiple sequences.

    Results are filtered using a FDR procedure.

    The procedure is as follows:

    The sequence is randomised x times and motif scores for
    all arrangements and locations are collected. From these scores an
    empirical P-Value distribution is derived. The number of samples
    can be small, as the number of positions is likely to be large.
    This randomisation ensures that sequences 

    Motif scores are computed for all arrangements and possible locations
    within the original sequence. Each motif score is assigned a pvalue based
    on the randomized positions.

    From the pvalues for motif scores for the original sequence, qvalues
    are computed using the method from Storey et al. (2002).

    Only values less than the q-value threshold are returned.
    '''

    def __init__(self, matrix, samples = 100, *args, **kwargs ):
        Matcher.__init__( self, matrix, *args, **kwargs )
        self.samples = samples
        self.binsize = 0.001
        
    def collectMatches( self, seq, arrangements ):
        '''compute matching score and pvalue for each site and arrangement.
        '''

        lalphabet, lmotif = self.matrix.shape
        sense_halfsite_scores, anti_halfsite_scores = [], []
        s_hs = computeScores( seq, self.matrix, self.weights ) \
               / self.max_score
        a_hs = computeScores( seq, self.antisense_matrix, self.antisense_weights ) \
               / self.max_score 

        results = []
         
        # iterate over arrangements
        append = results.append
        for arrangement in arrangements:

            s_cs, as_cs = combineScores( s_hs,
                                         a_hs,
                                         lmotif,
                                         arrangement )

            append( (arrangement, s_cs, as_cs ) )

        return results
        
    def run( self,
             sequences,
             arrangements,
             qvalue_threshold = 0.05,
             pvalue_threshold = 0.01 ):
        
        lalphabet, lmotif = self.matrix.shape
        assert lalphabet == 5

        # do the sampling per sequence
        all_pvalues = []
        
        # keep only matches below a certain threshold -
        # otherwise there are too many.
        putative_matches = []
        binsize = self.binsize
        append = all_pvalues.append
        cdef int start, l
        cdef float log_pvalue_threshold = numpy.log(pvalue_threshold)

        for seq_id, seq in enumerate(sequences):

            seq = mapSequence(seq)
            
            bg_pvalues, t, bg_mean, bg_std = self.computeBackground( seq, arrangements )

            # convert pvalues to logspace
            # minimum P-value: number of values
            bg_pvalues = numpy.log( bg_pvalues )
            min_p = 1.0 / t
            bg_pvalues[ numpy.isneginf(bg_pvalues) ] = numpy.log( min_p )
            # compute all matches
            matches = self.collectMatches( seq, arrangements )
        
            for arrangement, sense_scores, anti_scores in matches:
                l = len(sense_scores)
                
                # compute p-values (in log-space)
                sense_p = numpy.array( [ bg_pvalues[int(x/binsize)] for x in sense_scores ] )
                anti_p = numpy.array( [ bg_pvalues[int(x/binsize)] for x in anti_scores ] )
                
                all_pvalues.append( sense_p )
                all_pvalues.append( anti_p )
                
                # record putative matches:
                for 0 <= start < l:
                    if sense_p[start] < log_pvalue_threshold:
                        putative_matches.append( (arrangement, seq_id, start, "+", sense_scores[start], sense_p[start]) )
                    if anti_p[start] < log_pvalue_threshold:
                        putative_matches.append( (arrangement, seq_id, start, "-", anti_scores[start], anti_p[start]) )

        # there are too many p-values to compute q-values for all
        # instead - work within a histogram.
        #
        # compute empirical qvalues based on all pvalues
        all_pvalues = numpy.concatenate(all_pvalues)
        hist, bin_edges = numpy.histogram( all_pvalues,
                                           bins=10000,
                                           normed = False,
                                           new = True )

        # qvalue = pi0 * t * pvalue / {number of observations <= pvalue}
        t = sum(hist)
        pi0 = Stats.getPi0( numpy.exp(all_pvalues) )
        # factor in pi0 to nobservations
        nobservations = hist.cumsum() / float(t) / pi0

        results = []
        # filter by qvalue
        nq = len(nobservations) -1
        for 0 <= x < len(putative_matches):
            match = putative_matches[x]
            arrangement, seq_id, start, strand, score, log_pvalue = match
            pvalue = math.exp(log_pvalue)
            # number of observations <= pvalue
            nobs = nobservations[ min(bisect.bisect( bin_edges, log_pvalue), nq) ]
            qvalue = pvalue / nobs
            if qvalue > qvalue_threshold: continue
            
            orientation, distance = self.parseArrangement( arrangement )
            zscore = (score - bg_mean) / bg_std
            results.append( NubiscanMatch._make( (seq_id,
                                                  len(sequences[seq_id]),
                                                  start,
                                                  start + distance + 2*lmotif,
                                                  strand,
                                                  arrangement,
                                                  score,
                                                  zscore,
                                                  pvalue,
                                                  qvalue,
                                                  None
                                                  ) ) )

        # save results for calls
        self.bin_edges = bin_edges
        self.hist = hist
        self.nobservations = nobservations
        self.pi0 = pi0
        
        return results

