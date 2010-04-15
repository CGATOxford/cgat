'''utility methods for computing expression differences.'''

import Stats
import math
import numpy
import Experiment as E
import rpy
from rpy import r as R
import collections

#import rpy2, rpy2.robjects as ro
#import rpy2.robjects.numpy2ri

class GeneExpressionResult(object):
    pass

class WelchsTTest(object):
    '''base class for computing expression differences.
    '''

    def __call__(self, probesets, treatments, controls):
        
        assert len(probesets) == len(treatments[0])
        assert len(probesets) == len(controls[0])
        
        nskipped = 0
        results = []

        for probeset, treatment, control in zip( probesets, zip(*treatments), zip(*controls) ):

            nval1, nval2 = len(treatment), len(control)
            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)
            stddev1, stddev2 = numpy.std(treatment), numpy.std(control)
            
            try:
                s = Stats.doWelchsTTest( nval1, mean1, stddev1,
                                         nval2, mean2, stddev2,
                                         alpha = 0.05 )
            except ValueError:
                E.warn("expressionDifferences: standard deviations are 0 for probeset %s - skipped" % probeset )
                nskipped += 1
                continue

            s.mProbeset = probeset
            results.append( s )

        qvalues =  Stats.doFDR( [x.mPValue for x in results ] ).mQValues

        for s, qvalue in zip(results, qvalues ):
            s.mQValue = qvalue

        return results, nskipped

class SAM( object ):
    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the siggenes library. Note that there is also
    an rsam package at:

    http://rss.acs.unt.edu/Rdoc/library/samr/html/samr.html

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.
    
    .. note:: 
        SAM requires log2 scaled expression levels.
    '''
    
    def __call__(self, probesets, treatments, controls,
                 pattern = None,
                 fdr = 0.10,
                 ngenes = None,
                 npermutations = 1000,
                 ndelta=10,
                 method = "ttest" ):

        if ngenes and fdr:
            raise ValueError( "either supply ngenes or fdr, but not both.")

        R.library("siggenes")

        m = numpy.matrix( treatments + controls )
        m = numpy.transpose(m)
        labels = numpy.array([1] * len(treatments) + [0] * len(controls))
        ## 1000 permutations for P-Values of down to 0.0001. Setting this
        ## to a high value improved reproducibility of results.

        # the option B needs to be not set if wilc.stat is chosen
        kwargs = {}
        if method == "ttest":
            kwargs["method"] = R('''d.stat''')
            kwargs["B"] = npermutations            
        elif method == "wilc":
            kwargs["method"] = R('''wilc.stat''')
        elif metod == "cat":
            kwargs["method"] = R('''cat.stat''')
        else:
            raise ValueError("unknown statistic `%s`" % method )
        
        a = R.sam( numpy.array(m),
                   labels,
                   gene_names=probesets,
                   n_delta=ndelta,
                   **kwargs )
        
        R.assign( "a", a )

        fdr_data = collections.namedtuple( "sam_fdr", ("delta", "p0", "false", "called", "fdr", "cutlow","cutup", "j2","j1" ) )
        cutoff_data = collections.namedtuple( "sam_cutoff", ("delta", "called", "fdr"))
        gene_data = collections.namedtuple( "sam_fdr", ("row","dvalue","stddev","rawp","qvalue","rfold" ) )

        # how to extract the fdr values
        # fdr_values = [ fdr_data( *x ) for x in R('''a@mat.fdr''') ]

        # find d cutoff
        if fdr != None and fdr > 0:
            try:
                cutoffs = [ cutoff_data( *x ) for x in ( R('''findDelta( a, %f )''' % fdr) ) ]
                E.debug( "sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs) ) )
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug( "could not get cutoff" )
                cutoff = None
        elif ngenes:
            try:
                cutoffs = [ cutoff_data( *x ) for x in ( R('''findDelta( a, genes = %i )''' % ngenes) ) ]
                E.debug( "sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs) ) )
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug( "could not get cutoff" )
                cutoff = None
        else:
            raise ValueError("either supply ngenes or fdr")
        
        # collect (unadjusted) p-values and qvalues for all probesets
        pvalues = R('''a@p.value''')
        qvalues = R('''a@q.value''')
        
        siggenes = {}        
        if cutoff != None:
            E.debug( "using cutoff %s" % str(cutoff) )

            summary = R.summary( a, cutoff.delta )
            R.assign( "summary", summary )
            r_result = R('''summary@mat.sig''') 
            if len(r_result) > 0:

                assert len(r_result) == 6, "expected six columns from siggenes module, got: %s" % str(r_result.keys())

                try:
                    # note that the qvalue can be higher than the threshold (I do not know why)
                    # the qvalue is thus bounded by the threshold in order to get consistent data
                    for x in zip( *[r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ] ):
                        if x[4] > fdr:
                            E.warn( "WARNING: %s has qvalue larger than cutoff, but is called significant. Set to %f" % (str(x), fdr))
                            x = list(x)
                            x[4] = fdr
                        siggenes[probesets[int(x[0])-1]] = gene_data( *x )
                except TypeError:
                    # only a single value
                    x = [r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ]
                    if x[4] > fdr:
                        E.warn( "WARNING: %s has qvalue larger than cutoff, but is called significant. Set to %f" % (str(x), fdr))
                        x = list(x)
                        x[4] = fdr
                    siggenes[probesets[int(x[0])-1]] = gene_data( *x )
                
            if pattern:
                outfile = pattern % "sam.pdf"
                R.pdf(outfile)
                R.plot( a, cutoff.delta )
                R.plot( a )
                R.dev_off()
        else:
            E.debug( "no cutoff found - no significant genes." )
            
        genes = []
        for probeset, treatment, control in zip( probesets, zip(*treatments), zip(*controls) ):

            nval1, nval2 = len(treatment), len(control)
            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)
            stddev1, stddev2 = numpy.std(treatment), numpy.std(control)
            result = GeneExpressionResult()

            result.probeset = probeset
            result.difference = mean1 - mean2
            result.fold = math.pow(2,result.difference)

            result.mean1, result.mean2 = mean1, mean2
            result.stddev1, result.stddev2 = stddev1, stddev2

            if probeset in siggenes:
                s = siggenes[probeset]
                result.pvalue = s.rawp
                result.qvalue = s.qvalue
            else:
                result.pvalue = pvalues[probeset]
                result.qvalue = qvalues[probeset]

            genes.append( result )

        return genes, cutoff
