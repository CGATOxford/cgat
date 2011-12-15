'''utility methods for computing expression differences.'''

import Stats
import math
import numpy
import Experiment as E
import sys

#import rpy
#from rpy import r as R

from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

import collections


def buildProbeset2Gene( infile, 
                        outfile, 
                        database = "hgu133plus2.db",
                        mapping = "hgu133plus2ENSEMBL" ):
    
    '''build map relating a probeset to an ENSEMBL gene_id'''
    
    R.library( database )

    # map is a Bimap object
    m = R(mapping)
    
    result = R.toTable(m)

    outf = open( outfile, "w")
    outf.write( "probe_id\tgene_id\n" )
    for probeset_id, gene_id in zip(result["probe_id"], 
                                    result["ensembl_id"] ):
        outf.write( "%s\t%s\n" % (probeset_id, gene_id))
    outf.close()

    E.info( "written %i mappings to %s: probes=%i, genes=%i" % \
                (len(result),
                 outfile,
                 len(set(result["probe_id"])),
                 len(set(result["ensembl_id"])) ) )

class GeneExpressionResult(object):
    pass

class WelchsTTest(object):
    '''base class for computing expression differences.
    '''

    def __call__(self, 
                 probesets, 
                 treatments, 
                 controls):
        
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

class SAMR( object ):
    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the samr library.

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See 

    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. This version is parameterised to reproduce Excel SAM
    by setting:
       var.equal = TRUE
       med = TRUE

    .. note:: 
        SAM requires log2 scaled expression levels.
    '''
    
    def __call__(self, probesets, 
                 treatments, 
                 controls,
                 pattern = None,
                 fdr = 0.10,
                 ngenes = None,
                 npermutations = 1000,
                 ndelta=10,
                 method = "ttest" ):

        if ngenes and fdr:
            raise ValueError( "either supply ngenes or fdr, but not both.")

        R.library("samr")

        m = numpy.matrix( treatments + controls )
        m = numpy.transpose(m)
        labels = numpy.array([1] * len(treatments) + [2] * len(controls))

        R.assign("x", numpy.array(m))
        R.assign("y", labels)
        R.assign("probesets", probesets)

        data = R('''data=list( x=x, y=y, geneid=1:length(probesets), genenames=probesets, logged2=TRUE)''' )
        result = R('''samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)''')
        R('''plot(samr.obj, delta=.4)''')

class SAM( object ):
    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the siggenes library. Note that there is also
    an rsam package at:

    http://rss.acs.unt.edu/Rdoc/library/samr/html/samr.html

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See 

    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. To parameterize the FDR to excel sam, set the
    flag *use_excel_sam*.

    .. note:: 
        SAM requires log2 scaled expression levels.

    I ran into trouble using this library. I was not able to
    reproduce the same results from the original SAM study getting
    differences in d and in the fdr.

    '''
    
    def __call__(self, probesets, 
                 treatments, 
                 controls,
                 pattern = None,
                 fdr = 0.10,
                 ngenes = None,
                 npermutations = 1000,
                 ndelta=10,
                 method = "ttest",
                 use_excel_sam = False ):

        if ngenes and fdr:
            raise ValueError( "either supply ngenes or fdr, but not both.")
        
        R.library("siggenes")

        m = numpy.matrix( treatments + controls )
        m = numpy.transpose(m)

        E.debug( "build expression matrix: %i x %i" % m.shape )

        labels = numpy.array([1] * len(treatments) + [0] * len(controls))
        ## 1000 permutations for P-Values of down to 0.0001. Setting this
        ## to a high value improved reproducibility of results.

        kwargs = {}
        # kwargs set to replicate excel SAM
        if use_excel_sam:
            kwargs.update( { "control" : R('''samControl( lambda = 0.5, n.delta = %(ndelta)s ) ''' % locals()),
                             "med" : True,
                             "var.equal": True } )
        else:
            kwargs.update( { "control" : R('''samControl( n.delta = %(ndelta)s ) ''' % locals()) }, 
                           )

        # the option B needs to be not set if wilc.stat is chosen

        if method == "ttest":
            kwargs["method"] = R('''d.stat''')
            kwargs["B"] = npermutations            
        elif method == "wilc":
            kwargs["method"] = R('''wilc.stat''')
        elif metod == "cat":
            kwargs["method"] = R('''cat.stat''')
        else:
            raise ValueError("unknown statistic `%s`" % method )

        E.info( "running sam with the following options: %s" % str(kwargs) )
        
        a = R.sam( numpy.array(m),
                   labels,
                   gene_names=numpy.array(probesets),
                   **kwargs )
        
        # E.debug("%s" % str(a))

        R.assign( "a", a )

        fdr_data = collections.namedtuple( "sam_fdr", ("delta", "p0", "false", "called", "fdr", "cutlow","cutup", "j2","j1" ) )
        cutoff_data = collections.namedtuple( "sam_cutoff", ("delta", "called", "fdr"))
        gene_data = collections.namedtuple( "sam_fdr", ("row","dvalue","stddev","rawp","qvalue","rfold" ) )

        def _totable( robj ):
            '''convert robj to a row-wise table.'''
            s = numpy.matrix( robj )
            t = [ numpy.array(x).reshape(-1,) for x in s ]
            return t

        # extract the fdr values
        # returns R matrix
        t = _totable( a.do_slot('mat.fdr') )
        assert len(t[0]) == len(fdr_data._fields)
        for x in t:
            E.debug( "x=%s" % str(x))
        fdr_values = [ fdr_data( *x ) for x in t ]

        # find d cutoff
        if fdr != None and fdr > 0:
            s = numpy.matrix( R.findDelta( a, fdr ) )
            try:
                cutoffs = [ cutoff_data( *numpy.array(x).reshape(-1,) ) for x in s ]
                E.debug( "sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs) ) )
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug( "could not get cutoff" )
                cutoff = None
        elif ngenes:
            s = numpy.matrix( R.findDelta( a, ngenes ) )
            try:
                cutoffs = [ cutoff_data( *numpy.array(x).reshape(-1,) ) for x in s ]
                E.debug( "sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs) ) )
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug( "could not get cutoff" )
                cutoff = None
        else:
            raise ValueError("either supply ngenes or fdr")

        # collect (unadjusted) p-values and qvalues for all probesets
        pvalues = dict( zip( probesets, R('''a@p.value''') ) )
        qvalues = dict( zip( probesets, R('''a@q.value''') ) )
        
        siggenes = {}        
        called_genes = set()
        if cutoff != None:
            E.debug( "using cutoff %s" % str(cutoff) )
            
            summary = R.summary( a, cutoff.delta )
            R.assign( "summary", summary )

            called_genes = set( [probesets[int(x)-1] for x in R('''summary@row.sig.genes''')] )
            E.debug( "called genes=%s" % str(called_genes))
            
            r_result = zip(*_totable( summary.do_slot( 'mat.sig' ) ))
            
            if len(r_result) > 0:

                assert len(r_result[0]) == 6, "expected six columns from siggenes module, got: %s" % len(r_result[0])
                
                for x in r_result:
                    if x[4] > fdr:
                        E.warn( "%s has qvalue (%f) larger than cutoff, but is called significant." % (str(x), x[4]))
                            
                # except TypeError:
                #     # only a single value
                #     x = [r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ]
                #     if x[4] > fdr:
                #         E.warn( "%s has qvalue (%f) larger than cutoff, but is called significant." % (str(x), x[4]))

                siggenes[probesets[int(x[0])-1]] = gene_data( *x )                

            if pattern:
                outfile = pattern % "sam.pdf"
                R.pdf(outfile)
                R.plot( a, cutoff.delta )
                R.plot( a )
                R['dev.off']()

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

            result.called = probeset in called_genes

            genes.append( result )

        return genes, cutoff, fdr_values
