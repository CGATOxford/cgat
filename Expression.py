'''utility methods for computing expression differences.'''

import math
import numpy
import sys, os
import collections

from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

import Experiment as E
import Pipeline as P
import Database
import IOTools
import Stats
import sqlite3

try:
    PARAMS = P.getParameters()
except IOError:
    pass

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

GeneExpressionResult = collections.namedtuple( "GeneExpressionResult", \
                                               "test_id treatment_name treatment_mean treatment_std " \
                                                   " control_name control_mean control_std " \
                                                   " pvalue qvalue l2fold fold called status" )


def writeExpressionResults( outfile, result ):
    '''output expression results table.'''
    outfile.write( "%s\n" % "\t".join(GeneExpressionResult._fields))
    for x in result:
        outfile.write("%s\n" % "\t".join( map(str,x)))

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

    fold change is treatment / control.

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
            # E.debug( "called genes=%s" % str(called_genes))
            
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

            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)

            if probeset in siggenes:
                s = siggenes[probeset]
                pvalue = s.rawp
                qvalue = s.qvalue
            else:
                pvalue = pvalues[probeset]
                qvalue = qvalues[probeset]

            called = (0,1)[probeset in called_genes]

            genes.append( GeneExpressionResult._make( (probeset,
                                                       "treatment",
                                                       mean1,
                                                       numpy.std( treatment ),
                                                       "control",
                                                       mean2,
                                                       numpy.std( control ),
                                                       pvalue,
                                                       qvalue,
                                                       mean1 - mean2,
                                                       math.pow(2,mean1 - mean2),
                                                       called,
                                                       "OK" ) ) )

        return genes, cutoff, fdr_values



#########################################################################
#########################################################################
#########################################################################
def loadTagData( infile, design_file ):
    '''load tag data for deseq/edger analysis.
    
    *Infile* is a tab-separated file with counts.

    *design_file* is a tab-separated file with the
    experimental design with four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

    track
        name of track - should correspond to column header in *infile*
    include
        flag to indicate whether or not to include this data
    group
        group indicator - experimental group
    pair
        pair that sample belongs to (for paired tests)

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups
    pairs  : vector with pairs

    It returns (groups,pairs)
    '''

    R( '''counts_table = read.delim( '%(infile)s', header = TRUE,
                                                   row.names = 1,
                                                   stringsAsFactors = TRUE )''' % locals() )

    E.info( "read data: %i observations for %i samples" % tuple(R('''dim(counts_table)''')))

    # Load comparisons from file
    R('''pheno = read.delim( '%(design_file)s', header = TRUE, stringsAsFactors = TRUE )''' % locals() )

    # Make sample names R-like - substitute - for . and add the .prep suffix
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')

    # Ensure pheno rows match count columns
    R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''' )

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''conds <- pheno2$group[ includedSamples ]''')

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- pheno2$group[ includedSamples ]''')
    R('''pairs = factor(pheno2$pair[ includedSamples ])''')

    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')

    E.info( "filtered data: %i observations for %i samples" % tuple( R('''dim(countsTable)''') ) )

    return groups, pairs

def runEdgeR( infile, 
              design_file, 
              outfile, 
              outfile_prefix = "edger.",
              fdr = 0.1,
              prefix = ""
              ):
    '''run DESeq on.

    See loadTagData on the input form format for *infile* and
    *design_file*.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    '''
    
    # load library 
    R('''suppressMessages(library('edgeR'))''')
    R('''suppressMessages(library('limma'))''')
    
    to_cluster = True

    logf = IOTools.openFile( outfile + ".log", "w" )
    
    groups, pairs = Expression.loadTagData( infile, design_file )

    # build DGEList object
    R( '''countsTable = DGEList( countsTable, group = groups )''' )

    # calculate normalisation factors
    E.info( "calculating normalization factors" )
    R('''countsTable = calcNormFactors( countsTable )''' )
    E.info( "output")
    # logf.write( str(R('''countsTable''')) + "\n" )

    # Remove windows with few counts
    R( '''countsTable = countsTable[rowSums( 
             1e+06 * countsTable$counts / 
             expandAsMatrix ( countsTable$samples$lib.size, dim(countsTable)) > 1 ) >= 2, ]''')

    E.info( "trimmed data: %i observations for %i samples" % tuple( R('''dim(countsTable)''') ) )

    # output MDS plot
    R.png( '''%(outfile_prefix)smds.png''' % locals() )
    R('''plotMDS( countsTable )''')
    R['dev.off']()

    # build design matrix
    R('''design = model.matrix( ~pairs + countsTable$samples$group )''' )
    R('''rownames(design) = rownames( countsTable$samples )''')
    R('''colnames(design)[length(colnames(design))] = "CD4" ''' )
    
    # logf.write( R('''design''') + "\n" )

    # estimate common dispersion
    R('''countsTable = estimateGLMCommonDisp( countsTable, design )''')
    
    # fitting model to each tag
    R('''fit = glmFit( countsTable, design, dispersion = countsTable$common.dispersion )''')

    # perform LR test
    R('''lrt = glmLRT( countsTable, fit)''' )

    E.info("Generating output")

    # compute adjusted P-Values
    R('''padj = p.adjust( lrt$table$p.value, 'BH' )''' )

    outf = IOTools.openFile( outfile, "w" )
    isna = R["is.na"]

    outf.write( "Group1\tGroup2\tinterval_id\tlogConc\tlfold\tLR\tpvalue\tpadj\tstatus\tsignificant\n" )
    rtype = collections.namedtuple( "rtype", "logConc lfold LR pvalue" )
    
    # output differences between pairs
    R.png( '''%(outfile_prefix)smaplot.png''' % locals() )
    R('''plotSmear( countsTable, pair=c('%s') )''' % "','".join( groups) )
    R('''abline( h = c(-2,2), col = 'dodgerblue') ''' )
    R['dev.off']()

    # I am assuming that logFC is the base 2 logarithm foldchange.
    # Parse results and parse to file
    results = []
    counts = E.Counter()

    for interval, data, padj in zip( R('''rownames(lrt$table)'''),
                                     zip( *R('''lrt$table''')), 
                                     R('''padj''')) :
        d = rtype._make( data )
        counts.input += 1

        # set significant flag
        if padj <= fdr: 
            signif = 1
            counts.significant += 1
            if d.lfold > 0:
                counts.significant_over += 1
            else:
                counts.significant_under += 1
        else: 
            signif = 0
            counts.insignificant += 1

        if d.lfold > 0:
            counts.all_over += 1
        else:
            counts.all_under += 1
        
        if isna( d.pvalue ): status = "OK"
        else: status = "FAIL"

        counts[status] += 1
        
        results.append( GeneExpressionResult._make( ( \
                    interval,
                    groups[0],
                    d.logConc,
                    0,
                    groups[1],
                    0,
                    0,
                    d.pvalue,
                    padj,
                    d.lfold,
                    math.pow( 2.0, d.lfold ),
                    str(signif),
                    status) ) )
            
    outf.close()

    with IOTools.openFile( outfile, "w" ) as outf:
        writeExpressionResults( outf, results )

    outf = IOTools.openFile( "(outfile_prefix)ssummary.tsv", "w" )
    outf.write( "category\tcounts\n%s\n" % counts.asTable() )
    outf.close()

def runDESeq( infile, 
              design_file, 
              outfile, 
              outfile_prefix = "deseq.",
              fdr = 0.1,
              prefix = ""
              ):
    '''run DESeq on.

    See loadTagData on the input form format for *infile* and
    *design_file*.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    DESeq ignores any pair information in the design matrix.
    
    Various plots are generate - annotation is from the manual (version 1.4)


    SVCPlot:
       squared coefficient of variation. Ratio of variance at base level to the
       square of the base mean.

       Solid lines are for the raw variances (biological noise).

       On top of the variance, there is shot noise, i.e., the Poissonean variance inherent to the
       process of counting reads. The amount of shot noise depends on the size factor, and hence, for
       each sample, a dotted line in the colour of its condition is plotted above the solid line. The dotted
       line is the base variance, i.e., the full variance, scaled down to base level by the size factors. The
       vertical distance between solid and dotted lines is the shot noise.
       The solid black line is a density estimate of the base means: Only were there is an appreciable
       number of base mean values, the variance estimates can be expected to be accurate.
       It is instructive to observe at which count level the biological noise starts to dominate the shot
       noise. At low counts, where shot noise dominates, higher sequencing depth (larger library size)
       will improve the signal-to-noise ratio while for high counts, where the biological noise dominates,
       only additional biological replicates will help.

   fit.png

       One should check whether the base variance functions seem to follow the empirical variance
       well. To this end, two diagnostic functions are provided. The function varianceFitDiagnostics
       returns, for a speci?ed condition, a data frame with four columns: the mean base level for each
       gene, the base variance as estimated from the count values of this gene only, and the ?tted base
       variance, i.e., the predicted value from the local ?t through the base variance estimates from
       all genes. As one typically has few replicates, the single-gene estimate of the base variance can
       deviate wildly from the ?tted value. To see whether this might be too wild, the cumulative prob-
       ability for this ratio of single-gene estimate to ?tted value is calculated from the ?2 distribution,
       as explained in the paper.      

       We may now plot the per-gene estimates of the base variance against the base levels and draw
       a line with the ?t from the local regression

    residuals.png
       Another way to study the diagnostic data is to check whether the probabilities in the fourth
       column of the diagnostics data frame are uniform, as they should be. One may simply look at the
       histogram of diagForGB$pchisq but a more convenient way is the function residualsEcdfPlot,
       which show empirical cumulative density functions (ECDF) strati?ed by base level.
       
    '''

    # load library 
    R('''suppressMessages(library('DESeq'))''')

    groups, pairs = loadTagData( infile, design_file )

    # Remove windows with no data
    R( '''max_counts = apply(counts_table,1,max)''' )
    R( '''counts_table = counts_table[max_counts>0,]''')
    E.info( "removed %i empty rows" % tuple( R('''sum(max_counts == 0)''') ) )
    E.info( "trimmed data: %i observations for %i samples" % tuple( R('''dim(counts_table)''') ) )

    # Test if replicates exist
    min_reps = R('''min(table(groups)) ''')
    no_replicates = False
    if min_reps < 2:
        no_replicates = True

    ######## Run DESeq
    # Create Count data object
    E.info( "running DESeq" )
    R('''cds <-newCountDataSet( countsTable, groups) ''')

    # Estimate size factors
    R('''cds <- estimateSizeFactors( cds )''')

    # Estimate variance
    if no_replicates:
        R('''cds <- estimateVarianceFunctions( cds, method="blind" )''')
    else:
        R('''cds <- estimateVarianceFunctions( cds )''')

    # Plot size factors
    R.png( '''%(outfile_prefix)ssize_factors.png''' % locals() )
    R('''par(mar=c(8,4,4,2))''')
    R('''barplot( sizeFactors( cds ), main="size factors", las=2)''')
    R['dev.off']()

    # output size factors
    size_factors = R('''sizeFactors( cds )''')
    samples = R('''names(sizeFactors(cds))''')
    with IOTools.openFile( "%(outfile_prefix)ssize_factors.tsv" % locals(), "w" ) as outf:
        outf.write( "sample\tfactor\n" )
        for name, x in zip( samples, size_factors):
            outf.write( "%s\t%s\n" % (name, str(x)))

    # Plot pairs
    R.png( '''%(outfile_prefix)spairs.png''' % locals(), width=960, height=960 )
    R('''panel.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y))
            txt <- format(c(r, 0.123456789), digits=digits)[1]
            txt <- paste(prefix, txt, sep="")
            if(missing(cex.cor)) cex <- 0.6/strwidth(txt)
            x = 0.5;
            y = 0.5;
            if (par("xlog")) { x = 10^x }; 
            if (par("ylog")) { y = 10^y }; 
            text(x, y, txt, cex = cex);
            }
       ''')
    R('''pairs( counts_table, lower.panel = panel.pearson, pch=".", log="xy" )''')
    R['dev.off']()

    # in DESeq versions > 1.6 the following can be used
    # to output normalized data
    # R('''write.table( counts(cds, normalized=TRUE), file='%(outfile_prefix)scounts.tsv.gz', sep='\t') ''' % locals())
    # output counts
    R('''write.table( counts(cds), file='%(outfile_prefix)scounts.tsv.gz', sep='\t') ''' % locals())    

    R.png( '''%(outfile_prefix)sscvplot.png''' % locals() )
    R('''scvPlot( cds, ylim = c(0,3))''')
    R['dev.off']()

    # Generate heatmap of variance stabilised data
    R('''vsd <- getVarianceStabilizedData( cds )''' )
    R('''dists <- dist( t( vsd ) )''')
    R.png( '''%(outfile_prefix)sheatmap.png''' % locals() )
    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )
    R['dev.off']()

    for group in groups:
        if not no_replicates:
            R.png( '''%(outfile_prefix)s%(group)s_fit.png''' % locals() )
            R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % group )
            R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
            R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
            R['dev.off']()
            R.png( '''%(outfile_prefix)s%(group)s_residuals.png''' % locals()  )
            R('''residualsEcdfPlot( cds, "%s" )''' % group )
            R['dev.off']()

    # Differential expression
    E.info("calling differential expression")
    R('''res <- nbinomTest( cds, '%s', '%s' )''' % (groups[0],groups[1]))

    # Plot significance
    R.png( '''%(outfile_prefix)ssignificance.png''' % locals() )
    R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, 
                    col = ifelse( res$padj < %(fdr)s, "red", "black" ) )''' % locals() )
    R['dev.off']()

    outf = IOTools.openFile( "%(outfile_prefix)sall.txt", "w" )
    isna = R["is.na"]

    E.info("Generating output")
    # Get column names from output and edit
    names = None
    if not names:
        names = list(R['res'].names)
        m = dict( [ (x,x) for x in names ])
        m.update( dict(
                pval = "pvalue", 
                baseMeanA = "value1", 
                baseMeanB = "value2",
                id = "interval_id", 
                log2FoldChange = "lfold") )
        
        header = [ m[x] for x in names ] 
        outf.write( "Group1\tGroup2\t%s\tstatus\tsignificant\n" % "\t".join(header))
    else:
        if names != list(R['res'].names):
            raise ValueError( "different column headers in DESeq output: %s vs %s" % (names, list(R['res'].names)))

    # Parse results and parse to file
    rtype = collections.namedtuple( "rtype", names )
    counts = E.Counter()
    
    results = []

    for data in zip( *R['res']) :
        counts.input += 1
        d = rtype._make( data )
        outf.write( "%s\t%s\t" % (groups[0],groups[1]))
        # set significant flag
        if d.padj <= fdr: 
            signif = 1
            counts.significant += 1
            if d.log2FoldChange > 0:
                counts.significant_over += 1
            else:
                counts.significant_under += 1
        else: 
            signif = 0
            counts.insignificant += 1

        if d.log2FoldChange > 0:
            counts.all_over += 1
        else:
            counts.all_under += 1

        # set lfold change to 0 if both are not expressed
        if d.baseMeanA == 0.0 and d.baseMeanB == 0.0:
            d = d._replace( foldChange = 0, log2FoldChange = 0 )

        if isna( d.pval ): status = "OK"
        else: status = "FAIL"

        counts[status] += 1

        outf.write( "\t".join( map(str, d) ))
        outf.write("\t%s\t%s\n" % (status, str(signif)))
        counts.output += 1

        results.append( GeneExpressionResult._make( ( \
                    d.id,
                    groups[0],
                    d.baseMeanA,
                    0,
                    groups[1],
                    d.baseMeanB,
                    0,
                    d.pval,
                    d.padj,
                    d.log2FoldChange,
                    d.foldChange,
                    str(signif),
                    status) ) )
                    

    outf.close()

    E.info( counts )

    with IOTools.openFile( outfile, "w" ) as outf:
        writeExpressionResults( outf, results )

    outf = IOTools.openFile( "(outfile_prefix)ssummary.tsv", "w" )
    outf.write( "category\tcounts\n%s\n" % counts.asTable() )
    outf.close()

Design = collections.namedtuple( "Design", ("include", "group", "pair") )

def readDesignFile( design_file ):
    '''reads a design file.'''
    
    design = collections.OrderedDict()
    with IOTools.openFile( design_file ) as inf:
        for line in inf:
            design = readDesignFile( design_file )
            if line.startswith("track"): continue
            track, include, group, pair = line[:-1].split("\t")
            if track in design: raise ValueError( "duplicate track '%s'" % track )
            design[track] = Design._make( (int(include), group, pair))
    return design

#########################################################################
#########################################################################
#########################################################################
def loadCuffdiff( infile, outfile ):
    '''load results from differential expression analysis and produce
    summary plots.

    Note: converts from ln(fold change) to log2 fold change.
   
    The cuffdiff output is parsed. 

    Pairwise comparisons in which one gene is not expressed (fpkm < fpkm_silent)
    are set to status 'NOCALL'. These transcripts might nevertheless be significant.
    '''

    prefix = P.toTable( outfile )
    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    to_cluster = False
    dbhandle = sqlite3.connect( PARAMS["database"] )

    CuffdiffResult = collections.namedtuple("CuffdiffResult",
                                            "test_id gene_id gene  locus   sample_1        sample_2  "\
                                            " status  value_1 value_2 l2fold  "\
                                            "test_stat p_value q_value significant " )

    tmpname = P.getTempFilename()    
    min_fpkm = PARAMS["cuffdiff_fpkm_expressed"]

    # ignore promoters and splicing - no fold change column, but  sqrt(JS)
    for fn, level in ( ("cds_exp.diff", "cds"),
                       ("gene_exp.diff", "gene"),
                       ("isoform_exp.diff", "isoform"),
                       # ("promoters.diff", "promotor"),
                       # ("splicing.diff", "splice"), 
                       ("tss_group_exp.diff", "tss") ):
        
        tablename = prefix + "_" + level + "_diff"

        results = []
        for line in IOTools.openFile( os.path.join( indir, fn) ):
            if line.startswith("test_id"): continue
            data = CuffdiffResult._make( line[:-1].split("\t"))
            status = data.status
            significant = [0,1][data.significant == "yes"]
            if status == "OK" and (float(data.value_1) < min_fpkm or float(data.value_2) < min_fpkm):
                status = "NOCALL"
            try:
                fold = math.pow(2.0, float(data.l2fold))
            except OverflowError:
                fold = "na"

            results.append( GeneExpressionResult._make( (
                        data.test_id,
                        data.sample_1,
                        data.value_1,
                        0,
                        data.sample_2,
                        data.value_2,
                        0,
                        data.p_value,
                        data.q_value,
                        data.l2fold,
                        fold,
                        significant,
                        status ) ) )
                
        with IOTools.openFile( tmpname, "w" ) as outf:
            writeExpressionResults( outf, results )
            
        # max/minimum fold change seems to be (-)1.79769e+308
        # ln to log2: multiply by log2(e)
        statement = '''cat %(tmpname)s 
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=treatment_name
              --index=control_name
              --index=test_id
              --table=%(tablename)s 
         >> %(outfile)s.log
         '''
        
        P.run()

    for fn, level in ( ("cds.fpkm_tracking", "cds" ),
                       ("genes.fpkm_tracking", "gene"),
                       ("isoforms.fpkm_tracking", "isoform"),
                       ("tss_groups.fpkm_tracking", "tss") ):

        tablename = prefix + "_" + level + "_levels" 

        statement = '''cat %(indir)s/%(fn)s
        | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=tracking_id
              --table=%(tablename)s 
         >> %(outfile)s.log
         '''
        
        P.run()

    ## build convenience table with tracks
    tablename = prefix + "_isoform_levels"
    tracks = Database.getColumnNames( dbhandle, tablename )
    tracks = [ x[:-len("_FPKM")] for x in tracks if x.endswith("_FPKM") ]
    
    tmpfile = P.getTempFile()
    tmpfile.write( "track\n" )
    tmpfile.write("\n".join(tracks) + "\n" )
    tmpfile.close()
    
    statement = P.load( tmpfile.name, outfile )
    os.unlink( tmpfile.name )

def runCuffdiff( bamfiles, 
                 design_file,
                 geneset_file,
                 outfile,
                 cuffdiff_options = "",
                 threads = 4,
                 fdr = 0.1 ):
    '''estimate differential expression using cuffdiff.

    infiles
       bam files

    geneset_file
       geneset to use for the analysis

    design_file
       design file describing which differential expression to test

    Replicates within each track are grouped.
    '''

    design = readDesignFile( design_file )

    to_cluster = True

    outdir = outfile + ".dir" 
    try: os.mkdir( outdir )
    except OSError: pass

    job_options= "-pe dedicated %i -R y" % threads

    # replicates are separated by ","
    reps = collections.defaultdict( list )
    for bamfile in bamfiles:
        groups = collections.defaultdict()
        track = P.snip( os.path.basename( bamfile ), ".accepted.bam" )
        if track not in design:
            E.warn( "bamfile '%s' not part of design - skipped" % bamfile )
            continue
        
        d = design[track]
        if not d.include: continue
        reps[d.group].append( bamfile )
        
    groups = sorted(reps.keys())
    labels = ",".join( groups )

    reps = "   ".join( [ ",".join( reps[group] ) for group in groups ] )

    statement = '''date > %(outfile)s; hostname >> %(outfile)s.log;
    cuffdiff --output-dir %(outdir)s
             --verbose
             --num-threads %(threads)i
             --labels %(labels)s
             --FDR %(fdr)f
             %(cuffdiff_options)s
             <(gunzip < %(geneset_file)s )
             %(reps)s
    >> %(outfile)s.log 2>&1;
    date >> %(outfile)s.log;
    '''
    P.run()
