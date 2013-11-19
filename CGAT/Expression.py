################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
Expression.py - wrap various differential expression tools
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides tools for differential expression analysis 
for a variety of methods.

Methods implemented are:

   DESeq
   EdgeR
   cuffdiff

There is also a command line interface. Note that the module is incomplete
as a stand-alone script as it requires to be executed in the context of
an existing pipeline for parameterization.

Usage
-----

Documentation
-------------

Code
----

'''

import math
import numpy
import sys
import os
import collections
import itertools
import re

from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

try:
    import CGAT.Experiment as E
    import CGAT.Pipeline as P
    import CGAT.Database as Database
    import CGAT.IOTools as IOTools
    import CGAT.Stats as Stats
except ImportError:
    import Experiment as E
    import Pipeline as P
    import Database,IOTools,Stats

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
                                                   " pvalue qvalue l2fold fold significant status" )


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
    by setting::
       
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
                 use_excel_sam = False,
                 treatment_label = "treatment",
                 control_label = "control" ):
        
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

        fdr_data = collections.namedtuple( "sam_fdr", ("delta", "p0", "false", "significant", "fdr", "cutlow","cutup", "j2","j1" ) )
        cutoff_data = collections.namedtuple( "sam_cutoff", ("delta", "significant", "fdr"))
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
        # for x in t: E.debug( "x=%s" % str(x))
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

        if pattern:
            outfile = pattern % "sam.pdf"
            R.pdf(outfile)
            if cutoff:
                R.plot( a, cutoff.delta )
            else:
                R.plot( a )
            R['dev.off']()
        
        siggenes = {}        
        significant_genes = set()
        if cutoff != None:
            E.debug( "using cutoff %s" % str(cutoff) )
            
            summary = R('''summary( a, %f )''' % cutoff.delta)

            # summary = R.summary( a, cutoff.delta )
            R.assign( "summary", summary )

            significant_genes = set( [probesets[int(x)-1] for x in R('''summary@row.sig.genes''')] )
            # E.debug( "significant genes=%s" % str(significant_genes))
            
            r_result = zip(*_totable( summary.do_slot( 'mat.sig' ) ))
            
            if len(r_result) > 0:

                assert len(r_result[0]) == 6, "expected six columns from siggenes module, got: %s" % len(r_result[0])
                
                for x in r_result:
                    if x[4] > fdr:
                        E.warn( "%s has qvalue (%f) larger than cutoff, but is significant significant." % (str(x), x[4]))
                            
                # except TypeError:
                #     # only a single value
                #     x = [r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ]
                #     if x[4] > fdr:
                #         E.warn( "%s has qvalue (%f) larger than cutoff, but is called significant." % (str(x), x[4]))

                siggenes[probesets[int(x[0])-1]] = gene_data( *x )                

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

            significant = (0,1)[probeset in significant_genes]

            genes.append( GeneExpressionResult._make( (probeset,
                                                       treatment_label,
                                                       mean1,
                                                       numpy.std( treatment ),
                                                       control_label,
                                                       mean2,
                                                       numpy.std( control ),
                                                       pvalue,
                                                       qvalue,
                                                       mean1 - mean2,
                                                       math.pow(2,mean1 - mean2),
                                                       significant,
                                                       "OK" ) ) )

        return genes, cutoff, fdr_values



#########################################################################
#########################################################################
#########################################################################
def loadTagData( tags_filename, design_filename ):
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

    '''

    E.info( "loading tag data from %s" % tags_filename)

    R( '''counts_table = read.delim( '%(tags_filename)s', 
                                     header = TRUE,
                                     row.names = 1,
                                     stringsAsFactors = TRUE,
                                     comment.char = '#' )''' % locals() )

    E.info( "read data: %i observations for %i samples" % tuple(R('''dim(counts_table)''')))
    E.debug( "sample names: %s" % R('''colnames(counts_table)'''))

    # Load comparisons from file
    R('''pheno = read.delim( '%(design_filename)s', 
                             header = TRUE, 
                             stringsAsFactors = TRUE,
                             comment.char = '#')''' % locals() )

    # Make sample names R-like - substitute - for .
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')
    E.debug( "design names: %s" % R('''pheno[,1]'''))
    
    # Ensure pheno rows match count columns
    pheno = R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''' )
    missing = R('''colnames(counts_table)[is.na(pheno2)][1]''')
    if missing:
        E.warn( "missing samples from design file are ignored: %s" % missing)
        
    # Subset data & set conditions
    R('''includedSamples <- !(is.na(pheno2$include) | pheno2$include == '0') ''')
    E.debug( "included samples: %s" % R('''colnames(counts_table)[includedSamples]''') )
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- factor(pheno2$group[ includedSamples ])''')
    R('''conds <- pheno2$group[ includedSamples ]''')
    R('''pairs = factor(pheno2$pair[ includedSamples ])''')
 
    E.info( "filtered data: %i observations for %i samples" % tuple( R('''dim(countsTable)''') ) )

def filterTagData( filter_min_counts_per_row = 1, 
                   filter_min_counts_per_sample = 10,
                   filter_percentile_rowsums = 0):
    '''filter tag data.

    * remove rows with at least x number of counts

    * remove samples with a maximum of *min_sample_counts*

    * remove the lowest percentile of rows in the table, sorted
       by total tags per row
    '''
    
    # Remove windows with no data
    R( '''max_counts = apply(countsTable,1,max)''' )
    R( '''countsTable = countsTable[max_counts>%i,]''' % filter_min_counts_per_row )
    E.info( "removed %i empty rows" % tuple( R('''sum(max_counts == 0)''') ) )
    observations, samples = tuple( R('''dim(countsTable)'''))
    E.info( "trimmed data: %i observations for %i samples" % (observations, samples ))

    # remove samples without data
    R('''max_counts = apply(countsTable,2,max)''' )

    empty_samples = tuple(R('''max_counts < %i''' % filter_min_counts_per_sample))
    sample_names = R('''colnames(countsTable)''')
    nempty_samples = sum( empty_samples)

    if nempty_samples:
        E.warn( "%i empty samples are being removed: %s" % \
                    (nempty_samples, ",".join( [sample_names[x] for x,y in enumerate( empty_samples) if y]) ) )
        R('''countsTable <- countsTable[, max_counts >= %i]''' % filter_min_counts_per_sample)
        R('''groups <- groups[max_counts >= %i]''' % filter_min_counts_per_sample)
        R('''pairs <- pairs[max_counts >= %i]''' % filter_min_counts_per_sample)
        observations, samples = tuple( R('''dim(countsTable)'''))


    # percentile filtering
    if filter_percentile_rowsums > 0:
        percentile = float( filter_percentile_rowsums) / 100.0
        R( '''sum_counts = rowSums( countsTable )''')
        R( '''take = (sum_counts > quantile( sum_counts, probs = %(percentile)f))''' % locals() )
        discard, keep = R('''table( take )''')
        E.info( "percentile filtering at level %f: keep=%i, discard=%i" % (filter_percentile_rowsums,
                                                                           keep, discard ) )
        R('''countsTable = countsTable[take,]''')

    observations, samples = tuple( R('''dim(countsTable)'''))

    return observations, samples

def groupTagData(ref_group = None):
    '''compute groups and pairs from tag data table.'''

    # Relevel the groups so that the reference comes first
    if not ref_group is None:
        R('''groups <- relevel(groups, ref = "%s")''' % ref_group)

    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')
    
    # Test if replicates exist - at least to samples pre replicate
    min_per_group = R('''min(table(groups)) ''')[0]
    has_replicates = min_per_group >= 2

    # Test if pairs exist:
    npairs = R('''length(table(pairs)) ''')[0]
    has_pairs = npairs == 2

    # at least two samples per pair
    if has_pairs:
        min_per_pair = R('''min(table(pairs)) ''')[0]
        has_pairs = min_per_pair >= 2

    return groups, pairs, has_replicates, has_pairs
    
def plotHeatmap( method = "correlation" ):
    '''plot a heatmap from countsTable.
    '''

    if method == "correlation":
        R('''dists <- dist( (1 - cor(countsTable)) / 2 )''')
    else:
        R('''dists <- dist( t(as.matrix(countsTable)), method = '%s' )''' % method)

    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )


def plotPCA():
    '''plot a PCA plot from countsTable.'''
    
    R('''library( RColorBrewer )''')
    R('''library( lattice )''')
    R('''pca = prcomp(t(countsTable))''')
    R('''if (length(groups) >= 3) colours = brewer.pal(nlevels(groups), "Paired") else colours = c("green", "blue")''')
    R('''xyplot( PC2 ~ PC1, data = as.data.frame(pca$x ), 
                 groups=groups, col = colours, 
                 main = draw.key(key = list(rect = list(col = colours), text = list(levels(groups)))),
                 cex=2,
                 pch=16,
                 )''')

def runEdgeR( outfile, 
              outfile_prefix = "edger.",
              fdr = 0.1,
              prefix = "",
              dispersion = None,
              ref_group = None
              ):
    '''run EdgeR on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The dispersion is usually measuered from replicates. If there are no 
    replicates, you need to set the *dispersion* explicitely.

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your experience with similar data,
       and use that. Although subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data, dispersion=0.1 for data
       on genetically identical model organisms or dispersion=0.01 for technical replicates.

    '''
    
    # load library 
    R('''suppressMessages(library('edgeR'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group)
    
    # output heatmap plot
    R.png( '%(outfile_prefix)sheatmap.png' % locals() )
    plotHeatmap()
    R['dev.off']()    

    E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' % \
               (groups, pairs, has_replicates, has_pairs))
    
    if has_pairs:
        # output difference between groups
        R.png( '''%(outfile_prefix)sbalance_groups.png''' % locals() )
        first = True
        for g1, g2 in itertools.combinations(groups, 2 ):
            R('''a = rowSums( countsTable[groups == '%s'] ) ''' % g1 )
            R('''b = rowSums( countsTable[groups == '%s'] ) ''' % g2 )
            if first:
                R('''plot( cumsum( sort(a - b) ), type = 'l') ''' )
                first = False
            else:
                R('''lines( cumsum( sort(a - b) )) ''' )

        R['dev.off']()

        R('''library('ggplot2')''')
        R('''library('reshape')''')

        # output difference between pairs within groups
        first = True
        legend = []
        for pair in pairs:
            for g1, g2 in itertools.combinations(groups, 2 ):
                key = re.sub( "-", "_", "pair_%s_%s_vs_%s" % (pair, g1,g2))
                legend.append( key )
                #print R('''colnames( countsTable) ''')
                #print R(''' pairs=='%s' ''' % pair)
                #print R(''' groups=='%s' ''' % g1)
                R('''a = rowSums( countsTable[pairs == '%s' & groups == '%s'] ) ''' % (pair,g1) )
                R('''b = rowSums( countsTable[pairs == '%s' & groups == '%s'] ) ''' % (pair,g2) )
                R('''c = cumsum( sort(a - b) )''' )
                R('''c = c - min(c)''')
                if first:
                    data = R( '''d = data.frame( %s = c)''' % key)
                    first = False
                else:
                    R('''d$%s = c''' % key)

        # remove row names (gene idenitifiers)
        R('''row.names(d) = NULL''')
        # add numbers of genes (x-axis)
        R('''d$genes=1:nrow(d)''')

        # merge data for ggplot
        R('''d = melt( d, 'genes', variable_name = 'comparison' )''')

        # plot
        R('''gp = ggplot(d)''')
        R('''pp = gp + \
            geom_line(aes(x=genes,y=value,group=comparison,color=comparison))''')
                    
        R.ggsave( '''%(outfile_prefix)sbalance_pairs.png''' % locals() )
        R['dev.off']()

        
    # build DGEList object
    # ignore message: "Calculating library sizes from column totals"
    R( '''countsTable = suppressMessages(DGEList( countsTable, group = groups ))''' )

    # Relevel groups to make the results predictable - IMS
    if not ref_group is None:
        R('''countsTable$samples$group <- relevel(countsTable$samples$group, ref = "%s")''' % ref_group)
    else:
        #if no ref_group provided use first group in groups
        R('''countsTable$sample$group <- relevel(countsTable$samples$group, ref = "%s")''' % groups[0])

    # calculate normalisation factors
    E.info( "calculating normalization factors" )
    R('''countsTable = calcNormFactors( countsTable )''' )
    E.info( "output")

    # output MDS plot
    R.png( '''%(outfile_prefix)smds.png''' % locals() )
    try:
        R('''plotMDS( countsTable )''')
    except rpy2.rinterface.RRuntimeError, msg:
        E.warn( "can not plot mds: %s" % msg)
    R['dev.off']()

    # build design matrix
    if has_pairs:
        R('''design = model.matrix( ~pairs + countsTable$samples$group )''' )
    else:
        R('''design = model.matrix( ~countsTable$samples$group )''' )

    # R('''rownames(design) = rownames( countsTable$samples )''')
    # R('''colnames(design)[length(colnames(design))] = "CD4" ''' )
    
    # fitting model to each tag
    if has_replicates:
        # estimate common dispersion
        R('''countsTable = estimateGLMCommonDisp( countsTable, design )''')
        # estimate tagwise dispersion
        R('''countsTable = estimateGLMTagwiseDisp( countsTable, design )''')
        # fitting model to each tag
        R('''fit = glmFit( countsTable, design )''')
    else:
        # fitting model to each tag
        if dispersion == None:
            raise ValueError( "no replicates and no dispersion" )
        E.warn("no replicates - using a fixed dispersion value" )
        R('''fit = glmFit( countsTable, design, dispersion = %f )''' % dispersion )

    # perform LR test
    R('''lrt = glmLRT(fit)''' )

    E.info("Generating output")

    # output cpm table
    R('''library(reshape2)''')
    R('''countsTable.cpm <- cpm(countsTable,  normalized.lib.sizes=TRUE)''')
    R('''countsTable.cpm.melt <- melt(countsTable.cpm)''')
    R('''names(countsTable.cpm.melt) <- c("id","sample","ncpm")''')
    R('''gz = gzfile( "%(outfile_prefix)scpm.tsv.gz", "w" )''' % locals() )
    R('''write.table(countsTable.cpm.melt, file=gz, sep = "\t", row.names=FALSE, quote=FALSE)''' )
    R('''close( gz )''')

    # compute adjusted P-Values
    R('''padj = p.adjust( lrt$table$PValue, 'BH' )''' )

    isna = R["is.na"]

    rtype = collections.namedtuple( "rtype", "lfold logCPM LR pvalue" )

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

        try:
            fold = math.pow( 2.0, d.lfold )
        except OverflowError:
            E.warn( "%s: fold change out of range: lfold=%f" % (interval, d.lfold ))
            # if out of range set to 0
            fold = 0
            
        # fold change is determined by the alphabetical order of the factors. 
        # Is the following correct?
        results.append( GeneExpressionResult._make( ( \
                    interval,
                    groups[1],
                    d.logCPM,
                    0,
                    groups[0],
                    d.logCPM,
                    0,
                    d.pvalue,
                    padj,
                    d.lfold,
                    fold,
                    str(signif),
                    status) ) )
            
    if outfile == sys.stdout:
        writeExpressionResults( outfile, results )
    else:
        with IOTools.openFile( outfile, "w" ) as outf:
            writeExpressionResults( outf, results )

    outf = IOTools.openFile( "%(outfile_prefix)ssummary.tsv" % locals(), "w" )
    outf.write( "category\tcounts\n%s\n" % counts.asTable() )
    outf.close()

## needs to put into class
##
def deseqPlotSizeFactors(outfile):
    '''plot size factors - requires cds object.'''
    R.png( outfile )
    R('''par(mar=c(8,4,4,2))''')
    R('''barplot( sizeFactors( cds ), main="size factors", las=2)''')
    R['dev.off']()

def deseqOutputSizeFactors( outfile ):
    '''output size factors - requires cds object.'''
    size_factors = R('''sizeFactors( cds )''')
    samples = R('''names(sizeFactors(cds))''')
    with IOTools.openFile( outfile, "w" ) as outf:
        outf.write( "sample\tfactor\n" )
        for name, x in zip( samples, size_factors):
            outf.write( "%s\t%s\n" % (name, str(x)))

def deseqPlotHeatmap( outfile ):
    '''plot a heatmap

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    
    R('''dists <- dist( t( exprs(vsd) ) )''')
    R.png( outfile )
    R('''heatmap.2( as.matrix( dists ), trace='none', margin=c(10,10) )''' )
    R['dev.off']()

def deseqPlotPCA( outfile ):
    '''plot a PCA

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    R.png( outfile )    
    R('''plotPCA( vsd )''')
    R['dev.off']()

def deseqPlotPairs( outfile ):
    '''requires counts table'''
    # Plot pairs
    R.png( outfile, width=960, height=960 )
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
    R('''pairs( countsTable, lower.panel = panel.pearson, pch=".", log="xy" )''')
    R['dev.off']()

def deseqPlotPvaluesAgainstRowsums( outfile ):
    '''plot pvalues against row sum rank.

    This plot is useful to see if quantile filtering could
    be applied.
    '''
    
    R('''counts_sum = rowSums( countsTable )''')
    R.png( outfile )    
    R('''plot( rank( counts_sum)/length(counts_sum),
               -log10( res$pval),
               pch = 16,
               cex= 0.1)''')

    R('''abline( a=3, b=0, col='red')''')
    R['dev.off']()

def deseqParseResults( control_name, treatment_name, fdr, vsd = False):
    '''parse deseq output.

    retrieve deseq results from object 'res' in R namespace.

    The 'res' object is a dataframe with the following columns (from the DESeq manual):

    id: The ID of the observable, taken from the row names of the
          counts slots.

    baseMean: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for both conditions

    baseMeanA: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for condition A

    baseMeanB: The base mean for condition B

    foldChange: The ratio meanB/meanA

    log2FoldChange: The log2 of the fold change

    pval: The p value for rejecting the null hypothesis 'meanA==meanB'

    padj: The adjusted p values (adjusted with 'p.adjust( pval,
          method="BH")')

    Here, 'conditionA' is 'control' and 'conditionB' is 'treatment' such that
    a foldChange of 2 means that treatment is twice upregulated compared to control.

    Returns a list of results.

    If vsd is True, the log fold change will be computed from the variance
    stabilized data.

    '''

    results = []
    isna = R["is.na"]

    # Get column names from output and edit
    names = list(R['res'].names)
    m = dict( [ (x,x) for x in names ])
    m.update( dict(
            pval = "pvalue", 
            baseMeanA = "value1", 
            baseMeanB = "value2",
            id = "interval_id", 
            log2FoldChange = "lfold") )
    
    rtype = collections.namedtuple( "rtype", names )
    counts = E.Counter()

    for data in zip( *R['res']) :
        counts.input += 1

        d = rtype._make( data )

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

        counts.output += 1

        # check if our assumptions about the direction of fold change
        # are correct
        assert (d.foldChange > 1) == (d.baseMeanB > d.baseMeanA )

        # note that fold change is computed as second group (B) divided by first (A)
        results.append( GeneExpressionResult._make( ( \
                    d.id,
                    treatment_name,
                    d.baseMeanB,
                    0,
                    control_name,
                    d.baseMeanA,
                    0,
                    d.pval,
                    d.padj,
                    d.log2FoldChange,
                    d.foldChange,
                    str(signif),
                    status) ) )
                    
    return results, counts

def runDESeq( outfile,
              outfile_prefix = "deseq.",
              fdr = 0.1,
              prefix = "",
              fit_type = "parametric",
              dispersion_method = "pooled",
              sharing_mode = "maximum",
              ref_group = None,
              ):
    '''run DESeq on countsTable.

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

    The output is treatment and control. Fold change values are computed as treatment divided by control.
    '''

    # load library 
    R('''suppressMessages(library('DESeq'))''')
    R('''suppressMessages(library('gplots'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData( ref_group )
 
    ######## Run DESeq
    # Create Count data object
    E.info( "running DESeq: replicates=%s" % (has_replicates))
    R('''cds <-newCountDataSet( countsTable, groups) ''')

    # Estimate size factors
    R('''cds <- estimateSizeFactors( cds )''')
    
    no_size_factors = R('''is.na(sum(sizeFactors(cds)))''' )[0]
    if no_size_factors:
        E.warn( "no size factors - can not estimate - no output" )
        return

    # estimate variance
    if has_replicates:
        E.info("replicates - estimating variance from replicates" )
    else:
        E.info("no replicates - estimating variance with method='blind'" )
        dispersion_method = "blind"

    E.info( "Dispersion method = %s, fit type =%s" % (dispersion_method, fit_type ) )
    R('''cds <- estimateDispersions( cds, 
                                     method='%(dispersion_method)s',
                                     fitType='%(fit_type)s',
                                     sharingMode='%(sharing_mode)s' )''' % locals())

    # plot fit - if method == "pooled":
    if dispersion_method == "pooled":
        R.png( '''%(outfile_prefix)sdispersion_estimates_pooled.png''' % locals() )
        R('''plotDispEsts( cds )''')
        R['dev.off']()
    else:
        dispersions = R('''ls(cds@fitInfo)''')
        for dispersion in dispersions:
            R.png( '''%(outfile_prefix)sdispersion_estimates_%(dispersion)s.png''' % locals() )
        R('''plotDispEsts( cds, name = '%(dispersion)s' )''' % locals())
        R['dev.off']()

    # plot size factors
    deseqPlotSizeFactors( '%(outfile_prefix)ssize_factors.png' % locals() )

    # output size factors
    deseqOutputSizeFactors( "%(outfile_prefix)ssize_factors.tsv" % locals() ) 
    
    # plot scatter plots of pairs
    deseqPlotPairs('%(outfile_prefix)spairs.png' % locals()) 

    if dispersion_method not in ("blind", "pooled"): 
        # also do a blind/pooled dispersion estimate
        R('''cds_blind <- estimateDispersions( cds, 
                                         method='blind',
                                         fitType='%(fit_type)s',
                                         sharingMode='%(sharing_mode)s' )''' % locals())
    else:
        R('''cds_blind = cds''')

    # perform variance stabilization for log2 fold changes
    R('''vsd = varianceStabilizingTransformation( cds_blind )''')

    # in DESeq versions > 1.6 the following can be used
    # to output normalized data
    # R('''write.table( counts(cds, normalized=TRUE), file='%(outfile_prefix)scounts.tsv.gz', sep='\t') ''' % locals())
    # output counts
    R('''write.table( counts(cds), file=gzfile('%(outfile_prefix)scounts.tsv.gz'), sep='\t') ''' % locals())
    
    # plot heatmap
    deseqPlotHeatmap( '%(outfile_prefix)sheatmap.png' % locals() )

    # plot PCA
    deseqPlotPCA( '%(outfile_prefix)spca.png' % locals() )

    for group in groups:
        if has_replicates:
            #R.png( '''%(outfile_prefix)s%(group)s_fit.png''' % locals() )
            #R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % group )
            #R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
            #R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
            #R['dev.off']()
            #R.png( '''%(outfile_prefix)s%(group)s_residuals.png''' % locals()  )
            #R('''residualsEcdfPlot( cds, "%s" )''' % group )
            #R['dev.off']()
            pass

    # Call diffential expression for all pairings of groups included in the design
    
    all_results = []
    for combination in itertools.combinations(groups,2):
        control, treatment = combination
        gfix = "%s_vs_%s_" % (control,treatment)

        outfile_groups_prefix = outfile_prefix + gfix
        E.info("calling differential expression for control=%s vs treatment=%s" % (control, treatment))
        R('''res <- nbinomTest( cds, '%s', '%s' )''' % (control,treatment) )

        # Plot significance
        R.png( '''%(outfile_groups_prefix)ssignificance.png''' % locals() )
        R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, 
                        col = ifelse( res$padj < %(fdr)s, "red", "black" ) )''' % locals() )
        R['dev.off']()

        # Plot pvalues against rowsums
        deseqPlotPvaluesAgainstRowsums( '%(outfile_groups_prefix)spvalue_rowsums.png' % locals() )
        
        E.info("Generating output (%s vs %s)" % (control, treatment))

        # Get variance stabilized fold changes - note the reversal of treatment/control 
        R('''vsd_l2f = (rowMeans( exprs(vsd)[,conditions(cds) == '%s', drop=FALSE] ) 
                      - rowMeans( exprs(vsd)[,conditions(cds) == '%s', drop=FALSE] ))''' % (treatment,control))

        # Plot vsd correlation, see Figure 14 in the DESeq manual
        # if you also want to colour by expression level
        R.png( '''%(outfile_groups_prefix)sfold_transformation.png''' % locals() )
        R('''plot( res$log2FoldChange, vsd_l2f,
                        pch=20, cex=.1, 
                        col = ifelse( res$padj < %(fdr)s, "red", "black" ) )''' % locals() )
        R['dev.off']()

        # Plot pvalue histogram
        R.png( '''%(outfile_groups_prefix)spvalue_histogram.png''' % locals() )
        R('''pvalues = res$pval''')
        R('''hist(pvalues, breaks=50, col='skyblue' )''')
        R['dev.off']()

        # Plot diagnostic plots for FDR
        R.png( '''%(outfile_groups_prefix)sfdr.png''' % locals() )
        R('''orderInPlot = order(pvalues)''')
        R('''showInPlot = (pvalues[orderInPlot] < 0.08)''')
        R('''plot( seq( along=which(showInPlot)), 
                   pvalues[orderInPlot][showInPlot], 
                   pch='.',
                   xlab=expression( rank(p[i]) ), 
                   ylab=expression( p[i] ) )''')
        R('''abline( a=0,b=%(fdr)f/length(pvalues), col="red") ''' % locals() )
        R['dev.off']()

        # Substitute log2 fold with variance stabilized l2fold value
        R('''res$log2FoldChange = vsd_l2f''' )

        # Parse results and parse to file
        results, counts = deseqParseResults( control, 
                                             treatment, 
                                             fdr = fdr )

        all_results += results

        E.info( counts )
        outf = IOTools.openFile( "%(outfile_groups_prefix)ssummary.tsv" % locals(), "w" )
        outf.write( "category\tcounts\n%s\n" % counts.asTable() )
        outf.close()

    if outfile == sys.stdout:
        writeExpressionResults(outfile, all_results)
    else:
        with IOTools.openFile(outfile, "w") as outf:
            writeExpressionResults(outf, all_results)

Design = collections.namedtuple( "Design", ("include", "group", "pair") )

def readDesignFile( design_file ):
    '''reads a design file.'''
    
    design = collections.OrderedDict()
    with IOTools.openFile( design_file ) as inf:
        for line in inf:
            if line.startswith("track"): continue
            track, include, group, pair = line[:-1].split("\t")
            if track in design: raise ValueError( "duplicate track '%s'" % track )
            design[track] = Design._make( (int(include), group, pair))
    return design

#########################################################################
#########################################################################
#########################################################################
def plotTagStats( infile, design_file, outfile ):
    '''provide summary plots for tag data.'''

    loadTagData( infile, design_file )

    nobservations, nsamples = filterTagData()

    if nobservations == 0:
        E.warn( "no observations - no output" )
        return

    if nsamples == 0:
        E.warn( "no samples remain after filtering - no output" )
        return

    groups, pairs, has_replicates, has_pairs = groupTagData() 

    # import rpy2.robjects.lib.ggplot2 as ggplot2

    R('''library('ggplot2')''')
    R('''library('reshape')''')

    R('''d = melt( log10(countsTable + 1), variable_name = 'sample' )''')
    R('''gp = ggplot(d)''')
    R('''pp = gp + \
        geom_density(aes(x=value,group=sample,color=sample,fill=sample),alpha=I(1/3))''')
    
    R.ggsave( outfile + ".densities.png" )
    R['dev.off']()

    R('''gp = ggplot(d)''')
    R('''pp = gp + \
        geom_boxplot(aes(x=sample,y=value,color=sample,fill=sample),size=0.3,alpha=I(1/3)) + 
        opts( axis.text.x = theme_text( angle=90, hjust=1, size=8 ) )''')

    R.ggsave( outfile + ".boxplots.png" )
    R['dev.off']()

#########################################################################
#########################################################################
#########################################################################
def plotDETagStats( infile, outfile ):
    '''provide summary plots for tag data.

    Stratify boxplots and densities according to differential expression calls.
    '''

    # import rpy2.robjects.lib.ggplot2 as ggplot2

    R('''library('ggplot2')''')
    R('''library('grid')''')
    R('''data = read.table( '%s', header = TRUE, row.names=1 )''' % infile ) 

    R(''' gp = ggplot(data)''')
    R('''a = gp + 
        geom_density(aes(x=log10(treatment_mean+1),group=factor(significant),
                                         color='factor(significant)',fill='factor(significant)'),alpha=I(1/3))''')

    R('''b = gp + 
        geom_density(aes(x=log10(control_mean+1),group=factor(significant),
                                         color=factor(significant),fill=factor(significant)),alpha=I(1/3))''')
    

    fn = outfile + ".densities.png" 
    R.png( fn )
    try:
        R('''grid.newpage()''')
        R.pushViewport(R.viewport( layout = R('''grid.layout''')(2,1)))
        R('''print( a, vp = viewport( layout.pos.row = 1, layout.pos.col = 1 ) )''')
        R('''print( b, vp = viewport( layout.pos.row = 2, layout.pos.col = 1 ) )''')
    except rpy2.rinterface.RRuntimeError:
        E.warn( "could not create %s" % fn )
    R['dev.off']()

    
    R('''grid.newpage()''')
    R.pushViewport(R.viewport( layout = R('''grid.layout''')(2,1)))

    R(''' gp = ggplot(data)''')
    R('''a = gp + 
        geom_boxplot(aes(x=factor(significant), y=log10(treatment_mean+1),
                                         color=factor(significant),fill=factor(significant)),
                             size=0.3,
                             alpha=I(1/3))''') 

    R('''b = gp + 
      geom_boxplot(aes(x=factor(significant), 
                       y=log10(control_mean+1),
                         color=factor(significant),
                         fill=factor(significant)),
                         size=0.3,
                         alpha=I(1/3)) +\
        opts( axis.text.x = theme_text( angle=90, hjust=1, size=8 ) )''')

    fn = outfile + ".boxplots.png" 
    R.png( fn )
    try:
        R('''print( a, vp = viewport( layout.pos.row = 1, layout.pos.col = 1 ) )''')
        R('''print( b, vp = viewport( layout.pos.row = 2, layout.pos.col = 1 ) )''')
    except rpy2.rinterface.RRuntimeError:
        E.warn( "could not create %s" % fn )
    R['dev.off']()

def parseCuffdiff( infile):
    '''parse a cuffdiff .diff output file.'''
    min_fpkm = PARAMS["cuffdiff_fpkm_expressed"]

    CuffdiffResult = collections.namedtuple("CuffdiffResult",
                                            "test_id gene_id gene  locus   sample_1 sample_2  " 
                                            " status  value_1 value_2 l2fold  " 
                                            "test_stat p_value q_value significant " )
    
    results = []

    for line in IOTools.openFile( infile ):
        if line.startswith("test_id"): continue
        data = CuffdiffResult._make( line[:-1].split("\t"))
        status = data.status
        significant = [0,1][data.significant == "yes"]
        if status == "OK" and (float(data.value_1) < min_fpkm or float(data.value_2) < min_fpkm):
            status = "NOCALL"

        try: fold = math.pow(2.0, float(data.l2fold))
        except OverflowError: fold = "na"

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
                                            
    return results

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

    This requires the cummeRbund library to be present in R.
    '''

    prefix = P.toTable( outfile )
    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    # E.info( "building cummeRbund database" )
    # R('''library(cummeRbund)''')
    # cuff = R('''readCufflinks(dir = %(indir)s, dbfile=%(indir)s/csvdb)''' )
    # to be continued

    to_cluster = False
    dbhandle = sqlite3.connect( PARAMS["database"] )

    tmpname = P.getTempFilename()    

    # ignore promoters and splicing - no fold change column, but  sqrt(JS)
    for fn, level in ( ("cds_exp.diff", "cds"),
                       ("gene_exp.diff", "gene"),
                       ("isoform_exp.diff", "isoform"),
                       # ("promoters.diff", "promotor"),
                       # ("splicing.diff", "splice"), 
                       ("tss_group_exp.diff", "tss") ):
        
        tablename = prefix + "_" + level + "_diff"

        infile = os.path.join( indir, fn)                
        results = parseCuffdiff( infile )

        with IOTools.openFile( tmpname, "w" ) as outf:
            writeExpressionResults( outf, results )
            
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


    # Jethro - load tables of sample specific cuffdiff fpkm values into csvdb
    for fn, level in ( ("cds.read_group_tracking", "cds" ),
                       ("genes.read_group_tracking", "gene"),
                       ("isoforms.read_group_tracking", "isoform"),
                       ("tss_groups.read_group_tracking", "tss") ):

        tablename = prefix + "_" + level + "sample_fpkms"
        
        tmpf = P.getTempFilename(".")
        inf = IOTools.openFile( os.path.join( indir, fn ) ).readlines()
        outf = IOTools.openFile( tmpf, "w" )

        samples = []
        genes= {}

        x = 0
        for line in inf:
            if x == 0:
                x += 1
                continue
            line = line.split()  
            gene_id = line[0]
            condition = line[1]
            replicate = line[2]
            fpkm = line[6]
            status = line[8]
    
            sample_id = condition + "_" + replicate
            
            if sample_id not in samples:
                samples.append( sample_id )
                
                if gene_id not in genes:
                    genes[gene_id] = {}
                    genes[gene_id][sample_id] = fpkm
                else: 
                    if sample_id in genes[gene_id]:
                        raise ValueError(
                            'sample_id %s appears twice in file for gene_id %s'
                            % ( sample_id, gene_id ) )
                    else:
                        if status != "OK":
                            genes[gene_id][sample_id] = status
                        else:
                            genes[gene_id][sample_id] = fpkm 


        samples = sorted( samples )

        # IMS - CDS files might be empty if not cds has been calculated for the genes
        # in the long term need to add CDS annotation to denovo predicted genesets
        # in meantime just skip if cds tracking file is empty

        if len(samples) == 0:
            continue

        headers = "gene_id\t" + "\t".join( [ x for x in samples ] )
        outf.write( headers + "\n" )

        for gene in genes.iterkeys():
            outf.write( gene + "\t" )
            x = 0
            while x < len( samples ) - 1:
                outf.write( genes[gene][samples[x]] + "\t" )
                x += 1
                outf.write( genes[gene][samples[len(samples)-1]] + "\n" )

        outf.close()


        statement = ( "cat %(tmpf)s |"
                      " python %(scriptsdir)s/csv2db.py "
                      "  %(csv2db_options)s"
                      "  --allow-empty"
                      "  --index=gene_id"
                      "  --table=%(tablename)s"
                      " >> %(outfile)s.log" )
        P.run()

        os.unlink( tmpf )


    ## build convenience table with tracks
    tablename = prefix + "_isoform_levels"
    tracks = Database.getColumnNames( dbhandle, tablename )
    tracks = [ x[:-len("_FPKM")] for x in tracks if x.endswith("_FPKM") ]
    
    tmpfile = P.getTempFile(dir=".")
    tmpfile.write( "track\n" )
    tmpfile.write("\n".join(tracks) + "\n" )
    tmpfile.close()
    
    statement = P.load( tmpfile.name, outfile)
    os.unlink( tmpfile.name )

def runCuffdiff( bamfiles, 
                 design_file,
                 geneset_file,
                 outfile,
                 cuffdiff_options = "",
                 threads = 4,
                 fdr = 0.1,
                 mask_file = None ):
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
        # .accepted.bam kept for legacy reasons (see rnaseq pipeline)
        track = P.snip( os.path.basename( bamfile ), ".bam", ".accepted.bam" )
        if track not in design:
            E.warn( "bamfile '%s' not part of design - skipped" % bamfile )
            continue
        
        d = design[track]
        if not d.include: continue
        reps[d.group].append( bamfile )
        
    groups = sorted(reps.keys())
    labels = ",".join( groups )
    reps = "   ".join( [ ",".join( reps[group] ) for group in groups ] )

    # Nick - add mask gtf to not assess rRNA and ChrM
    extra_options = []

    if mask_file:
        extra_options.append( " -M %s" % os.path.abspath( mask_file ) )

    extra_options = " ".join( extra_options )
    
    #IMS added a checkpoint to catch cuffdiff errors
    statement = '''date > %(outfile)s.log; hostname >> %(outfile)s.log;
    cuffdiff --output-dir %(outdir)s
             --verbose
             --num-threads %(threads)i
             --labels %(labels)s
             --FDR %(fdr)f
             %(extra_options)s
             %(cuffdiff_options)s
             <(gunzip < %(geneset_file)s )
             %(reps)s
    >> %(outfile)s.log 2>&1;
    checkpoint;
    date >> %(outfile)s.log;
    '''
    P.run()

    results = parseCuffdiff( os.path.join( outdir, "gene_exp.diff") )
    
    if outfile == sys.stdout:
        writeExpressionResults( outfile, results )
    else:
        with IOTools.openFile( outfile, "w" ) as outf:
            writeExpressionResults( outf, results )

def runMockAnalysis( outfile, 
                     outfile_prefix, 
                     ref_group = None,
                     pseudo_counts = 0):
    '''run a mock analysis on a count table.

    compute fold enrichment values, but do not normalize or
    perform any test.
    '''
    
    groups, pairs, has_replicates, has_pairs = groupTagData( ref_group )

    all_results = []
    for combination in itertools.combinations(groups,2):
        control, treatment = combination
        
        R('''control_counts = rowSums( countsTable[groups == '%s'] )''' % control )
        R('''treatment_counts = rowSums( countsTable[groups == '%s'] )''' % treatment )
        
        # add pseudocounts to enable analysis of regions 
        # that are absent/present
        if pseudo_counts: 
            R('''control_counts = control_counts + %f''' % pseudo_counts )
            R('''treatment_counts = treatment_counts + %f''' % pseudo_counts )
            
        R('''fc = treatment_counts / control_counts''')

        results = []

        for identifier, treatment_count, control_count, foldchange in \
                zip( R('''rownames( countsTable)'''),
                     R('''treatment_counts'''),
                     R('''control_counts'''),
                     R('''fc''')):
            try:
                log2fold = math.log( foldchange )
            except ValueError:
                log2fold = "Inf"
            
            results.append( GeneExpressionResult._make( ( \
                        identifier,
                        treatment,
                        treatment_count,
                        0,
                        control,
                        control_count,
                        0,
                        1,
                        1,
                        log2fold,
                        foldchange,
                        "0",
                        "OK") ) )
            
            
        all_results.extend( results )

    if outfile == sys.stdout:
        writeExpressionResults(outfile, all_results)
    else:
        with IOTools.openFile(outfile, "w") as outf:
            writeExpressionResults(outf, all_results)

        
def outputTagSummary( filename_tags, 
                      outfile, output_filename_pattern,
                      filename_design = None):
    '''output summary values for a count table.'''

    E.info( "loading tag data from %s" % filename_tags)

    if filename_design != None:
        # load all tag data
        loadTagData( filename_tags, filename_design )

        # filter
        nobservations, nsamples = filterTagData()
        
    else:
        # read complete table
        R( '''countsTable = read.delim( '%(filename_tags)s', 
                                         header = TRUE,
                                         row.names = 1,
                                         stringsAsFactors = TRUE,
                                         comment.char = '#' )''' % locals() )

        nobservations, nsamples = tuple(R('''dim(countsTable)'''))
        E.info( "read data: %i observations for %i samples" % (nobservations,nsamples))
        E.debug( "sample names: %s" % R('''colnames(countsTable)'''))
        R('''groups = factor(colnames( countsTable ))''')

    nrows, ncolumns = tuple(R('''dim(countsTable)'''))

    outfile.write( "metric\tvalue\tpercent\n" )
    outfile.write( "number of observations\t%i\t100\n" % nobservations )
    outfile.write( "number of samples\t%i\t100\n" % nsamples )
    
    # Count windows with no data
    R( '''max_counts = apply(countsTable,1,max)''' )
    
    # output distribution of maximum number of counts per window
    outfilename = output_filename_pattern + "max_counts.tsv.gz"
    E.info( "outputting maximum counts per window to %s" % outfilename )
    R( '''write.table( table(max_counts), file='%(outfilename)s', sep="\t", row.names=FALSE, quote=FALSE)''' % locals())

    # removing empty rows
    E.info( "removing rows with no counts in any sample" )
    R( '''countsTable = countsTable[max_counts>0,]''')

    for x in range( 0,20):
        nempty = tuple( R('''sum(max_counts <= %i)''' % x))[0]
        outfile.write( "max per row<=%i\t%i\t%f\n" % (x, nempty, 100.0 * nempty / nrows ) )
                       
    E.info( "removed %i empty rows" % tuple( R('''sum(max_counts == 0)''') ) )
    observations, samples = tuple( R('''dim(countsTable)'''))
    E.info( "trimmed data: %i observations for %i samples" % (observations, samples ))

    # build correlation
    R('''correlations = cor(countsTable)''')
    outfilename = output_filename_pattern + "correlation.tsv"
    E.info( "outputting sample correlations to %s" % outfilename )
    R('''write.table( correlations, file='%(outfilename)s', sep="\t", 
                      row.names=TRUE,
                      col.names=NA, 
                      quote=FALSE)''' % locals())

    # output heatmap based on correlations
    outfilename = output_filename_pattern + "heatmap.svg"
    R.svg( outfilename )
    plotHeatmap( method = "correlation" )
    R['dev.off']()    

    # output PCA
    outfilename = output_filename_pattern + "pca.svg"
    R.svg( outfilename )
    plotPCA()
    R['dev.off']()    
    
    # output an MDS plot
    R('''suppressMessages(library('limma'))''')
    outfilename = output_filename_pattern + "mds.svg"
    R.svg( outfilename )
    R('''plotMDS( countsTable )''')
    R['dev.off']()    

def dumpTagData( filename_tags, filename_design, outfile ):
    '''output filtered tag table.'''

    if outfile == sys.stdout:
        outfilename = ""
    else:
        outfilename = outfile.name

    # load all tag data
    loadTagData( filename_tags, filename_design )

    # filter
    nobservations, nsamples = filterTagData()

    # output
    R('''write.table( countsTable, 
                      file='%(outfilename)s',
                      sep='\t',
                      quote=FALSE)''' % locals() )
    

