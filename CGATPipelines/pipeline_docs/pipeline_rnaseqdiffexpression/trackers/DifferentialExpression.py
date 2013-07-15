import os, sys, re, types, itertools, math, sqlite3

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from collections import OrderedDict as odict

from RnaseqDiffExpressionReport import *

##############################################################
##############################################################
##############################################################
class TrackerDESeqFit( Tracker ):
    method = "deseq"
    tracks = [ x.asFile() for x in DESIGNS ]
    slices = [ x.asFile() for x in GENESETS ]

    def __call__(self, track, slice = None ):
        design = track
        level = "gene"
        geneset = slice
        method = self.method

        filename = "%(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_fit_%(track)s.png" % locals()

        # fitting information will not exist if there are no replicates
        if not os.path.exists( filename ): return None

        rst_text = '''
%(level)s %(track)s %(geneset)s
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: %(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_fit_%(track)s.png

.. figure:: %(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_residuals_%(track)s.png

''' % locals()

        return odict( (("text", rst_text),) )

##############################################################
##############################################################
##############################################################
class TrackerDESummaryDESeq( RnaseqTracker, SingleTableTrackerRows ):
    table = "deseq_stats"
    fields = ("level", "geneset", "treatment_name", "control_name", "design" )

class TrackerDESummaryEdgeR( RnaseqTracker, SingleTableTrackerRows ):
    table = "edger_stats"
    fields = ("level", "geneset", "treatment_name", "control_name" ,"design")

class TrackerDESummaryCuffdiff( RnaseqTracker, SingleTableTrackerRows ):
    table = "cuffdiff_stats"
    fields = ("level", "geneset", "treatment_name", "control_name", "design" )

##############################################################
##############################################################
##############################################################
class TrackerDESummaryPlots( Tracker ):

    tracks = [ x.asFile() for x in DESIGNS ]
    slices = [ x.asFile() for x in GENESETS ]

    def __call__(self, track, slice = None ):
        design = track
        geneset = slice
        method = self.method

        rst_text = []

        for level in self.levels:
            for fn in (
                "%(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_heatmap.png",
                "%(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_scvplot.png",
                "%(method)s.dir/%(design)s_%(geneset)s_%(method)s_%(level)s_pvalue_vs_length.png" ):
                f = fn % locals()
                if not os.path.exists(f): continue
                rst_text.append( ".. figure:: %(f)s" % locals() )
                
        if rst_text:
            rst_text = '''
%(design)s %(geneset)s
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

''' % locals() + "\n\n".join( rst_text )
        else:
            rst_text = ""

        return odict( (("text", rst_text),) )

class TrackerDESummaryPlotsDESeq(TrackerDESummaryPlots ):
    method = "deseq"
    levels = ("gene",)

class TrackerDESummaryPlotsEdgeR(TrackerDESummaryPlots ):
    method = "edger"
    levels = ("gene",)

class TrackerDESummaryPlotsCuffdiff(TrackerDESummaryPlots ):
    method = "cuffdiff"
    levels = CUFFDIFF_LEVELS

################################################################
################################################################
################################################################
class TrackerDifferentialExpression( Tracker ):

    tracks = [ x.asFile() for x in GENESETS ]
    slices = [ x.asFile() for x in DESIGNS ]

    def __call__(self, track, slice = None ):

        method = self.method
        design = slice
        rst_text = []
        geneset = track

        for level in self.levels:
            for x,y in itertools.combinations( EXPERIMENTS, 2 ):
                filename = "%(method)s.dir/%(design)s_%(geneset)s_%(level)s_%(x)s_vs_%(y)s_significance.png" % locals()
                if not os.path.exists( filename ): continue

                rst_text.append('''
%(design)s %(geneset)s %(level)s %(x)s vs %(y)s 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. figure:: %(filename)s

''' % locals())
                            
        return odict( (("text", "\n".join( rst_text)),) )

class TrackerDEPairwiseDESeq( TrackerDifferentialExpression ):
    method = "deseq"
    levels = ("gene",)
    
class TrackerDEPairwiseCuffdiff( TrackerDifferentialExpression ):
    method = "cuffdiff"
    levels = CUFFDIFF_LEVELS


#############################################################
#############################################################
#############################################################
##
#############################################################
class DifferentialExpressionComparison( RnaseqTracker ):

    tracks = list( itertools.combinations( ("deseq", "cuffdiff", "edger"), 3 ))
 
    slices = [ "%s_%s" % (y,x.asFile()) for x,y in itertools.product(GENESETS,DESIGNS) ]

class DifferentialExpressionOverlap( DifferentialExpressionComparison ):

    def __call__(self, track, slice = None ):
       
        pair1, pair2, pair3 = track
        

        a = self.get('''SELECT test_id FROM %(slice)s_%(pair1)s_gene_diff WHERE significant = 1''')
        b = self.get('''SELECT test_id FROM %(slice)s_%(pair2)s_gene_diff WHERE significant = 1''')
        c = self.get('''SELECT test_id FROM %(slice)s_%(pair3)s_gene_diff WHERE significant = 1''')

        a = set(map(str,a))
        b = set(map(str,b))
        c = set(map(str,c))

        return odict( ( (pair1, len(a)),
                        (pair2, len(b)),
                        (pair3, len(c)),
                        ("shared", len(a.intersection(b).intersection(c)) ) ) )
        

class DifferentialExpressionCorrelationPValueCuffdiffDeseq( DifferentialExpressionComparison ):
    '''fold change estimates per gene set.
    '''

    def getPairs(self, track):
        self.pair1, self.pair2 = track[0], track[1]   
        return self.pair1, self.pair2


    def __call__(self, track, slice = None ):
        
        pair1, pair2 = self.getPairs(track)
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        pvalues = {pair1:[], pair2: []}
        for pvals in cc.execute("""
                   SELECT a.pvalue, b.pvalue
                          FROM %s_%s_gene_diff AS a, 
                          %s_%s_gene_diff AS b
                          WHERE a.test_id = b.test_id
                          AND ABS( a.l2fold ) != 10
                          AND ABS( b.l2fold ) != 10
                          AND a.pvalue IS NOT NULL
                          AND b.pvalue IS NOT NULL
                          AND a.status == 'OK' and b.status == 'OK'
                           """ % (slice, pair1, slice, pair2)).fetchall():
            pvalues[pair1].append(-math.log10(pvals[0] + 0.0000001))
            pvalues[pair2].append(-math.log10(pvals[1] + 0.0000001))

        return pvalues

class DifferentialExpressionCorrelationPValueCuffdiffEdger( DifferentialExpressionCorrelationPValueCuffdiffDeseq ):
    '''fold change estimates per gene set.'''
    
    def getPairs(self, track):
        self.pair1, self.pair2 = track[1], track[2]   
        return self.pair1, self.pair2

class DifferentialExpressionCorrelationPValueDeseqEdger(DifferentialExpressionCorrelationPValueCuffdiffDeseq ):
    '''fold change estimates per gene set.'''
    
    def getPairs(self, track):
        self.pair1, self.pair2 = track[0], track[2]   
        return self.pair1, self.pair2


class DifferentialExpressionCorrelationFoldChangeCuffdiffDeseq( DifferentialExpressionComparison ):
    '''fold change estimates per gene set.'''
    
    def getPairs(self, track):
        self.pair1, self.pair2 = track[0], track[1]   
        return self.pair1, self.pair2

    def __call__(self, track, slice = None ):
        
        pair1, pair2 = self.getPairs(track)
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        fold_changes = {pair1:[], pair2: []}
        for folds in cc.execute("""
                   SELECT a.l2fold, b.l2fold
                          FROM design_%s_%s_gene_diff AS a, 
                          %s_%s_gene_diff AS b
                          WHERE a.test_id = b.test_id
                          AND ABS( a.l2fold ) < 10
                          AND ABS( b.l2fold ) < 10
                          AND a.pvalue IS NOT NULL
                          AND b.pvalue IS NOT NULL
                          AND a.status == 'OK' and b.status == 'OK'
                           """ % (slice, pair1, slice, pair2)).fetchall():
            fold_changes[pair1].append(folds[0])
            fold_changes[pair2].append(folds[1])

        return fold_changes

class DifferentialExpressionCorrelationFoldChangeCuffdiffEdger( DifferentialExpressionCorrelationFoldChangeCuffdiffDeseq ):
    '''fold change estimates per gene set.'''
    
    def getPairs(self, track):
        self.pair1, self.pair2 = track[1], track[2]   
        return self.pair1, self.pair2

class DifferentialExpressionCorrelationFoldChangeDeseqEdger( DifferentialExpressionCorrelationFoldChangeCuffdiffDeseq ):
    '''fold change estimates per gene set.'''
    
    def getPairs(self, track):
        self.pair1, self.pair2 = track[0], track[2]   
        return self.pair1, self.pair2

#############################################################
#############################################################
#############################################################
##
#############################################################
class VolcanoTracker( RnaseqTracker ):
    '''insert volcano plots.'''
    tracks = [ x.asFile() for x in GENESETS ]
    
    def __call__(self, track, slice = None ):

        # disabled - takes potentially a lot of time
        return None

        data = self.getAll( """SELECT l2fold,
                                    pvalue, 
                                    CASE WHEN value1 < value2 THEN value2 ELSE value1 END AS max_fpkm 
                                    FROM %(track)s_%(method)s_%(slice)s_diff 
                                    WHERE pvalue IS NOT NULL AND pvalue != 'nan' AND 
                                          l2fold IS NOT NULL AND 
                                          max_fpkm > 0""" )
        if data:
            data["pvalue"] = [ -math.log10( x + 0.00001 ) for x in data["pvalue"] ]
            data["max_fpkm"] = [ math.log10( x + 0.00001 ) for x in data["max_fpkm"] ]
        return data
                   
class VolcanoPlotCuffdiff( VolcanoTracker ):
    method = "cuffdiff"
    slices = CUFFDIFF_LEVELS

class VolcanoPlotDESeq( VolcanoTracker ):
    method = "deseq"
    slices = ("gene",)

#############################################################
#############################################################
#############################################################
##
#############################################################
class ExonCounts( RnaseqTracker ):
    '''get unnormalized read counts in the exons for a gene.
    
    The gene name is given as the slice.'''

    pattern = "(.*)_exon_counts"

    def __call__(self, track, options = None ):
        
        if not options: raise ValueError( 'tracker requires gene_id' )

        data = self.getAll('''SELECT gene_id, * FROM 
                       %(track)s_exon_counts
                       WHERE gene_id = '%(options)s' ''' )
        
        if data: del data["gene_id"]
        return data

    
