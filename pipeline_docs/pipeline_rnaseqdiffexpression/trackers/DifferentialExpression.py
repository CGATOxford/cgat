import os, sys, re, types, itertools, math

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

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
    fields = ("level", "geneset", "treatment_name", "control_name" )

class TrackerDESummaryEdgeR( RnaseqTracker, SingleTableTrackerRows ):
    table = "edger_stats"
    fields = ("level", "geneset", "treatment_name", "control_name" )

class TrackerDESummaryCuffdiff( RnaseqTracker, SingleTableTrackerRows ):
    table = "cuffdiff_stats"
    fields = ("level", "geneset", "treatment_name", "control_name" )

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

    tracks = list( itertools.combinations( ("deseq", "cuffdiff"), 2 ))
    slices = [ x.asFile() for x in GENESETS ]

class DifferentialExpressionOverlap( DifferentialExpressionComparison ):

    def __call__(self, track, slice = None ):
        
        pair1, pair2 = track
        
        a = self.get('''SELECT test_id, treatment_name, control_name FROM %(slice)s_%(pair1)s_gene_diff WHERE significant''')
        b = self.get('''SELECT test_id, treatment_name, control_name FROM %(slice)s_%(pair2)s_gene_diff WHERE significant''')

        a = set(map(str,a))
        b = set(map(str,b))
        
        return odict( ( (pair1, len(a)),
                        (pair2, len(b)),
                        ("shared", len(a.intersection(b)) ) ) )
        

class DifferentialExpressionCorrelationPValue( DifferentialExpressionComparison ):
    '''fold change estimates per gene set.'''        
    def __call__(self, track, slice = None ):
        
        pair1, pair2 = track

        data = self.getAll( '''
                   SELECT a.pvalue as %(pair1)s, b.pvalue as %(pair2)s
                          FROM %(slice)s_%(pair1)s_gene_diff AS a, 
                               %(slice)s_%(pair2)s_gene_diff AS b 
                   WHERE a.test_id = b.test_id AND a.treatment_name = b.treatment_name AND a.control_name = b.control_name
                         AND ABS( a.l2fold ) != 10
                         AND ABS( b.l2fold ) != 10
                         AND a.pvalue IS NOT NULL
                         AND b.pvalue IS NOT NULL
                         AND a.status == 'OK' and b.status == 'OK'
                         ''' )

        for k in data.keys():
            data[k] = [ math.log10( x + 0.0000001 ) for x in data[k]  ]

        return data

class DifferentialExpressionCorrelationFoldChange( DifferentialExpressionComparison ):
    '''fold change estimates per gene set.'''
    def __call__(self, track, slice = None ):
        
        pair1, pair2 = track

        data = self.getAll( '''
                   SELECT a.l2fold as %(pair1)s, b.l2fold as %(pair2)s
                          FROM %(slice)s_%(pair1)s_gene_diff AS a, 
                               %(slice)s_%(pair2)s_gene_diff AS b 
                   WHERE a.test_id = b.test_id AND a.treatment_name = b.treatment_name AND a.control_name = b.control_name
                         AND ABS( a.l2fold ) != 10
                         AND ABS( b.l2fold ) != 10''' )


        return data

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

    
