import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

from RnaseqReport import *

class TrackerExpressionGeneset( TrackerSQL ):
    min_fpkm = 1
    @property
    def tracks( self ):
        tables = set(self.getTableNames())
        d = [ "%s_cuffdiff_%s_levels" % ( x.asTable(),y) for x,y in 
              itertools.product( GENESET_TRACKS,
                                 ("cds", "gene", "isoform", "tss") ) ]
        return [ x for x in d if x in tables ]

    slices = [ x.asTable() for x in EXPERIMENTS ]

class TrackerExpressionAbInitio( TrackerSQL ):
    mPattern = "_genes$"

class TrackerExpressionReference( TrackerSQL ):
    mPattern = "_ref_genes$"
    
class ExpressionFPKM( TrackerExpressionGeneset ):
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        if c not in self.getColumns( track ): return None
        statement = '''SELECT %(slice)s_fpkm FROM %(track)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        return odict( (("fpkm", data ),) )

class ExpressionNormalizedFPKM( TrackerExpressionGeneset ):
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        if c not in self.getColumns( track ): return None
        max_fpkm = float(self.getValue( '''SELECT max(%(slice)s_fpkm) FROM %(track)s'''))
        statement = '''SELECT CAST( %(slice)s_fpkm AS FLOAT) / %(max_fpkm)f FROM %(track)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        d = float(len(data))
        return odict( (("percent of max(fpkm)", data),) )

class ExpressionFPKMConfidence( TrackerExpressionGeneset ):

    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        if c not in self.getColumns( track ): return None
        statement = '''SELECT (%(slice)s_conf_hi - %(slice)s_conf_lo ) / %(slice)s_fpkm 
                       FROM %(track)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        return odict( (("relative_error", data ),) )

class ExpressionFPKMConfidenceCorrelation( TrackerExpressionGeneset ):

    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        if c not in self.getColumns( track ): return None
        statement = '''SELECT %(slice)s_fpkm AS fpkm, 
                              (%(slice)s_conf_hi - %(slice)s_conf_lo ) AS confidence
                       FROM %(track)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getAll( statement )
        return data

class ExpressionHighestExpressedGenes( TrackerExpressionGeneset ):
    '''output ten highest expressed genes.'''
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        if c not in self.getColumns( track ): return None
        statement = '''SELECT gene_short_name, locus, class_code, %(slice)s_fpkm AS fpkm FROM %(track)s
                              ORDER BY %(slice)s_fpkm DESC LIMIT 10'''
        data = self.getAll( statement )
        data["locus"] = [ linkToUCSC( *splitLocus( x ) ) for x in data["locus"] ]

        return data
