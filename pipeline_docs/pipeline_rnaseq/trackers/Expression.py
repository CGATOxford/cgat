import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

from RnaseqReport import *

class TrackerExpressionGeneset( RnaseqTracker ):
    min_fpkm = 1
    @property
    def tracks( self ):
        tables = set(self.getTableNames())
        d = [ "%s_cuffdiff_%s" % ( x.asTable(),y) for x,y in 
              itertools.product( GENESET_TRACKS,
                                 ("cds", "gene", "isoform", "tss") ) ]
        return [ x for x in d if ("%s_levels" % x) in tables ]

    slices = [ x.asTable() for x in EXPERIMENTS ]

class TrackerExpressionAbInitio( RnaseqTracker ):
    pattern = "(.*)_genes$"

class TrackerExpressionReference( RnaseqTracker ):
    pattern = "(.*)_ref_genes$"
    
class ExpressionFPKM( TrackerExpressionGeneset ):
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns( table ): return None
        statement = '''SELECT %(slice)s_fpkm FROM %(table)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        return odict( (("fpkm", data ),) )

class ExpressionNormalizedFPKM( TrackerExpressionGeneset ):
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns( table ): return None
        max_fpkm = float(self.getValue( '''SELECT max(%(slice)s_fpkm) FROM %(table)s'''))
        statement = '''SELECT CAST( %(slice)s_fpkm AS FLOAT) / %(max_fpkm)f FROM %(table)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        return odict( (("percent of max(fpkm)", data),) )

class ExpressionFPKMConfidence( TrackerExpressionGeneset ):

    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns( table ): return None
        statement = '''SELECT (%(slice)s_conf_hi - %(slice)s_conf_lo ) / %(slice)s_fpkm 
                       FROM %(table)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getValues( statement )
        return odict( (("relative_error", data ),) )

class ExpressionFPKMConfidenceCorrelation( TrackerExpressionGeneset ):

    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns( table ): return None
        statement = '''SELECT %(slice)s_fpkm AS fpkm, 
                              (%(slice)s_conf_hi - %(slice)s_conf_lo ) AS confidence
                       FROM %(table)s WHERE %(slice)s_fpkm > %(min_fpkm)f'''
        data = self.getAll( statement )
        return data

class ExpressionHighestExpressedGenes( TrackerExpressionGeneset ):
    '''output ten highest expressed genes.'''
    limit = 10
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        if c not in self.getColumns( table ): return None
        statement = '''SELECT tracking_id, gene_short_name, locus, class_code, %(slice)s_fpkm AS fpkm FROM %(table)s
                              ORDER BY %(slice)s_fpkm DESC LIMIT %(limit)i'''
        data = self.getAll( statement )
        data["locus"] = [ linkToUCSC( *splitLocus( x ) ) for x in data["locus"] ]

        return data

class ExpressionHighestExpressedGenesDetailed( TrackerExpressionGeneset ):
    '''output ten highest expressed genes.'''
    limit = 10
    def __call__(self, track, slice = None ):
        c = "%s_FPKM" % slice
        table = track + "_levels"
        geneset = track[:track.index("_")]
        if c not in self.getColumns( table ): return None
        statement = '''SELECT tracking_id, gene_short_name, locus, source, class, sense, %(slice)s_fpkm AS fpkm FROM %(table)s,
                              %(geneset)s_class AS class
                              WHERE tracking_id = class.transcript_id
                              ORDER BY %(slice)s_fpkm DESC LIMIT %(limit)i'''
        data = self.getAll( statement )
        data["locus"] = [ linkToUCSC( *splitLocus( x ) ) for x in data["locus"] ]

        return data


