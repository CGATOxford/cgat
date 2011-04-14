import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *

class TrackerCoverage( TrackerSQL ):
    @property
    def tracks( self ):
        tables = set(self.getTableNames())
        d = [ "%s_cuffdiff_%s_levels" % ( x.asTable(),y) for x,y in 
              itertools.product( GENESET_TRACKS,
                                 ("cds", "gene", "isoform", "tss") ) ]
        return [ x for x in d if x in tables ]

    slices = [ x.asTable() for x in EXPERIMENTS ]


