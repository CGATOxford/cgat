import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict
from exomeReport import *

class VariantSummary( ExomeTracker, SingleTableTrackerRows ):
    table = "vcf_stats"

class SnpSummary( ExomeTracker, SingleTableTrackerRows ):
    table = "snp_stats"

class IndelSummary( ExomeTracker ):

    @property
    def tracks( self ):
        d = self.get( "SELECT DISTINCT track FROM indel_stats")
        return tuple( [x[0] for x in d ] )

    def __call__(self, track, slice = None ):
        statement = '''SELECT indel_length, indel_count FROM indel_stats where track='%(track)s'  order by indel_length;'''
        return self.getAll( statement )

class SharedSummary( ExomeTracker ):

    @property
    def tracks( self ):
        d = self.get( "SELECT DISTINCT track FROM vcf_shared_stats")
        return tuple( [x[0] for x in d ] )

    def __call__(self, track, slice = None ):
        statement = '''SELECT no_samples, var_count FROM vcf_shared_stats where track='%(track)s'  order by no_samples;'''
        return self.getAll( statement )


