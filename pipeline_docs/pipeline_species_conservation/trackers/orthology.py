import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
class trackSummary( SingleTableTrackerRows ):
    table = "genelist_stats"

##################################################################################
class orthologyGroupCounts( TrackerSQL ):

    mPattern = "ortholog_groups_with_feature"

    def __call__(self, track, slice = None ):
        statement = '''SELECT species_count , count(set_id) as genes
                       FROM ortholog_groups_with_feature
                       GROUP BY species_count
                       ORDER BY species_count desc;'''
        return self.getAll( statement )

##################################################################################
class conservedGenesAllSpecies( TrackerSQL ):

    mPattern = "ortholog_groups_with_feature"

    def __call__(self, track, slice = None ):
        statement = '''SELECT set_id, gene_names
                       FROM ortholog_groups_with_feature
                       WHERE species_count=4;'''
        return self.getAll( statement )


