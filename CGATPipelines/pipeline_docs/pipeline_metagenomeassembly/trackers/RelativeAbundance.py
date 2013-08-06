from SphinxReport.Tracker import *
import sqlite3
import collections

class RelativeAbundance(TrackerSQL):
    '''
    summarises the relative abundance at
    different taxonomic levels - results
    from metaphlan
    '''
    pattern = "(.*)_relab"

    def __call__(self, track, slice = None):

        '''
        only displays those above 1% relative abundance
        '''
        result = {"phylum":{}, "class":{}, "order":{}, "family":{}, "genus":{}, "species":{}}
        for taxon in result.keys():
            statement = """SELECT taxon, rel_abundance FROM %s_relab
                           WHERE taxon_level == '%s' AND rel_abundance > 1""" % (track, taxon)
            for tax, rel in self.execute(statement).fetchall():
                result[taxon][tax] = rel
        return result


class TotalSpecies(TrackerSQL):
    '''
    Summarises the total number of groups at
    each level that were detected in the samples
    '''
    pattern = "(.*)_relab"
    
    def __call__(self, track, slice = None):

        '''
        only displays those above 1% relative abundance
        '''
        taxon_levels = ["phylum", "class", "order", "family", "genus", "species"]
        result = collections.defaultdict(int)
        for taxon in taxon_levels:
            statement = """SELECT taxon FROM %s_relab
                           WHERE taxon_level == '%s'""" % (track, taxon)
            for tax in self.execute(statement).fetchall():
                result[taxon] += 1
        return result
