from SphinxReport.Tracker import *
import sqlite3
import collections


class LcaContributingReads(TrackerSQL):

    '''
    returns the proportion of reads 
    that contribute to the relative abundance
    estimations. 
    '''

    def __call__(self, track, slice=None):
        statement = """SELECT * FROM test_control_R1_taxa_readcounts"""
        return self.getAll(statement)


class LcaSpeciesAbundanceDistribution(TrackerSQL):

    '''
    summarise the relative abundance distribution of 
    species
    '''
    pattern = "(.*)_taxa_count"

    def __call__(self, track, slice=None):
        '''
        return species relative abundance distribution
        '''
        species = []
        statement = """SELECT proportion FROM %s_taxa_count
                   WHERE level == 'species'""" % track
        for rel in self.execute(statement).fetchall():
            species.append(float(rel[0]) * 100)
        species.sort()
        species.reverse()
        return {"species": species}


class LcaRelativeAbundance(TrackerSQL):

    '''
    summarises the relative abundance at
    different taxonomic levels - results
    from lcamapper
    '''
    pattern = "(.*)_taxa_count"

    def __call__(self, track, slice=None):
        '''
        only displays those above 1% relative abundance
        '''
        result = collections.OrderedDict()
        taxon_levels = [
            "phylum", "class", "order", "family", "genus", "species"]
        for taxon in taxon_levels:
            result[taxon] = {}
        for taxon in result.keys():
            statement = """SELECT taxa, proportion FROM %s_taxa_count
                           WHERE level == '%s'""" % (track, taxon)
            for tax, rel in self.execute(statement).fetchall():
                if float(rel) > 0.01:
                    result[taxon][tax] = float(rel) * 100
        return result


class LcaTotalTaxa(TrackerSQL):

    '''
    Summarises the total number of groups at
    each level that were detected in the samples
    '''
    pattern = "(.*)_level_count"

    def __call__(self, track, slice=None):
        '''
        number of taxa found at each level
        '''
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        result = collections.OrderedDict()
        statement = """SELECT nphylum, nclass, norder, nfamily, ngenus, nspecies FROM %s_level_count""" % track
        for data in cc.execute(statement).fetchall():
            phylum, _class, order, family, genus, species = map(
                int, [data[0], data[1], data[2], data[3], data[4], data[5]])
            result["phylum"], result["class"], result["order"], result["family"], result[
                "genus"], result["species"] = phylum, _class, order, family, genus, species
        return result
