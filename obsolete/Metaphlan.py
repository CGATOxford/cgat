from SphinxReport.Tracker import *
import sqlite3
import collections


class ContributingReads(TrackerSQL):

    '''
    returns the proportion of reads 
    that contribute to the relative abundance
    estimations. Metaphlan appears to ommit low
    abundance reads. This tracker will count
    all reads that have an assignment regardless of
    whether they contributed to the relab estimate.
    It is therfore an overestimation of contributing 
    reads.
    '''
    pattern = "(.*)_readmap"

    def __call__(self, track, slice=None):

        if len(track.split("_")) == 4:
            dtrack = track.split("_")
            dtrack = dtrack[0] + "-" + dtrack[1] + \
                "_" + dtrack[2] + "-" + dtrack[3]
        else:
            dtrack = track.replace("_", "-")
        total_stmt = """SELECT total_reads FROM reads_summary WHERE track == '%s'""" % dtrack
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        total = cc.execute(total_stmt).fetchone()[0]
        # returns the number at the phylum level i.e max possible
        statement = """SELECT count(*) FROM %s_readmap""" % track

        return {"pct_reads": (float(self.execute(statement).fetchone()[0]) / int(total)) * 100,
                "n_reads": float(self.execute(statement).fetchone()[0])}


class SpeciesAbundanceDistribution(TrackerSQL):

    '''
    plot the distribution of SPECIES relative abundances
    '''
    pattern = "(.*)_relab"

    def __call__(self, track, slice=None):
        '''
        display the rel abundance distribution
        '''
        result = {"species": []}
        statement = """SELECT rel_abundance FROM %s_relab
                       WHERE taxon_level == 'species'""" % track
        for rel in self.execute(statement).fetchall():
            result["species"].append(float(rel[0]))
        return result


class RelativeAbundance(TrackerSQL):

    '''
    summarises the relative abundance at
    different taxonomic levels - results
    from metaphlan
    '''
    pattern = "(.*)_relab"

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
            statement = """SELECT taxon, rel_abundance FROM %s_relab
                           WHERE taxon_level == '%s' AND rel_abundance > 1""" % (track, taxon)
            for tax, rel in self.execute(statement).fetchall():
                result[taxon][tax] = float(rel)
        return result


class TotalTaxa(TrackerSQL):

    '''
    Summarises the total number of groups at
    each level that were detected in the samples
    '''
    pattern = "(.*)_relab"

    def __call__(self, track, slice=None):
        '''
        only displays those above 1% relative abundance
        '''
        taxon_levels = [
            "phylum", "class", "order", "family", "genus", "species"]
        result = collections.OrderedDict()
        for taxon in taxon_levels:
            result[taxon] = 0
            statement = """SELECT taxon FROM %s_relab
                           WHERE taxon_level == '%s'""" % (track, taxon)
            for tax in self.execute(statement).fetchall():
                result[taxon] += 1
        return result
