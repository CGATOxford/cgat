from GeneSetsReport import *
import itertools
import scipy.stats


class BiasDistribution(GeneSetsTracker):
    '''Average values of bias attribute columns in gene list
    grouped by pathway.

    See GOSeq manuscript (:pmid:`20132535`)

    '''

    def getTracks(self):
        tables = self.getTables()
        genelists = [x for x in tables if x.startswith("genelist_")]
        pathways = [x for x in tables if x.startswith("pathways_")]
        return ["%s:%s" % (x, y) for x, y in
                itertools.product(genelists, pathways)]

    def getSlices(self):
        # extract bias fields from all options that end in "_bias"
        take = [y.split(",") for x, y in P.items() if x.endswith("_bias")]
        return list(set(itertools.chain(*take)))

    def __call__(self, track, slice):

        genelist, pathway = track.split(":")

        return self.getValues("""
        SELECT avg(%(slice)s)
        FROM
        %(genelist)s as g,
        %(pathway)s as p
        WHERE g.gene_id = p.gene_id
        GROUP by go_id""")


class BiasTest(BiasDistribution):
    '''test for bias in pathways.

    See GOSeq manuscript (:pmid:`20132535`)

    The test applies a Wilcoxon signed rank test
    for each pathway to see if it is different
    from the bulk of all pathways.

    .. note::

       The test should only be run on pathways with at
       least 20 genes.

    '''

    def __call__(self, track, slice):

        genelist, pathway = track.split(":")

        dataframe = self.getDataFrame("""
        SELECT go_id, %(slice)s
        FROM
        %(genelist)s as g,
        %(pathway)s as p
        WHERE g.gene_id = p.gene_id""")

        # drop all na values
        dataframe = dataframe.dropna(axis=0)

        # get all data
        all_data = dataframe[slice]

        results = []
        for name, group in dataframe.groupby('go_id'):
            if len(group) < 20:
                continue
            t, pvalue = scipy.stats.ranksums(
                group[slice],
                all_data)
            results.append(pvalue)

        return results
