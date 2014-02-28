import os
import sys
import re
import types
from VariantsReport import *

#####################################################
#####################################################
#####################################################


class TrackerAlleles(VariantsTracker):

    '''default tracker for allele analysis.'''
    pattern = "(.*)_alleles_genes$"


class AlleleCounts(TrackerAlleles):

    '''return counts over alleles.'''

    def __call__(self, track, slice=None):

        data = odict()
        data["total"] = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_alleles" % locals())
        data["wildtype"] = self.getValue(
            "SELECT SUM(is_wildtype) FROM %(track)s_alleles" % locals())
        data["stop"] = self.getValue(
            "SELECT SUM(is_stop_truncated) FROM %(track)s_alleles" % locals())
        data["nmd"] = self.getValue(
            "SELECT SUM(is_nmd_knockout) FROM %(track)s_alleles" % locals())
        data["splice"] = self.getValue(
            "SELECT SUM(is_splice_truncated) FROM %(track)s_alleles" % locals())
        data["variants"] = 2 * data["total"] - sum(data.values())

        return data


class WildtypeCountsPerTranscript(TrackerAlleles):

    '''return number of wildtype alleles per transcipt.'''

    def __call__(self, track, slice=None):

        data = odict()

        for x in range(0, 3, 1):
            data[str(x)] = len(self.getValues("""SELECT COUNT(*)
                                            FROM %(track)s_alleles 
                                            GROUP BY transcript_id 
                                            HAVING SUM(is_wildtype) = %(x)i
                                            """ % locals()))
        return data


class TranscriptsCounts(TrackerAlleles):

    '''return counts over transcripts. Counting is strict,
    i.e. a transcript is knocked out only, if both its
    alleles are knocked out.
    '''
    mPattern = "_alleles_transcripts$"
    suffix = "_alleles_transcripts"

    def __call__(self, track, slice=None):

        table = track + self.suffix

        data = odict()
        data["total"] = self.getValue(
            "SELECT COUNT(*) FROM %(table)s" % locals())
        data["wildtype"] = self.getValue("""SELECT COUNT(*) 
                                            FROM %(table)s
                                            WHERE is_wildtype""" % locals())
        data["truncated"] = self.getValue("""SELECT COUNT(*) 
                                            FROM %(table)s
                                            WHERE is_truncated""" % locals())
        data["knockout"] = self.getValue("""SELECT COUNT(*) 
                                            FROM %(table)s
                                            WHERE is_knockout""" % locals())
        data["partially disrupted"] = self.getValue("""SELECT COUNT(*) 
                                            FROM %(table)s
                                            WHERE is_affected""" % locals())
        data["variants"] = 2 * data["total"] - sum(data.values())
        return data


class GenesCounts(TranscriptsCounts):

    '''return counts over transcripts. Here,
    a transcript is counted, if flags are set in
    both alleles.
    '''
    mPattern = "_alleles_genes$"
    suffix = "_alleles_genes"


class GenesNMDKnockouts(TrackerAlleles):

    '''return a list of all genes that have been knocked out
    due to NMD
    '''

    def __call__(self, track, slice=None):

        headers = ("gene_id",
                   "gene_name",
                   "ntranscripts",
                   "contig",
                   "strand",
                   "stops-start",
                   "stops-end")

        data = self.get( '''SELECT DISTINCT g.gene_id, 
                                               i.gene_name, 
                                               g.ntranscripts,
                                               g.contig, g.strand,
                                               g.stop_codons_start,
                                               g.stop_codons_end
                               FROM %(track)s_alleles_genes as g,
                                       annotations.transcript_info AS i
                               WHERE g.gene_id = i.gene_id AND g.is_nmd_knockout''' % self.members(locals()))

        return odict(zip(headers, zip(*data)))


class GenesNMDKnockoutsOverview(VariantsTracker):

    '''return a list of all genes that have been knocked out
    due to NMD. Do not list genes without knockouts.
    '''
    mPattern = "summary_genes_alleles_"

    def __call__(self, track, slice=None):

        columns = set(self.getColumns("summary_alleles_genes_%s" % track))
        columns.remove("gene_id")
        columns.remove("total")
        columns = list(columns)
        columns.sort()
        columns = ["gene_id", "total"] + columns
        cols = ",".join(columns)
        data = self.get( '''SELECT %(cols)s 
                               FROM summary_alleles_genes_%(track)s WHERE total > 0''' % self.members(locals()) )
        d = zip(columns, zip(*data))
        return odict(d)


class GenesSingleExonKnockouts(TrackerAlleles):

    '''return a list of all single genes that have been 
    knocked out. 

    A single exon gene is reported, if it is truncated
    at less that 40 amino acids of tranlated sequence
    (the minimum domain size) are left.
    '''

    # maximum length (in bases) on truncated genes
    min_length = 120

    def __call__(self, track, slice=None):

        headers = ("gene_id",
                   "gene_name",
                   "ntranscripts",
                   "cds_len"
                   "original_cds_len",
                   "contig",
                   "strand",
                   "stops-start",
                   "stops-end")

        data = self.get( '''SELECT DISTINCT g.gene_id, 
                                               i.gene_name, 
                                               g.ntranscripts,
                                               a.cds_len, 
                                               a.cds_original_len,
                                               g.contig, g.strand,
                                               g.stop_codons_start,
                                               g.stop_codons_end
                               FROM %(track)s_alleles_genes as g,
                                       %(track)s_alleles as a,
                                       gene_stats as s,
                                       annotations.transcript_info AS i
                               WHERE g.gene_id = i.gene_id AND 
                                     g.gene_id = s.gene_id AND
                                     g.gene_id = a.gene_id AND
                                     a.cds_len < %(min_length)i AND
                                     s.nval = 1 AND
                                     g.is_truncated''' % self.members(locals()))

        return odict(zip(headers, zip(*data)))


class GenesNMDKnockoutsSummary(VariantsTracker):

    '''summarize genes knocked out.
    '''
    mPattern = "_alleles_genes$"

    def getTracks(self, subset=None):
        return ["view_genes"]

    def __call__(self, track, slice=None):

        tracks = VariantsTracker.getTracks(self)

        columns = ["%s_nmd_knockout" % t for t in tracks]
        info = []

        fields = ",".join(columns + info)

        statement = '''
        SELECT DISTINCT gene_id, gene_name, nmd_knockout_total, %(fields)s FROM %(track)s WHERE nmd_knockout_total > 0
        ''' % locals()

        data = self.get(statement)

        if len(data) == 0:
            return None

        d = odict(
            zip(["gene_id", "gene_name", "nmd_knockout_total"] + columns + info, zip(*data)))

        c = []
        for x in range(len(d["gene_id"])):
            s = []
            for t in tracks:
                if d["%s_nmd_knockout" % t][x] != 0:
                    s.append(t)
            c.append(",".join(s))

        for t in tracks:
            del d["%s_nmd_knockout" % t]
        d["tracks"] = c

        return d


class GenesNMDKnockoutsWithOMIM(VariantsTracker):

    '''return a list of all genes that have been knocked out
    due to NMD. Do not list genes without knockouts.
    '''
    mPattern = "_alleles_genes$"
    url_omim = "http://www.ncbi.nlm.nih.gov/omim/%s"

    def getTracks(self, subset=None):
        return ["view_genes", ]

    def getSlices(self, subset=None):
        if subset:
            return subset
        return ["omim_gene", "omim_phenotype"]

    def __call__(self, track, slice=None):

        tracks = VariantsTracker.getTracks(self)

        columns = ["%s_nmd_knockout" % t for t in tracks]

        info = ["hs_gene_id", "hs_ds",
                "omim_gene_id", "omim_morbid_id", "omim_description"]

        fields = ",".join(columns + info)

        if slice == "omim_gene":
            statement = '''
        SELECT gene_id, gene_name, nmd_knockout_total, %(fields)s FROM %(track)s WHERE nmd_knockout_total > 0 AND omim_gene_id != ''
        ''' % locals()
        elif slice == "omim_phenotype":
            statement = '''
        SELECT gene_id, gene_name, nmd_knockout_total, %(fields)s FROM %(track)s WHERE nmd_knockout_total > 0 AND omim_morbid_id != ''
        ''' % locals()

        data = self.get(statement)
        if len(data) == 0:
            return None

        d = odict(
            zip(["gene_id", "gene_name", "nmd_knockout_total"] + columns + info, zip(*data)))

        d["omim_gene_id"] = ["`%s <%s>`_" %
                             (x, self.url_omim % x) for x in d["omim_gene_id"]]
        d["omim_morbid_id"] = ["`%s <%s>`_" %
                               (x, self.url_omim % x) for x in d["omim_morbid_id"]]
        c = []
        for x in range(len(d["omim_gene_id"])):
            s = []
            for t in tracks:
                if d["%s_nmd_knockout" % t][x] != 0:
                    s.append(t)
            c.append(",".join(s))
        for t in tracks:
            del d["%s_nmd_knockout" % t]
        d["tracks"] = c

        return d


class GenesSpliceTruncated(TrackerAlleles):

    '''return a list of all genes that have been truncated
    due to the abrogation of a splice site.
    '''

    def __call__(self, track, slice=None):

        headers = ("gene_id",
                   "gene_name",
                   "ntranscripts",
                   "contig",
                   "strand")

        data = self.get( '''SELECT DISTINCT g.gene_id, 
                                               i.gene_name, 
                                               g.ntranscripts,
                                               g.contig, 
                                               g.strand
                               FROM %(track)s_alleles_genes as g,
                                       annotations.transcript_info AS i
                               WHERE g.gene_id = i.gene_id AND g.is_splice_truncated''' % self.members(locals()))

        return odict(zip(headers, zip(*data)))
