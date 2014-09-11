import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import scipy.stats
import numpy.ma
import Stats
import Histogram

from CGATReport.Tracker import *
from cpgReport import *

##########################################################################


class IntervalList(cpgTracker):

    '''list of intervals.'''

    nresults = 20
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold")
    mPattern = "_macs_merged_intervals$"

    def getSQLStatement(self, track, slice=None):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY i.avgval DESC''' % locals()

        if self.nresults:
            statement += " LIMIT %i" % self.nresults

        return statement

    def __call__(self, track, slice=None):

        statement = self.getSQLStatement(track, slice)
        data = self.get(statement)
        ucsc_genome = UCSC_GENOME
        n = odict()
        for d in data:
            id, contig, start, end, length = d[:5]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_genome)s&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[str(id)] = odict(
                zip(self.mColumnsFixed + self.mColumnsVariable, (pos, length,) + d[5:]))

        return n

##########################################################################


class IntervalListFull(cpgTracker):

    '''list of all intervals. Table for export. '''

    nresults = None
    mPattern = "_macs_merged_intervals$"

    def __call__(self, track, slice=None):

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY i.peakval DESC''' % locals()

        data = self.get(statement)
        return odict(zip(("contig", "start", "end", "peakval", "avgval"),  zip(*data)))

##########################################################################


class IntervalListPeakval(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY i.peakval DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListAvgval(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY i.avgval DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListFoldChange(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY fold DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListLength(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold
                       FROM %(track)s_macs_merged_intervals AS i
                       ORDER BY length DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class LongIntervals(IntervalList):

    '''list of intervals >10kb in length sorted by fold change'''

    nresults = 1000
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = (
        "avgval", "fold", "gene_id", "gene_name", "genes_pover1", "genes_pover2")
    mPattern = "_replicated_intervals$"

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT distinct i.interval_id, i.contig, i.start, i.end, i.length, i.avgval, i.fold, t.gene_id, t.gene_name, o.genes_pover1, o.genes_pover2
                       FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss s, %(track)s_replicated_ensembl_gene_overlap o, annotations.transcript_info t
                       WHERE (substr(s.closest_id,1,18)=t.transcript_id
                       or substr(s.closest_id,20,18)=t.transcript_id
                       or substr(s.closest_id,39,18)=t.transcript_id
                       or substr(s.closest_id,58,18)=t.transcript_id
                       or substr(s.closest_id,77,18)=t.transcript_id)
                       AND i.interval_id=s.gene_id
                       AND o.gene_id=i.interval_id
                       AND s.is_overlap > 0
                       AND t.gene_biotype='protein_coding'
                       AND i.length > 3000
                       AND o.genes_pover2 > 20
                       ORDER BY i.length desc
                       LIMIT 500''' % locals()
        return statement

##########################################################################


class IntervalListLowGC(IntervalList):

    '''list of intervals with pGC <0.5 and pCpG <0.6 sorted by fold change'''

    mPattern = "_replicated_intervals$"
    nresults = 100
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold", "pGC", "CpG_ObsExp")

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, round(i.avgval,2), i.fold, 
                       round(c.pGC,2), round(c.CpG_ObsExp,2)
                       FROM %(track)s_replicated_intervals AS i, %(track)s_replicated_composition AS c
                       WHERE i.interval_id=c.gene_id
                       AND c.CpG_ObsExp < 0.6
                       AND c.pGC < 0.5
                       ORDER BY i.fold DESC
                       LIMIT %(nresults)s''' % locals()
        return statement

##########################################################################


class IntervalListCDS(IntervalList):

    '''list of intervals overlapping CDS.'''
    nresults = 100
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval", "fold", "nover_CDS",
                        "pover1_CDS", "pover2_CDS", "closest_id", "gene_id", "gene_name")
    mPattern = "-?_macs_merged_intervals$"

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults

        statement = '''SELECT distinct i.interval_id, i.contig, i.start, i.end, i.length, i.peakval, i.avgval, i.fold, 
                       a.nover_CDS, a.pover1_CDS, a.pover2_CDS, 
                       tss.closest_id, tr.gene_id, tr.gene_name
                       FROM %(track)s_macs_merged_intervals AS i, %(track)s_merged_annotations AS a,
                       %(track)s_merged_tss AS tss, annotations.transcript_info AS tr
                       WHERE i.interval_id=a.gene_id 
                       AND i.interval_id=tss.gene_id
                       AND tr.transcript_id=tss.closest_id
                       AND a.nover_CDS>0
                       ORDER BY a.pover2_CDS DESC, a.pover1_CDS DESC
                       LIMIT %(nresults)s''' % locals()
        return statement
