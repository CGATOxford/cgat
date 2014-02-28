import os
import sys
import re
import types
import itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
from cpgReport import *
from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##########################################################################
##########################################################################


class genesCoveredByNMI(cpgTracker):

    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_genes_capseq_overlap$"

    def __call__(self, track, slice=None):
        query = '''select distinct o.gene_id, i.gene_name, capseq_nover, length, capseq_pover1, capseq_pover2 
                   from %(track)s_replicated_%(ANNOTATIONS_NAME)s_genes_capseq_overlap o, annotations.transcript_info i
                   where capseq_pover1>90
                   and o.gene_id=i.gene_id
                   and o.length > 1000
                   order by length desc '''
        print query
        data = self.getAll(query)
        return data

##########################################################################


class genesWithNMI(cpgTracker):

    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_genes_capseq_overlap$"

    def __call__(self, track, slice=None):
        query = '''select distinct o.gene_id, i.gene_name, capseq_nover, length, capseq_pover1, capseq_pover2 
                   from %(track)s_replicated_%(ANNOTATIONS_NAME)s_genes_capseq_overlap o, annotations.transcript_info i
                   where capseq_pover1 <10
                   and capseq_pover1 >0
                   and o.gene_id=i.gene_id
                   and o.length > 1000
                   and o.length < 15000
                   order by length desc '''
        data = self.getAll(query)
        return data

##########################################################################


class overlappedGenesGOAnalysisBP(cpgTracker):

    '''GO analysis biological process'''
    mPattern = "_overlapped_genes_go_biol_process$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, bcount as genes_with_term, spercent as percent_of_list, ratio as Enrichment, fdr
                   from %(track)s_overlapped_genes_go_biol_process
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class overlappedGenesGOAnalysisCL(cpgTracker):

    '''GO analysis '''
    mPattern = "_overlapped_genes_go_cell_location$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, bcount as genes_with_term, spercent as percent_of_list, ratio as Enrichment, fdr 
                   from %(track)s_overlapped_genes_go_cell_location
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class overlappedGenesGOAnalysisMF(cpgTracker):

    '''GO analysis '''
    mPattern = "_overlapped_genes_go_mol_function$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, bcount as genes_with_term, spercent as percent_of_list, ratio as Enrichment, fdr 
                   from %(track)s_overlapped_genes_go_mol_function
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class overlappedGenesGOSlimAnalysisBP(cpgTracker):

    '''GO slim analysis '''
    mPattern = "_overlapped_genes_goslim_biol_process$"

    def __call__(self, track, slice=None):
        query = '''select distinct goid, description, scount as genes, bcount as genes_with_term, spercent as percent_of_list, ratio as Enrichment, fdr
                   from %(track)s_overlapped_genes_goslim_biol_process
                   where fdr < 0.05
                   order by fdr asc, scount desc '''
        data = self.getAll(query)
        return data

##########################################################################


class overlappedGenesCapseqProfile(TrackerImages):

    """CAPseq profile per gene """

##########################################################################


class overlappedGenesH3K27Profile(TrackerImages):

    """Chromatin profile per gene"""

##########################################################################


class overlappedGenesH3K4Profile(TrackerImages):

    """Chromatin profile per gene"""

##########################################################################


class overlappedGenesH3K27Venn(TrackerImages):

    """intersection of overlapped genes and H3K27Me3 intervals"""

##########################################################################


class overlappedGenesTissueVenn(TrackerImages):

    """Conservation of overlapped genes across tissues"""

##########################################################################


class polycombGAT(cpgTracker):

    """genomic assocation of H3K27Me3 intervals and genes overlapped >90% by NMIs"""
    mPattern = "overlapped_genes_gat_results$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT track, annotation, round(expected,0) as expected, observed, round(fold,1) as fold, pvalue FROM overlapped_genes_gat_results ")
        return odict(zip(("Dataset1", "Dataset2", "Expected overlap", "Observed overlap", "Fold Enrichment", "P-value"), zip(*data)))

##########################################################################


class polycombIntersection(cpgTracker):

    """Intersection of H3K27Me3 intervals and genes overlapped >90% by NMIs"""
    mPattern = "overlapped_genes_h3k27me3_venn$"

    def __call__(self, track, slice=None):
        query = '''select track, chromatin_track, total_merged_intervals, track_and_chromatin_track, track_only, chromatin_track_only,
                   round((0.0+track_and_chromatin_track)/(track_and_chromatin_track+track_only+0.0)*100,2) as percent_track,
                   round((0.0+track_and_chromatin_track)/(track_and_chromatin_track+chromatin_track_only+0.0)*100,2) as percent_chromatin_track
                   from overlapped_genes_h3k27me3_venn'''
        data = self.getAll(query)
        return data
