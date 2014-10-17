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
from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict

##########################################################################
##########################################################################


class genesWithNMItranscript(cpgTracker):

    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_transcript_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT a.nmi_genes, c.cgi_genes, b.total_genes, 
                   round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi, 
                   round((c.cgi_genes+0.0)/b.total_genes, 2) as fraction_cgi
                   FROM 
                   (SELECT count(distinct t.gene_id) as nmi_genes
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_transcript_tss_distance d,
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_transcript_mapping m,
                   annotations.transcript_info t
                   WHERE t.gene_biotype='protein_coding'
                   AND d.closest_dist < 1000
                   AND t.transcript_id=m.transcript_id
                   AND m.interval_id=d.gene_id) a,
                   (SELECT count(distinct gene_id) as total_genes 
                   FROM annotations.transcript_info
                   WHERE gene_biotype='protein_coding') b,
                    (SELECT count(distinct t.gene_id) as cgi_genes
                   FROM cgi_%(ANNOTATIONS_NAME)s_transcript_tss_distance d,
                   cgi_%(ANNOTATIONS_NAME)s_interval_transcript_mapping m,
                   annotations.transcript_info t
                   WHERE t.gene_biotype='protein_coding'
                   AND d.closest_dist < 1000
                   AND t.transcript_id=m.transcript_id
                   AND m.interval_id=d.gene_id) c'''
        data = self.getAll(query)
        return data

##########################################################################


class genesWithNMIGene(cpgTracker):

    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_" + ANNOTATIONS_NAME + "_gene_tss_distance$"

    def __call__(self, track, slice=None):
        query = '''SELECT a.nmi_genes, b.total_genes, round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi
                   FROM (
                   SELECT count(distinct t.gene_id) as nmi_genes
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_gene_tss_distance d,
                   %(track)s_replicated_%(ANNOTATIONS_NAME)s_interval_gene_mapping m,
                   annotations.transcript_info t
                   WHERE t.gene_biotype='protein_coding'
                   AND d.closest_dist < 1000
                   AND t.gene_id=m.gene_id
                   AND m.interval_id=d.gene_id) a,
                   (SELECT count(distinct gene_id) as total_genes 
                   FROM annotations.transcript_info
                   WHERE gene_biotype='protein_coding') b'''
        data = self.getAll(query)
        return data
