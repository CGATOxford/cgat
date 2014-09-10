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
import cpgReport

from CGATReport.Tracker import *
from CGATReport.odict import OrderedDict as odict

##########################################################################
##########################################################################


class genesWithNMItranscript(cpgReport.cpgTracker):

    ''''''
    mPattern = "_replicated_tss$"

    def __call__(self, track, slice=None):
        query = '''SELECT a.nmi_genes, b.total_genes, round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi
                   FROM (SELECT count(distinct a.gene_id) as nmi_genes
                   FROM %(track)s_replicated_tss t, annotations.transcript_info a
                   WHERE t.closest_dist < 1000
                   AND substr(t.closest_id,1,18)=a.transcript_id
                   or substr(t.closest_id,20,18)=a.transcript_id
                   or substr(t.closest_id,39,18)=a.transcript_id
                   or substr(t.closest_id,58,18)=a.transcript_id
                   or substr(t.closest_id,77,18)=a.transcript_id)
                   AND a.gene_biotype='protein_coding') a,
                   (SELECT count(distinct gene_id) as total_genes 
                   FROM annotations.transcript_info
                   WHERE gene_biotype='protein_coding') b'''
        data = self.getAll(query)
        return data

##########################################################################


class genesWithNMIGene(cpgReport.cpgTracker):

    ''''''
    mPattern = "_replicated_gene_tss$"

    def __call__(self, track, slice=None):
        query = '''SELECT a.nmi_genes, b.total_genes, round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi
                   FROM (SELECT count(distinct t.closest_id) as nmi_genes
                   FROM %(track)s_replicated_gene_tss t, annotations.transcript_info a
                   WHERE t.closest_dist < 1000
                   AND (substr(t.closest_id,1,18)=a.gene_id
                   or substr(t.closest_id,20,18)=a.gene_id
                   or substr(t.closest_id,39,18)=a.gene_id
                   or substr(t.closest_id,58,18)=a.gene_id
                   or substr(t.closest_id,77,18)=a.gene_id)
                   AND a.gene_biotype='protein_coding') a,
                   (SELECT count(distinct gene_id) as total_genes 
                   FROM annotations.transcript_info
                   WHERE gene_biotype='protein_coding') b'''
        data = self.getAll(query)
        return data
