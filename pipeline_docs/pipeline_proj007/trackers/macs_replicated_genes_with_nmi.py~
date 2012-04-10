import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
from cpgReport import *
from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
##################################################################################
class genesWithNMItranscript(cpgTracker):
    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_transcript_tss_distance$"

    def __call__(self, track, slice = None ):
        query = '''SELECT a.nmi_genes, b.total_genes, round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi
                   FROM (SELECT count(distinct a.gene_id) as nmi_genes
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_transcript_tss_distance t, annotations.transcript_info a
                   WHERE t.closest_dist < 1000
                   AND (substr(t.closest_id,1,18)=a.transcript_id
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

##################################################################################
class genesWithNMIGene(cpgTracker):
    ''''''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance$"

    def __call__(self, track, slice = None ):
        query = '''SELECT a.nmi_genes, b.total_genes, round((a.nmi_genes+0.0)/b.total_genes, 2) as fraction_nmi
                   FROM (SELECT count(distinct t.closest_id) as nmi_genes
                   FROM %(track)s_replicated_%(ANNOTATIONS_NAME)s_gene_tss_distance t, annotations.transcript_info a
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
        
