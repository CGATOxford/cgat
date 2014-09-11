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

from CGATReport.Tracker import *
from cpgReport import *
from CGATReport.odict import OrderedDict as odict

##########################################################################


class cgiEnsemblTranscriptOverlap(featureOverlap):
    mPattern = "cgi_ensembl_transcript_overlap$"
    mTable = "cgi_ensembl_transcript_overlap"
    mWhere = "tss_transcript_extended_pover1"

##########################################################################


class cgiEnsemblGeneOverlap(featureOverlap):
    mPattern = "cgi_ensembl_gene_overlap$"
    mTable = "cgi_ensembl_gene_overlap"
    mWhere = "tss_gene_extended_pover1"

##########################################################################


class cgiEnsemblTranscriptTSSOverlap(cpgTracker):

    """overlap of predicted CGIs with TSS """
    mPattern = "cgi_ensembl_transcript_tss_venn"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT track, intervals from cgi_ensembl_transcript_tss_venn")
        return data

##########################################################################


class cgiEnsemblGeneTSSOverlap(cpgTracker):

    """overlap of predicted CGIs with TSS """
    mPattern = "cgi_ensembl_gene_tss_venn"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT track, intervals from cgi_ensembl_gene_tss_venn")
        return data

##########################################################################


class CGI_CpGObsExp(cpgTracker):

    """CpG Obs/Exp of predicted CGI and TSS intervals """
    mPattern = "_comp$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT CpG_ObsExp FROM %(track)s_comp" % locals())
        return data

##########################################################################


class CGI_GCContent(cpgTracker):

    """GC content of proedicted CGI and TSS intervals """
    mPattern = "_comp$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT pGC FROM %(track)s_comp" % locals())
        return data
