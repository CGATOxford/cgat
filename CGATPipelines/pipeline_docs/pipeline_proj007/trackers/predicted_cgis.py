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


class cgiTranscriptOverlap(featureOverlap):
    mPattern = "cgi_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "cgi_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_transcript_extended_pover1"

##########################################################################


class cgiGeneOverlap(featureOverlap):
    mPattern = "cgi_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "cgi_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_gene_extended_pover1"

##########################################################################


class cgiTranscriptTSSOverlap(cpgTracker):

    """overlap of predicted CGIs with TSS """
    mPattern = "cgi_" + ANNOTATIONS_NAME + "_transcript_tss_venn"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT track, intervals from cgi_" + ANNOTATIONS_NAME + "_transcript_tss_venn")
        return data

##########################################################################


class cgiGeneTSSOverlap(cpgTracker):

    """overlap of predicted CGIs with TSS """
    mPattern = "cgi_" + ANNOTATIONS_NAME + "_gene_tss_venn"

    def __call__(self, track, slice=None):
        data = self.getAll(
            "SELECT track, intervals from cgi_" + ANNOTATIONS_NAME + "_gene_tss_venn")
        return data

##########################################################################


class CGI_CpGObsExp(cpgTracker):

    """CpG Obs/Exp of predicted CGI and TSS intervals """
    mPattern = "cgi_comp$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT CpG_ObsExp FROM cgi_comp" % locals())
        return data

##########################################################################


class CGI_GCContent(cpgTracker):

    """GC content of proedicted CGI and TSS intervals """
    mPattern = "cgi_comp$"

    def __call__(self, track, slice=None):
        data = self.getAll("SELECT pGC FROM cgi_comp" % locals())
        return data
