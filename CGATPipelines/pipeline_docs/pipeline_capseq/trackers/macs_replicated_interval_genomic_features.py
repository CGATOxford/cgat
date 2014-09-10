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


class replicatedIntervalEnsemblTranscriptOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding transcripts """
    mPattern = "_ensembl_transcript_overlap$"
    mTable = "_ensembl_transcript_overlap"
    mWhere = "tss_transcript_extended_pover1"

##########################################################################


class replicatedIntervalEnsemblGeneOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding genes """
    mPattern = "_ensembl_gene_overlap$"
    mTable = "_ensembl_gene_overlap"
    mWhere = "tss_gene_extended_pover1"
