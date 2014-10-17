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


class replicatedIntervalTranscriptOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding transcripts """
    mPattern = "_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_transcript_extended_pover1"

##########################################################################


class replicatedIntervalGeneOverlap(featureOverlap):

    """return overlap of interval with Ensembl protein-coding genes """
    mPattern = "_" + ANNOTATIONS_NAME + "_overlap$"
    mTable = "_" + ANNOTATIONS_NAME + "_overlap"
    mWhere = "tss_gene_extended_pover1"
