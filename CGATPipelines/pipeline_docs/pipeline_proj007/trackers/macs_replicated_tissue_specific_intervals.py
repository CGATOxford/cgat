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


class liverVsTestes(cpgReport.cpgTracker):

    '''Liver versus tetstes replicated capseq intervals'''
    mPattern = "liver_testes_venn$"

    def __call__(self, track, slice=None):
        query = '''SELECT category, intervals FROM liver_testes_venn'''
        data = self.getAll(query)
        return data

##########################################################################


class liverVsTestesTSS(cpgReport.cpgTracker):

    '''Liver versus tetstes replicated capseq intervals'''
    mPattern = "liver_testes_tss_venn$"

    def __call__(self, track, slice=None):
        query = '''SELECT category, intervals FROM liver_testes_tss_venn'''
        data = self.getAll(query)
        return data

##########################################################################


class liverVsTestesIntergenic(cpgReport.cpgTracker):

    '''Liver versus tetstes replicated capseq intervals'''
    mPattern = "liver_testes_intergenic_venn$"

    def __call__(self, track, slice=None):
        query = '''SELECT category, intervals FROM liver_testes_intergenic_venn'''
        data = self.getAll(query)
        return data
