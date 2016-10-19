#########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
Expression.py - wrap various differential expression tools
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides tools for differential expression analysis
for a variety of methods.

Methods implemented are:

   DESeq
   EdgeR
   ttest

The aim of this module is to run these individual tools and
output a table in a common format.

Usage
-----

Documentation
-------------

Requirements:

* DESeq >= 1.17
* DESeq2 >= 1.5.62
* edgeR >= 3.7.16
* gplots >= 2.14.2
* ggplot2 >= 1.0.0
* reshape >= 0.8.5
* RColorBrewer >= 1.0.5
* grid >= 3.1.1
* limma >= 3.21.18
* samr >= 2.0 (optional)
* siggenes >= 1.39.0 (optional)

Code
----

To do:
--check contrasts against design model

'''

import math
import numpy
import sys
import collections
import itertools
import re
import pandas
import ggplot
import copy
import numpy as np
from scipy.stats import ttest_ind
import matplotlib
import matplotlib.pyplot as plt
import rpy2
from rpy2.robjects import r as R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import os

try:
    import CGAT.Experiment as E
    import CGAT.IOTools as IOTools
    import CGAT.Stats as Stats
except ImportError:
    import Experiment as E
    import IOTools
    import Stats

import CGAT.Expression as Expression
import CGAT.Counts as Counts

# activate pandas/rpy conversion
pandas2ri.activate()


class buildDifferentialExpressionRscript(object):
    # note other things to be included - thresholds for IHW and FDR 
                                        # plots 
                                            # dispersion
    # ma plot
                                            # volcano plot
                                            # counts
                                            # heatmaps 
    ''' '''

    def __init__(self, counts_infile, design_infile, full_model, raw_outfile,
                 contrasts, relevel=None, ihw=False):
        self.counts_inf = counts_infile
        self.design_inf = design_infile
        self.full_model = full_model
        self.relevel = relevel
        self.IHW = ihw
        self.raw_outfile = raw_outfile
        self.contrasts = contrasts

    def validate(self):
        design = Expression.ExperimentalDesign(self.design_inf)
        counts = Counts.Counts(self.counts_inf)
        design.validate(counts, self.full_model)

    def loadLibraries(self):
        return ""

    def readDesign(self):
        return ""

    def readCounts(self):
        return ""

    def differentialExpression(self):
        return ""

    def plotResults(self):
        #def plotMA(self):
            #return ""

        #def plotVolcano(self):
            #return ""

        #def plotPValHistogram(self):
            #return ""

        #def plotAdjPValHistogram(self):
            #return ""
        return  ''


    def buildRscript(self):
        Rstatement = []

        Rstatement.append(self.loadLibraries())

        Rstatement.append(self.readCounts())

        Rstatement.append(self.readDesign())

        Rstatement.append(self.differentialExpression())

        Rstatement.append(self.plotResults())

        Rstatement = "\n".join(Rstatement)

        return Rstatement


class buildDESeq2Rscript(buildDifferentialExpressionRscript):

    def loadLibraries(self):
        load_statement = ""
        libraries = ["DESeq2", "ggplot2", "reshape2", "IHW"]
        for library in libraries:
            load_statement += "suppressMessages(library(%s))\n" % library
        load_statement += "print('finished loading')"
        return load_statement

    def readDesign(self):
        statement = []
        statement.append('''
        #' read in design table

        design = read.table(
        "%s", sep="\\t", header=T, row.names=1)\n''' % self.design_inf)

        statement.append('''
        #' make all design columns factors

        for(column in colnames(design)) {
            design[[column]] = factor(design[[column]])
        }''')

        return "\n".join(statement)

    def readCounts(self):
        statement = []
        statement.append('''
        #' read in counts table

        counts = read.table(
        "%s", sep="\\t", header=T, row.names=1)\n''' % self.counts_inf)

        return "\n".join(statement)

    def differentialExpression(self):
        statement = []

        statement.append('''
        #' creating the DESeq2 object

        dds <- suppressMessages(DESeqDataSetFromMatrix(
        countData= counts,
        colData = design,
        design = %s))
        ''' % self.full_model)

        if self.relevel:
            statement.append("#' set the base group for comparisons")
            # how to make this P3 & P2 compatible?
            for key, value in self.relevel.iteritems():
                statement.append('''dds$%s <- relevel(dds$%s, ref="%s")
                ''' % (key, key, value))

        statement.append("\n#' run differential expression testing")
        statement.append("#' note this is a simple DE pairwise comparision")
        statement.append("#' using wald test - you might want to try a LRT if
        statement.append("#' you have a model with many variables and/or your main 
        statement.append("#' variable has more then 2 factors")
        statement.append("dds <- DESeq(dds)")

        if self.IHW:
            statement.append("\n#' calculate IHW")
            statement.append('''
            results = as.data.frame(results(dds, contrast=list(%s),
            filterFun=ihw))''' % ('"' + '",'.join(self.contrasts) + '"'))

        else:
            statement.append("\n#' these results were generated without IHW")
            statement.append('''# see vignette to see how to get IHW results and
            # why you might want to try it
            results = as.data.frame(results(dds, contrast=list(%s)
            ))''' % ('"' + '",'.join(self.contrasts) + '"'))

        statement.append("\n#' write out raw differential expression results")
        statement.append('''
        write.table(results, file="%s", sep="\\t")''' % self.raw_outfile)

        return "\n".join(statement)


    def plotResults(self):

       statement = []
       statement.append("\n#' ## create deseq2 dispersion plot")
       statement.append("plotDispEsts(dds)")

       statement.append("\n#' ##histogram of p values")
       statement.append('\nhist(results$pvalue[results$baseMean > 1],breaks=0:20/20, col="grey50",border="white")')
       statement.append('\nhist(results$padj[results$baseMean > 1],breaks=0:20/20, col="grey50", border="white")')

       return "\n".join(statement)
