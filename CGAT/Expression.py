'''
Expression.py - wrap various differential expression tools
===========================================================

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
from rpy2.robjects import r
from rpy2.robjects import r as R
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from rpy2.rinterface import RRuntimeError
import os

try:
    import CGAT.Experiment as E
    import CGAT.IOTools as IOTools
    import CGAT.Stats as Stats
except ImportError:
    import Experiment as E
    import IOTools
    import Stats


# activate pandas/rpy conversion
pandas2ri.activate()

# AH: Only do this on demand, module might not be
#     be able to be imported if there are any issues.
# grdevices = importr('grDevices')


def runDETest(raw_DataFrame,
              design_file,
              outfile,
              de_caller,
              **kwargs):
    ''' provide higher level API to run tools with default setting '''
    if de_caller.lower() == "deseq":
        pass
    else:
        raise ValueError("Unknown caller")


def splitModel(model):
    '''returns the terms in the model'''
    return [x for x in
            re.split("[\.:,~+\s*]", re.sub("~(\s*)?", "", model)) if
            len(x) > 0]


def adjustPvalues(p_values):
    '''return a list of BH adjusted pvalues'''
    # import r stats module to adjust pvalues
    stats = importr('stats')
    adj_pvalues = list(stats.p_adjust(FloatVector(p_values), method='BH'))
    return adj_pvalues


def pvaluesToSignficant(p_values, fdr):
    '''return a list of bools for significance'''
    return [int(x < fdr) for x in p_values]


class ExperimentalDesign(object):
    """Objects representing experimental designs.

    This class takes an experimental design in tabular
    form and exports several convenience functions and
    attributes.

    `filename_or_table` can be a filename of a tab-separated table
    with the following columns.

    track
       the sample name

    include
       whether or not this sample should be included
       in the design

    groups
       a label grouping several samples into a group

    pair
       for paired tests, a numeric identifier linking
       samples that are paired.

    An example of an experimtal design with two groups and paired
    samples is below::

        track   include group     pair
        sample1 1       treatment 1
        sample2 1       treatment 2
        sample3 1       control   1
        sample4 1       control   2

    When providing `filename_or_table`, the `include` column is used
    to directly filter the design to remove any non-included samples.

    Additional columns will be added as factors to the design.

    Alternatively, `filename_or_table` can be a pandas DataFrame with
    sample names as row index and the appropriate columns.

    Attributes
    -----------

    table : pandas DataFrame
       dataframe object describing the design
    groups : list
       list of groups in the design
    conditions : list
       group for each sample
    pairs : list
       pair status for each sample
    samples: list
       sample names
    factors: list
       factors for each sample
    has_replicates : bool
       True if at least one group has multiple samples
    has_pairs : bool
       True if design is a paired design

    """

    def __init__(self, filename_or_table):
        # read in table in the constructor for ExpDesign
        # e.g design = ExpDesign(pd.read_csv(...))

        if isinstance(filename_or_table, str):
            self.table = pandas.read_csv(filename_or_table, sep="\t",
                                         index_col=0)
        elif isinstance(filename_or_table, pandas.core.frame.DataFrame):
            self.table = filename_or_table
        else:
            raise ValueError("Type needs to be string or pandas data frame."
                             "Type = %s", type(filename_or_table))

        assert self.table.shape, "design table is empty"

        # parse the design table. Users probably expect this
        # to happen once data is uploaded.
        self._update()

    def _update(self):
        """parse design file and fill class attributes.

        Call this function whenever self.table changes.
        """

        # remove all entries that should not be included
        self.table = self.table[self.table["include"] != 0]

        # define attributes
        self.conditions = self.table['group'].tolist()
        self.pairs = self.table['pair'].tolist()
        # TS - use OrderedDict to retain order in unique
        self.groups = (list(collections.OrderedDict.fromkeys(
            self.conditions)))
        self.samples = self.table.index.tolist()

        # Test if replicates exist, i.e at least one group has multiple samples
        # TS - does this need to be extended to check whether replicates exist
        # for each group?
        max_per_group = max([self.conditions.count(x) for x in self.groups])
        self.has_replicates = max_per_group >= 2

        # Test if pairs exist:
        npairs = len(set(self.pairs))
        has_pairs = npairs == 2

        # ..if so, at least two samples are required per pair
        if has_pairs:
            min_per_pair = min([self.pairs.count(x) for x in set(self.pairs)])
            self.has_pairs = min_per_pair >= 2
        else:
            self.has_pairs = False

        # all columns except "include" may be considered as factors
        self.factors = self.table.drop(["include"], axis=1)
        # remove "pair" from factor if design does not include pairs
        if not self.has_pairs:
            self.factors.drop("pair", inplace=True, axis=1)

    def validate(self, counts=None, model=None):

        if counts is not None:
            missing = set(self.samples).difference(set(counts.table.columns))
            if len(missing) > 0:
                raise ValueError(
                    "following samples in design table are missing"
                    " from counts table: %s" % ", ".join(missing))

        if model is not None:
            # check all model terms exist

            model_terms = splitModel(model)

            missing = set(model_terms).difference(
                set(self.table.columns.tolist()))

            if len(missing) > 0:
                raise ValueError("following terms in the model are missing"
                                 " from the design table: %s" %
                                 ", ".join(missing))

            # check there are at least two values for each level
            for term in model_terms:
                levels = set(self.table.ix[:, term])
                if len(levels) < 2:
                    raise ValueError("term '%s' in the model has less "
                                     "than two levels (%s) in the "
                                     " design table" %
                                     (term, ", ".join(levels)))

    def restrict(self, counts):
        ''' return design with samples not in counts table removed '''

        self.table = self.table.ix[counts.table.columns, :]

    def revalidate(self, counts, model=None):
        ''' re-validate, i.e post filtering of counts table '''

        if len(set(self.samples).symmetric_difference(
                set(counts.table.columns))) > 0:
            self.restrict(counts)
            self._update()
            self.validate(counts, model)

        else:
            pass

    def firstPairOnly(self):
        '''restrict the design table to the first pair only.

        If unpaired will retain whole design table
        '''

        if not self.pairs:
            self.pairs = self.table['pair'].tolist()
        self.table = self.table.ix[self.table['pair'] == min(self.pairs), ]

    def getSamplesInGroup(self, group):
        """return list of sample names belonging to group."""
        if group not in self.groups:
            raise KeyError("group '%s' not present")
        return self.table[self.table["group"] == group].index.tolist()

    def getGroupForSample(self, sample):
        """return group a sample belongs to"""
        return self.table.loc[sample]["group"]

    def getGroups2Samples(self):
        """return a dictionary mapping a group to samples within the group.

        Returns
        -------
        dict
            with groups as keys and list of samples within a group as values.

        """
        groups_to_tracks = {}

        for group in self.groups:
            match_group = (self.table['group'] == group).tolist()
            subset = self.table.iloc[match_group, ]
            groups_to_tracks[group] = subset.index.tolist()
        return groups_to_tracks

    def mapGroupsSuffix(self, shuffle_suffix, keep_suffix):
        '''use suffixes supplied to extract groups from the
        design table and return dictionaries mapping each group to tracks
        for keeping with tracks which should be shuffled
        '''

        groups_to_keep_tracks = {}
        groups_to_spike_tracks = {}

        keep_suffix = keep_suffix.split(",")

        for group in self.groups:
            match_group = (self.table['group'] == group).tolist()
            tmp_design = self.table.iloc[match_group, ]

            groups_to_spike_tracks[group] = [
                x + shuffle_suffix for x in tmp_design.index.tolist()]

            groups_to_keep_tracks[group] = copy.copy(
                groups_to_spike_tracks[group])
            groups_to_keep_tracks[group].extend(
                [x + y for x in tmp_design.index.tolist() for y in keep_suffix])

        return groups_to_keep_tracks, groups_to_spike_tracks


class DEExperiment(object):
    ''' base clase for DE experiments '''

    def __init__(self):
        pass

    def run(self):
        ''' Custom DE functions
        return a DEResult object'''


class DEResult(object):
    ''' base class for DE result '''

    def __init__(self, testTable=None):
        self.table = testTable

    def getResults(self):
        ''' post-process results into generic output
        columns are:
        - contrast
        - treatment_name
        - control_name
        - test_id
        - control_mean
        - treatment_mean
        - control_std
        - treatment_std
        - p_value
        - p_value_adj
        - significant
        - l2fold
        - transformed_l2fold
        - fold
        - status
        '''

    def calculateIHW(self, alpha=0.1):
        ''' Use the Independent Hypothesis Weighting method from
        IGNATIADIS et al (2016) to perform weighted FDR'''

        if not ('control_mean' in self.table.columns and
                'treatment_mean' in self.table.columns and
                'p_value' in self.table.columns):
            E.error("IHW requires control_mean, treatment_mean and p_value "
                    "columns, have you first run the getResults method?")

        runIHW = r('''function(df){
        library(IHW)
        mean_expression = (df$control_mean + df$treatment_mean)/2
        ihw_res = ihw(df$p_value ~ mean_expression,  alpha = %(alpha)s)
        df$p_value_adj = adj_pvalues(ihw_res)
        return(df)
        }''' % locals())

        self.table = pandas2ri.ri2py(runIHW(pandas2ri.py2ri(self.table)))
        self.table["significant"] = pvaluesToSignficant(
            self.table["p_value_adj"], alpha)

    def summariseDEResults(self):
        ''' summarise DE results. Counts instances of possible outcomes'''

        # TS: the summarising is now split by the comparison being made and a
        # dict returned with keys=comparisons, value=E.Counter per comparison

        self.Summary = {}

        control_names = set(self.table['control_name'])
        treatment_names = set(self.table['treatment_name'])

        for control, treatment in itertools.product(control_names,
                                                    treatment_names):

            tmp_table = self.table[self.table['control_name'] == control]
            tmp_table = tmp_table[tmp_table['treatment_name'] == treatment]
            tmp_table.reset_index(inplace=True)

            # check control, treatment combination exists
            n_rows = tmp_table.shape[0]
            if n_rows > 0:

                if control != treatment:
                    label = control + "_" + treatment
                else:
                    label = control
                label = re.sub(":", "_int_", label)

                counts = E.Counter()

                counts.signficant = sum(tmp_table['significant'])
                counts.insignficant = (len(tmp_table['significant']) -
                                       counts.signficant)
                counts.all_over = sum([x > 0 for x in tmp_table['l2fold']])
                counts.all_under = sum([x < 0 for x in tmp_table['l2fold']])
                counts.signficant_over = sum(
                    [tmp_table['significant'][x] == 1 and
                     tmp_table['l2fold'][x] > 0 for x in range(0, n_rows)])
                counts.signficant_under = sum(
                    [tmp_table['significant'][x] == 1 and
                     tmp_table['l2fold'][x] < 0 for x in range(0, n_rows)])

                self.Summary[label] = counts

    def plotMA(self, contrast=None, outfile_prefix=None,
               point_alpha=1, point_size=1, R=None):
        ''' base function for making a MA plot '''

        if not R:
            R = rpy2.robjects.r

        ro.globalenv['tmp_df'] = pandas2ri.py2ri(self.table)

        R('''

        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))

        l_txt = element_text(size=20)

        tmp_df = tmp_df[tmp_df$contrast=="%(contrast)s",]
        tmp_df = tmp_df[order(-tmp_df$p_value_adj),]
        p = ggplot(tmp_df, aes(log((control_mean+treatment_mean)/2,2),
                            transformed_l2fold,
                            colour=as.factor(significant))) +

        geom_point(size=%(point_size)f, alpha=%(point_alpha)f) +
        xlab("log2 mean expression") + ylab("log2 fold change")+
        ggtitle("%(contrast)s") +
        scale_colour_manual(name="Significant", values=c("black", "red")) +
        guides(colour = guide_legend(override.aes = list(size=10)))+
        theme_bw() +
        theme(axis.text.x = l_txt, axis.text.y = l_txt,
              axis.title.x = l_txt, axis.title.y = l_txt,
              legend.title = l_txt, legend.text = l_txt,
              title=l_txt, legend.key.size=unit(1, "cm"),
              aspect.ratio=1)

        suppressMessages(
          ggsave(file="%(outfile_prefix)s_%(contrast)s_MA_plot.png",
          width=10, height=10))''' % locals())

    def plotVolcano(self, contrast=None, outfile_prefix=None, R=None):
        ''' base function for Volcano plotting'''

        if not R:
            R = rpy2.robjects.r

        ro.globalenv['tmp_df'] = pandas2ri.py2ri(self.table)

        R('''
        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))

        l_txt = element_text(size=20)

        tmp_df = tmp_df[tmp_df$contrast=="%(contrast)s",]

        p = ggplot(tmp_df, aes(transformed_l2fold, -log(p_value,10),
                   colour=as.factor(significant))) +
        geom_point() + xlab("log2 fold change") + ylab("p-value (-log10)") +
        ggtitle("%(contrast)s") +
        scale_colour_manual(name="Significant", values=c("black", "#619CFF")) +
        guides(colour = guide_legend(override.aes = list(size=10))) +
        theme_bw() +
        theme(axis.text.x = l_txt, axis.text.y = l_txt,
              axis.title.x = l_txt, axis.title.y = l_txt,
              legend.title = l_txt, legend.text = l_txt,
              title=l_txt, legend.key.size=unit(1, "cm"))

        suppressMessages(
          ggsave(file="%(outfile_prefix)s_%(contrast)s_volcano_plot.png",
          width=10, height=10))''' % locals())

    def plotPvalueHist(self, contrast=None, outfile_prefix=None, R=None):
        ''' base function for Volcano plotting'''

        if not R:
            R = rpy2.robjects.r

        ro.globalenv['tmp_df'] = pandas2ri.py2ri(self.table)

        R('''
        suppressMessages(library(ggplot2))
        suppressMessages(library(grid))

        l_txt = element_text(size=20)

        tmp_df = tmp_df[tmp_df$contrast=="%(contrast)s",]

        p = ggplot(tmp_df, aes(p_value)) +
        geom_histogram(fill="dodgerblue4") +
        xlab("p-value") + ylab("count") +
        ggtitle("p-value histogram - %(contrast)s") +
        theme_bw() +
        theme(axis.text.x = l_txt, axis.text.y = l_txt,
              axis.title.x = l_txt, axis.title.y = l_txt,
              title=l_txt)

        suppressMessages(
          ggsave(file="%(outfile_prefix)s_%(contrast)s_p_value_histogram.png",
          width=10, height=10))''' % locals())

    def plotPvalueQQ(self, contrast=None, outfile_prefix=None, R=None):
        ''' base function for Volcano plotting'''

        if not R:
            R = rpy2.robjects.r

        ro.globalenv['tmp_df'] = pandas2ri.py2ri(self.table)

        R('''
        log_obs_pvalues = sort(-log10(tmp_df[['p_value']]))
        uni_pvalues=runif(length(log_obs_pvalues))
        log_uni_pvalues= -log10(uni_pvalues)
        log_uni_pvalues = sort(log_uni_pvalues)

        png(file="%(outfile_prefix)s_%(contrast)s_p_value_qq_plot.png")
        plot(log_uni_pvalues,log_obs_pvalues,
             xlab=expression(Theoretical~~-log[10](italic(p))),
             ylab=expression(Observed~~-log[10](italic(p))),
             main="P-value QQ-plot",
             pch=20)
        abline(0,1)''' % locals())


class DEExperiment_TTest(DEExperiment):
    '''DECaller object to run TTest on counts data'''

    # TS: to do: deal with genes/regions with zero counts

    def run(self, counts, design, normalise=True,
            normalise_method="deseq-size-factors"):

        # TS: normalisation performed here rather than earlier as
        # the method of normalisation is dependent upon the DE test
        if normalise is True:
            counts.normalise(method=normalise_method)

        df_dict = collections.defaultdict(list)

        for combination in itertools.combinations(design.groups, 2):

            control, treatment = combination
            n_rows = counts.table.shape[0]
            df_dict["control_name"].extend((control,)*n_rows)
            df_dict["treatment_name"].extend((treatment,)*n_rows)
            df_dict["test_id"].extend(counts.table.index.tolist())

            # set all status values to "OK"
            df_dict["status"].extend(("OK",)*n_rows)

            # subset counts table for each combination
            c_keep = [x == control for
                      x in design.conditions]
            control_counts = counts.table.iloc[:, c_keep]
            t_keep = [x == treatment for
                      x in design.conditions]
            treatment_counts = counts.table.iloc[:, t_keep]

            c_mean = control_counts.mean(axis=1)
            df_dict["control_mean"].extend(c_mean)
            df_dict["control_std"].extend(control_counts.std(axis=1))

            t_mean = treatment_counts.mean(axis=1)
            df_dict["treatment_mean"].extend(t_mean)
            df_dict["treatment_std"].extend(treatment_counts.std(axis=1))

            t, prob = ttest_ind(control_counts, treatment_counts, axis=1)
            df_dict["p_value"].extend(prob)

        result = DEResult_TTest(testTable=pandas.DataFrame(df_dict))
        result.table.set_index("test_id", inplace=True)

        return result


class DEResult_TTest(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output '''

        # TS - what about zero values?!
        self.table["fold"] = (
            self.table["treatment_mean"] / self.table["control_mean"])

        self.table["p_value_adj"] = adjustPvalues(self.table["p_value"])

        self.table["significant"] = pvaluesToSignficant(
            self.table["p_value_adj"], fdr)

        self.table["l2fold"] = list(numpy.log2(self.table["fold"]))

        # note: the transformed log2 fold change is not transformed for TTest
        self.table["transformed_l2fold"] = self.table["l2fold"]
        self.table["contrast"] = "_vs_".join((self.table['control_name'],
                                              self.table['treatment_name']))


class DEExperiment_edgeR(DEExperiment):
    '''DEExperiment object to run edgeR on counts data

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.
    '''

    def run(self,
            counts,
            design,
            model=None,
            contrast=None,
            outfile_prefix=None,
            ref_group=None,
            fdr=0.1,
            dispersion=None):

        if not design.has_replicates and dispersion is None:
            raise ValueError("no replicates and no dispersion")

        # create r objects
        r_counts = pandas2ri.py2ri(counts.table)
        r_groups = ro.StrVector(design.conditions)
        r_pairs = ro.StrVector(design.pairs)
        r_has_pairs = ro.BoolVector([design.has_pairs])
        r_has_replicates = ro.BoolVector([design.has_replicates])

        if model is not None:
            r_factors_df = pandas2ri.py2ri(design.factors)
        else:
            r_factors_df = ro.default_py2ri(False)

        E.info('running edgeR: groups=%s, replicates=%s, pairs=%s, '
               'additional_factors:%s' %
               (design.groups, design.has_replicates, design.has_pairs,
                design.factors))

        levels = set(design.table[contrast])
        if len(levels) > 2:

            E.warn(
                "There are more than 2 levels for the "
                "contrast specified" "(%s:%s). The log2fold changes in the "
                "results table and MA plots will be for the first two "
                "levels in the contrast. The p-value will be the p-value "
                "for the overall significance of the contrast. Hence, some "
                "genes will have a signficant p-value but 0-fold change "
                "between the first two levels" % (contrast, levels))

        # build DGEList object
        buildDGEList = r('''
        suppressMessages(library('edgeR'))

        function(counts){

        countsTable = DGEList(counts)

        countsTable = calcNormFactors(countsTable)

        return(countsTable)}''' % locals())

        r_countsTable = buildDGEList(r_counts)

        # build design matrix
        buildDesign = r('''
        function(factors_df){

        for (level in colnames(factors_df)){
           factors_df[[level]] <- factor(factors_df[[level]])
        }

        factors_df$%(contrast)s <- relevel(
           factors_df$%(contrast)s, ref="%(ref_group)s")

        design <- model.matrix(%(model)s, data=factors_df)

        return(design)}''' % locals())

        r_design = buildDesign(r_factors_df)

        # fit model
        fitModel = r('''
        function(countsTable, design, has_replicates){

        if (has_replicates[1] == TRUE) {

            # estimate common dispersion
            countsTable = estimateGLMCommonDisp( countsTable, design )

            # estimate trended dispersion
            countsTable <- estimateGLMTrendedDisp( countsTable, design)

            # estimate tagwise dispersion
            countsTable = estimateGLMTagwiseDisp( countsTable, design )

            # fitting model to each tag
            fit = glmFit( countsTable, design ) }

        else {
            # fitting model to each tag
            fit = glmFit(countsTable, design, dispersion=%(dispersion)s) }

        return(fit)}''' % locals())

        r_fit = fitModel(r_countsTable, r_design, r_has_replicates)

        E.info("Conducting likelihood ratio tests")

        lrtTest = r('''
        function(fit, design, factors_df, countsTable){
        suppressMessages(library(reshape2))

        lrt = glmLRT(fit)

        lrt_table = as.data.frame(lrt$table)

        lrt_table$contrast <- "%(contrast)s"

        for (level in colnames(factors_df)){
           factors_df[[level]] <- factor(factors_df[[level]])
        }

        factors_df$%(contrast)s <- relevel(
           factors_df$%(contrast)s, ref="%(ref_group)s")

        contrast_levels = as.vector(levels(factor(factors_df[["%(contrast)s"]])))

        lrt_table$control_name <- contrast_levels[1]
        lrt_table$treatment_name <- contrast_levels[2]

        dt <- decideTestsDGE(lrt, adjust.method="BH", p.value=%(fdr)s)
        isDE <- as.logical(dt)
        DEnames <- rownames(fit)[isDE]

        png(paste0(c("%(outfile_prefix)s", "MA.png"), collapse="_"))
        plotSmear(lrt, de.tags=DEnames, cex=0.35, main="%(contrast)s")
        abline(h=c(-1,1), col="blue")
        dev.off()

        return(lrt_table)}''' % locals())

        r_lrt_table = lrtTest(r_fit, r_design, r_factors_df, r_countsTable)

        result = DEResult_edgeR(testTable=pandas2ri.ri2py(r_lrt_table))

        return result


class DEResult_edgeR(DEResult):

    def getResults(self, fdr, DEtype="GLM"):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = self.table['treatment_name']
        df_dict["control_name"] = self.table['control_name']
        df_dict["contrast"] = self.table['contrast']
        df_dict["test_id"] = self.table.index
        df_dict["control_mean"] = self.table['logCPM']
        df_dict["treatment_mean"] = self.table['logCPM']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['PValue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['PValue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = (self.table['logFC'])

        # TS: the transformed log2 fold change is not transformed!
        df_dict["transformed_l2fold"] = df_dict["l2fold"]

        # TS: check what happens when no fold change is available
        # TS: may need an if/else in list comprehension. Raise E.warn too?
        df_dict["fold"] = [math.pow(2, float(x)) for
                           x in self.table['logFC']]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)


class DEExperiment_DESeq2(DEExperiment):
    '''DEExperiment object to run DESeq2 on counts data'''

    def run(self,
            counts,
            design,
            fdr=0.1,
            fit_type="parametric",
            model=None,
            outfile_prefix=None,
            ref_group=None,
            contrast=None,
            DEtest="Wald",
            R=None):

        if not R:
            R = rpy2.robjects.r

        pandas2ri.activate()

        # R will replace any "-" with "." in rownames.
        # Here, we make sure the design and counts samples are the same
        design.table.index = [x.replace("-", ".") for x in design.table.index]
        design.factors.index = [x.replace("-", ".") for x in design.factors.index]
        counts.table.columns = [x.replace("-", ".") for x in counts.table.columns]

        # create r objects
        ro.globalenv['counts'] = pandas2ri.py2ri(counts.table)
        ro.globalenv['design'] = pandas2ri.py2ri(design.table)
        ro.globalenv['factors_df'] = pandas2ri.py2ri(design.factors)

        model_terms = [x for x in re.split("[\+~ ]+", model)[1:]
                       if x != "0"]

        E.info('running DESeq2: groups=%s, replicates=%s, pairs=%s, '
               'DE test: %s, additional_factors:%s, ' %
               (design.groups, design.has_replicates, design.has_pairs,
                DEtest, design.factors))

        # load DESeq
        R('''suppressMessages(library('DESeq2'))''')

        # build DESeq2 Datasets (dds)
        assert contrast, ("must supply a contrast for wald or LRT "
                          "(for LRT, contrast is used to derive reduced model")
        if DEtest == "wald":
            assert ref_group, "Must supply a ref_group to perform Wald test"

        if ref_group:
            R('''

            for(column in colnames(factors_df)){
                factors_df[[column]] = factor(factors_df[[column]])
            }

            full_model <- formula("%(model)s")

            factors_df$%(contrast)s <- relevel(
               factors_df$%(contrast)s, ref="%(ref_group)s")

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= counts,
                     colData = factors_df,
                     design = full_model))
            ''' % locals())

        else:
            R('''
            for(column in colnames(factors_df)){
              factors_df[[column]] = factor(factors_df[[column]])
            }

            full_model <- formula("%(model)s")

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= counts,
                     colData = factors_df,
                     design = full_model))
            ''' % locals())

        if DEtest == "wald":

            levels = set(design.table[contrast])
            if len(levels) > 2:
                E.warn('''Using Wald test for factor with more than 2
                levels (%s:%s), Consider LRT''' % (contrast, levels))

            contrast = model_terms[-1]
            contrast_levels = set(design.factors[contrast])

            # performDifferentialTesting
            R('''

            dds = suppressMessages(
               DESeq(dds, test="Wald", fitType="%(fit_type)s"))


            contrast_levels = as.vector(levels(dds@colData$%(contrast)s))

            png("%(outfile_prefix)s_dispersion.png")
            plotDispEsts(dds)
            dev.off()

            res = suppressMessages(results(dds))
            png(paste0(c("%(outfile_prefix)s", "MA.png"), collapse="_"))
            plotMA(res, alpha=%(fdr)s)
            dev.off()

            res = as.data.frame(res)

            c = counts(dds, normalized = TRUE)

            res$contrast = "%(contrast)s"
            contrast_levels = levels(dds@colData$%(contrast)s)

            res$control = contrast_levels[1]
            res$treatment = contrast_levels[2]

            res['test_id'] = rownames(res)

            ''' % locals())

            results = pandas2ri.ri2py(ro.globalenv['res'])

        # liklihood ratio test
        # Note that if there are more than 3 levels for the contrast,
        # the results table will include a log2-fold change from the
        # first two levels only, however, MA plots will be generated
        # for each combination of levels
        elif DEtest == "lrt":

            levels = set(design.table[contrast])
            if len(levels) > 2:

                E.warn('''There are more than 2 levels for the
                contrast specified" "(%s:%s). The log2fold changes in the
                results table and MA plots will be for the first two
                levels in the contrast. The p-value will be the p-value
                for the overall significance of the contrast. Hence, some
                genes may have a signficant p-value but ~0-fold change
                between the first two levels''' % (contrast, levels))

            n = 0

            reduced_model = [x for x in model_terms if x != contrast]

            if len(reduced_model) > 0:
                reduced_model = "~" + "+".join(reduced_model)
            else:
                reduced_model = "~1"

            print('''
            ddsLRT <- suppressMessages(
            DESeq(dds, test="LRT", reduced=formula("%(reduced_model)s"),
                  betaPrior=TRUE, fitType="%(fit_type)s"))

            png("%(outfile_prefix)s_dispersion.png")
            plotDispEsts(ddsLRT)
            dev.off()

            contrast_levels = as.vector(levels(dds@colData$%(contrast)s))

            res = suppressMessages(results(ddsLRT, addMLE=TRUE,
                                   contrast=c("%(contrast)s",
                                   contrast_levels[2], contrast_levels[1])))

            png(paste0(c("%(outfile_prefix)s", "MA.png"), collapse="_"))
            plotMA(res, alpha=%(fdr)s)
            dev.off()

            res = as.data.frame(res)
            res$contrast = "%(contrast)s"

            if(length(contrast_levels)==2){
              res$control = contrast_levels[1]
              res$treatment = contrast_levels[2]
              }
            else{
              res$control = "%(contrast)s"
              res$treatment = "%(contrast)s"
              }

            res['test_id'] = rownames(res)
            ''' % locals())

            R('''
            ddsLRT <- suppressMessages(
            DESeq(dds, test="LRT", reduced=formula("%(reduced_model)s"),
                  betaPrior=TRUE, fitType="%(fit_type)s"))

            png("%(outfile_prefix)s_dispersion.png")
            plotDispEsts(ddsLRT)
            dev.off()

            contrast_levels = as.vector(levels(dds@colData$%(contrast)s))

            res = suppressMessages(results(ddsLRT, addMLE=TRUE,
                                   contrast=c("%(contrast)s",
                                   contrast_levels[2], contrast_levels[1])))

            png(paste0(c("%(outfile_prefix)s", "MA.png"), collapse="_"))
            plotMA(res, alpha=%(fdr)s)
            dev.off()

            res = as.data.frame(res)
            res$contrast = "%(contrast)s"

            if(length(contrast_levels)==2) {
               res$control = contrast_levels[1]
               res$treatment = contrast_levels[2]
            } else {
            res$control = "%(contrast)s"
            res$treatment = "%(contrast)s"
            }

            res['test_id'] = rownames(res)
            ''' % locals())

            results = pandas2ri.ri2py(ro.globalenv['res'])

        else:
            raise ValueError("DEtest must be 'wald' or 'lrt'")

        final_result = DEResult_DESeq2(testTable=results)

        return final_result


class DEResult_DESeq2(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = self.table['treatment']
        df_dict["control_name"] = self.table['control']
        df_dict["test_id"] = self.table['test_id']
        df_dict["contrast"] = self.table['contrast']
        df_dict["control_mean"] = self.table['baseMean']
        df_dict["treatment_mean"] = self.table['baseMean']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['pvalue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['pvalue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = self.table['log2FoldChange']

        # Transformed l2fold is the shrunken values
        df_dict["transformed_l2fold"] = self.table['log2FoldChange']

        # TS: check what happens when no fold change is available
        # TS: may need an if/else in list comprehension. Raise E.warn too?
        df_dict["fold"] = [math.pow(2, float(x)) for
                           x in df_dict["l2fold"]]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        # causes errors if multiple instance of same test_id exist, for example
        # if multiple constrasts have been tested
        # self.table.set_index("test_id", inplace=True)


class DEExperiment_DEXSeq(DEExperiment):
    '''DEExperiment object to run DEXSeq on counts data'''

    def run(self,
            design,
            base_dir,
            model=None,
            flattenedfile=None,
            outfile_prefix=None,
            ref_group=None,
            contrast=None,
            fdr=0.1):

        pandas2ri.activate()

        # create r objects
        E.info('running DEXSeq: groups=%s, pairs=%s, replicates=%s, pairs=%s,'
               ' additional_factors:' %
               (design.groups, design.pairs, design.has_replicates,
                design.has_pairs))

        # load DEXSeq
        R('''suppressMessages(library('DEXSeq'))''')

        sampleTable = design.table

        allfiles = [file for file in os.listdir(base_dir)]
        countfiles = []
        for item in list(design.table.index):
            countfiles += [base_dir+"/"+x for x in allfiles if item in x]

        E.info("Processing Samples. Sample table:")
        E.info("%s" % sampleTable)
        buildCountDataSet = R('''
        function(countFiles, gff, sampleTable, model){

        full_model <- formula("%(model)s")

        dxd <- suppressMessages(DEXSeqDataSetFromHTSeq(
                     countFiles,
                     sampleData=sampleTable,
                     flattenedfile=gff,
                     design=full_model))

        contrast_levels = as.vector(levels(dxd@colData$%(contrast)s))

        dxd = estimateSizeFactors(dxd)
        dxd = estimateDispersions(dxd)

        png("%(outfile_prefix)s_dispersion.png")
        plotDispEsts(dxd)
        dev.off()

        dxd = testForDEU(dxd)
        dxd = estimateExonFoldChanges( dxd, fitExpToVar="%(contrast)s")

        result = DEXSeqResults(dxd)
        result = as.data.frame(result)

        result$contrast = "%(contrast)s"
        result$log2FoldChange = result$log2fold

        if(length(contrast_levels)==2) {
           result$control = contrast_levels[1]
           result$treatment = contrast_levels[2]
        } else {
            result$control = "%(contrast)s"
            result$treatment = "%(contrast)s"
        }
        return(result)
        }''' % locals())
        result = pandas2ri.ri2py(
            buildCountDataSet(countfiles, flattenedfile, sampleTable, model))

        result['test_id'] = result.index
        result['contrast'] = contrast
        final_result = DEResult_DEXSeq(result)

        return final_result


class DEResult_DEXSeq(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output '''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = self.table['treatment']
        df_dict["control_name"] = self.table['control']
        df_dict["test_id"] = self.table['test_id']
        df_dict["contrast"] = self.table['contrast']
        df_dict["control_mean"] = self.table['exonBaseMean']
        df_dict["treatment_mean"] = self.table['exonBaseMean']
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['pvalue']
        df_dict["p_value_adj"] = adjustPvalues(self.table['pvalue'])
        df_dict["significant"] = pvaluesToSignficant(
            df_dict["p_value_adj"], fdr)
        df_dict["l2fold"] = ("NA",)*n_rows

        # Transformed l2fold is the shrunken values
        df_dict["transformed_l2fold"] = self.table['log2FoldChange']
        df_dict["fold"] = ("NA",)*n_rows
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        # causes errors if multiple instance of same test_id exist, for example
        # if multiple constrasts have been tested
        # self.table.set_index("test_id", inplace=True)

    def plotMAplot(self, design, outfile_prefix):
        # need to implement DEXSeq specific MA plot
        raise ValueError("MA plotting is not yet implemented for DESeq")


class DEExperiment_Sleuth(DEExperiment):
    '''DEExperiment object to run sleuth on kallisto bootstrap files
    Unlike the other DEExperiment instances, this does not operate on
    a Counts.Counts object but instead reads the bootstrap hd5 files
    from kallisto into memory in R and then performs the differential
    testing

    The run method expects all kallisto abundance.h5 files to be under
    a single directory with a subdirectory for each sample

    Note: LRT does not generate fold change estimates (see DEResult_Sleuth)

    use dummy_run = True if you don't want to perform differential
    testing but want the counts/tpm outfiles

    '''

    def run(self,
            design,
            base_dir,
            model=None,
            contrast=None,
            outfile_prefix=None,
            counts=None,
            tpm=None,
            fdr=0.1,
            DE_test="wald",
            reduced_model=None,
            dummy_run=False,
            genewise=False,
            gene_biomart=None,
            ref_group=None):

        if DE_test == "lrt":
            E.info("Note: LRT will not generate fold changes")
            assert reduced_model is not None, ("need to provide a reduced "
                                               "model to use LRT")

        # Design table needs a "sample" column
        design.table['sample'] = design.table.index
        r_design_df = pandas2ri.py2ri(design.table)

        E.info('running sleuth: groups=%s, pairs=%s, replicates=%s, pairs=%s, '
               'additional_factors:' %
               (design.groups, design.pairs, design.has_replicates,
                design.has_pairs))

        # load sleuth
        r('''suppressMessages(library('sleuth'))''')

        # make variates string to ensure all model terms are in the
        # design dataframe for sleuth
        model_terms = [x for x in re.split("[\+~ ]+", model)[1:]
                       if x != "0"]

        variates = "c(%s)" % ",".join(model_terms)

        # need to code in option to not use a reference group (e.g for LRT)
        if genewise:
            assert gene_biomart, ("for genewise analysis, "
                                  "must provide a 'gene_biomart'")

            createSleuthObject = r('''
            function(design_df){

            library(biomaRt)

            sample_id = design_df$sample
            kal_dirs <- sapply(sample_id,
                function(id) file.path('%(base_dir)s', id))

            design_df <- dplyr::select(design_df, sample = sample,
                                       %(variates)s)
            design_df <- dplyr::mutate(design_df, path = kal_dirs)

            %(contrast)s <- factor(design_df$%(contrast)s)
            %(contrast)s <- relevel(%(contrast)s,ref='%(ref_group)s')
            md <- model.matrix(%(model)s, design_df)
            colnames(md)[grep("%(contrast)s", colnames(md))] <- '%(contrast)s%(ref_group)s'

            mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
            #dataset = "hsapiens_gene_ensembl",
            dataset = "%(gene_biomart)s",
            host="www.ensembl.org")

            t2g <- biomaRt::getBM(
                attributes = c("ensembl_transcript_id","ensembl_gene_id",
                               "external_gene_name"), mart = mart)
            t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                                 ens_gene = ensembl_gene_id,
                                 ext_gene = external_gene_name)

            so <- sleuth_prep(design_df, md,
                              target_mapping = t2g, aggregation_column = 'ens_gene')
            so <- suppressMessages(sleuth_fit(so))

            return(so)
            }''' % locals())

        else:
            createSleuthObject = r('''
            function(design_df){

            sample_id = design_df$sample
            kal_dirs <- sapply(sample_id,
                function(id) file.path('%(base_dir)s', id))

            design_df <- dplyr::select(design_df, sample = sample,
                                       %(variates)s)
            design_df <- dplyr::mutate(design_df, path = kal_dirs)

            %(contrast)s <- factor(design_df$%(contrast)s)
            %(contrast)s <- relevel(%(contrast)s,ref='%(ref_group)s')
            md <- model.matrix(%(model)s, design_df)
            colnames(md)[grep("%(contrast)s", colnames(md))] <- '%(contrast)s%(ref_group)s'

            so <- sleuth_prep(design_df, md)
            so <- sleuth_fit(so)

            return(so)
            }''' % locals())

        so = createSleuthObject(r_design_df)

        # write out counts and tpm tables if required
        if counts:

            makeCountsTable = r('''
            function(so){

            library('reshape')

            df = cast(so$obs_raw, target_id~sample, value = "est_counts")
            colnames(df)[1] <- "transcript_id"
            write.table(df, "%(counts)s", sep="\t", row.names=F, quote=F)

            }''' % locals())

            makeCountsTable(so)

        if tpm:

            makeTPMTable = r('''
            function(so){

            library('reshape')

            df = cast(so$obs_raw, target_id~sample, value = "tpm")
            colnames(df)[1] <- "transcript_id"
            write.table(df, "%(tpm)s", sep="\t", row.names=F, quote=F)

            }''' % locals())

            makeTPMTable(so)

        if dummy_run:
            return None

        if DE_test == "lrt":
            differentialTesting = r('''
            function(so){
            so <- suppressMessages(sleuth_fit(so, formula = %(reduced_model)s,
                                   fit_name = "reduced"))
            so <- suppressMessages(sleuth_lrt(so, "reduced", "full"))

            results_table <- sleuth_results(so, test = 'reduced:full',
                                            test_type = 'lrt')
            return(results_table)

            } ''' % locals())

            final_result = DEResult_Sleuth(pandas2ri.ri2py(
                differentialTesting(so)))

        elif DE_test == "wald":

            differentialTesting = r('''
            function(so){

            so <- sleuth_wt(so, which_beta = '%(contrast)s%(ref_group)s')

            p_ma = plot_ma(so, '%(contrast)s%(ref_group)s')
            ggsave("%(outfile_prefix)s_%(contrast)s_sleuth_ma.png",
                width=15, height=15, units="cm")

            p_vars = plot_vars(so, '%(contrast)s%(ref_group)s')
            ggsave("%(outfile_prefix)s_%(contrast)s_sleuth_vars.png",
                width=15, height=15, units="cm")

            p_mean_var = plot_mean_var(so)
            ggsave("%(outfile_prefix)s_%(contrast)s_sleuth_mean_var.png",
            width=15, height=15, units="cm")

            results_table <- sleuth_results(so, test = '%(contrast)s%(ref_group)s')

            return(results_table)
            } ''' % locals())

            results = pandas2ri.ri2py(differentialTesting(so))
            results['contrast'] = contrast

        else:
            raise ValueError("unknown DE test type, use 'wald' or 'lrt'")

        final_result = DEResult_Sleuth(results)

        return final_result


class DEResult_Sleuth(DEResult):

    def getResults(self, fdr):
        ''' post-process test results table into generic results output
        expression and fold changes from Sleuth are natural logs'''

        E.info("Generating output - results table")

        df_dict = collections.defaultdict()

        n_rows = self.table.shape[0]

        df_dict["treatment_name"] = ("NA",)*n_rows
        df_dict["control_name"] = ("NA",)*n_rows
        df_dict["test_id"] = self.table['target_id']
        df_dict["contrast"] = self.table['contrast']
        df_dict["control_mean"] = [math.exp(float(x)) for
                                   x in self.table['mean_obs']]
        df_dict["treatment_mean"] = df_dict["control_mean"]
        df_dict["control_std"] = (0,)*n_rows
        df_dict["treatment_std"] = (0,)*n_rows
        df_dict["p_value"] = self.table['pval']
        df_dict["p_value_adj"] = adjustPvalues(self.table['pval'])
        df_dict["significant"] = pvaluesToSignficant(df_dict["p_value_adj"],
                                                     fdr)
        df_dict["fold"] = [math.exp(float(x)) for
                           x in self.table['b']]
        df_dict["l2fold"] = [math.log(float(x), 2) for x in df_dict['fold']]
        df_dict["transformed_l2fold"] = df_dict["l2fold"]

        # set all status values to "OK"
        # TS: again, may need an if/else to check...
        df_dict["status"] = ("OK",)*n_rows

        self.table = pandas.DataFrame(df_dict)
        # causes errors if multiple instance of same test_id exist, for example
        # if multiple constrasts have been tested
        # self.table.set_index("test_id", inplace=True)


###############################################################################


def buildProbeset2Gene(infile,
                       outfile,
                       database="hgu133plus2.db",
                       mapping="hgu133plus2ENSEMBL"):
    '''build map relating a probeset to an ENSEMBL gene_id'''

    R.library(database)

    # map is a Bimap object
    m = r(mapping)

    result = R.toTable(m)

    outf = open(outfile, "w")
    outf.write("probe_id\tgene_id\n")
    for probeset_id, gene_id in zip(result["probe_id"],
                                    result["ensembl_id"]):
        outf.write("%s\t%s\n" % (probeset_id, gene_id))
    outf.close()

    E.info("written %i mappings to %s: probes=%i, genes=%i" %
           (len(result),
            outfile,
            len(set(result["probe_id"])),
            len(set(result["ensembl_id"]))))

GeneExpressionResult = collections.namedtuple(
    "GeneExpressionResult",
    "test_id treatment_name treatment_mean treatment_std "
    "control_name control_mean control_std "
    "pvalue qvalue l2fold fold transformed_l2fold "
    "significant status")


def writeExpressionResults(outfile, result):
    '''output expression results table.'''
    if outfile == sys.stdout:
        outf = outfile
    else:
        outf = IOTools.openFile(outfile, "w")

    outf.write("%s\n" % "\t".join(GeneExpressionResult._fields))
    for x in sorted(result):
        outf.write("%s\n" % "\t".join(map(str, x)))

    if outf != sys.stdout:
        outf.close()


class WelchsTTest(object):

    '''base class for computing expression differences.
    '''

    def __call__(self,
                 probesets,
                 treatments,
                 controls):

        assert len(probesets) == len(treatments[0])
        assert len(probesets) == len(controls[0])

        nskipped = 0
        results = []

        for probeset, treatment, control in zip(
                probesets, zip(*treatments), zip(*controls)):

            nval1, nval2 = len(treatment), len(control)
            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)
            stddev1, stddev2 = numpy.std(treatment), numpy.std(control)

            try:
                s = Stats.doWelchsTTest(nval1, mean1, stddev1,
                                        nval2, mean2, stddev2,
                                        alpha=0.05)
            except ValueError:
                E.warn(
                    "expressionDifferences: standard deviations are 0 for "
                    "probeset %s - skipped" % probeset)
                nskipped += 1
                continue

            s.mProbeset = probeset
            results.append(s)

        qvalues = Stats.doFDR([x.mPValue for x in results]).mQValues

        for s, qvalue in zip(results, qvalues):
            s.mQValue = qvalue

        return results, nskipped


class SAMR(object):

    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the samr library.

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See

    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. This version is parameterised to reproduce Excel SAM
    by setting::

       var.equal = TRUE
       med = TRUE

    .. note::
        SAM requires log2 scaled expression levels.
    '''

    def __call__(self, probesets,
                 treatments,
                 controls,
                 pattern=None,
                 fdr=0.10,
                 ngenes=None,
                 npermutations=1000,
                 ndelta=10,
                 method="ttest"):

        if ngenes and fdr:
            raise ValueError("either supply ngenes or fdr, but not both.")

        R.library("samr")

        m = numpy.matrix(treatments + controls)
        m = numpy.transpose(m)
        labels = numpy.array([1] * len(treatments) + [2] * len(controls))

        R.assign("x", numpy.array(m))
        R.assign("y", labels)
        R.assign("probesets", probesets)

        data = r(
            '''data=list( x=x, y=y, geneid=1:length(probesets), genenames=probesets, logged2=TRUE)''')
        result = r(
            '''samr.obj<-samr(data,  resp.type="Two class unpaired", nperms=100)''')
        r('''plot(samr.obj, delta=.4)''')


class SAM(object):

    '''SAM analysis of microarray data.

    Use the Two-Class Unpaired Case Assuming Unequal Variances.

    This uses the siggenes library. Note that there is also
    an rsam package at:

    http://rss.acs.unt.edu/Rdoc/library/samr/html/samr.html

    Significant genes are either called at *fdr* or the
    top *ngenes* are returned.

    *treatments* and *control* are arrays of
    arrays of expression values.

    See
    https://stat.ethz.ch/pipermail/bioconductor/2008-July/023251.html

    for an explanation of the differences between siggens SAM
    and Excel SAM. To parameterize the FDR to excel sam, set the
    flag *use_excel_sam*.

    .. note::
        SAM requires log2 scaled expression levels.

    I ran into trouble using this library. I was not able to
    reproduce the same results from the original SAM study getting
    differences in d and in the fdr.

    fold change is treatment / control.

    '''

    def __call__(self, probesets,
                 treatments,
                 controls,
                 pattern=None,
                 fdr=0.10,
                 ngenes=None,
                 npermutations=1000,
                 ndelta=10,
                 method="ttest",
                 use_excel_sam=False,
                 treatment_label="treatment",
                 control_label="control"):

        if ngenes and fdr:
            raise ValueError("either supply ngenes or fdr, but not both.")

        R.library("siggenes")

        m = numpy.matrix(treatments + controls)
        m = numpy.transpose(m)

        E.debug("build expression matrix: %i x %i" % m.shape)

        labels = numpy.array([1] * len(treatments) + [0] * len(controls))
        # 1000 permutations for P-Values of down to 0.0001. Setting this
        # to a high value improved reproducibility of results.

        kwargs = {}
        # kwargs set to replicate excel SAM
        if use_excel_sam:
            kwargs.update(
                {"control":
                 r('''samControl( lambda = 0.5, n.delta = %(ndelta)s) ''' %
                   locals()),
                 "med": True,
                 "var.equal": True})
        else:
            kwargs.update({"control":
                           r('''samControl( n.delta = %(ndelta)s ) ''' %
                             locals())},)

        # the option B needs to be not set if wilc.stat is chosen

        if method == "ttest":
            kwargs["method"] = r('''d.stat''')
            kwargs["B"] = npermutations
        elif method == "wilc":
            kwargs["method"] = r('''wilc.stat''')
        elif method == "cat":
            kwargs["method"] = r('''cat.stat''')
        else:
            raise ValueError("unknown statistic `%s`" % method)

        E.info("running sam with the following options: %s" % str(kwargs))

        a = R.sam(numpy.array(m),
                  labels,
                  gene_names=numpy.array(probesets),
                  **kwargs)

        # E.debug("%s" % str(a))

        R.assign("a", a)

        fdr_data = collections.namedtuple("sam_fdr", (
            "delta", "p0", "false", "significant", "fdr", "cutlow",
            "cutup", "j2", "j1"))
        cutoff_data = collections.namedtuple(
            "sam_cutoff", ("delta", "significant", "fdr"))
        gene_data = collections.namedtuple(
            "sam_fdr", ("row", "dvalue", "stddev", "rawp", "qvalue", "rfold"))

        def _totable(robj):
            '''convert robj to a row-wise table.'''
            s = numpy.matrix(robj)
            t = [numpy.array(x).reshape(-1,) for x in s]
            return t

        # extract the fdr values
        # returns R matrix
        t = _totable(a.do_slot('mat.fdr'))
        assert len(t[0]) == len(fdr_data._fields)
        # for x in t: E.debug( "x=%s" % str(x))
        fdr_values = [fdr_data(*x) for x in t]

        # find d cutoff
        if fdr is not None and fdr > 0:
            s = numpy.matrix(R.findDelta(a, fdr))
            try:
                cutoffs = [cutoff_data(*numpy.array(x).reshape(-1,))
                           for x in s]
                E.debug("sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs)))
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug("could not get cutoff")
                cutoff = None
        elif ngenes:
            s = numpy.matrix(R.findDelta(a, ngenes))
            try:
                cutoffs = [cutoff_data(*numpy.array(x).reshape(-1,))
                           for x in s]
                E.debug("sam cutoffs for fdr %f: %s" % (fdr, str(cutoffs)))
                cutoff = cutoffs[-1]
            except TypeError:
                E.debug("could not get cutoff")
                cutoff = None
        else:
            raise ValueError("either supply ngenes or fdr")

        # collect (unadjusted) p-values and qvalues for all probesets
        pvalues = dict(zip(probesets, r('''a@p.value''')))
        qvalues = dict(zip(probesets, r('''a@q.value''')))

        if pattern:
            outfile = pattern % "sam.pdf"
            R.pdf(outfile)
            if cutoff:
                R.plot(a, cutoff.delta)
            else:
                R.plot(a)
            r['dev.off']()

        siggenes = {}
        significant_genes = set()
        if cutoff is not None:
            E.debug("using cutoff %s" % str(cutoff))

            summary = r('''summary( a, %f )''' % cutoff.delta)

            # summary = R.summary( a, cutoff.delta )
            R.assign("summary", summary)

            significant_genes = set(
                [probesets[int(x) - 1] for x
                 in r('''summary@row.sig.genes''')])
            # E.debug( "significant genes=%s" % str(significant_genes))

            r_result = zip(*_totable(summary.do_slot('mat.sig')))

            if len(r_result) > 0:

                assert len(r_result[0]) == 6, \
                    "expected six columns from siggenes module, got: %s" % \
                    len(r_result[0])

                for x in r_result:
                    if x[4] > fdr:
                        E.warn("%s has qvalue (%f) larger than cutoff, but "
                               "is significant significant." % (str(x), x[4]))

                # except TypeError:
                # only a single value
                #     x = [r_result[y] for y in ("Row", "d.value", "stdev", "rawp", "q.value", "R.fold") ]
                #     if x[4] > fdr:
                #         E.warn( "%s has qvalue (%f) larger than cutoff, but is called significant." % (str(x), x[4]))

                siggenes[probesets[int(x[0]) - 1]] = gene_data(*x)

        else:
            E.debug("no cutoff found - no significant genes.")

        genes = []
        for probeset, treatment, control in zip(
                probesets, zip(*treatments), zip(*controls)):

            mean1, mean2 = numpy.mean(treatment), numpy.mean(control)

            if probeset in siggenes:
                s = siggenes[probeset]
                pvalue = s.rawp
                qvalue = s.qvalue
            else:
                pvalue = pvalues[probeset]
                qvalue = qvalues[probeset]

            significant = (0, 1)[probeset in significant_genes]

            genes.append(GeneExpressionResult._make((probeset,
                                                     treatment_label,
                                                     mean1,
                                                     numpy.std(treatment),
                                                     control_label,
                                                     mean2,
                                                     numpy.std(control),
                                                     pvalue,
                                                     qvalue,
                                                     mean1 - mean2,
                                                     math.pow(
                                                         2, mean1 - mean2),
                                                     math.pow(
                                                         2, mean1 - mean2),
                                                     significant,
                                                     "OK")))

        return genes, cutoff, fdr_values


#########################################################################
#########################################################################
#########################################################################
def loadTagData(tags_filename, design_filename):
    '''load tag data for deseq/edger analysis.

    *tags_file* is a tab-separated file with counts.

    *design_file* is a tab-separated file with the
    experimental design with a minimum of four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

    track
        name of track - should correspond to column header in *infile*
    include
        flag to indicate whether or not to include this data
    group
        group indicator - experimental group
    pair
        pair that sample belongs to (for paired tests)

    Additional columns in design file are taken to contain levels for
    additional factors and may be included for tests that allow multi-factor
    model designs.

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups
    pairs  : vector with pairs
    factors : df of additional factors for more complex model designs

    '''

    # Load counts table
    E.info("loading tag data from %s" % tags_filename)

    r('''counts_table = read.table('%(tags_filename)s',
    header=TRUE,
    row.names=1,
    stringsAsFactors=TRUE,
    comment.char='#')''' % locals())

    E.info("read data: %i observations for %i samples" %
           tuple(r('''dim(counts_table)''')))
    E.debug("sample names: %s" % r('''colnames(counts_table)'''))

    # Load comparisons from file
    r('''pheno = read.delim(
    '%(design_filename)s',
    header=TRUE,
    stringsAsFactors=TRUE,
    comment.char='#')''' % locals())

    # Make sample names R-like - substitute - for .
    r('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')
    E.debug("design names: %s" % r('''pheno[,1]'''))

    # Ensure pheno rows match count columns
    pheno = r(
        '''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''')
    missing = r('''colnames(counts_table)[is.na(pheno2)][1]''')
    if missing:
        E.warn("missing samples from design file are ignored: %s" %
               missing)

    # Subset data & set conditions
    r('''includedSamples <- !(is.na(pheno2$include) | pheno2$include == '0') ''')
    E.debug("included samples: %s" %
            r('''colnames(counts_table)[includedSamples]'''))
    r('''countsTable <- counts_table[ , includedSamples ]''')
    r('''groups <- factor(pheno2$group[ includedSamples ])''')
    r('''conds <- pheno2$group[ includedSamples ]''')
    r('''pairs <- factor(pheno2$pair[ includedSamples ])''')

    # JJ if additional columns present, pass to 'factors'
    r('''if (length(names(pheno2)) > 4) {
           factors <- data.frame(pheno2[includedSamples,5:length(names(pheno2))])
         } else {
           factors <- NA
         }''')

    E.info("filtered data: %i observations for %i samples" %
           tuple(r('''dim(countsTable)''')))


def filterTagData(filter_min_counts_per_row=1,
                  filter_min_counts_per_sample=10,
                  filter_percentile_rowsums=0):
    '''filter tag data.

    * remove rows with fewer than x counts in most highly expressed sample

    * remove samples with fewer than x counts in most highly expressed row

    * remove the lowest percentile of rows in the table, sorted
       by total tags per row
    '''

    # Remove windows with no data
    r('''max_counts = apply(countsTable,1,max)''')
    r('''countsTable = countsTable[max_counts>=%i,]''' %
      filter_min_counts_per_row)
    E.info("removed %i empty rows" %
           tuple(r('''sum(max_counts == 0)''')))
    observations, samples = tuple(r('''dim(countsTable)'''))
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # remove samples without data
    r('''max_counts = apply(countsTable,2,max)''')

    empty_samples = tuple(
        r('''max_counts < %i''' % filter_min_counts_per_sample))
    sample_names = r('''colnames(countsTable)''')
    nempty_samples = sum(empty_samples)

    if nempty_samples:
        E.warn("%i empty samples are being removed: %s" %
               (nempty_samples,
                ",".join([sample_names[x]
                          for x, y in enumerate(empty_samples) if y])))
        r('''countsTable <- countsTable[, max_counts >= %i]''' %
          filter_min_counts_per_sample)
        r('''groups <- groups[max_counts >= %i]''' %
          filter_min_counts_per_sample)
        r('''pairs <- pairs[max_counts >= %i]''' %
          filter_min_counts_per_sample)
        r('''if (!is.na(factors)) {factors <- factors[max_counts >= %i,]}''' %
          filter_min_counts_per_sample)
        observations, samples = tuple(r('''dim(countsTable)'''))

    # percentile filtering
    if filter_percentile_rowsums > 0:
        percentile = float(filter_percentile_rowsums) / 100.0
        r('''sum_counts = rowSums( countsTable )''')
        r('''take = (sum_counts > quantile(sum_counts, probs = %(percentile)f))''' %
          locals())
        discard, keep = r('''table( take )''')
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (filter_percentile_rowsums,
                keep, discard))
        r('''countsTable = countsTable[take,]''')

    observations, samples = tuple(r('''dim(countsTable)'''))

    return observations, samples


def groupTagData(ref_group=None, ref_regex=None):
    '''compute groups and pairs from tag data table.'''

    if ref_regex is not None and ref_group is None:
        groups = r('''levels(groups)''')
        for g in groups:
            rx = re.compile(ref_regex)
            if rx.search(g):
                ref_group = g

    # Relevel the groups so that the reference comes first
    if ref_group is not None:
        E.info("reference group (control) is '%s'" % ref_group)
        r('''groups <- relevel(groups, ref="%s")''' % ref_group)

    groups = r('''levels(groups)''')
    pairs = r('''levels(pairs)''')
    factors = r('''factors''')

    # JJ - check whether there are additional factors in design file...
    # warning... isintance(df, rpy.robjects.vectors.Vector) returns True
    if isinstance(factors, rpy2.robjects.vectors.DataFrame):
        E.warn("There are additional factors in design file that are ignored"
               " by groupTagData: %s" % factors.r_repr())
    # AH: uncommented as it causes an issue with rpy2 2.7.7
    # else:
    #     # Hack... must be a better way to evaluate r NA instance in python?
    #     import pdb; pdb.set_trace()
    #     assert len(list(factors)) == 1 and bool(list(factors)[0]) is False, \
    #         "factors must either be DataFrame or NA in R global namespace"

    # Test if replicates exist - at least one group must have multiple samples
    max_per_group = r('''max(table(groups)) ''')[0]

    has_replicates = max_per_group >= 2

    # Test if pairs exist:
    npairs = r('''length(table(pairs)) ''')[0]
    has_pairs = npairs == 2

    # at least two samples per pair
    if has_pairs:
        min_per_pair = r('''min(table(pairs)) ''')[0]
        has_pairs = min_per_pair >= 2

    return groups, pairs, has_replicates, has_pairs


def plotCorrelationHeatmap(method="correlation"):
    '''plot a heatmap of correlations derived from
    countsTable.
    '''

    if method == "correlation":
        r('''dists <- dist( (1 - cor(countsTable)) / 2 )''')
    else:
        r('''dists <- dist( t(as.matrix(countsTable)), method = '%s' )''' %
          method)

    r('''heatmap( as.matrix( dists ), symm=TRUE )''')


def plotPairs():
    '''requires counts table'''
    # Plot pairs
    r('''panel.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y))
            txt <- format(c(r, 0.123456789), digits=digits)[1]
            txt <- paste(prefix, txt, sep="")
            if(missing(cex.cor)) cex <- 0.6/strwidth(txt)
            x = 0.5;
            y = 0.5;
            if (par("xlog")) { x = 10^x };
            if (par("ylog")) { y = 10^y };
            text(x, y, txt, cex = cex);
            }
       ''')
    try:
        r('''pairs(countsTable,
        lower.panel = panel.pearson,
        pch=".",
        labels=colnames(countsTable),
        log="xy")''')
    except RRuntimeError:
        E.warn("can not plot pairwise scatter plot")


def plotPCA(groups=True):
    '''plot a PCA plot from countsTable using ggplot.

    If groups is *True*, the variable ``groups`` is
    used for colouring. If *False*, the groups are
    determined by sample labels.
    '''
    r('''suppressMessages(library(ggplot2))''')
    r('''pca = prcomp(t(countsTable))''')
    # Build factor groups by splitting labels at "."
    r('''colour=groups''')
    r('''shape=0''')
    r('''size=1''')
    if groups is False:
        r('''mm = matrix(
        unlist(sapply(colnames(countsTable),strsplit,'[.]')),
        nrow=length(colnames(countsTable)),
        byrow=T)''')
        nrows, nlevels = r('''dim(mm)''')
        if nlevels > 1:
            r('''colour=mm[,1]''')
        if nlevels > 2:
            r('''shape=mm[,2]''')

    try:
        r('''p1 = ggplot(
        as.data.frame(pca$x),
        aes(x=PC1, y=PC2,
        colour=colour,
        shape=shape,
        label=rownames(pca$x))) \
        + geom_text(size=4, vjust=1) \
        + geom_point()''')
        r('''p2 = qplot(x=PC1, y=PC3,
        data = as.data.frame(pca$x),
        label=rownames(pca$x),
        shape=shape,
        colour=colour)''')
        r('''p3 = qplot(x=PC2, y=PC3,
        data = as.data.frame(pca$x),
        label=rownames(pca$x),
        shape=shape,
        colour=colour)''')
        # TODO: plot all in a multi-plot with proper scale
        # the following squishes the plots
        # r('''source('%s')''' %
        #   os.path.join(os.path.dirname(E.__file__),
        #                "../R",
        #                "multiplot.R"))
        # r('''multiplot(p1, p2, p3, cols=2)''')
        r('''plot(p1)''')
    except RRuntimeError as msg:
        E.warn("could not plot in plotPCA(): %s" % msg)


def runEdgeR(outfile,
             outfile_prefix="edger.",
             fdr=0.1,
             prefix="",
             dispersion=None,
             ref_group=None,
             ref_regex=None,
             ):
    '''run EdgeR on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The dispersion is usually measuered from replicates. If there are no
    replicates, you need to set the *dispersion* explicitely.

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.

    '''

    # load library
    r('''suppressMessages(library('edgeR'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group,
                                                            ref_regex)

    # output heatmap plot
    fn = '%(outfile_prefix)sheatmap.png' % locals()
    E.info("outputing heatmap to {}".format(fn))
    R.png(fn)
    plotCorrelationHeatmap()
    r['dev.off']()

    E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
           (groups, pairs, has_replicates, has_pairs))

    if has_pairs:
        # output difference between groups
        R.png('''%(outfile_prefix)sbalance_groups.png''' % locals())
        first = True
        for g1, g2 in itertools.combinations(groups, 2):
            r('''a = rowSums( countsTable[groups == '%s'] ) ''' % g1)
            r('''b = rowSums( countsTable[groups == '%s'] ) ''' % g2)
            if first:
                r('''plot(cumsum(sort(a - b)), type = 'l') ''')
                first = False
            else:
                r('''lines(cumsum(sort(a - b))) ''')

        r['dev.off']()

        r('''suppressMessages(library('ggplot2'))''')
        r('''suppressMessages(library('reshape'))''')

        # output difference between pairs within groups
        first = True
        legend = []
        for pair in pairs:
            for g1, g2 in itertools.combinations(groups, 2):
                key = re.sub("-", "_", "pair_%s_%s_vs_%s" % (pair, g1, g2))
                legend.append(key)
                # print r('''colnames( countsTable) ''')
                # print r(''' pairs=='%s' ''' % pair)
                # print r(''' groups=='%s' ''' % g1)
                r('''a = rowSums( countsTable[pairs == '%s' & groups == '%s'])''' % (
                    pair, g1))
                r('''b = rowSums( countsTable[pairs == '%s' & groups == '%s'])''' % (
                    pair, g2))
                r('''c = cumsum( sort(a - b) )''')
                r('''c = c - min(c)''')
                if first:
                    data = r('''d = data.frame( %s = c)''' % key)
                    first = False
                else:
                    r('''d$%s = c''' % key)

        # remove row names (gene idenitifiers)
        r('''row.names(d) = NULL''')
        # add numbers of genes (x-axis)
        r('''d$genes=1:nrow(d)''')

        # merge data for ggplot
        r('''d = melt( d, 'genes', variable_name = 'comparison' )''')

        # plot
        r('''gp = ggplot(d)''')
        r('''pp = gp + \
            geom_line(aes(x=genes,y=value,group=comparison,color=comparison))''')

        try:
            R.ggsave('''%(outfile_prefix)sbalance_pairs.png''' % locals())
            r['dev.off']()
        except RRuntimeError:
            E.warn("could not plot")

    # build DGEList object
    # ignore message: "Calculating library sizes from column totals"
    r('''countsTable = suppressMessages(DGEList(countsTable, group=groups))''')

    # Relevel groups to make the results predictable - IMS
    if ref_group is not None:
        r('''countsTable$samples$group <- relevel(countsTable$samples$group,
        ref = "%s")''' % ref_group)
    else:
        # if no ref_group provided use first group in groups
        r('''countsTable$sample$group <- relevel(countsTable$samples$group,
        ref = "%s")''' % groups[0])

    # calculate normalisation factors
    E.info("calculating normalization factors")
    r('''countsTable = calcNormFactors( countsTable )''')
    E.info("output")

    # output MDS plot
    R.png('''%(outfile_prefix)smds.png''' % locals())
    try:
        r('''plotMDS( countsTable )''')
    except RRuntimeError:
        E.warn("can not plot mds")
    r['dev.off']()

    # build design matrix
    if has_pairs:
        r('''design = model.matrix(~pairs + countsTable$samples$group)''')
    else:
        r('''design = model.matrix(~countsTable$samples$group)''')

    # r('''rownames(design) = rownames( countsTable$samples )''')
    # r('''colnames(design)[length(colnames(design))] = "CD4" ''' )

    # fitting model to each tag
    if has_replicates:
        # estimate common dispersion
        r('''countsTable = estimateGLMCommonDisp(countsTable, design)''')
        # estimate tagwise dispersion
        r('''countsTable = estimateGLMTagwiseDisp(countsTable, design)''')
        # fitting model to each tag
        r('''fit = glmFit(countsTable, design)''')
    else:
        # fitting model to each tag
        if dispersion is None:
            raise ValueError("no replicates and no dispersion")
        E.warn("no replicates - using a fixed dispersion value")
        r('''fit = glmFit(countsTable, design, dispersion=%f)''' %
          dispersion)

    # perform LR test
    r('''lrt = glmLRT(fit)''')

    E.info("Generating output")

    # output cpm table
    r('''suppressMessages(library(reshape2))''')
    r('''countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)''')
    r('''melted <- melt(countsTable.cpm)''')
    r('''names(melted) <- c("test_id", "sample", "ncpm")''')
    # melt columns are factors - convert to string for sorting
    r('''melted$test_id = levels(melted$test_id)[as.numeric(melted$test_id)]''')
    r('''melted$sample = levels(melted$sample)[as.numeric(melted$sample)]''')
    # sort cpm table by test_id and sample
    r('''sorted = melted[with(melted, order(test_id, sample)),]''')
    r('''gz = gzfile("%(outfile_prefix)scpm.tsv.gz", "w" )''' % locals())
    r('''write.table(sorted, file=gz, sep = "\t",
                     row.names=FALSE, quote=FALSE)''')
    r('''close(gz)''')

    # compute adjusted P-Values
    r('''padj = p.adjust(lrt$table$PValue, 'BH')''')

    rtype = collections.namedtuple("rtype", "lfold logCPM LR pvalue")

    # output differences between pairs
    if len(groups) == 2:
        R.png('''%(outfile_prefix)smaplot.png''' % locals())
        r('''plotSmear(countsTable, pair=c('%s'))''' % "','".join(groups))
        r('''abline(h=c(-2, 2), col='dodgerblue') ''')
        r['dev.off']()

    # I am assuming that logFC is the base 2 logarithm foldchange.
    # Parse results and parse to file
    results = []
    counts = E.Counter()

    df = r('''lrt$table''')
    df["padj"] = r('''padj''')

    counts.significant = sum(df.padj <= fdr)
    counts.insignificant = sum(df.padj > fdr)
    counts.significant_over = sum((df.padj <= fdr) & (df.logFC > 0))
    counts.significant_under = sum((df.padj <= fdr) & (df.logFC < 0))
    counts.input = len(df)
    counts.all_over = sum(df.logFC > 0)
    counts.all_under = sum(df.logFC < 0)
    counts.fail = sum(df.PValue.isnull())
    counts.ok = counts.input - counts.fail

    df["fold"] = df.logFC.pow(2.0)
    df["significant"] = df.padj <= fdr

    # TODO: use pandas throughout
    for interval, d in df.iterrows():
        # fold change is determined by the alphabetical order of the factors.
        # Is the following correct?
        results.append(GeneExpressionResult._make((
            interval,
            groups[1],
            d.logCPM,
            0,
            groups[0],
            d.logCPM,
            0,
            d.PValue,
            d.padj,
            d.logFC,
            d.fold,
            d.logFC,  # no transform of lfold
            str(int(d.significant)),
            "OK")))

    writeExpressionResults(outfile, results)

    outf = IOTools.openFile("%(outfile_prefix)ssummary.tsv" % locals(), "w")
    outf.write("category\tcounts\n%s\n" % counts.asTable())
    outf.close()

# needs to put into class
##


def deseqPlotSizeFactors(outfile):
    '''plot size factors - requires cds object.'''
    R.png(outfile)
    r('''par(mar=c(8,4,4,2))''')
    r('''barplot( sizeFactors( cds ), main="size factors", las=2)''')
    r['dev.off']()


def deseqOutputSizeFactors(outfile):
    '''output size factors - requires cds object.'''
    size_factors = r('''sizeFactors( cds )''')
    samples = r('''names(sizeFactors(cds))''')
    with IOTools.openFile(outfile, "w") as outf:
        outf.write("sample\tfactor\n")
        for name, x in zip(samples, size_factors):
            outf.write("%s\t%s\n" % (name, str(x)))


def deseqPlotCorrelationHeatmap(outfile, vsd):
    '''plot a heatmap

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''

    # rpy2.4.2 - passing of arrays seems to be broken - do it in R
    # dists = r['as.matrix'](R.dist(R.t(R.exprs(vsd))))
    dists = r('''as.matrix(dist(t(exprs(vsd))))''')
    R.png(outfile)
    r['heatmap.2'](
        dists,
        trace='none',
        margin=ro.IntVector((10, 10)))
    r['dev.off']()


def deseqPlotGeneHeatmap(outfile,
                         data,
                         Rowv=False,
                         Colv=False):
    '''plot a heatmap of all genes

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    if len(data) == 0:
        return

    # do not print if not enough values in one
    # direction (single row or column)
    if min(R.dim(data)) < 2:
        return

    R.png(outfile, width=500, height=2000)
    hmcol = R.colorRampPalette(r['brewer.pal'](9, "GnBu"))(100)

    r['heatmap.2'](
        data,
        col=hmcol,
        trace="none",
        dendrogram="none",
        Rowv=Rowv,
        Colv=Colv,
        labRow=False,
        margin=ro.IntVector((5, 5)),
        lhei=ro.IntVector((1, 10)),
        key=False)

    r['dev.off']()


def deseqPlotPCA(outfile, vsd, max_genes=500):
    '''plot a PCA

    Use variance stabilized data in object vsd.
    Should be 'blind', as then the transform is
    not informed by the experimental design.
    '''
    R.png(outfile)
    #  if there are more than 500 genes (after filtering)
    #  use the 500 most variable in the PCA
    #  else use the number of genes
    r('''ntop = ifelse(as.integer(dim(vsd))[1] >= %(max_genes)i,
    %(max_genes)i,
    as.integer(dim(vsd))[1])''' % locals())
    try:
        r('''plotPCA(vsd)''')
    except RRuntimeError as msg:
        E.warn("can not plot PCA: %s" % msg)
    r['dev.off']()


def deseqPlotPairs(outfile):
    '''requires counts table'''
    # Plot pairs
    R.png(outfile, width=960, height=960)
    plotPairs()
    r['dev.off']()


def deseqPlotPvaluesAgainstRowsums(outfile):
    '''plot pvalues against row sum rank.

    This plot is useful to see if quantile filtering could
    be applied.
    '''

    r('''counts_sum = rowSums( countsTable )''')
    R.png(outfile)
    r('''plot( rank( counts_sum)/length(counts_sum),
               -log10( res$pval),
               pch = 16,
               cex= 0.1)''')

    r('''abline( a=3, b=0, col='red')''')
    r['dev.off']()


def deseqParseResults(control_name, treatment_name, fdr, vsd=False):
    '''parse deseq output.

    retrieve deseq results from object 'res' in R namespace.

    The 'res' object is a dataframe with the following columns (from the
    DESeq manual):

    id: The ID of the observable, taken from the row names of the
          counts slots.

    baseMean: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for both conditions

    baseMeanA: The base mean (i.e., mean of the counts divided by the size
          factors) for the counts for condition A

    baseMeanB: The base mean for condition B

    foldChange: The ratio meanB/meanA

    log2FoldChange: The log2 of the fold change

    pval: The p value for rejecting the null hypothesis 'meanA==meanB'

    padj: The adjusted p values (adjusted with 'p.adjust( pval,
          method="BH")')

    vsd_log2FoldChange: The log2 fold change after variance stabilization.
          This data field is not part of DESeq proper, but has been added
          in this module in the runDESeq() method.

    Here, 'conditionA' is 'control' and 'conditionB' is 'treatment'
    such that a foldChange of 2 means that treatment is twice
    upregulated compared to control.

    Returns a list of results.

    If vsd is True, the log fold change will be computed from the variance
    stabilized data.

    '''

    results = []
    isna = r["is.na"]

    counts = E.Counter()

    for index, data in r['res'].iterrows():
        counts.input += 1

        # set significant flag
        if data['padj'] <= fdr:
            signif = 1
            counts.significant += 1
            if data['log2FoldChange'] > 0:
                counts.significant_over += 1
            else:
                counts.significant_under += 1
        else:
            signif = 0
            counts.insignificant += 1

        if data['log2FoldChange'] > 0:
            counts.all_over += 1
        else:
            counts.all_under += 1

        # set lfold change to 0 if both are not expressed
        if data['baseMeanA'] == 0.0 and data['baseMeanB'] == 0.0:
            data['foldChange'] = 0
            data['log2FoldChange'] = 0

        if data['pval'] != r('''NA'''):
            status = "OK"
        else:
            status = "FAIL"

        counts[status] += 1

        counts.output += 1

        # check if our assumptions about the direction of fold change
        # are correct
        assert (data['foldChange'] > 1) == (data['baseMeanB'] > data['baseMeanA'])

        # note that fold change is computed as second group (B) divided by
        # first (A)
        results.append(GeneExpressionResult._make((
            data['id'],
            treatment_name,
            data['baseMeanB'],
            0,
            control_name,
            data['baseMeanA'],
            0,
            data['pval'],
            data['padj'],
            data['log2FoldChange'],
            data['foldChange'],
            data['transformed_log2FoldChange'],
            str(signif),
            status)))

    return results, counts


def deseq2ParseResults(treatment_name, control_name, fdr, vsd=False):
    '''
    Standardises the output format from deseq2.

    Deseq2 has the following output columns:
    baseMean log2FoldChange lfcSE stat pvalue padj
    described in
    https://bioconductor.org/packages/release/bioc/
                 vignettes/DESeq2/inst/doc/DESeq2.pdf

    Standardised columns are generated from this output as follows:

    test_id - the gene or region tested, row names from raw output

    treatment_name - the first group in this differential expression
        comparison (from the design file)

    treatment_mean - the mean expression value for this treatment from the
    normalised count table generated by deseq2

    treatment_std - the standard deviation of experssion for this treatment
    from the normalised count table generated by deseq2

    control_name - the second group in this differential expression
        comparison (from the design file)

    control_mean - the mean expression value for this treatment from the
    normalised count table generated by deseq2

    control_std - the standard deviation of experssion for this treatment
    from the normalised count table generated by deseq2

    pvalue - the pvalue generated by deseq2 (from the pvalue column)

    qvalue - the adjusted pvalue generated by deseq2 (from the padj column)

    l2fold - the log2fold change between normalised counts generated by
    deseq2 (log2FoldChange column). If betaPrior is set to TRUE, this is the
    shrunken log2 fold change. If set to FALSE, no shrinkage.

    fold - control mean / treatment mean

    transformed_l2fold - not applicable, set to 0 (see deseq2 manual,
    "The variance stabilizing and rlog transformations are provided
    for applications other than differential testing,
    for example clustering of samples or other machine learning applications.
    For differential testing we recommend the DESeqfunction applied to raw
    counts"

    signif = True if the qvalue is less than the FDR set in `term`PARAMS.

    status = OK if a pvalue has been generated, else FAIL
    '''
    r(''' fdr=%s ''' % fdr)
    #  assign standard column names
    r('''cols = c("test_id",
                  "treatment_name",
                  "treatment_mean",
                  "treatment_std",
                  "control_name",
                  "control_mean",
                  "control_std",
                  "pvalue",
                  "qvalue",
                  "l2fold",
                  "fold",
                  "transformed_l2fold",
                  "signif",
                  "status") ''')

    #  extract normalised counts
    r('''normalcounts = counts(dds, normalized=T)''')
    #  build empty dataframe
    r('''res2 = data.frame(matrix(nrow=nrow(res), ncol=length(cols)))''')
    r('''colnames(res2) = cols''')

    #  fill columns with values described above
    r('''res2['test_id'] = rownames(res)''')
    r('''g = unique(groups[groups == "%s" | groups == "%s"])''' % (treatment_name, control_name))
    r('''g1 = which(groups == g[1])''')
    r('''g2 = which(groups == g[2])''')
    r('''res2['treatment_name'] = g[1]''')
    r('''res2['treatment_mean'] = rowMeans(normalcounts[,g1])''')
    r('''res2['treatment_std'] = apply(normalcounts[,g1], 1, sd)''')
    r('''res2['control_name'] = g[2]''')
    r('''res2['control_mean'] = rowMeans(normalcounts[,g2])''')
    r('''res2['control_std'] = apply(normalcounts[,g2], 1, sd)''')
    r('''res2['pvalue'] = res$pvalue''')
    r('''res2['qvalue'] = res$padj''')
    r('''res2['l2fold'] = res$log2FoldChange''')

    # Fold change here does not reflect the shrinkage applied to
    # log2fold changes
    r('''res2['fold'] = res2$control_mean / res2$treatment_mean''')
    r('''res2['signif'] = as.integer(res2$qvalue <= fdr)''')
    r('''res2['status'] = ifelse(is.na(res2$pvalue), "FAIL", "OK")''')

    #  replace l2fold change and fold for expression levels of 0 in treatment
    #  and control with 0
    r('''z1  = which(res2$treatment_mean == 0)''')
    r('''z2 = which(res2$control_mean == 0)''')
    r('''zeros = intersect(z1, z2)''')
    r('''res2$l2fold[zeros] = 0''')
    r('''res2$fold[zeros] = 0''')

    #  occupy transformed l2fold with 0s
    r('''res2$transformed_l2fold = 0''')

    D = r('res2')
    D.index = D['test_id']
    D = D.drop('test_id', 1)
    return D


def runDESeq(outfile,
             outfile_prefix="deseq.",
             fdr=0.1,
             prefix="",
             fit_type="parametric",
             dispersion_method="pooled",
             sharing_mode="maximum",
             ref_group=None,
             ref_regex=None,
             ):
    '''run DESeq on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The current analysis follows the analysis as outlined in version
    1.14.0

    DESeq ignores any pair information in the design matrix.

    The output is treatment and control. Fold change values are
    computed as treatment divided by control.

    '''

    # load library
    r('''suppressMessages(library('DESeq'))''')
    r('''suppressMessages(library('gplots'))''')
    r('''suppressMessages(library('RColorBrewer'))''')

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group,
                                                            ref_regex)

    # Run DESeq
    # Create Count data object
    E.info("running DESeq: replicates=%s" % (has_replicates))
    r('''cds <-newCountDataSet( countsTable, groups) ''')

    # Estimate size factors
    r('''cds <- estimateSizeFactors( cds )''')

    no_size_factors = r('''is.na(sum(sizeFactors(cds)))''')[0]
    if no_size_factors:
        E.warn("no size factors - can not estimate - no output")
        return

    # estimate variance
    if has_replicates:
        E.info("replicates - estimating variance from replicates")
    else:
        E.info("no replicates - estimating variance with method='blind'")
        dispersion_method = "blind"

    E.info("dispersion_method=%s, fit_type=%s, sharing_mode=%s" %
           (dispersion_method, fit_type, sharing_mode))
    r('''cds <- estimateDispersions( cds,
    method='%(dispersion_method)s',
    fitType='%(fit_type)s',
    sharingMode='%(sharing_mode)s')''' % locals())

    # bring into python namespace
    cds = r('''cds''')

    # plot fit - if method == "pooled":
    if dispersion_method == "pooled":
        R.png('%sdispersion_estimates_pooled.png' %
              outfile_prefix)
        R.plotDispEsts(cds)
        r['dev.off']()
    elif not has_replicates:
        # without replicates the following error appears
        # in the rpy2.py2ri conversion:
        #   'dims' cannot be of length 0
        pass
    else:
        dispersions = r('''ls(cds@fitInfo)''')
        for dispersion in dispersions:
            R.png('%sdispersion_estimates_%s.png' %
                  (outfile_prefix, dispersion))
        R.plotDispEsts(cds, name=dispersion)
        r['dev.off']()

    # plot size factors
    deseqPlotSizeFactors('%(outfile_prefix)ssize_factors.png' % locals())

    # output size factors
    deseqOutputSizeFactors("%(outfile_prefix)ssize_factors.tsv" % locals())

    # plot scatter plots of pairs
    deseqPlotPairs('%(outfile_prefix)spairs.png' % locals())

    if dispersion_method not in ("blind",):
        # also do a blind dispersion estimate for
        # a variance stabilizing transform
        r('''cds_blind <- estimateDispersions( cds,
        method='blind',
        fitType='%(fit_type)s',
        sharingMode='%(sharing_mode)s')''' % locals())
    else:
        r('''cds_blind = cds''')

    # perform variance stabilization for log2 fold changes
    vsd = r('''vsd = varianceStabilizingTransformation(cds_blind)''')
    # output normalized counts (in order)
    # gzfile does not work with rpy 2.4.2 in python namespace
    # using R.gzfile, so do it in R-space

    r('''t = counts(cds, normalized=TRUE);
    write.table(t[order(rownames(t)),],
    file=gzfile('%(outfile_prefix)scounts.tsv.gz'),
    row.names=TRUE,
    col.names=NA,
    quote=FALSE,
    sep='\t') ''' % locals())

    # output variance stabilized counts (in order)
    r('''t = exprs(vsd);
    write.table(t[order(rownames(t)),],
    file=gzfile('%(outfile_prefix)svsd.tsv.gz'),
    row.names=TRUE,
    col.names=NA,
    quote=FALSE,
    sep='\t')
    ''' % locals())

    # plot correlation heatmap of variance stabilized data
    deseqPlotCorrelationHeatmap(
        '%scorrelation_heatmap.png' % outfile_prefix,
        vsd)

    # plot PCA
    deseqPlotPCA('%spca.png' % outfile_prefix,
                 vsd)

    # plot gene heatmap for all genes - order by average expression
    # subtract one to get numpy indices
    select = R.order(R.rowMeans(R.counts(cds)), decreasing=True)
    # the following uses R-based indexing
    deseqPlotGeneHeatmap(
        '%sgene_heatmap.png' % outfile_prefix,
        r['as.matrix'](R.exprs(vsd).rx(select)))

    # plot heatmap of top 200 expressed genes
    deseqPlotGeneHeatmap(
        '%sgene_heatmap_top200.png' % outfile_prefix,
        r['as.matrix'](R.exprs(vsd).rx(select[:200])))

    # Call diffential expression for all pairings of groups included in the
    # design

    all_results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination
        gfix = "%s_vs_%s_" % (control, treatment)

        outfile_groups_prefix = outfile_prefix + gfix
        E.info(("calling differential expression for "
                "control=%s vs treatment=%s") %
               (control, treatment))
        res = r('''res = nbinomTest(cds, '%s', '%s')''' % (control, treatment))

        # plot significance
        R.png('''%(outfile_groups_prefix)ssignificance.png''' % locals())
        r('''plot(
        res$baseMean,
        res$log2FoldChange,
        log="x",
        pch=20, cex=.1,
        col = ifelse( res$padj < %(fdr)s, "red", "black"))''' % locals())
        r['dev.off']()

        # plot pvalues against rowsums
        deseqPlotPvaluesAgainstRowsums(
            '%(outfile_groups_prefix)spvalue_rowsums.png' % locals())

        E.info("Generating output (%s vs %s)" % (control, treatment))

        # get variance stabilized fold changes - note the reversal of
        # treatment/control
        r('''vsd_l2f =
        (rowMeans(exprs(vsd)[,conditions(cds) == '%s', drop=FALSE])
        - rowMeans( exprs(vsd)[,conditions(cds) == '%s', drop=FALSE]))''' %
          (treatment, control))

        # plot vsd correlation, see Figure 14 in the DESeq manual
        # if you also want to colour by expression level
        R.png('''%(outfile_groups_prefix)sfold_transformation.png''' %
              locals())
        r('''plot(
        res$log2FoldChange, vsd_l2f,
        pch=20, cex=.1,
        col = ifelse( res$padj < %(fdr)s, "red", "black" ) )''' % locals())
        r['dev.off']()

        # plot heatmap of differentially expressed genes
        # plot gene heatmap for all genes - order by average expression
        select = r('''select = res['padj'] < %f''' % fdr)

        if r('''sum(select)''')[0] > 0:
            E.info('%s vs %s: plotting %i genes in heatmap' %
                   (treatment, control, len(select)))
            data = R.exprs(vsd).rx(select)

            if not isinstance(data, rpy2.robjects.vectors.FloatVector):
                order = R.order(R.rowMeans(data), decreasing=True)
                deseqPlotGeneHeatmap(
                    '%sgene_heatmap.png' % outfile_groups_prefix,
                    r['as.matrix'](data[order]),
                    Colv=False,
                    Rowv=True)
            else:
                E.warn('can not plot differentially expressed genes')
        else:
            E.warn('no differentially expressed genes at fdr %f' % fdr)

        # Plot pvalue histogram
        R.png('''%(outfile_groups_prefix)spvalue_histogram.png''' % locals())
        r('''pvalues = res$pval''')
        r('''hist(pvalues, breaks=50, col='skyblue' )''')
        r['dev.off']()

        # Plot diagnostic plots for FDR
        if has_replicates:
            r('''orderInPlot = order(pvalues)''')
            r('''showInPlot = (pvalues[orderInPlot] < 0.08)''')
            # Jethro - previously plotting x =
            # pvalues[orderInPlot][showInPlot]
            # pvalues[orderInPlot][showInPlot] contains all NA values
            # from pvalues which(showInPlot) doesn't... removing NA
            # values
            r('''true.pvalues <- pvalues[orderInPlot][showInPlot]''')
            r('''true.pvalues <- true.pvalues[is.finite(true.pvalues)]''')
            if r('''sum(showInPlot)''')[0] > 0:
                R.png('''%(outfile_groups_prefix)sfdr.png''' % locals())
                # failure when no replicates:
                # rpy2.rinterface.RRuntimeError:
                # Error in plot.window(...) : need finite 'xlim' values
                r('''plot( seq( along=which(showInPlot)),
                           true.pvalues,
                           pch='.',
                           xlab=expression(rank(p[i])),
                           ylab=expression(p[i]))''')
                r('''abline(a = 0, b = %(fdr)f / length(pvalues), col = "red")
                ''' % locals())
                r['dev.off']()
            else:
                E.warn('no p-values < 0.08')

        # Add log2 fold with variance stabilized l2fold value
        r('''res$transformed_log2FoldChange = vsd_l2f''')

        # Parse results and parse to file
        results, counts = deseqParseResults(control,
                                            treatment,
                                            fdr=fdr)

        all_results += results

        E.info(counts)
        outf = IOTools.openFile(
            "%(outfile_groups_prefix)ssummary.tsv" % locals(), "w")
        outf.write("category\tcounts\n%s\n" % counts.asTable())
        outf.close()

    writeExpressionResults(outfile, all_results)


def runDESeq2(outfile,
              outfile_prefix="deseq2",
              fdr=0.1,
              ref_group=None,
              model=None,
              contrasts=None,
              plot=1,
              ):
    """
    Run DESeq2 on counts table.

    If no model is passed, then defaults to the group column in design file

    Does not make use of group tag data bc function doesn't accomodate
    multi-factor designs

    To Do: Parse results into standard output format.
    KB: I have done this but I'm not sure if it is compatible with complex
    design tables
    Fix fact that plotMA is hardcoded.
    """

    # load libraries
    r('''suppressMessages(library('DESeq2'))''')

    # Create metadata... this will eventually be a pandas dataframe
    if isinstance(r('''factors'''), rpy2.robjects.vectors.DataFrame):
        E.info("DESeq2: Merging additional factors in design file to"
               "create metadata table")
        r('''mdata <- cbind(groups, factors)''')
        mdata = tuple(r('''names(mdata)'''))
    else:
        r('''mdata <- data.frame(group=groups)''')
        mdata = "group"
    E.info("DESeq2 colData headers are: %s" % mdata)

    # Check for model and that model terms are in metadata table
    if model:
        assert contrasts, "Must specifiy contrasts if model design provided"
        terms = set([x for x in re.split("\W", model) if x != ''])
        assert terms.issubset(mdata), \
            "DESeq2: design formula has terms not present in colData"
    else:
        if mdata != "group":
            E.warn("DESeq2 model specified, with no metadata in design file")
        terms = ["group", ]
        model = "~ group"
    E.info("DESeq2 design formula is: %s" % model)

    # Create DESeqDataSet, using countsTable, mdata, model
    r('''suppressMessages(dds <- DESeqDataSetFromMatrix(countData=countsTable,
                                       colData=mdata,
                                       design=%(model)s))''' % locals())
    # WARNING: This is not done automatically... I don't know why?
    r('''colnames(dds) <- colnames(countsTable)''')
    E.info("Combined colata, design formula and counts table to create"
           " DESeqDataSet instance")

    # Rlog transform
    r('''suppressMessages(rld <- rlog(dds))''')

    if plot == 1:
        # Plot PCA of rlog transformed count data for top 500
        for factor in terms:
            outf = outfile_prefix + factor + "_PCAplot500.tiff"
            E.info("Creating PCA plot for factor: %s" % outf)
            r('''x <- plotPCA(rld, intgroup="%(factor)s")''' % locals())
            # r('''saveRDS(x, '%(outf)s')''' % locals())
            r('''tiff("%(outf)s")''' % locals())
            r('''plot(x)''')
            r('''dev.off()''')

    # Extract rlog transformed count data...

    r('''rlogtab = as.data.frame(assay(rld))''')
    r('''rlogtab$test_id = rownames(rlogtab)''')
    r('''rlogtab = rlogtab[, c(ncol(rlogtab), 1:ncol(rlogtab)-1)]''')
    r('''rlogtab = as.data.frame(rlogtab)''')
    R.data('rlogtab')
    rlog_out = r('rlogtab')
    rlogoutf = outfile_prefix + "rlog.tsv"

    rlog_out.to_csv(rlogoutf, sep="\t", index=False)
    os.system("gzip %s" % rlogoutf)

    # Run DESeq2
    r('''suppressMessages(dds <- DESeq(dds))''')
    E.info("Completed DESeq2 differential expression analysis")

    # Extract contrasts...
    if contrasts:
        contrasts = (x.split(":") for x in contrasts.split(","))
    else:
        # created by loadTagData...
        groups = r('''levels(groups)''')
        contrasts = (("group",) + x for x in itertools.combinations(groups, 2))

    df_final = pandas.DataFrame()
    raw_final = pandas.DataFrame()
    all_results = []
    for combination in contrasts:
        variable, control, treatment = combination

        # Fetch results
        gfix = "%s_%s_vs_%s" % (variable, control, treatment)
        outfile_groups_prefix = outfile_prefix + gfix + "_MAplot.png"
        r('''res <- results(dds, contrast=c("%(variable)s",
                                            "%(treatment)s",
                                            "%(control)s"))''' % locals())
        E.info("Extracting contrast for levels %s (treatment) vs %s (control)"
               " for factor %s" % (treatment, control, variable))

        # plot MA plot
        if plot == 1:
            r('''png("%(outfile_groups_prefix)s")''' % locals())
            r('''plotMA(res, alpha=%f)''' % fdr)
            r('''dev.off()''')
            E.info("Plotted MA plot for levels %s (treatment) vs %s (control)"
                   " for factor %s" % (treatment, control, variable))
        r('''res_df <- as.data.frame(res)''')
        r('''res_df$test_id = rownames(res_df)''')
        r('''res_df = res_df[, c(ncol(res_df), 1:ncol(res_df)-1)]''')
        R.data('res_df')
        raw_out = r('res_df')

        #  manipulate output into standard format
        df_out = deseq2ParseResults(treatment, control, fdr, vsd=False)

        #  label the deseq2 raw output file and append it to the raw output tab
        raw_out["treatment"] = [treatment, ]*len(df_out.index)
        raw_out["control"] = [control, ]*len(df_out.index)
        raw_out["variable"] = [variable, ]*len(df_out.index)
        raw_final = raw_final.append(raw_out, ignore_index=True)

        #  write the standardised output table
        df_out.to_csv(IOTools.openFile(outfile_prefix + gfix + ".tsv.gz", "w"),
                      sep="\t",
                      index_label="test_id")

        E.info("Extracted results table for contrast  '%s' (treatment) vs '%s'"
               " (control) for factor '%s'" % (treatment, control, variable))

        # append to final dataframe
        df_out.reset_index(inplace=True)
        df_out.rename(columns={"index": "test_id"}, inplace=True)
        df_final = df_final.append(df_out, ignore_index=True)

    results = df_final.values.tolist()

    # write final dataframe into standard format
    writeExpressionResults(outfile, results)

    rawoutf = outfile_prefix + "raw.tsv"
    raw_final.to_csv(rawoutf, sep="\t", index=False)
    os.system("gzip %s" % rawoutf)

Design = collections.namedtuple("Design", ("include", "group", "pair"))


def readDesignFile(design_file):
    '''reads a design file.'''

    design = collections.OrderedDict()
    with IOTools.openFile(design_file) as inf:
        for line in inf:
            if line.startswith("track"):
                continue
            track, include, group, pair = line.split("\t")[:4]
            if track in design:
                raise ValueError("duplicate track '%s'" % track)
            design[track] = Design._make((int(include), group, pair))
    return design


def plotTagStats(infile, design_file, outfile_prefix):
    '''provide summary plots for tag data.'''

    loadTagData(infile, design_file)

    nobservations, nsamples = filterTagData()

    if nobservations == 0:
        E.warn("no observations - no output")
        return

    if nsamples == 0:
        E.warn("no samples remain after filtering - no output")
        return

    groups, pairs, has_replicates, has_pairs = groupTagData()

    # import rpy2.robjects.lib.ggplot2 as ggplot2

    r('''suppressMessages(library('ggplot2'))''')
    r('''suppressMessages(library('reshape'))''')

    r('''d = melt( log10(countsTable + 1), variable_name = 'sample' )''')

    # Note that ggsave does not work if there is
    # X display.
    R.png(outfile_prefix + ".densities.png")
    r('''gp = ggplot(d)''')
    r('''pp = gp + geom_density(aes(x=value, group=sample,
    color=sample, fill=sample), alpha=I(1/3))''')
    r('''plot(pp)''')
    r['dev.off']()

    R.png(outfile_prefix + ".boxplots.png")
    r('''gp = ggplot(d)''')
    r('''pp = gp +
    geom_boxplot(aes(x=sample,y=value,color=sample,fill=sample),
    size=0.3,
    alpha=I(1/3)) +
    theme(axis.text.x = element_text( angle=90, hjust=1, size=8 ) )''')
    r('''plot(pp)''')
    r['dev.off']()


def plotDETagStats(infile, outfile_prefix,
                   additional_file=None,
                   join_columns=None,
                   additional_columns=None):
    '''provide summary plots for tag data.

    Stratify boxplots and densities according to differential
    expression calls.

    The input file is the output of any of the DE
    tools, see GeneExpressionResults for column names.

    Additional file will be joined with infile and any additional
    columns will be output as well.
    '''

    table = pandas.read_csv(IOTools.openFile(infile),
                            sep="\t")

    if additional_file is not None:
        additional_table = pandas.read_csv(
            IOTools.openFile(additional_file),
            sep="\t")
        table = pandas.merge(table,
                             additional_table,
                             on=join_columns,
                             how="left",
                             sort=False)

    # remove index. If it is numbered starting from 1, there is a bug
    # in ggplot, see https://github.com/yhat/ggplot/pull/384
    table.reset_index(inplace=True)

    # add log-transformed count data
    table['log10_treatment_mean'] = numpy.log10(table['treatment_mean'] + 1)
    table['log10_control_mean'] = numpy.log10(table['control_mean'] + 1)
    table['dmr'] = numpy.array(["insignicant"] * len(table))
    table.loc[
        (table["l2fold"] > 0) & (table["significant"] == 1), "dmr"] = "up"
    table.loc[
        (table["l2fold"] < 0) & (table["significant"] == 1), "dmr"] = "down"

    def _dplot(table, outfile, column):

        plot = ggplot.ggplot(
            ggplot.aes(column,
                       colour='dmr',
                       fill='dmr'),
            data=table) + \
            ggplot.geom_density(alpha=0.5)

        try:
            plot.save(filename=outfile)
        except Exception as msg:
            E.warn("no plot for %s: %s" % (column, msg))

    def _bplot(table, outfile, column):

        plot = ggplot.ggplot(
            ggplot.aes(x='dmr', y=column),
            data=table) + \
            ggplot.geom_boxplot()
        try:
            plot.save(filename=outfile)
        except ValueError as msg:
            # boxplot fails if all values are the same
            # see https://github.com/yhat/ggplot/issues/393
            E.warn(msg)

    _dplot(table,
           outfile_prefix + ".densities_tags_control.png",
           "log10_control_mean")
    _dplot(table,
           outfile_prefix + ".densities_tags_treatment.png",
           "log10_treatment_mean")
    _bplot(table,
           outfile_prefix + ".boxplot_tags_control.png",
           "log10_control_mean")
    _bplot(table,
           outfile_prefix + ".boxplot_tags_treatment.png",
           "log10_treatment_mean")

    if additional_columns:
        for column in additional_columns:
            _dplot(table,
                   outfile_prefix + ".densities_%s.png" % column,
                   column)
            _bplot(table,
                   outfile_prefix + ".boxplot_%s.png" % column,
                   column)
    return


def runMockAnalysis(outfile,
                    outfile_prefix,
                    ref_group=None,
                    ref_regex=None,
                    pseudo_counts=0):
    '''run a mock analysis on a count table.

    compute fold enrichment values, but do not normalize or
    perform any test.
    '''

    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group,
                                                            ref_regex)

    all_results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination

        r('''control_counts = rowSums( countsTable[groups == '%s'] )''' %
          control)
        r('''treatment_counts = rowSums( countsTable[groups == '%s'] )''' %
          treatment)

        # add pseudocounts to enable analysis of regions
        # that are absent/present
        if pseudo_counts:
            r('''control_counts = control_counts + %f''' % pseudo_counts)
            r('''treatment_counts = treatment_counts + %f''' % pseudo_counts)

        r('''fc = treatment_counts / control_counts''')

        results = []

        for identifier, treatment_count, control_count, foldchange in \
                zip(r('''rownames( countsTable)'''),
                    r('''treatment_counts'''),
                    r('''control_counts'''),
                    r('''fc''')):
            try:
                log2fold = math.log(foldchange)
            except ValueError:
                log2fold = "Inf"

            results.append(GeneExpressionResult._make((
                identifier,
                treatment,
                treatment_count,
                0,
                control,
                control_count,
                0,
                1,
                1,
                log2fold,
                foldchange,
                log2fold,
                "0",
                "OK")))

        all_results.extend(results)

    writeExpressionResults(outfile, all_results)


def outputTagSummary(filename_tags,
                     outfile,
                     output_filename_pattern,
                     filename_design=None):
    '''output summary values for a count table.'''

    E.info("loading tag data from %s" % filename_tags)

    if filename_design is not None:
        # load all tag data
        loadTagData(filename_tags, filename_design)

        # filter
        nobservations, nsamples = filterTagData()

    else:
        # read complete table
        r('''countsTable = read.delim('%(filename_tags)s',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = TRUE,
        comment.char = '#')''' % locals())

        nobservations, nsamples = tuple(r('''dim(countsTable)'''))
        E.info("read data: %i observations for %i samples" %
               (nobservations, nsamples))
        # remove samples without data
        r('''max_counts = apply(countsTable,2,max)''')

        filter_min_counts_per_sample = 1

        empty_samples = tuple(
            r('''max_counts < %i''' % filter_min_counts_per_sample))
        sample_names = r('''colnames(countsTable)''')
        nempty_samples = sum(empty_samples)

        if nempty_samples:
            E.warn("%i empty samples are being removed: %s" %
                   (nempty_samples,
                    ",".join([sample_names[x]
                              for x, y in enumerate(empty_samples) if y])))
            r('''countsTable <- countsTable[, max_counts >= %i]''' %
              filter_min_counts_per_sample)
            nobservations, nsamples = tuple(r('''dim(countsTable)'''))

        r('''groups = factor(colnames( countsTable ))''')
        E.debug("sample names: %s" % r('''colnames(countsTable)'''))

    nrows, ncolumns = tuple(r('''dim(countsTable)'''))

    outfile.write("metric\tvalue\tpercent\n")
    outfile.write("number of observations\t%i\t100\n" % nobservations)
    outfile.write("number of samples\t%i\t100\n" % nsamples)

    # Count windows with no data
    r('''max_counts = apply(countsTable,1,max)''')

    # output distribution of maximum number of counts per window
    outfilename = output_filename_pattern + "max_counts.tsv.gz"
    E.info("outputting maximum counts per window to %s" % outfilename)
    r('''write.table(table(max_counts),
    file=gzfile('%(outfilename)s'),
    sep="\t",
    row.names=FALSE,
    quote=FALSE)''' %
      locals())

    # removing empty rows
    E.info("removing rows with no counts in any sample")
    r('''countsTable = countsTable[max_counts>0,]''')

    for x in range(0, 20):
        nempty = tuple(r('''sum(max_counts <= %i)''' % x))[0]
        outfile.write("max per row<=%i\t%i\t%f\n" %
                      (x, nempty, 100.0 * nempty / nrows))

    E.info("removed %i empty rows" % tuple(r('''sum(max_counts == 0)''')))
    observations, samples = tuple(r('''dim(countsTable)'''))
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # build correlation
    r('''correlations = cor(countsTable)''')
    outfilename = output_filename_pattern + "correlation.tsv"
    E.info("outputting sample correlations table to %s" % outfilename)
    r('''write.table(correlations, file='%(outfilename)s',
    sep="\t",
    row.names=TRUE,
    col.names=NA,
    quote=FALSE)''' % locals())

    # output scatter plots
    outfilename = output_filename_pattern + "scatter.png"
    E.info("outputting scatter plots to %s" % outfilename)
    R.png(outfilename, width=960, height=960)
    plotPairs()
    r['dev.off']()

    # output heatmap based on correlations
    outfilename = output_filename_pattern + "heatmap.svg"
    E.info("outputting correlation heatmap to %s" % outfilename)
    R.svg(outfilename)
    plotCorrelationHeatmap(method="correlation")
    r['dev.off']()

    # output PCA
    outfilename = output_filename_pattern + "pca.svg"
    E.info("outputting PCA plot to %s" % outfilename)
    R.svg(outfilename)
    plotPCA(groups=False)
    r['dev.off']()

    # output an MDS plot
    r('''suppressMessages(library('limma'))''')
    outfilename = output_filename_pattern + "mds.svg"
    E.info("outputting mds plot to %s" % outfilename)
    R.svg(outfilename)
    try:
        r('''plotMDS(countsTable)''')
    except RRuntimeError:
        E.warn("can not plot mds")
    r['dev.off']()


def dumpTagData(filename_tags, filename_design, outfile):
    '''output filtered tag table.'''

    if outfile == sys.stdout:
        outfilename = ""
    else:
        outfilename = outfile.name

    # load all tag data
    loadTagData(filename_tags, filename_design)

    # filter
    nobservations, nsamples = filterTagData()

    # output
    r('''write.table( countsTable,
                      file='%(outfilename)s',
                      sep='\t',
                      quote=FALSE)''' % locals())


def runTTest(outfile,
             outfile_prefix,
             fdr=0.1,
             ref_group=None,
             ref_regex=None):
    '''apply a ttest on the data.

    For the T-test it is best to use FPKM values as
    this method does not perform any library normalization.
    '''
    groups, pairs, has_replicates, has_pairs = groupTagData(ref_group,
                                                            ref_regex)

    results = []
    for combination in itertools.combinations(groups, 2):
        control, treatment = combination
        r = r('''r = apply(countsTable, 1,
        function(x) { t.test(
        x[groups == '%(treatment)s'],
        x[groups == '%(control)s']) } )
        ''' % locals())

        for test_id, ttest in zip(r.names, r):
            # TS, swapped order below as assignment was incorrect
            treatment_mean, control_mean = tuple(ttest.rx2('estimate'))
            fold_change = treatment_mean / control_mean
            pvalue = tuple(ttest.rx2('p.value'))[0]
            significant = (0, 1)[pvalue < fdr]
            results.append(GeneExpressionResult._make((test_id,
                                                       treatment,
                                                       treatment_mean,
                                                       0,
                                                       control,
                                                       control_mean,
                                                       0,
                                                       pvalue,
                                                       pvalue,
                                                       numpy.log2(fold_change),
                                                       fold_change,
                                                       numpy.log2(fold_change),
                                                       significant,
                                                       "OK")))

    writeExpressionResults(outfile, results)


#####################################################################
# Pandas-based functions and matplotlib-based plotting functions ####
#####################################################################

def loadTagDataPandas(tags_filename, design_filename):
    '''load tag data for deseq/edger analysis.

    *Infile* is a tab-separated file with counts.

    *design_file* is a tab-separated file with the
    experimental design with four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

    track
        name of track - should correspond to column header in *infile*
    include
        flag to indicate whether or not to include this data
    group
        group indicator - experimental group
    pair
        pair that sample belongs to (for paired tests)

    This method creates various R objects:

    countsTable : data frame with counts.
    groups : vector with groups
    pairs  : vector with pairs

    '''

    E.info("loading tag data from %s" % tags_filename)

    inf = IOTools.openFile(tags_filename)
    counts_table = pandas.read_csv(inf,
                                   sep="\t",
                                   index_col=0,
                                   comment="#")
    inf.close()

    E.info("read data: %i observations for %i samples" %
           counts_table.shape)

    E.debug("sample names: %s" % list(counts_table.columns))

    inf = IOTools.openFile(design_filename)
    design_table = pandas.read_csv(inf, sep="\t", index_col=0)
    inf.close()

    E.debug("design names: %s" % list(design_table.index))

    missing = set(counts_table.columns).difference(design_table.index)

    if missing:
        E.warn("missing samples from design file are ignored: %s" % missing)

    # remove unnecessary samples
    design_table = design_table[design_table["include"] != 0]
    E.debug("included samples: %s" % list(design_table.index))

    counts_table = counts_table[list(design_table.index)]
    E.info("filtered data: %i observations for %i samples" %
           counts_table.shape)

    return counts_table, design_table


def filterTagDataPandas(counts_table,
                        design_table,
                        filter_min_counts_per_row=1,
                        filter_min_counts_per_sample=10,
                        filter_percentile_rowsums=0):
    '''filter tag data.

    * remove rows with at least x number of counts

    * remove samples with a maximum of *min_sample_counts*

    * remove the lowest percentile of rows in the table, sorted
       by total tags per row
    '''

    # Remove windows with no data
    max_counts_per_row = counts_table.max(1)
    counts_table = counts_table[
        max_counts_per_row >= filter_min_counts_per_row]
    observations, samples = counts_table.shape
    E.info("trimmed data: %i observations for %i samples" %
           (observations, samples))

    # remove samples without data
    max_counts_per_sample = counts_table.max()

    empty_samples = max_counts_per_sample < filter_min_counts_per_sample
    sample_names = counts_table.columns
    nempty_samples = sum(empty_samples)

    if nempty_samples:
        E.warn("%i empty samples are being removed: %s" %
               (nempty_samples,
                ",".join([sample_names[x] for x, y in

                          enumerate(empty_samples) if y])))
        raise NotImplementedError("removing empty samples needs to be done")
        # r('''countsTable <- countsTable[, max_counts >= %i]''' % filter_min_counts_per_sample)
        # r('''groups <- groups[max_counts >= %i]''' % filter_min_counts_per_sample)
        # r('''pairs <- pairs[max_counts >= %i]''' % filter_min_counts_per_sample)
        # observations, samples = tuple( r('''dim(countsTable)'''))

    # percentile filtering
    if filter_percentile_rowsums > 0:
        percentile = float(filter_percentile_rowsums) / 100.0
        sum_counts = counts_table.sum(1)
        take = sum_counts > sum_counts.quantile(percentile)
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (filter_percentile_rowsums,
                sum(take),
                len(take) - sum(take)))
        counts_table = counts_table[take]

    return counts_table


def identifyVariablesPandas(design_table):

    # design table should have been processed by loadTagDataPandas already
    # just in case, re-filter for not included samples here

    design_table = design_table[design_table["include"] != 0]
    conds = design_table['group'].tolist()
    pairs = design_table['pair'].tolist()

    # TS, adapted from JJ code for DESeq2 design tables:
    # if additional columns present, pass to 'factors'
    if len(design_table.columns) > 3:
        factors = design_table.iloc[:, 3:]
    else:
        factors = None

    return conds, pairs, factors


def checkTagGroupsPandas(design_table, ref_group=None):
    '''compute groups and pairs from tag data table.'''

    conds, pairs, factors = identifyVariablesPandas(design_table)
    groups = list(set(conds))

    # Relevel the groups so that the reference comes first
    # how to do this in python?
    # if ref_group is not None:
    #    r('''groups <- relevel(groups, ref = "%s")''' % ref_group)

    # check this works, will need to make factors from normal df
    # TS adapted from JJ code for DESeq2 -
    # check whether there are additional factors in design file...
    if factors:
        E.warn("There are additional factors in design file that are ignored"
               " by groupTagData: ", factors)
    else:
        pass

    # Test if replicates exist - at least one group must have multiple samples
    max_per_group = max([conds.count(x) for x in groups])

    has_replicates = max_per_group >= 2

    # Test if pairs exist:
    npairs = len(set(pairs))
    has_pairs = npairs == 2

    # ..if so, at least two samples are required per pair
    if has_pairs:
        min_per_pair = min([pairs.count(x) for x in set(pairs)])
        has_pairs = min_per_pair >= 2

    return groups, pairs, conds, factors, has_replicates, has_pairs


ResultColumns = ["test_id", "treatment_name", "treatment_mean",
                 "treatment_std", "control_name", "control_mean",
                 "control_std", "p_value", "p_value_adj", "l2fold", "fold",
                 "transformed_l2fold", "significant", "status"]

ResultColumns_dtype = {"test_id": object, "treatment_name": object,
                       "treatment_mean": float, "treatment_std":
                       float, "control_name": object, "control_mean":
                       float, "control_std": float, "p_value": float,
                       "p_value_adj": float, "l2fold": float, "fold":
                       float, "transformed_l2fold": float,
                       "significant": int, "status": object}


def makeEmptyDataFrameDict():
    return {key: [] for key in ResultColumns}


def runTTestPandas(counts_table,
                   design_table,
                   outfile,
                   outfile_prefix,
                   fdr,
                   ref_group=None):
    '''apply a ttest on the data.

    For the T-test it is best to use FPKM values as
    this method does not perform any library normalization.
    Alternatively, perform normalisation on counts table using Counts.py
    '''
    stats = importr('stats')

    (groups, pairs, conds, factors, has_replicates,
     has_pairs) = checkTagGroupsPandas(design_table, ref_group)

    df_dict = makeEmptyDataFrameDict()

    for combination in itertools.combinations(groups, 2):
        # as each combination may have different numbers of samples in control
        # and treatment, calculations have to be performed on a per
        # combination basis

        control, treatment = combination
        n_rows = counts_table.shape[0]
        df_dict["control_name"].extend((control,)*n_rows)
        df_dict["treatment_name"].extend((treatment,)*n_rows)
        df_dict["test_id"].extend(counts_table.index.tolist())

        # subset counts table for each combination
        c_keep = [x == control for x in conds]
        control_counts = counts_table.iloc[:, c_keep]
        t_keep = [x == treatment for x in conds]
        treatment_counts = counts_table.iloc[:, t_keep]

        c_mean = control_counts.mean(axis=1)
        df_dict["control_mean"].extend(c_mean)
        df_dict["control_std"].extend(control_counts.std(axis=1))

        t_mean = treatment_counts.mean(axis=1)
        df_dict["treatment_mean"].extend(t_mean)
        df_dict["treatment_std"].extend(treatment_counts.std(axis=1))

        t, prob = ttest_ind(control_counts, treatment_counts, axis=1)
        df_dict["p_value"].extend(prob)
        # what about zero values?!
        df_dict["fold"].extend(t_mean / c_mean)

    df_dict["p_value_adj"].extend(
        list(stats.p_adjust(FloatVector(df_dict["p_value"]), method='BH')))

    df_dict["significant"].extend(
        [int(x < fdr) for x in df_dict["p_value_adj"]])

    df_dict["l2fold"].extend(list(numpy.log2(df_dict["fold"])))
    # note: the transformed log2 fold change is not transformed!
    df_dict["transformed_l2fold"].extend(list(numpy.log2(df_dict["fold"])))

    # set all status values to "OK"
    df_dict["status"].extend(("OK",)*n_rows)

    results = pandas.DataFrame(df_dict)
    results.set_index("test_id", inplace=True)
    results.to_csv(outfile, sep="\t", header=True, index=True)


def plotCorrelationHeatmapMatplot(counts, outfile, method="correlation",
                                  cor_method="pearson"):
    '''plot a heatmap of correlations derived from
    countsTable.
    '''
    # to do: add other methods?

    # define outside function? - Will we reuse?
    heatmap_cdict_b_to_y = {
        'red': ((0.0, 0.4, .4), (0.01, .4, .4), (1., .95, .95)),
        'green': ((0.0, 0.4, 0.4), (0.01, .4, .4), (1., .95, .95)),
        'blue': ((0.0, .9, .9), (0.01, .9, .9), (1., 0.4, 0.4))}

    cm = matplotlib.colors.LinearSegmentedColormap(
        '', heatmap_cdict_b_to_y, 256)

    df = counts.corr(method=cor_method)
    plt.pcolor(np.array(df), cmap=cm)
    plt.colorbar()
    plt.title("%(cor_method)s correlation heatmap" % locals())
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    plt.tight_layout()
    plt.savefig(outfile)


def runEdgeRPandas(counts,
                   design_table,
                   outfile,
                   outfile_prefix="edger.",
                   fdr=0.1,
                   prefix="",
                   dispersion=None,
                   ref_group=None):
    '''run EdgeR on countsTable.

    Results are stored in *outfile* and files prefixed by *outfile_prefix*.

    The dispersion is usually measuered from replicates. If there are no
    replicates, you need to set the *dispersion* explicitely.

    See page 13 of the EdgeR user guide::

       2. Simply pick a reasonable dispersion value, based on your
       experience with similar data, and use that. Although
       subjective, this is still more defensible than assuming Poisson
       variation. Typical values are dispersion=0.4 for human data,
       dispersion=0.1 for data on genetically identical model
       organisms or dispersion=0.01 for technical replicates.

    '''

    # load library
    r('''suppressMessages(library('edgeR'))''')

    (groups, pairs, conds, factors, has_replicates,
     has_pairs) = checkTagGroupsPandas(design_table, ref_group)

    if not has_replicates and dispersion is None:
        raise ValueError("no replicates and no dispersion")

    # output heatmap plot
    plotCorrelationHeatmapMatplot(counts,
                                  '%(outfile_prefix)sheatmap.png' % locals(),
                                  cor_method="spearman")

    E.info('running EdgeR: groups=%s, pairs=%s, replicates=%s, pairs=%s' %
           (groups, pairs, has_replicates, has_pairs))

    r_counts = pandas2ri.py2ri(counts)

    passDFtoRGlobalEnvironment = r('''function(df){
    countsTable <<- df}''')
    passDFtoRGlobalEnvironment(r_counts)

    if has_pairs:
        # output difference between groups
        # TS #####
        # this is performed on non-normalised data
        # should we use Counts.py to normalise first?
        # also, this isn't edgeR specific, should this be
        # moved to a seperate summary function?
        # also move the MDS plotting?
        # #####
        first = True
        pairs_df = pandas.DataFrame()
        nrows = len(counts.index)
        n = 0
        for g1, g2 in itertools.combinations(groups, 2):

            keep_a = [x == g1 for x in conds]
            counts_a = counts.iloc[:, keep_a]
            keep_b = [x == g2 for x in conds]
            counts_b = counts.iloc[:, keep_b]
            index = range(n, n+nrows)
            n += nrows
            a = counts_a.sum(axis=1)
            b = counts_b.sum(axis=1)
            diff = a-b
            diff.sort()
            temp_df = pandas.DataFrame({"cumsum": np.cumsum(diff).tolist(),
                                        "comb": "_vs_".join([g1, g2]),
                                        "id": range(0, nrows)},
                                       index=index)
            pairs_df = pairs_df.append(temp_df)
        plot_pairs = r('''function(df, outfile){
        suppressMessages(library('ggplot2'))
        p = ggplot(df, aes(y=cumsum, x=id)) +
        geom_line(aes(col=factor(comb))) +
        scale_color_discrete(name="Comparison") +
        xlab("index") + ylab("Cumulative sum")
        ggsave("%(outfile_prefix)sbalance_groups.png", plot = p)}
        ''' % locals())

        r_pairs_df = pandas2ri.py2ri(pairs_df)
        plot_pairs(r_pairs_df)

        # output difference between pairs within groups
        first = True
        legend = []
        n = 0
        pairs_in_groups_df = pandas.DataFrame()
        for pair in set(pairs):
            for g1, g2 in itertools.combinations(groups, 2):

                key = "pair-%s-%s-vs-%s" % (pair, g1, g2)
                legend.append(key)
                keep_a = [x == pair and y == g1 for x, y in zip(pairs, conds)]
                keep_b = [x == pair and y == g2 for x, y in zip(pairs, conds)]

                # check pair and group combination actually exists
                if sum(keep_a) > 0 and sum(keep_b) > 0:
                    counts_a = counts.iloc[:, keep_a]
                    counts_b = counts.iloc[:, keep_b]
                    index = range(n, n+nrows)
                    n += nrows
                    a = counts_a.sum(axis=1)
                    b = counts_b.sum(axis=1)
                    diff = a-b
                    diff.sort()
                    comparison = "pair-%s-%s-vs-%s" % (pair, g1, g2)
                    temp_df = pandas.DataFrame({"cumsum": np.cumsum(diff).tolist(),
                                                "comb": comparison,
                                                "id": range(0, nrows)},
                                               index=index)
                    pairs_in_groups_df = pairs_in_groups_df.append(temp_df)

        plot_pairs_in_groups = r('''function(df, outfile){
        suppressMessages(library('ggplot2'))
        p = ggplot(df, aes(y=cumsum, x=id)) +
        geom_line(aes(col=factor(comb))) +
        scale_color_discrete(name="Comparison") +
        xlab("index") + ylab("Cumulative sum")
        ggsave("%(outfile_prefix)sbalance_pairs.png", plot = p)}
        ''' % locals())

        r_pairs_in_groups_df = pandas2ri.py2ri(pairs_in_groups_df)
        plot_pairs_in_groups(r_pairs_in_groups_df)

    # create r objects
    r_counts = pandas2ri.py2ri(counts)
    r_groups = ro.StrVector(conds)
    r_pairs = ro.StrVector(pairs)
    r_has_pairs = ro.default_py2ri(has_pairs)
    r_has_replicates = ro.default_py2ri(has_replicates)

    if dispersion is not None:
        r_dispersion = ro.default_py2ri(dispersion)
    else:
        r_dispersion = ro.default_py2ri(False)

    if ref_group is not None:
        r_ref_group = ro.default_py2ri(ref_group)
    else:
        r_ref_group = ro.default_py2ri(groups[0])

    # build DGEList object
    buildDGEList = r('''
    suppressMessages(library('edgeR'))
    function(counts, groups, ref_group){
    countsTable = DGEList(counts, group=groups)
    countsTable$samples$group <- relevel(countsTable$samples$group,
    ref = ref_group)
    countsTable = calcNormFactors(countsTable)
    return(countsTable)
    }''')

    r_countsTable = buildDGEList(r_counts, r_groups, r_ref_group)

    # output MDS plot
    # TT - should this be performed on the normalised counts table?(see above)
    R.png('''%(outfile_prefix)smds.png''' % locals())
    try:
        MDSplot = r('''function(counts){
        plotMDS(counts)}''')
        MDSplot(r_counts)
    except RRuntimeError:
        E.warn("can not plot mds")
    r['dev.off']()

    # build design matrix
    buildDesign = r('''function(countsTable, has_pairs, pairs){
    if (has_pairs==TRUE) {
        design <- model.matrix( ~pairs + countsTable$samples$group ) }
    else {
        design <- model.matrix( ~countsTable$samples$group ) }
    return(design)
    }''')

    r_design = buildDesign(r_countsTable, r_has_pairs, r_pairs)

    fitModel = r('''function(countsTable, design, has_replicates, dispersion){
    if (has_replicates == TRUE) {
        # estimate common dispersion
        countsTable = estimateGLMCommonDisp( countsTable, design )
        # estimate tagwise dispersion
        countsTable = estimateGLMTagwiseDisp( countsTable, design )
        # fitting model to each tag
        fit = glmFit( countsTable, design ) }
    else {
        # fitting model to each tag
        fit = glmFit(countsTable, design, dispersion=dispersion) }
    return(fit)
    }''')

    r_fit = fitModel(r_countsTable, r_design, r_has_replicates, r_dispersion)

    E.info("Generating output")

    # perform LR test
    lrtTest = r('''function(fit, prefix){
    lrt = glmLRT(fit)
    # save image for access to the whole of the lrt object
    save.image(paste0(prefix,"lrt.RData"))
    return(lrt)
    }''')
    r_lrt = lrtTest(r_fit, outfile_prefix)

    # return statistics table - must be a better way to do this?
    extractTable = r('''function(lrt){ return(lrt$table)}''')
    r_lrt_table = extractTable(r_lrt)

    # output cpm table
    outputCPMTable = r('''function(countsTable, outfile_prefix, lrt){
    suppressMessages(library(reshape2))
    countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)
    melted <- melt(countsTable.cpm)
    names(melted) <- c("test_id", "sample", "ncpm")

    # melt columns are factors - convert to string for sorting
    melted$test_id = levels(melted$test_id)[as.numeric(melted$test_id)]
    melted$sample = levels(melted$sample)[as.numeric(melted$sample)]

    # sort cpm table by test_id and sample
    sorted <- melted[with(melted, order(test_id, sample)),]
    gz <- gzfile(paste0(outfile_prefix,"cpm.tsv.gz"), "w" )
    write.table(sorted, file=gz, sep = "\t", row.names=FALSE, quote=FALSE)
    close(gz)}''')

    outputCPMTable(r_countsTable, outfile_prefix, r_lrt)

    # output differences between pairs
    if len(groups) == 2:
        plotSmear = r('''function(countsTable, outfile_prefix){
        png(paste0(outfile_prefix,"smaplot.png"
        plotSmear(countsTable, pair=c('%s'))
        abline(h=c(-2, 2), col='dodgerblue')
        dev.off()}''' % "','".join(groups))

    lrt_table = pandas2ri.ri2py(r_lrt_table)

    n_rows = lrt_table.shape[0]
    df_dict = makeEmptyDataFrameDict()

    # import r stats module for BH adjustment
    stats = importr('stats')

    df_dict["control_name"].extend((groups[0],)*n_rows)
    df_dict["treatment_name"].extend((groups[1],)*n_rows)
    df_dict["test_id"].extend(lrt_table.index)
    df_dict["control_mean"].extend(lrt_table['logCPM'])
    df_dict["treatment_mean"].extend(lrt_table['logCPM'])
    df_dict["control_std"].extend((0,)*n_rows)
    df_dict["treatment_std"].extend((0,)*n_rows)
    df_dict["p_value"].extend(lrt_table['PValue'])
    df_dict["p_value_adj"].extend(
        list(stats.p_adjust(FloatVector(df_dict["p_value"]), method='BH')))
    df_dict["significant"].extend(
        [int(float(x) < fdr) for x in df_dict["p_value_adj"]])
    df_dict["l2fold"].extend(list(numpy.log2(lrt_table['logFC'])))

    # TS -note: the transformed log2 fold change is not transformed!
    df_dict["transformed_l2fold"].extend(list(numpy.log2(lrt_table['logFC'])))

    # TS -note: check what happens when no fold change is available
    # TS -may need an if/else in list comprehension. Raise E.warn too?
    df_dict["fold"].extend([math.pow(2, float(x)) for x in lrt_table['logFC']])

    # set all status values to "OK"
    # TS - again, may need an if/else, check...
    df_dict["status"].extend(("OK",)*n_rows)

    results = pandas.DataFrame(df_dict)
    results.set_index("test_id", inplace=True)
    results.to_csv(outfile, sep="\t", header=True, index=True)

    counts = E.Counter()
    counts.signficant = sum(results['significant'])
    counts.insignficant = (len(results['significant']) - counts.signficant)
    counts.all_over = sum([x > 0 for x in results['l2fold']])
    counts.all_under = sum([x < 0 for x in results['l2fold']])
    counts.signficant_over = sum([results['significant'][x] == 1 and
                                 results['l2fold'][x] > 1 for
                                 x in range(0, n_rows)])
    counts.signficant_under = sum([results['significant'][x] == 1 and
                                   results['l2fold'][x] < 1 for
                                   x in range(0, n_rows)])

    outf = IOTools.openFile("%(outfile_prefix)ssummary.tsv" % locals(), "w")
    outf.write("category\tcounts\n%s\n" % counts.asTable())
    outf.close()


def outputSpikeIns(filename_tags,
                   outfile,
                   output_filename_pattern,
                   filename_design=None,
                   foldchange_max=10.0,
                   expression_max=5.0,
                   max_counts_per_bin=100,
                   expression_bin_width=0.5,
                   foldchange_bin_width=0.5,
                   iterations=1):

    E.info("loading tag data from %s" % filename_tags)

    if filename_design is not None:
        # load all tag data
        counts_table, design_table = loadTagDataPandas(
            filename_tags, filename_design)

        # filter
        counts_table = filterTagDataPandas(counts_table, design_table)

    else:
        raise NotImplementedError("reading full table not implemented")

    nobservations, nsamples = counts_table.shape

    groups = list(set(design_table["group"]))
    if len(groups) != 2:
        raise NotImplementedError(
            "spike in only implemented for one pairwise comparison")

    # select group data
    group1 = counts_table[design_table[design_table.group == groups[0]].index]
    group2 = counts_table[design_table[design_table.group == groups[1]].index]

    outfile.write("interval_id\t%s\t%s\n" %
                  ("\t".join(group1.columns), "\t".join(group2.columns)))
    outf_info = IOTools.openFile(output_filename_pattern + "info.tsv.gz", "w")
    outf_info.write("interval_id\tl10_expression\tl2fold\n")

    # important: convert to matrixes, otherwise there will be a per-name lookup
    # when l10average or l2fold are computed.
    group1 = group1.as_matrix()
    group2 = group2.as_matrix()

    # compute bins
    expression_bins = numpy.arange(0, expression_max, expression_bin_width)
    fold_change_bins = numpy.arange(-foldchange_max,
                                    foldchange_max, foldchange_bin_width)

    E.info("l10expression bins=%s" % expression_bins)
    E.info("l2fold change bins=%s" % fold_change_bins)

    # output values
    output_counts = numpy.zeros(
        (len(expression_bins) + 1, len(fold_change_bins) + 1))

    for iteration in range(iterations):

        # randomize order
        group1 = numpy.random.permutation(group1)
        group2 = numpy.random.permutation(group2)

        # compute means and foldchanges
        group1_mean = group1.mean(1)
        group2_mean = group2.mean(1)

        # compute average expression by means
        l10average = numpy.log((group1_mean + group2_mean) / 2.0)

        # compute a fold change with pseudo count of 1
        l2fold = numpy.log2((group1_mean + 1.0) / (group2_mean + 1.0))

        # digitize arrays with bins
        l10average_idx = numpy.digitize(l10average, expression_bins)
        l2fold_idx = numpy.digitize(l2fold, fold_change_bins)

        interval_id = 0
        for idx, coord in enumerate(zip(l10average_idx, l2fold_idx)):
            # assert expression_bins[coord[0]-1] <= l10average[idx] < expression_bins[coord[0]]
            # assert fold_change_bins[coord[1]-1] <= l2fold[idx] <  expression_bins[coord[1]]

            if output_counts[coord] >= max_counts_per_bin:
                continue

            output_counts[coord] += 1
            outf_info.write("spike%i\t%f\t%f\n" % (interval_id,
                                                   l10average[idx],
                                                   l2fold[idx]))

            outfile.write("spike%i\t%s\t%s\n" %
                          (interval_id,
                           "\t".join(
                               map(str, list(group1[idx]))),
                           "\t".join(map(str, list(group2[idx])))))
            interval_id += 1

    outf_info.close()

    df = pandas.DataFrame(output_counts)
    df.index = list(expression_bins) + [expression_max]
    df.columns = list(fold_change_bins) + [foldchange_max]

    df.to_csv(IOTools.openFile(output_filename_pattern + "counts.tsv.gz", "w"),
              sep="\t")

    E.info("output %i spike in intervals" % interval_id)
