###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
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
###############################################################################
"""
==============================
Pipeline timeseries_clusters2metrics.py
================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline for analysing RNAseq timeseries gene expression data

Overview
========

A pipeline to perform time-dependent differential expression and
hierarchical clustering of timeseries RNA-seq data.
Input files are bam files and reference transcriptome
GTF(s).  This pipeline generates, via a replicate resampling
iterative clustering approach or direct computation on sample replicates,
clusters of genes using a consensus clustering.  Consensus clustering metrics
are also produced.  Cluster eigengenes are derived as well
as loading plots for each gene/cluster.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

targets are:

diffExpression - perform time-point and condition differential expression
analysis

clustering - perform hierarchical clustering on time series, per condition

full - perform both differential expression  and clustering analysis


Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----
Input are *.bam files in the format:

    sample-time-replicate.bam

and a gtf/gff file of transcripts/genes/features over which to
count an analyse.

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------+-----------+---------------------------+
|*Program*     |*Version*  |*Purpose*                  |
+--------------+-----------+---------------------------+
|featureCounts |           |count mapped reads over    |
|              |           |target gene models         |
+--------------+-----------+---------------------------+
|R             | 2.15/3.0  |PCA and data transformation|
+--------------+-----------+---------------------------+
|featureCounts |           |read counting over genes   |
+--------------+-----------+---------------------------+


Pipeline output
===============

The major output are a set of file with gene:cluster assignments and
clustering metrics

Requirements:

* featurecounts >= 1.4.6
* deseq
* deseq2
* masigpro
* wgcna
* rcolorbrewer
* ggplot2
* gplots
* reshape2
* venndiagram
* dtw
* flashclust

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import glob
import gzip
import os
import itertools
import re
import math
import types
import collections
import time
import optparse
import shutil
import sqlite3
import random
import itertools
import tempfile
import numpy as np
import pandas as pd
import pandas.rpy.common as com
from pandas.io import sql
import rpy2.rinterface as rinterface
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.GTF as GTF
import PipelineTimeseries

###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# sample = PipelineTracks.AutoSample

# assign tracks
GENESETS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("*.gtf.gz"),
    "(\S+).gtf.gz")
TRACKS3 = PipelineTracks.Tracks(PipelineTracks.Sample3)
TRACKS = TRACKS3.loadFromDirectory(glob.glob("*.bam"), "(\S+).bam")

REPLICATE = PipelineTracks.Aggregate(TRACKS, labels=("replicate", ))
TIME = PipelineTracks.Aggregate(TRACKS, labels=("condition", "tissue"))

###################################################################
###################################################################
###################################################################


def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh
###################################################################
###################################################################
###################################################################


@follows(connect, mkdir("plots.dir"))
@transform("reference.gtf.gz",
           suffix("reference.gtf.gz"),
           "refcoding.gtf.gz")
def buildCodingGeneSet(infile, outfile):
    '''build a gene set with only protein coding transcripts.

    Genes are selected via their gene biotype in the GTF file.
    Note that this set will contain all transcripts of protein
    coding genes, including processed transcripts.

    This set includes UTR and CDS.
    '''

    statement = '''
    zcat %(infile)s | awk '$2 == "protein_coding"' | gzip > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################


@follows(mkdir("feature_counts.dir"))
@files([(("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile()),
         ("feature_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(), y.asFile())))
        for x, y in itertools.product(TRACKS, GENESETS)])
def buildFeatureCounts(infiles, outfile):
    '''counts reads falling into "features", which by default are genes.

    A read overlaps if at least one bp overlaps.

    Pairs and strandedness can be used to resolve reads falling into
    more than one feature. Reads that cannot be resolved to a single
    feature are ignored.

    '''

    infile, annotations = infiles

    # featureCounts cannot handle gzipped in or out files
    outfile = P.snip(outfile, ".gz")
    annotations_tmp = P.getTempFilename()

    # -p -B specifies count fragments rather than reads, and both
    # reads must map to the feature
    if PARAMS['featurecounts_paired'] == "1":
        paired = "-p -B"
    else:
        paired = ""

    job_options = "-pe dedicated %i" % PARAMS['featurecounts_threads']

    statement = '''
    zcat %(annotations)s > %(annotations_tmp)s;
    checkpoint;
    featureCounts %(featurecounts_options)s
    -T %(featurecounts_threads)s
    -s %(featurecounts_strand)s
    -b
    -a %(annotations_tmp)s
    -o %(outfile)s
    %(infile)s
    > %(outfile)s.log;
    checkpoint;
    gzip %(outfile)s;
    checkpoint;
    rm %(annotations_tmp)s '''

    P.run()

###################################################################
###################################################################
###################################################################


@collate(buildFeatureCounts,
         regex("feature_counts.dir/(.+)-(.+)-(.+)_vs_(.+).tsv.gz"),
         r"feature_counts.dir/\1-\4-feature_counts.tsv.gz")
def aggregateFeatureCounts(infiles, outfile):
    ''' build a matrix of counts with genes and tracks dimensions.
    '''

    # Use column 7 as counts. This is a possible source of bugs, the
    # column position has changed before.

    infiles = " ".join(infiles)
    statement = '''
    python %(scriptsdir)s/combine_tables.py
    --columns=1
    --take=7
    --use-file-prefix
    --regex-filename='(.+)_vs.+.tsv.gz'
    --log=%(outfile)s.log
    %(infiles)s
    | sed 's/Geneid/gene_id/'
    | sed 's/\-/\./g'
    | tee %(outfile)s.table.tsv
    | gzip > %(outfile)s '''

    P.run()
###################################################################
###################################################################
###################################################################


@transform(aggregateFeatureCounts,
           suffix(".tsv.gz"),
           ".load")
def loadFeatureCounts(infile, outfile):
    P.load(infile, outfile, "--add-index=gene_id")
###################################################################
###################################################################
###################################################################


@follows(mkdir("combined_analysis.dir"), aggregateFeatureCounts)
@collate(aggregateFeatureCounts,
         regex("feature_counts.dir/(.+)-(.+)-feature_counts.tsv.gz"),
         r"combined_analysis.dir/\2-combined.tsv.gz")
def buildCombinedExpression(infiles, outfile):
    '''
    aggregate together all of the datasets for a combined
    all-vs-all analysis
    '''
    infiles = " ".join(infiles)
    statement = '''
    python %(scriptsdir)s/combine_tables.py
    --columns=1
    --log=%(outfile)s.log
    %(infiles)s
    | sed 's/Geneid/gene_id/'
    | sed 's/\-/\./g'
    | tee %(outfile)s.table.tsv
    | gzip > %(outfile)s
    '''

    P.run()
###################################################################
###################################################################
###################################################################


@transform(buildCombinedExpression,
           suffix("combined.tsv.gz"),
           "combined.load")
def loadCombinedExpression(infile, outfile):
    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################


@follows(mkdir("deseq.dir"), loadFeatureCounts)
@transform(aggregateFeatureCounts,
           regex(r"feature_counts.dir/(.+)-(.+)-feature_counts.tsv.gz"),
           r"deseq.dir/\1-\2-vst.tsv")
def DESeqNormalize(infile, outfile):
    ''' Use DESeq normalization and variance
    stabilizing transformation on all data.
    Use `blind` dispersion method and `fit-only`
    sharingMode.
    '''

    time_agg = TIME.__dict__['track2groups'].keys()
    time_points = [int(str(x).split("-")[1]) for x in time_agg]
    time_points = set(time_points)
    time_points = list(time_points)
    time_points.sort()
    time_points = [str(x) for x in time_points]
    rep_agg = REPLICATE.__dict__['track2groups'].keys()
    replicates = [str(x).split("-")[2] for x in rep_agg]

    time_points = ",".join(time_points)
    replicates = ",".join(replicates)

    statement = '''
    python %(scriptsdir)s/expression2expression.py
    --task=deseq
    --log=%(outfile)s.log
    --replicates=%(replicates)s
    --time=%(time_points)s
    %(infile)s
    > %(outfile)s '''
    P.run()

###################################################################
###################################################################
###################################################################


@transform(DESeqNormalize,
           suffix("-vst.tsv"),
           "-vst.load")
def loadDESeqNormalize(infile, outfile):
    P.load(infile, outfile, transpose=True)
##################################################################
###################################################################
###################################################################

if len([PARAMS['refs']]) > 1:

    @follows(buildCombinedExpression)
    @collate(buildCombinedExpression,
             regex(r"combined_analysis.dir/(.+)-combined.tsv.gz"),
             r"combined_analysis.dir/merged-combined.tsv.gz")
    def mergeExpressionTables(infile, outfile):
        '''
        Merge refcoding and lncRNA count tables
        '''
        file1 = infile[0]
        file2 = infile[1]

        tmpfile = P.getTempFilename(shared=True)

        df1 = pd.read_table(file1,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression="gzip")

        df2 = pd.read_table(file2,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression="gzip")

        out_frame = df1.append(df2)

        out_frame.to_csv(tmpfile, sep="\t")

        statement = '''cat %(tmpfile)s | gzip > %(outfile)s; rm -rf %(tmpfile)s'''

        P.run()

    @follows(aggregateFeatureCounts)
    @collate(aggregateFeatureCounts,
             regex(r"feature_counts.dir/(.+)-(.+)-feature_counts.tsv.gz"),
             r"combined_analysis.dir/\1-merged.tsv.gz")
    def mergeSingleExpressionTables(infile, outfile):
        '''
        Merge refcoding and lncRNA count tables from a single condition
        if there are separate input reference gtfs.
        '''

        file1 = infile[0]
        file2 = infile[1]

        tmpfile = P.getTempFilename(shared=True)

        df1 = pd.read_table(file1,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression="gzip")

        df2 = pd.read_table(file2,
                            sep="\t",
                            index_col=0,
                            header=0,
                            compression="gzip")

        out_frame = df1.append(df2)

        out_frame.to_csv(tmpfile, sep="\t")

        statement = '''cat %(tmpfile)s | gzip > %(outfile)s; rm -rf %(tmpfile)s'''

        P.run()

else:
    @follows(buildCombinedExpression)
    @collate(buildCombinedExpression,
             regex(r"combined_analysis.dir/(.+)-combined.tsv.gz"),
             r"combined_analysis.dir/merged-combined.tsv.gz")
    def mergeExpressionTables(infile, outfile):
        '''
        Only a single reference gtf, copy combined expression
        table
        '''

        infile = infile[0]

        statement = '''zcat %(infile)s | gzip > %(outfile)s'''

        P.run()

    @follows(aggregateFeatureCounts)
    @transform(aggregateFeatureCounts,
               regex("feature_counts.dir/(.+)-(.+)-feature_counts.tsv.gz"),
               r"combined_analysis.dir/\1-merged.tsv.gz")
    def mergeSingleExpressionTables(infile, outfile):
        '''
        Only a single reference gtf, copy expression table, no need to merge
        '''

        statement = ''' zcat %(infile)s | gzip > %(outfile)s'''

        P.run()
##################################################################
###################################################################
###################################################################


@follows(loadCombinedExpression,
         mergeExpressionTables,
         mkdir("diff_condition.dir"))
@transform([buildCombinedExpression, mergeExpressionTables],
           suffix("combined.tsv.gz"),
           "condition-diff.tsv.gz")
def conditionDiffExpression(infile, outfile):
    '''
    Call DEGs showing statistically significantly
    different expression based on interaction terms between condition
    and time point.  Uses DESeq2.
    '''

    statement = '''
    zcat %(infile)s |
    python %(scriptsdir)s/timeseries2diffgenes.py
    --log=%(outfile)s.log
    --method=condition
    --alpha=%(deseq_alpha)s
    --results-directory=diff_condition.dir
    '''

    P.run()

    P.touch(outfile)
##################################################################
###################################################################
###################################################################


@follows(conditionDiffExpression)
@transform("diff_condition.dir/*.tsv",
           regex(r"diff_condition.dir/(.+).tsv"),
           r"diff_condition.dir/\1.load")
def loadConditionDiffExpression(infile, outfile):
    P.load(infile, outfile)
##################################################################
###################################################################
###################################################################


@follows(mergeSingleExpressionTables,
         mkdir("diff_timepoints.dir"))
@transform("combined_analysis.dir/*-merged.tsv.gz",
           suffix("merged.tsv.gz"),
           "diff-time.tsv.gz")
def timePointDiffExpression(infile, outfile):
    '''
    Within each condition test for differentially expressed
    genes against the baseline time point.  Uses DESeq2.
    '''

    statement = '''
    python %(scriptsdir)s/timeseries2diffgenes.py
    --log=%(outfile)s.log
    --method=timepoint
    --alpha=%(deseq_alpha)s
    --results-directory=diff_timepoints.dir
    %(infile)s
    '''

    P.run()

    P.touch(outfile)
##################################################################
###################################################################
###################################################################


@follows(timePointDiffExpression)
@transform("diff_timepoints.dir/*.tsv",
           regex(r"diff_timepoints.dir/(.+).tsv"),
           r"diff_timepoints.dir/\1.load")
def loadTimePointDiffExpression(infile, outfile):
    P.load(infile, outfile)
##################################################################
###################################################################
###################################################################


@follows(loadConditionDiffExpression)
@collate(r"diff_condition.dir/*.tsv",
         regex(r"diff_condition.dir/(.+).(.+)_(.+).tsv"),
         r"images.dir/\1-venn.png")
def drawConditionVennDiagram(infiles, outfile):
    '''
    Generates a Venn Diagram for the overlap of differentially
    expressed genes and lncRNAs for up to the first 5 time points
    for the time:condition interaction analysis.
    '''

    # select up to 5 time points to plot
    select = []
    time_points = PARAMS['venn_timepoints'].split(",")
    if len(infiles) >= len(time_points):
        for te in time_points:
            fle = [x for x in infiles if re.search(r"0_%s" % te, x)]
            if fle:
                select.append(fle[0])
            else:
                pass
    else:
        select = infiles

    select = ",".join(select)

    statement = '''
    python %(scriptsdir)s/diffgene2venn.py
    --alpha=%(deseq_alpha)s
    --log=condition-venn.log
    --file-list=%(select)s
    --output-directory=images.dir
    '''

    P.run()
##################################################################
###################################################################
###################################################################


@follows(loadConditionDiffExpression)
@collate(r"diff_timepoints.dir/*.tsv",
         regex(r"diff_timepoints.dir/(.+)_(.+)_(.+)-time.tsv"),
         r"images.dir/\1-time-venn.png")
def drawTimeVennDiagram(infiles, outfile):
    '''
    Generates a Venn Diagram for the overlap of differentially
    expressed genes and lncRNAs for up to the first 5 time points
    from the time point differential analysis.
    '''

    # select up to 5 time points to plot
    select = []
    time_points = PARAMS['venn_timepoints'].split(",")
    if len(infiles) >= len(time_points):
        for te in time_points:
            fle = [x for x in infiles if re.search(r"%s-time" % te, x)]
            if fle:
                select.append(fle[0])
            else:
                pass
    else:
        select = infiles
    select = ",".join(select)

    statement = '''
    python %(scriptsdir)sdiffgene2venn.py
    --alpha=%(deseq_alpha)s
    --log=condition-venn.log
    --file-list=%(select)s
    --output-directory=images.dir
    '''

    P.run()
##################################################################
###################################################################
###################################################################


@follows(loadDESeqNormalize)
@transform(DESeqNormalize,
           regex("deseq.dir/(.+)-(.+)-vst.tsv"),
           r"deseq.dir/\1-\2-filtered-vst.tsv")
def sumCovarFilter(infile, outfile):
    '''
    Filter gene list based on the distribution of the
    sums of the covariance of each gene.  This is highly
    recommended to reduce the total number of genes used
    in the dynamic time warping clustering to reduce the
    computational time.  The threshold is placed at the
    intersection of the expected and observed value
    for the given quantile.
    '''

    time_agg = TIME.__dict__['track2groups'].keys()
    time_points = [int(str(x).split("-")[1]) for x in time_agg]
    time_points = set(time_points)
    time_points = list(time_points)
    time_points.sort()
    time_points = [str(x) for x in time_points]
    rep_agg = REPLICATE.__dict__['track2groups'].keys()
    replicates = [str(x).split("-")[2] for x in rep_agg]

    time_points = ",".join(time_points)
    replicates = ",".join(replicates)

    statement = '''
    python %(scriptsdir)s/expression2expression.py
    --log=%(outfile)s.log
    --task=sumcovar
    --time=%(time_points)s
    --replicates=%(replicates)s
    --quantile=%(filtering_quantile)s
    %(infile)s
    > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@transform(sumCovarFilter,
           suffix("-filtered-vst.tsv"),
           "-filtered-vst.load")
def loadFilteredData(infile, outfile):
    P.load(infile, outfile)

###################################################################
###################################################################
###################################################################
ANALYSIS = PARAMS['resampling_analysis_type']
if ANALYSIS == 'replicates':
    @follows(sumCovarFilter,
             mkdir("clustering.dir"))
    @transform(sumCovarFilter,
               regex("deseq.dir/(.+)-(.+)-filtered-vst.tsv"),
               r"clustering.dir/\1-\2-replicates.tsv")
    def genReplicateData(infile, outfile):
        '''
        Split each replicate into a separate file for clustering
        within each replicate.  Relies on each replicate being the
        same across the whole time series.
        '''

        outdir = outfile.split("/")[0]
        PipelineTimeseries.splitReplicates(infile=infile,
                                           axis="column",
                                           group_var="replicates",
                                           outdir=outdir)

        P.touch(outfile)
        ###################################################################
        ###################################################################
        ###################################################################

    if PARAMS['resampling_parallel']:
        @follows(genReplicateData,
                 mkdir("parallel_files.dir"))
        @subdivide("clustering.dir/*-expression.tsv",
                   regex("clustering.dir/(.+)-(.+)-(.+)-expression.tsv"),
                   r"parallel_files.dir/\1-\2-\3-split.sentinel")
        def splitFiles(infile, outfile):
            '''
            Arbitrarily split files into chunks for parallelisation
            '''

            PipelineTimeseries.splitFiles(infile=infile,
                                          nchunks=PARAMS['resampling_chunks'],
                                          out_dir="parallel_files.dir")
            P.touch(outfile)
        ###################################################################
        ###################################################################
        ###################################################################

        @follows(splitFiles)
        @transform("parallel_files.dir/*-split.tsv",
                   regex("parallel_files.dir/(.+)-(.+)-(.+)-(.+)-split.tsv"),
                   add_inputs(r"clustering.dir/\1-\2-\3-expression.tsv"),
                   r"parallel_files.dir/\1-\2-\3-\4-distance.tsv")
        def splitDistance(infiles, outfile):
            '''
            Calculate distances on split files
            '''

            if PARAMS['clustering_lag']:
                clustering_options = " --lag=%s " % PARAMS['clustering_lag']
            else:
                clustering_options = " "

            if PARAMS['clustering_k']:
                clustering_options = " --k=%s " % PARAMS['clustering_k']
            else:
                clustering_options = " "

            infile = infiles[0]
            expression_file = infiles[1]
            statement = '''
            python %(scriptsdir)s/expression2distance.py
            --log=%(outfile)s.log
            --parallel
            --distance-metric=%(clustering_metric)s
            --expression-file=%(expression_file)s
            %(clustering_options)s
            --out=%(outfile)s
            %(infile)s
            '''

            P.run()
        ###################################################################
        ###################################################################
        ###################################################################

        @follows(splitDistance)
        @collate(splitDistance,
                 regex("parallel_files.dir/(.+)-(.+)-(.+)-(.+)-distance.tsv"),
                 r"clustering.dir/\1-\2-\3-distance.tsv")
        def distanceCalculation(infiles, outfile):
            '''
            Merge split files in the correct order, defined by their file
            names.
            '''

            infiles = ",".join(infiles)

            statement = '''
            python /ifs/devel/projects/proj036/pipeline_timeseries/src/scripts/distance2merge.py
            --log=%(outfile)s.log
            --outfile=%(outfile)s
            %(infiles)s
            '''
            P.run()

        ###################################################################
        ###################################################################
        ###################################################################

    else:
        @follows(genReplicateData)
        @transform("clustering.dir/*-expression.tsv",
                   regex("clustering.dir/(.+)-(.+)-(.+)-expression.tsv"),
                   r"clustering.dir/\1-\2-\3-distance.tsv")
        def distanceCalculation(infile, outfile):
            '''
            Calls the dtw script for each resampled file.
            Calculates the dynamic time-warping distance for each
            pairwise gene combination.
            '''

            if PARAMS['clustering_lag']:
                clustering_options = " --lag=%s " % PARAMS['clustering_lag']
            else:
                clustering_options = " "

            if PARAMS['clustering_k']:
                clustering_options = " --k=%s " % PARAMS['clustering_k']
            else:
                clustering_options = " "

            statement = '''
            python %(scriptsdir)/expression2distance.py
            --distance-metric=%(clustering_metric)s
            --log=%(outfile)s.log
            %(clustering_options)s
            --out=%(outfile)s
            %(infile)s
            '''

            P.run()
    ###################################################################
    ###################################################################
    ###################################################################

    @follows(mkdir("tmp.dir/"),
             distanceCalculation)
    @transform(distanceCalculation,
               suffix("-distance.tsv"),
               "-clusters.tsv")
    def clusterCut(infile, outfile):
        '''
        Use dynamic tree cutting to derive clusters for each
        resampled distance matrix
        '''

        condition = (str(infile).split("-")[0]).lstrip("clustering.dir/")
        gtf_source = str(infile).split("-")[1]
        resample_prefix = str(infile).rstrip("-distance.tsv")
        expression_file = "deseq.dir/%s-%s-expression.tsv" % (condition,
                                                              gtf_source)
        wgcna_out = "clustering.dir/WGCNA.out"
        cluster_file = "%s-clusterlabels.tsv" % (resample_prefix)

        job_options = "-l mem_free=6G"

        if PARAMS['clustering_deepsplit']:
            options = " --split-clusters "
        else:
            options = " "

        statement = '''
        python %(scriptsdir)s/distance2clusters.py
        --log=%(outfile)s.log
        --task=cluster
        --cluster-algorithm=%(clustering_algorithm)s
        --cluster-file=%(cluster_file)s
        --expression-file=%(expression_file)s
        %(options)s
        %(infile)s
        > %(outfile)s'''

        P.run()

    ###################################################################
    ###################################################################
    ##################################################################

    @follows(clusterCut,
             mkdir("consensus_clustering.dir"))
    @collate(distanceCalculation,
             regex("clustering.dir/(.+)-(.+)-(.+)-distance.tsv"),
             r"consensus_clustering.dir/\1-\2-cocluster.tsv")
    def clusterAgree(infiles, outfile):
        '''
        Calculate average distance matrix over replicates
        '''

        infiles = ",".join(infiles)

        statement = '''
        python %(scriptsdir)s/distance2clusters.py
        --task=clustagree
        --method=replicate
        --log=%(outfile)s.log
        %(infiles)s
        > %(outfile)s
        '''

        P.run()
    ###################################################################
    ###################################################################
    ###################################################################

    @follows(mkdir("plots.dir"))
    @transform(clusterAgree,
               suffix("-cocluster.tsv"),
               "-consensus.tsv")
    def consensusClustering(infile, outfile):
        '''
        hierachical clustering based on clustering correlation
        across all resampled data sets
        '''

        if PARAMS['clustering_deepsplit']:
            options = " --split-clusters "
        else:
            options = " "

        statement = '''
        python %(scriptsdir)s/distance2clusters.py
        --log=%(outfile)s.log
        --task=consensus-cluster
        --cut-height=%(clustering_cut)s
        --cluster-algorithm=%(clustering_consensus_algorithm)s
        --cluster-size=%(clustering_min_size)s
        %(options)s
        %(infile)s
        > %(outfile)s '''

        P.run()
###################################################################
###################################################################
###################################################################


elif PARAMS["resampling_analysis_type"] == 'resample':
    @follows(mkdir("clustering.dir"),
             loadFilteredData)
    @transform(sumCovarFilter,
               regex("deseq.dir/(.+)-(.+)-filtered-vst.tsv"),
               r"clustering.dir/\1-\2-resampled.tsv")
    def genResampleData(infile, outfile):
        '''
        Resample the data n-times with replacement - generates
        n flat files which are then propagated at later stages.
        Files are generally small though.
        '''

        time_agg = TIME.__dict__['track2groups'].keys()
        time_points = [int(str(x).split("-")[1]) for x in time_agg]
        time_points.sort()
        time_points = list(set(time_points))
        rep_agg = REPLICATE.__dict__['track2groups'].keys()
        replicates = [str(x).split("-")[2] for x in rep_agg]
        time_rep_comb = [x for x in itertools.product(time_points, replicates)]
        time_cond = ro.StrVector([x[0] for x in time_rep_comb])
        rep_cond = ro.StrVector([x[1] for x in time_rep_comb])
        ref_gtf = str(infile).split("-")[1]
        condition = (str(infile).split("-")[0]).strip("deseq.dir/")

        time_points = ",".join([str(i) for i in time_points])
        replicates = ",".join(replicates)

        statement = '''
        python %(scriptsdir)s/data2resamples.py
        --log=%(outfile)s.log
        --time=%(time_points)s
        --replicates=%(replicates)s
        --condition=%(condition)s
        --resamples=%(resampling_resample)s
        --input-gtf=%(ref_gtf)s
        --output-file-directory=clustering.dir
        --seed=%(resampling_seed)s
        %(infile)s
        '''
        P.run()

        P.touch(outfile)
    ###################################################################
    ###################################################################
    ###################################################################
    # for randomly resampled data files
    # see below for single-replicate data files
    # using parallelisation will split jobs into no more than n x 500 genes
    # per job.  This may mean > number specified in resampling_chunks
    # parameter
    if PARAMS['resampling_parallel']:
        @follows(genResampleData,
                 mkdir("parallel_files.dir"))
        @subdivide("clustering.dir/*-expression.tsv",
                   regex("clustering.dir/(.+)-(.+)-(.+)-expression.tsv"),
                   r"parallel_files.dir/\1-\2-\3-split.sentinel")
        def splitFiles(infile, outfile):
            '''
            Arbitrarily split files into chunks for parallelisation
            '''

            PipelineTimeseries.splitFiles(infile=infile,
                                          nchunks=PARAMS['resampling_chunks'],
                                          out_dir="parallel_files.dir")
            P.touch(outfile)
        ###################################################################
        ###################################################################
        ###################################################################

        @follows(splitFiles)
        @transform("parallel_files.dir/*-split.tsv",
                   regex("parallel_files.dir/(.+)-(.+)-(.+)-(.+)-split.tsv"),
                   add_inputs(r"clustering.dir/\1-\2-\3-expression.tsv"),
                   r"parallel_files.dir/\1-\2-\3-\4-distance.tsv")
        def splitDistance(infiles, outfile):
            '''
            Calculate distances on split files
            '''

            if PARAMS['clustering_lag']:
                clustering_options = " --lag=%s " % PARAMS['clustering_lag']
            else:
                clustering_options = " "

            if PARAMS['clustering_k']:
                clustering_options = " --k=%s " % PARAMS['clustering_k']
            else:
                clustering_options = " "

            infile = infiles[0]
            expression_file = infiles[1]
            statement = '''
            python %(scriptsdir)s/expression2distance.py
            --log=%(outfile)s.log
            --parallel
            --distance-metric=%(clustering_metric)s
            --expression-file=%(expression_file)s
            %(clustering_options)s
            --out=%(outfile)s
            %(infile)s
            '''

            P.run()
        ###################################################################
        ###################################################################
        ###################################################################

        @follows(splitDistance)
        @collate(splitDistance,
                 regex("parallel_files.dir/(.+)-(.+)-(.+)-(.+)-distance.tsv"),
                 r"clustering.dir/\1-\2-\3-distance.tsv")
        def mergeSplitFiles(infiles, outfile):
            '''
            Merge split files in the correct order, defined by their file
            names.
            '''

            infiles = ",".join(infiles)

            statement = '''
            python %(scriptsdir)s/distance2merge.py
            --log=%(outfile)s.log
            --outfile=%(outfile)s
            %(infiles)s
            '''
            P.run()

        ###################################################################
        ###################################################################
        ###################################################################

        @follows(mkdir("tmp.dir/"))
        @transform(mergeSplitFiles,
                   suffix("-distance.tsv"),
                   "-clusters.tsv")
        def clusterCut(infile, outfile):
            '''
            Use dynamic tree cutting to derive clusters for each
            resampled distance matrix
            '''

            condition = (str(infile).split("-")[0]).lstrip("clustering.dir/")
            gtf_source = str(infile).split("-")[1]
            resample_prefix = str(infile).rstrip("-distance.tsv")
            expression_file = "deseq.dir/%s-%s-diffgenes.tsv" % (condition,
                                                                 gtf_source)
            wgcna_out = "clustering.dir/WGCNA.out"
            cluster_file = "%s-clusterlabels.tsv" % (resample_prefix)

            job_options = "-l mem_free=6G"
            statement = '''
            python %(scriptsdir)s/distance2clusters.py
            --log=%(outfile)s.log
            --task=cluster
            --cluster-algorithm=%(clustering_algorithm)s
            --cluster-file=%(cluster_file)s
            --expression-file=%(expression_file)s
            %(infile)s
            > %(outfile)s'''

            P.run()

    else:
        @follows(genResampleData)
        @transform("clustering.dir/*-expression.tsv",
                   regex("clustering.dir/(.+)-(.+)-(.+)-expression.tsv"),
                   r"clustering.dir/\1-\2-\3-distance.tsv")
        def distanceCalculation(infile, outfile):
            '''
            Calls the dtw script for each resampled file.
            Calculates the dynamic time-warping distance for each
            pairwise gene combination.
            '''

            if PARAMS['clustering_lag']:
                clustering_options = " --lag=%s " % PARAMS['clustering_lag']
            else:
                clustering_options = " "

            if PARAMS['clustering_k']:
                clustering_options = " --k=%s " % PARAMS['clustering_k']
            else:
                clustering_options = " "

            statement = '''
            python %(scriptsdir)s/expression2distance.py
            --distance-metric=%(clustering_metric)s
            %(clustering_options)s
            --log=%(outfile)s.log
            --out=%(outfile)s
            %(infile)s
            '''

            P.run()
        ###################################################################
        ###################################################################
        ###################################################################

        @follows(mkdir("tmp.dir/"))
        @transform(distanceCalculation,
                   suffix("-distance.tsv"),
                   "-clusters.tsv")
        def clusterCut(infile, outfile):
            '''
            Use dynamic tree cutting to derive clusters for each
            resampled distance matrix
            '''

            condition = (str(infile).split("-")[0]).lstrip("clustering.dir/")
            gtf_source = str(infile).split("-")[1]
            resample_prefix = str(infile).rstrip("-distance.tsv")
            expression_file = "deseq.dir/%s-%s-diffgenes.tsv" % (condition,
                                                                 gtf_source)
            wgcna_out = "clustering.dir/WGCNA.out"
            cluster_file = "%s-clusterlabels.tsv" % (resample_prefix)

            job_options = "-l mem_free=6G"

            if PARAMS['clustering_deepsplit']:
                options = " --split-clusters "
            else:
                options = " "

            statement = '''
            python %(scriptsdir)s/distance2clusters.py
            --log=%(outfile)s.log
            --task=cluster
            --cluster-algorithm=%(clustering_algorithm)s
            --cluster-file=%(cluster_file)s
            --expression-file=%(expression_file)s
            %(options)s
            %(infile)s
            > %(outfile)s'''

            P.run()
        ###################################################################
        ###################################################################
        ###################################################################

    @transform(clusterCut,
               suffix("-clusters.tsv"),
               "-clean_clusters.tsv")
    def sanitizeClusters(infile, outfile):
        '''remove extraneous quotations from files
        and add in tab-delimiter'''

        statement = ''' sed '/cluster_matched.cluster/d' %(infile)s |
        awk '{printf("%%s\\t%%s\\n", $1,$2)}'
        > %(outfile)s'''

        P.run()
    ###################################################################
    ###################################################################
    ###################################################################

    @collate(sanitizeClusters,
             regex("clustering.dir/(.+)-(.+)-(.+)-clean_clusters.tsv"),
             r"clustering.dir/\1-\2-match_cluster.tsv")
    def aggregateClusters(infiles, outfile):
        '''
        Pull together all clusterings for each condition and reference
        gtf file
        '''

        infiles = " ".join(infiles)
        statement = '''
        python %(scriptsdir)s/combine_tables.py
        --columns=1
        --log=%(outfile)s.log
        %(infiles)s
        > %(outfile)s'''

        P.run()

    ###################################################################
    ###################################################################
    ###################################################################

    @follows(aggregateClusters)
    @transform("clustering.dir/*-expression.tsv",
               regex("clustering.dir/(.+)-(.+)-(.+)-expression.tsv"),
               r"clustering.dir/\1-\2-\3-clean_expression.tsv")
    def sanitizeExpression(infile, outfile):
        '''
        strip timepoints from header of expression files to match
        against clustering
        '''

        statement = '''awk '{if(NR>1){print $0}}' %(infile)s > %(outfile)s'''

        P.run()

    ###################################################################
    ###################################################################
    ###################################################################

    @follows(sanitizeExpression)
    @transform(sanitizeClusters,
               regex("clustering.dir/(.+)-(.+)-(.+)-(.+).tsv"),
               add_inputs(r"clustering.dir/\1-\2-\3-clean_expression.tsv"),
               r"clustering.dir/\1-\2-\3-matched.tsv")
    def matchClusterExpression(infiles, outfile):
        '''
        merge cluster assignemnt and expression data
        '''

        all_files = " ".join(infiles)

        statement = '''python %(scriptsdir)s/combine_tables.py
        --columns=1
        --log=%(outfile)s.log
        < %(all_files)s
        | sed '/times/d'
        | sed '/replicates/d'
        > %(outfile)s'''

        P.run()
    ###################################################################
    ###################################################################
    ###################################################################

    @follows(matchClusterExpression)
    @transform(matchClusterExpression,
               suffix("-matched.tsv"),
               "-matched.load")
    def loadMatchClusterExpression(infile, outfile):
        P.load(infile, outfile)

    ###################################################################
    ###################################################################
    ###################################################################

    @follows(mkdir("consensus_cluster.dir/"))
    @transform(aggregateClusters,
               regex("clustering.dir/(.+)-(.+)-match_cluster.tsv"),
               r"consensus_cluster.dir/\1-\2-cocluster.tsv")
    def clusterAgree(infile, outfile):
        '''
        Generate consensus clustering across all resampled data sets
        '''

        statement = '''
        python %(scriptsdir)s/distance2clusters.py
        --log=%(outfile)s.log
        --task=clustagree
        --method=resample
        %(infile)s
        > %(outfile)s'''

        P.run()

    ###################################################################
    ###################################################################
    ###################################################################

    @follows(mkdir("plots.dir"))
    @transform(clusterAgree,
               suffix("-cocluster.tsv"),
               "-consensus.tsv")
    def consensusClustering(infile, outfile):
        '''
        hierachical clustering based on clustering correlation
        across all resampled data sets
        '''

        if PARAMS['clustering_deepsplit']:
            options = " --split-clusters "
        else:
            options = " "

        statement = '''
        python %(scriptsdir)s/distance2clusters.py
        --log=%(outfile)s.log
        --task=consensus-cluster
        --cut-height=%(clustering_cut)s
        --cluster-algorithm=%(clustering_consensus_algorithm)s
        --cluster-size=%(clustering_min_size)s
        %(options)s
        %(infile)s
        > %(outfile)s '''

        P.run()
###################################################################
###################################################################
###################################################################


@follows(consensusClustering,
         mkdir("images.dir"),
         mkdir("eigengenes.dir"))
@transform(consensusClustering,
           regex("consensus_clustering.dir/(.+)-(.+)-consensus.tsv"),
           add_inputs(r"deseq.dir/\1-\2-filtered-vst.tsv"),
           r"eigengenes.dir/\1-PCA-\2-cluster_eigengenes.tsv")
def clusterEigengenePCA(infile, outfile):
    '''
    Using the consensus clustering gene modules
    to derive eigengene expression profiles
    and loadings - uses R prcomp
    This function is complementary to `moduleEigenes`
    '''

    infile = ",".join(infile)

    image_dir = os.path.join(os.getcwd(),
                             "images.dir")

    statement = '''
    python %(scriptsdir)s/distance2clusters.py
    --log=%(outfile)s.log
    --task=pca
    --image-dir=%(image_dir)s
    %(infile)s
    > %(outfile)s
    '''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(clusterEigengenePCA)
@transform(clusterEigengenePCA,
           suffix(".tsv"),
           ".load")
def loadClusterEigengenes(infile, outfile):
    '''
    Load module eignengene expression profiles in to DB
    '''

    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################


@transform(consensusClustering,
           suffix("-consensus.tsv"),
           "-consensus.load")
def loadConsensusClustering(infile, outfile):
    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################


@follows(loadConsensusClustering)
@transform("clustering.dir/*-match_cluster.tsv",
           regex("clustering.dir/(.+)-(.+)-match_cluster.tsv"),
           r"consensus_clustering.dir/\1-\2-consensus.metrics.tsv")
def consensusMetrics(infile, outfile):
    '''
    Calculate clustering consensus metrics
    across resampled data sets
    '''

    job_options = "-l mem_free=14G"

    statement = '''
    python %(scriptsdir)s/clusters2metrics.py
    --log=%(outfile)s.log
    %(infile)s
    > %(outfile)s'''

    P.run()
###################################################################
###################################################################
###################################################################


@follows(consensusMetrics)
@transform(consensusMetrics,
           suffix(".tsv"),
           ".load")
def loadConsensusMetrics(infile, outfile):
    P.load(infile, outfile)
###################################################################
###################################################################
###################################################################
###################################################################
# Targets
###################################################################


@follows(clusterAgree,
         loadClusterEigengenes,
         loadConsensusClustering,
         loadConsensusMetrics)
def clustering():
    pass


@follows(mergeExpressionTables,
         loadTimePointDiffExpression,
         loadConditionDiffExpression,
         drawTimeVennDiagram,
         drawConditionVennDiagram)
def diff_expression():
    pass


@follows(clustering,
         diff_expression)
def full():
    pass
###################################################################
# primary targets
###################################################################


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit(P.main(sys.argv))
