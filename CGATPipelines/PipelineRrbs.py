'''
PipelineRrbs.py - Utility functions for analysing rrbs sequence data
==============================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra` files)
of different data (single end, paired end) with different mapping
algorithms (bowtie, tophat, stampy). This module provides utility
functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra` archive
might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

The module currently is able to deal with:

   * tophat mapping against genome
   * bowtie mapping against transcriptome, genome and junctions
   * bwa against genome
   * stampy against genome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''
import random
import os
import shutil
import glob
import collections
import re
import gzip
import itertools
import copy
import CGAT.Pipeline as P
import logging as L
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Fastq as Fastq
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Fasta as Fa
import pysam
import pandas as pd
from CGAT.Pipeline import cluster_runnable
import numpy as np
import pandas.rpy.common as com
from rpy2.robjects import r
import rpy2.robjects as robj
from rpy2.robjects.packages import importr


@cluster_runnable
def plotReadBias(infile, outfile):
    '''take the m-bias outfile from bismark and plot read position biases'''

    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    def dfFromTsvTable(lines, line, context):
        row = 2
        rows = []
        row_values = lines[line + row].strip().split("\t")
        header = [re.sub(" ", "_", re.sub("%", "perc", x))
                  for x in row_values]
        row += 1
        while lines[line + row].strip():
            print lines[line+row].strip()
            row_values = lines[line + row].strip().split("\t")
            row_values = map(num, row_values)
            rows.append(row_values)
            row += 1
        df = pd.DataFrame(rows, columns=header)
        df['context'] = context
        return(df)

    with open(infile, "r") as f:
        lines = f.readlines()
        line = 0
        while line < len(lines):
            if lines[line].strip() == "CpG context":
                df = dfFromTsvTable(lines, line, "CpG")
                print df
            elif lines[line].strip() == "CHG context":
                df2 = dfFromTsvTable(lines, line, "CHG")
                print df2
            elif lines[line].strip() == "CHH context":
                df3 = dfFromTsvTable(lines, line, "CHH")
                print df3
            line += 1
    final_df = pd.concat([df, df2, df3])
    final_df.index = range(len(final_df['context']))
    r_dataframe = com.convert_to_r_dataframe(final_df)

    plot1_out = (P.snip(outfile, ".read_position.tsv") +
                 "_read_position_methylation_bias.png")
    plot2_out = re.sub("_methylation_bias", "_count", plot1_out)

    title = re.sub(".fastq.*", "", os.path.basename(infile))
    plotMbiasReadPositionSummary = r("""
    library(ggplot2)
    library(reshape)
    require(scales)
    function(df){
    read_length = max(df$position)
    p = ggplot(df, aes(x=position, y=perc_methylation, col=context)) +
    geom_line() + geom_point(size=1) +
    facet_grid(context ~ ., scales="free") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5)) +
    ggtitle("%(title)s") +
    scale_colour_discrete(guide = FALSE) +
    xlab("Read position") + ylab("Methylation (percentage)")
    ggsave('%(plot1_out)s', plot = p, width = 5, height = 5)

    p2 = ggplot(df, aes(x=position,
    y=count_methylated+count_unmethylated, col=context)) +
    geom_line() + geom_point(size=1) +
    facet_grid(context ~ ., scales="free") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5)) +
    ggtitle("%(title)s") +
    scale_colour_discrete(guide = FALSE) +
    scale_y_continuous(labels = comma) +
    xlab("Read position") + ylab("Number of cytosines")
    ggsave('%(plot2_out)s', plot = p2, width = 5, height = 5)

    write.table(df, file="%(outfile)s", row.names = F, sep="\t")
    }""" % locals())

    plotMbiasReadPositionSummary(r_dataframe)


@cluster_runnable
def pandasMerge(infile1, infile2, outfile, merge_type, left, right,
                delim="\t", com="#"):
    '''function to merge two tsv files using pandas
    left and right are the columns to merge on'''

    def pandasRead(infile):
        return pd.io.parsers.read_csv(infile, sep=delim, comment=com)

    df1 = pandasRead(infile1)
    df2 = pandasRead(infile2)
    merged = df1.merge(df2, how=merge_type, left_on=left,
                       right_on=right, sort=False)
    merged.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def fasta2CpG(infile, outfile):
    '''perform an insilico digest at MspI sites and return all CpGs
    thorughout the genome, whether they are in a MspI fragment
    and if so, what their read position is'''
    # to do: paramterise (digestion site, read length, PE/SE)

    FA = Fa.Fasta(open(infile, "r"))
    temp_contig = FA.FetchOne()

    outfile = open(outfile, "w")
    outfile.write("contig\tposition\tstrand\tread_position\n")

    while len(temp_contig[1]) > 1:
        contig, contig_seq = temp_contig
        # find MspI sites
        iter_pos = re.finditer("[cC][cC][gG][gG]", contig_seq)
        pos_start = [x.start(0) for x in iter_pos]
        # loop through pairs of MspI sites and add positions to dictionary
        MspI_cpg_dict = {}
        for n in range(0, len(pos_start)-1):
            frag_length = pos_start[n+1] - pos_start[n]
            if frag_length > 25 and frag_length < 300:
                if frag_length >= 51:
                    # need to include additional base at the end of the read
                    # to check for final CG
                    read_f = contig_seq[pos_start[n]+1:pos_start[n]+53]
                    read_r = contig_seq[pos_start[n+1]-49:pos_start[n+1]+3]
                else:
                    read_f = contig_seq[pos_start[n]+1:pos_start[n+1]+1]
                    read_r = contig_seq[pos_start[n]+3:pos_start[n+1]+3]
                # find positions of CGs in in silico reads
                f_cgs = re.finditer("[cC][gG]", read_f)
                r_cgs = re.finditer("[cC][gG]", read_r)
                for x in f_cgs:
                    # 1-based read position
                    read_pos = x.start(0) + 1
                    # contig position 1-based
                    contig_position = pos_start[n] + 1 + read_pos
                    MspI_cpg_dict[contig_position] = read_pos
                for x in r_cgs:
                    # 1-based read position
                    read_pos = len(read_r)-(x.start(0)+1)
                    # contig position 1-based
                    contig_position = pos_start[n+1] + 4 - read_pos
                    MspI_cpg_dict[contig_position] = read_pos

        # find location of ALL cpgs
        all_cpgs = re.finditer("[cC][gG]", temp_contig[1])
        for cpg in all_cpgs:
            # 2 Cs for each CpG (one on each strand)
            c_pos, g_pos = (cpg.start(0) + 1, cpg.start(0) + 2)
            read_pos = "NA"
            # check if position of C is in the dictionary of MspI CpGs
            if c_pos in MspI_cpg_dict:
                read_pos = MspI_cpg_dict[c_pos]
            outfile.write("%s\t%s\t%s\t%s\n" % (contig, c_pos,
                                                "+", read_pos))
            read_pos = "NA"
            if g_pos in MspI_cpg_dict:
                read_pos = MspI_cpg_dict[g_pos]
            outfile.write("%s\t%s\t%s\t%s\n" % (contig, g_pos,
                                                "-", read_pos))
        temp_contig = FA.FetchOne()
    outfile.close()


@cluster_runnable
def mergeAndDrop(cpgs, infiles, outfile):
    # prepare dataframe with cpg locations
    # do I need to index on location column?
    all_cpgs = pd.io.parsers.read_csv(cpgs, sep='\t', comment='#')
    all_cpgs['location'] = (all_cpgs['contig'] + ":" +
                            all_cpgs['position'].apply(str))
    # all_cpgs.drop(['contig', 'position'], axis=1, inplace=True)

    # for loop through files and merge each time with previous merged df
    for infile in infiles:
        sample_name = re.sub("\.fastq.*bismark", "",
                             P.snip(os.path.basename(infile), ".cov"))

        # prepare temporary dataframe with meth values
        column_names = ["contig", "start", "stop", "meth_perc",
                        "count_meth", "count_unmeth"]
        temp_meth_df = pd.io.parsers.read_csv(infile, sep='\t',
                                              comment='#', header=None,
                                              names=column_names)
        temp_meth_df['location'] = (temp_meth_df['contig'] + ":" +
                                    temp_meth_df['start'].apply(str))
        temp_meth_df.drop(['contig', 'start', 'stop', 'meth_perc'],
                          axis=1, inplace=True)
        temp_meth_df.columns = [sample_name + "_meth", sample_name + "_unmeth",
                                'location']

        # merge into temp dataframe
        merged = temp_meth_df.merge(all_cpgs, how="outer",
                                    left_on='location', right_on='location',
                                    sort=False)
        all_cpgs = merged.copy(deep=True)

    merged.drop('location', axis=1, inplace=True)
    merged.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def subsetToCovered(infile, outfile, cov_threshold=10):
    '''take a flat file containing methylation calls for all cpgs
    and output a flat file containing only sites above coverage threshold'''

    raw = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    # only need the number of sample rows
    samples = len(identifyExperimentalFactors(raw)[0])

    def minCov(row, number_sample_rows):
        column_range = range(0, number_sample_rows)
        # take every other column (e.g meth count), add next column
        # (e.g unmeth count) and find min for row
        return min([np.nansum([row[i],  row[i+1]])
                    for i in column_range[::2]])

    raw['min_cov'] = raw.apply(minCov, axis=1, args=(samples,))
    high_cov = raw[raw['min_cov'] >= cov_threshold]
    high_cov['cpgi'].fillna('Non-CpGIsland', inplace=True)
    # note: check setting NA to above works
    high_cov.to_csv(outfile, sep="\t", index=False, na_rep="NA")


def identifyExperimentalFactors(df):
    '''takes a data frame and returns the samples, treatment groups,
    replicate numbers and tissues based on tissue-treatment-replicate
    nomenclature'''

    columns = []
    for col in list(df.columns.values):
        m = re.match("\S*-\S*-\S*", col)
        if m is not None:
            columns.append(col)

    samples = [re.sub("_\S*", "", x) for x in columns]
    treatments = list(set([x.split("-")[1] for x in samples]))
    tissues = list(set([x.split("-")[0] for x in samples]))
    replicates = list(set([x.split("-")[2] for x in samples]))

    return samples, treatments, tissues, replicates


@cluster_runnable
def addTreatmentMean(infile, outfile):
    '''take a flat file containing methylation calls and calculate
    treatment means based on column names'''

    def getMethylationFrequency(row):
        ''' this function assumes meth and unmeth columns are next to
        one another - generalise'''
        return row[0]/(row[0]+row[1])

    df = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    samples, treatments, tissues, replicates = identifyExperimentalFactors(df)

    for tissue in tissues:
        for treatment in treatments:
            treatment_columns = []
            for replicate in replicates:
                m_col_name = "-".join([tissue, treatment,
                                       replicate]) + "_meth"
                um_col_name = "-".join([tissue, treatment,
                                        replicate]) + "_unmeth"
                temp_df = df[[m_col_name, um_col_name]]
                meth_freq_col_name = tissue+"_"+treatment+"_"+replicate
                df[meth_freq_col_name] = temp_df.apply(getMethylationFrequency,
                                                       axis=1)
                treatment_columns.append(meth_freq_col_name)
            temp_df = df[treatment_columns]

            df[tissue+"_"+treatment+"_mean"] = temp_df.apply(np.mean, axis=1)
            df[tissue+"_"+treatment+"_var"] = temp_df.apply(np.var, axis=1)

    df.to_csv(outfile, sep="\t", index=False, na_rep="NA")


@cluster_runnable
def summaryPlots(infile, outfile):
    '''take a flat file containing average methylation values and
    produce various summary plots

    Plots:
    1. Read position vs. methylation level
    2. Distributions of average methylation, split by group and CpGI/Non-CpGI
    3. Methylation distributions by sample, split by CpGI/Non-CpGI
    4. Heatmap of Spearman's rho split by CpGI/Non-CpGI
    5. Correlation between mean value pairs for each treatment pair,
       split by CpGI/Non-CpGI
    6. Read position vs. methylation difference for each treatment pair

    to do:
    1. Add all vs. all correlations plots
    test!
    '''

    df = pd.io.parsers.read_csv(infile, sep='\t',  comment='#')
    samples, treatments, tissues, replicates = identifyExperimentalFactors(df)
    r_dataframe = com.convert_to_r_dataframe(df)

    base = P.snip(outfile, ".log")

    gr = importr('grDevices')

    out = open(outfile, "w")

    plotFuncReadPositionMethylation = r("""
    library(ggplot2)
    library(reshape)
    function(df){
    read_length = max(df$read_position[is.finite(df$read_position)])
    df2 = df[, grepl("mean", names(df))]
    df_agg = aggregate(df2, by=list(df$read_position),FUN=mean, na.rm=TRUE)
    df_agg_m = melt(df_agg,id="Group.1")
    p = ggplot(df_agg_m, aes(x=Group.1,y = value)) +
    geom_line(aes(col=variable)) +
    geom_point(aes(col=variable), size=1.5) +
    ggtitle("Read Position vs. Methylation") +
    scale_colour_discrete(name = "Sample") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5))+
    xlab("Read position") +
    ylab("Methylation by read position")

    ggsave(paste0('%(base)s','-read_position_meth.png'),
    plot = p, width = 5, height = 5)
    }""" % locals())

    plotFuncReadPositionMethylation(r_dataframe)
    out.write("plot saved at: %(base)s-read_position_meth.png\n" % locals())

    plotMethBySampleMean = r("""
    library(ggplot2)
    library(reshape)
    function(df){
    df2 = df[df$cpgi != "CpGIsland",]
    df3 = df2[, grepl("mean", names(df))]
    df4 = melt(df3)

    p=ggplot(df4,
    aes(y=as.numeric(value), x=variable)) +
    geom_violin(aes(col=variable)) +
    scale_colour_discrete(name = "Sample", guide = FALSE) +
    ylab("Fraction Methylated") + xlab("") +
    scale_x_discrete(labels = gsub("_", "-", levels(df4$variable))) +
    theme(legend.position="bottom",
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5)) +
    ggtitle("Methylation by treatment group: Non-CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_type_non_cpgi.png'),
    plot = p, width = 5, height = 5)

    df5 = df[df$cpgi == "CpGIsland",]
    df6 = df5[, grepl("mean", names(df))]
    df7 = melt(df6)

    p2=ggplot(df7,
    aes(y=as.numeric(value), x=variable)) +
    geom_violin(aes(col=variable)) +
    scale_colour_discrete(name = "Sample", guide = FALSE) +
    ylab("Fraction Methylated") + xlab("") +
    theme(legend.position="bottom",
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5)) +
    ggtitle("Methylation by treatment group: CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_type_cpgi.png'),
    plot = p2, width = 5, height = 5)
    }""" % locals())

    plotMethBySampleMean(r_dataframe)
    out.write("plot saved at: \
    %(base)s-methylation_by_sample_type_non_cpgi.png\n" % locals())
    out.write("plot saved at: \
    %(base)s-methylation_by_sample_type_cpgi.png\n" % locals())

    meth_freq_columns = ["_".join((x, y, z))
                         for x in tissues
                         for y in treatments
                         for z in replicates]
    meth_freq_columns.append("cpgi")

    df_meth_freq = df[meth_freq_columns]
    r_dataframe_meth_freq = com.convert_to_r_dataframe(df_meth_freq)

    # re-write with subfunctions in r function to reduce length of code
    # each plot is repeated twice!
    drop_cpgi = '!(colnames(df) %in% c("cpgi"))'
    plotMethBySampleAndHeatmap = r("""
    library(ggplot2)
    library(reshape)
    library(amap)
    function(df){
    df2 = df[, %(drop_cpgi)s]
    var = apply(df2,1,var)
    var[is.na(var)] <- 0
    mean = apply(df2,1,mean)
    df3 = df[(var > 0 & mean >= 0.1 & mean <= 0.9),]
    df4 = df3[df3$cpgi != "CpGIsland", %(drop_cpgi)s]
    df4_m = melt(df4)
    samples=data.frame(matrix(unlist(strsplit(as.character(df4_m$variable),"_")),
    ncol=3,byrow=TRUE))
    colnames(samples) = c("tissue","treatment","replicate")
    df4_m = cbind(samples, df4_m)

    p=ggplot(df4_m,
    aes(x=as.numeric(value), group=variable)) +
    geom_density(aes(col=paste0(tissue,"-",treatment), linetype=replicate)) +
    scale_colour_discrete(name = "Treatment group") +
    scale_linetype_discrete(name = "Replicate") +
    xlab("Fraction Methylation") +
    ggtitle("Methylation by sample: Non-CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_non_cpgi.png'),
    plot = p, width = 5, height = 5)

    df5 = as.matrix(sapply(df4, as.numeric))
    df6 = cor(df4, method="spearman")
    df6_m = melt(df6)
    p = ggplot(df6_m, aes(x=X1, y=X2, fill=value))+
    geom_tile() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_continuous(name = "rho") +
    ggtitle("non-CpG Islands")

    ggsave(paste0('%(base)s','_spearmans_heatmap_non_cpgi.png'),
    plot = p, width = 5, height = 5)

    df7 = df3[df3$cpgi == "CpGIsland", %(drop_cpgi)s]
    df7_m = melt(df7)
    samples=data.frame(matrix(unlist(strsplit(as.character(df7_m$variable),"_")),
    ncol=3,byrow=TRUE))
    colnames(samples) = c("tissue","treatment","replicate")
    df7_m = cbind(samples, df7_m)

    p = ggplot(df7_m,
    aes(x=as.numeric(value), group=variable)) +
    geom_density(aes(col=paste0(tissue,"-",treatment), linetype=replicate)) +
    scale_colour_discrete(name = "Treatment group") +
    scale_linetype_discrete(name = "Replicate") +
    xlab("Fraction Methylation") +
    ggtitle("Methylation by sample: CpG Islands")

    ggsave(paste0('%(base)s','_methylation_by_sample_cpgi.png'),
    plot = p, width = 5, height = 5)

    df8 = as.matrix(sapply(df7, as.numeric))
    df9 = cor(df8, method="spearman")
    df9_m = melt(df9)
    p2 = ggplot(df9_m, aes(x=X1, y=X2, fill=value))+
    theme(axis.text.x = element_text(angle = 90)) +
    geom_tile()  + xlab("") + ylab("") +
    scale_fill_continuous(name = "rho")  +
    ggtitle("CpG Islands")

    ggsave(paste0('%(base)s','_spearmans_heatmap_cpgi.png'),
    plot = p2, width = 5, height = 5)


    library(pvclust)
    # redefine distance function to use spearman's rank correlation
    dist.pvclust <- function(x){
    res <- as.dist(1 - cor(x, method="spearman", use="pairwise.complete.obs"))
    attr(res,"method") <- "correlation"
    return(res)
    x <- as.matrix(x)
    P <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)}

    result <- pvclust(as.matrix(df5), method.hclust="ward", nboot=100)
    png(paste0('%(base)s','_hierarchical_clustering_non_cpgi.png'))
    par(mar=c(0, 4, 4, 2))
    plot(result, xlab="", sub="", print.num=F, cex.pv=0.7,
    col.pv=c(0,3,0), main="Non-CpG islands")
    dev.off()

    result2 <- pvclust(as.matrix(df8), method.hclust="ward", nboot=100)
    png(paste0('%(base)s','_hierarchical_clustering_cpgi.png'))
    par(mar=c(0, 4, 4, 2))
    plot(result2, xlab="", sub="", print.num=F, cex.pv=0.7,
    col.pv=c(0,3,0), main="CpG islands")
    dev.off()

    }""" % locals())

    plotMethBySampleAndHeatmap(r_dataframe_meth_freq)
    out.write("plot saved at: \
    %(base)s-_spearmans_heatmap_non_cpgi.png\n" % locals())
    out.write("plot saved at: \
    %(base)s-_spearmans_heatmap_cpgi.png\n" % locals())

    for tissue in tissues:
        for treatment_pair in itertools.combinations(treatments, 2):
            treatment1, treatment2 = map(str, treatment_pair)
            tissue = str(tissue)
            out.write("plotting tissue: %s, treatment: %s vs. %s\n"
                      % (tissue, treatment1, treatment2))

            label_base = "-".join([tissue, treatment1, treatment2])

            title = "Global Correlation: \n\
            %(tissue)s: %(treatment1)s vs. %(treatment2)s" % locals()

            # re-write this function to include r plotting function and then
            # pass two dfs in turn split by cpgi so that only one python
            # function is required
            plotFuncGlobalCorrelation = r("""
            library(ggplot2)
            function(df){
            df_r1rm = df[df$read_position!=1 | !(is.finite(df$read_position)),]
            p = ggplot(df_r1rm, aes(x=%(tissue)s_%(treatment1)s_mean,
            y=%(tissue)s_%(treatment2)s_mean)) +
            geom_point(size=0.5, alpha = 0.2) +
            xlab("%(treatment1)s") + ylab("%(treatment2)s") +
            theme(axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5)) +
            ggtitle('%(title)s') + facet_grid(.~cpgi)
            ggsave(paste0('%(base)s','-global_correlation',
            '%(label_base)s','.png'), plot = p, width = 5, height = 5)
            }""" % locals())

            plotFuncGlobalCorrelation(r_dataframe)
            out.write("plot saved at: \
            %(base)s-global_correlation_cpgi_%(label_base)s.png\n" % locals())

            title = "Read Position vs. Difference in Methylation:\n\
            %(tissue)s: %(treatment1)s vs. %(treatment2)s" % locals()
            plotFuncReadPositionDifference = r("""
            library(ggplot2)
            function(df){

            f <- function(x) {
            r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
            names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
            r}
            quantiles = quantile(df$%(tissue)s_%(treatment1)s_mean-
            df$%(tissue)s_%(treatment2)s_mean, probs = c(0.01,0.99))

            df$read_position = factor(df$read_position)
            p = ggplot(df, aes(x=factor(read_position),
            y = (%(tissue)s_%(treatment1)s_mean-\
            %(tissue)s_%(treatment2)s_mean))) +
            stat_summary(fun.data=f, geom="boxplot") + ggtitle('%(title)s') +
            theme(axis.text.x = element_text(size=5, angle=90)) +
            xlab("Read position") + ylim(quantiles) +
            ylab("Difference in methylation:\
            %(treatment1)s vs. %(treatment2)s")
            ggsave(paste0('%(base)s','-read_position_meth_diff_',
            '%(label_base)s','.png'), plot = p, width = 5, height = 5)
            }""" % locals())
            # note:need to check whether final plot function works:pep8 changes
            plotFuncReadPositionDifference(r_dataframe)
            out.write("plot saved at: \
            %(base)s-read_position_meth_diff_%(label_base)s.png\n" % locals())

    out.close()


@cluster_runnable
def spikeInClusters(infiles, outfile):
    '''explain purpose'''
    class CpG():
        def __init__(self):
            self.contig = ""
            self.pos = ""
            self.meth = []
            self.unmeth = []
            self.perc = []

        def __str__(self):
            return "\t".join(map(str, [self.contig, self.pos,
                                       self.meth, self.unmeth]))

        def row_values(self):
            row = [self.contig]
            row.append(self.pos)
            row.extend([x for x in self.perc])
            row.extend([x for x in self.meth])
            row.extend([x for x in self.unmeth])
            return row

    class cluster():
        def __init__(self, cpgs, contig):
            self.cpgs = cpgs
            self.contig = contig
            self.df = ""

        def __str__(self):
            return "\t".join(map(str, self.cpgs))

    class cluster_region():
        ''' class function to hold information on randomly selected
        regions from clusters prior to region swapp'''
        def __init__(self, cluster, avr_perc, s, e):
            self.cluster = cluster
            self.avr_perc = avr_perc
            self.s, self.e = (s, e)

        def __str__(self):
            return "\t".join(map(str, (self.cluster, self.avr_perc,
                                       self.s, self.e)))

    class spike_in_cluster():
        def __init__(self, region_initial, region_change, region_size, df):
            self.intial = region_initial
            self.change = region_change
            self.size = region_size
            self.df = df

        def __str__(self):
            return "\t".join(map(str, [self.change, self.size]))

    outfile = open(outfile, "w")
    outfile.write("reached checkpoint 1\n")

    samples = []
    CpGs = {}
    for infile in infiles:
        sample = re.sub(".fastq.gz_bismark.*", "", os.path.basename(infile))
        outfile.write(sample)
        sample_values = sample.split("-")
        sample_values[0] = sample_values[0][0]
        sample_values[1] = sample_values[1][0]
        samples.append("".join(sample_values))
        with open(infile, 'r') as f:
            for line in f.readlines():
                contig, start, end, perc, meth, unmeth = line.strip().split()
                if meth + unmeth >= 10:
                    start = str(start)
                    if start not in CpGs:
                        CpGs[start] = CpG()
                        CpGs[start].contig = contig
                        CpGs[start].pos = start
                    CpGs[start].perc.append(float(perc))
                    CpGs[start].meth.append(int(meth))
                    CpGs[start].unmeth.append(int(unmeth))
    # header = ["position", "unmeth", "meth"]
    meth_header = []
    meth_header.extend(samples)

    outfile.write("reached checkpoint 2\n")

    rows = {}
    thres = len(meth_header)
    for x in CpGs:
        if len(CpGs[x].meth) == thres:
            rows[x] = CpGs[x].row_values()

    df = pd.DataFrame(data=rows)
    df = df.transpose()
    columns = ["contig", "position"]
    columns.extend([x+"_perc" for x in meth_header])
    columns.extend([x+"_meth" for x in meth_header])
    columns.extend([x+"_unmeth" for x in meth_header])
    df.columns = columns
    df_con = df.convert_objects(convert_numeric=True)
    df_sort = df_con.sort(columns=["contig", "position"])

    d = 100  # max distance
    n = 10  # min number of CpGs

    positions = df_sort['position']
    contigs = df_sort['contig']

    clusters = {}
    c_current = positions[0]
    current_cpgs = []
    current_contig = contigs[0]
    for cpg in range(1, len(df_sort.index)):
        next_cpg = positions[cpg]
        next_contig = contigs[cpg]
        if ((next_cpg < c_current + d) &
            (next_contig == current_contig)):  # catch change of contig
            current_cpgs.append(next_cpg)
        else:
            if len(current_cpgs) >= n:
                start, end = map(str, [current_cpgs[0], current_cpgs[-1]])
                label = start + ":" + end
                clusters[label] = cluster(current_cpgs, current_contig)
                clusters[label].df = df_sort.loc[start:end]
            current_cpgs = []
        c_current = next_cpg
        current_contig = next_contig

    outfile.write("reached checkpoint 3\n")

    colnames = df_sort.columns.values
    swap_regex = "\SD\d+.*"
    keep_perc_regex = "\SS\d+_perc"
    swap_perc_regex = "\SD\d+_perc"
    swap_cols = [col for col in colnames if re.search(swap_regex, col)]
    keep_cols = [col for col in colnames if not re.search(swap_regex, col)]
    keep_perc_cols = [col for col in colnames if
                      re.search(keep_perc_regex, col)]
    swap_perc_cols = [col for col in colnames if
                      re.search(swap_perc_regex, col)]
    outfile.write("swap_cols:%s\n" % ",".join(swap_cols))
    outfile.write("keep_cols:%s\n" % ",".join(keep_cols))
    outfile.write("swap_perc_cols:%s\n" % ",".join(swap_perc_cols))
    outfile.write("keep_perc_cols:%s\n" % ",".join(keep_perc_cols))
    # need to bin clusters by average methylation so that a
    # range of intitial methylation levels can be sampled
    initial_cluster_regions = {}
    initial_step = 10  # bin size for initial methylation levels
    max_bin = 100
    min_in_bin = 100  # minimal number of clusters for bin to be included
    initial_meth_bins = range(0, max_bin+initial_step, initial_step)
    for meth_perc in initial_meth_bins:
        initial_cluster_regions[meth_perc] = []
    for x in clusters:
        old_cluster = clusters[x].df.copy()
        avr_perc = np.mean(np.mean(old_cluster.ix[:, keep_perc_cols]))
        bin_m = ((np.digitize([avr_perc], initial_meth_bins))[0]*initial_step)
        initial_cluster_regions[bin_m].append(x)
    for x in initial_meth_bins:
        outfile.write("for bin %s, there are %s regions" %
               (x, len(initial_cluster_regions[x])))
    initial_meth_bins_covered = [x for x in initial_meth_bins if
                                 len(initial_cluster_regions[x]) > min_in_bin]

    outfile.write("reached checkpoint 4\n")

    min_bin = 10  # minimum mean methylation in region to be shuffled
    max_bin = 120  # maximum mean methylation in region to be shuffled
    step = 10
    cluster_regions_dict = {}
    for size in range(1, 11, 2):  # number of cpgs in region for shuffling
        cluster_regions = {}
        meth_bins = range(min_bin, max_bin, step)
        for meth_perc in meth_bins:
            cluster_regions[meth_perc] = []
        for x in clusters:
            old_cluster = clusters[x].df.copy()
            '''only use clusters with average methylation of at least 10%
            otherwise most of the compute is wasted looking for regions of
            >10% methylation in clusters with very little methylation!'''
            if np.mean(np.mean(old_cluster.ix[:, keep_perc_cols])) > 10:
                cluster_size = len(old_cluster.index)
                success = 0
                for y in range(0, 200):
                    # try 200 times to find 20 regions with at least 10%
                    # methylation, then give up and move onto the next cluster
                    if success < 20:
                        s = random.randint(0, cluster_size-size)
                        e = s + size
                        avr_perc = np.mean(np.mean(
                            old_cluster.ix[s:e, keep_perc_cols]))
                        if avr_perc > min_bin:
                            temp_region = cluster_region(x, avr_perc, s, e)
                            bin_m = ((np.digitize([avr_perc], meth_bins))[0] *
                                     step)+min_bin
                            cluster_regions[bin_m].append(temp_region)
                            success += 1
        cluster_regions_dict[size] = cluster_regions
    outfile.write("reached checkpoint 5\n")

    outfile.write("size = %s\n" % size)
    for x in meth_bins:
        outfile.write("for bin %s, there are %s regions\n" %
               (x, len(cluster_regions_dict[size][x])))


    sample_col_dict = {}
    sample_dataframe_dict = {}
    for x in samples:
        sample_columns = ["contig", "position", "position"]
        regex = x
        sample_columns.extend([col for col in columns if
                               re.search(regex, col)])
        sample_col_dict[x] = sample_columns
        sample_dataframe_dict[x] = pd.DataFrame(columns=sample_columns)

    for size in range(1, 10, 2):
        meth_bins_covered = [x for x in meth_bins if
                             len(cluster_regions_dict[size][x]) > 49]
        for repeat in range(0, 11):
            rand_initial_bin = random.randint(0, len(
                initial_meth_bins_covered)-1)
            initial_meth_bin = initial_meth_bins_covered[rand_initial_bin]
            rand_x = random.randint(0, len(initial_cluster_regions
                                           [initial_meth_bin])-1)
            temp_df = clusters[initial_cluster_regions[initial_meth_bin]
                               [rand_x]].df.copy()
            rand_bin = random.randint(0, len(meth_bins_covered)-1)
            meth_bin = meth_bins_covered[rand_bin]
            rand_y = random.randint(0, len(cluster_regions_dict
                                           [size][meth_bin])-1)
            y = cluster_regions_dict[size][meth_bin][rand_y]
            temp_y, start, end = clusters[y.cluster], y.s, y.e
            temp_df_swap = temp_y.df.ix[start:end, swap_cols]

            rand_s = random.randint(0, len(temp_df.index)-size)
            rand_e = rand_s + size
            temp_df_swap.set_index(temp_df[rand_s:rand_e].index,
                                   drop=True,  inplace=True)
            initial_meth = np.mean(np.mean(temp_df.ix[rand_s:rand_e,
                                                      keep_perc_cols]))
            change_meth = np.mean(np.mean(temp_df_swap.ix
                                          [:, swap_perc_cols])) - initial_meth
            temp_df.ix[rand_s:rand_e, swap_cols] = temp_df_swap

            change_step = 5
            change_min_bin = -100
            change_max_bin = 100
            change_bins = range(change_min_bin, change_max_bin, change_step)

            bins = [(x*change_step) + change_min_bin for x in
                    np.digitize([initial_meth, change_meth], change_bins)]
            initial, change = map(str, bins)

            temp_df['contig'] = ["_".join([x, str(size), str(rand_s),
                                           str(rand_e), initial,
                                           re.sub("-", "m", change)])
                                 for x in temp_df['contig']]

            for x in samples:
                temp_df['position'] = temp_df['position'].astype(int)
                sample_temp_df = temp_df.ix[:, sample_col_dict[x]]
                sample_dataframe_dict[x] = sample_dataframe_dict[x].append(
                    sample_temp_df)
    outfile.write("reached checkpoint 6")

    # generate null clusters
    perc_regex = "\S\S\d+_perc"
    perc_cols = [x for x in range(0, len(columns))
                 if re.search(perc_regex, columns[x])]
    unmeth_regex = "\S\S\d+_unmeth"
    unmeth_cols = [x for x in range(0, len(columns))
                   if re.search(unmeth_regex, columns[x])]
    meth_regex = "\S\S\d+_meth"
    meth_cols = [x for x in range(0, len(columns))
                 if re.search(meth_regex, columns[x])]

    for repeat in range(0, 11):  # number of null clusters
        rand_initial_bin = random.randint(0, len(initial_meth_bins_covered)-1)
        initial_meth_bin = initial_meth_bins_covered[rand_initial_bin]
        rand_x = random.randint(0, len(initial_cluster_regions[
            initial_meth_bin])-1)
        temp_df = clusters[initial_cluster_regions[initial_meth_bin]
                           [rand_x]].df.copy()
        initial_meth = np.mean(np.mean(temp_df.ix[:, keep_perc_cols]))
        initial = str((np.digitize([initial_meth], change_bins)
                       [0]*change_step)+change_min_bin)

        randomise_samples = range(0, len(samples))
        random.shuffle(randomise_samples)

        columns_randomised = ["contig", "position"]
        columns_randomised.extend([columns[perc_cols[y]]
                                   for y in randomise_samples])
        columns_randomised.extend([columns[meth_cols[y]]
                                   for y in randomise_samples])
        columns_randomised.extend([columns[unmeth_cols[y]]
                                   for y in randomise_samples])
        temp_df.columns = columns_randomised
        temp_df['contig'] = ["_".join([x, "Null", "NA", "NA", initial, "NA"])
                             for x in temp_df['contig']]
        for x in samples:
            temp_df['position'] = temp_df['position'].astype(int)
            sample_temp_df = temp_df.ix[:, sample_col_dict[x]]
            sample_dataframe_dict[x] = sample_dataframe_dict[x].append(
                sample_temp_df)
    outfile.write("reached checkpoint 7\n")

    out_prefix = "/ifs/projects/proj034/spike_in/"
    for x in sample_dataframe_dict:
        outfile_cov = str(out_prefix + x + "_10_pipeline_spike_in.cov")
        temp_df = sample_dataframe_dict[x].copy()
        temp_df['position'] = temp_df['position'].astype(int)
        temp_df['%s_meth' % x] = temp_df['%s_meth' % x].astype(int)
        temp_df['%s_unmeth' % x] = temp_df['%s_unmeth' % x].astype(int)
        temp_df.to_csv(outfile_cov, index=False, header=False, sep="\t",
                       dtype={'position': int})
        outfile.write(outfile_cov)
    outfile.write("reached checkpoint 8\n")
    outfile.close()

@cluster_runnable
def spikeInClustersAnalysis(infile, outfile):

    infiles = open(infile, "r").readlines()
    samples = [re.sub(".fastq.gz_bismark*", "", os.path.basename(x))
               for x in infiles]

    samples_abv = [re.sub("-", "", x) for x in samples]
    samples_abv = '","'.join(samples_abv)

    samples_treatment = [x.split("-")[1] for x in samples]
    samples_treatment = '","'.join(samples_treatment)

    infiles = '","'.join(infiles)

    M3D = r('''library(BiSeq)
            library(M3D)
            files <- c("%(infiles)s")
            samples <- data.frame(row.names=c("%(samples_abv)s"),
                       group=c("%(samples_treatment)s"))
            rawData <-readBismark(files,samples)
            clust.unlim <- clusterSites(object = rawData, perc.samples = 6/6,
                           min.sites = 10, max.dist = 100)
            clust.unlim=clusterSitesToGR(clust.unlim)
            overlaps<-findOverlaps(clust.unlim,rowData(rawData))
            MMDlist<-M3D_Wrapper(rawData, overlaps)
            M3Dstat<-MMDlist$Full-MMDlist$Coverage

            ranges_df <- data.frame(seqnames=seqnames(clust.unlim),
                         start=start(clust.unlim)-1, end=end(clust.unlim))
            M3Dstat_df <- as.data.frame(M3Dstat)
            adjusted_colnames = gsub("  ", "_", colnames(M3Dstat_df))
            colnames(M3Dstat_df) = adjusted_colnames
            complete_df = cbind(ranges_df, M3Dstat_df)
            write.table(complete_df,file="%(outfile)s",sep="\t")'''
            % locals())
    M3D()


@cluster_runnable
def spikeInClustersPlot(infile, outfile, groups):

    analysis_df = pd.read_csv(infile, sep="\t")
    cluster_characteristics = ["contig", "size", "s", "e", "initial", "change"]
    seqnames_split = analysis_df['seqnames'].apply(lambda x: x.split("_"))
    cluster_values = pd.DataFrame(data=seqnames_split.to_dict())
    cluster_values_df = cluster_values.transpose()
    cluster_values_df.columns = cluster_characteristics
    concat_df = pd.concat([analysis_df, cluster_values_df], axis=1)

    groups = [x[0] for x in groups]

    within = copy.copy(cluster_characteristics)
    M3D_within_col = []
    for x in groups:
        regex = "\S%s\d_vs_\S%s\d" % (x, x)
        M3D_within_col.extend([col for col in concat_df.columns
                               if re.search(regex, col)])
    within.extend(M3D_within_col)

    between = copy.copy(cluster_characteristics)
    pairs = [x for x in itertools.combinations(groups, 2)]
    M3D_between_col = []
    for x in pairs:
        regex = "\S%s\d_vs_\S%s\d" % (x[0], x[1])
        M3D_between_col.extend([col for col in concat_df.columns
                                if re.search(regex, col)])
    between.extend(M3D_between_col)

    between_df = concat_df.ix[:, between]
    within_df = concat_df.ix[:, within]

    r_df_between = com.convert_to_r_dataframe(between_df)
    r_df_within = com.convert_to_r_dataframe(within_df)
    base = P.snip(infile, ".out")
    PowerPlot = r('''library(ggplot2)
                  function(df_between, df_within){
                  m_between_df=melt(between_df,id.var=cluster_characteristics)
                  m_between_df['between'] = 1
                  m_between_df$change = gsub("m", "-", m_between_df$change)
                  m_within_df=melt(within_df,id.var=cluster_characteristics)
                  m_within_df['between'] = 0
                  m_cat_df = rbind(m_between_df,m_within_df)
                  null_df = m_within_df[m_within_df$size=="Null",]
                  non_null_within_df = m_within_df[m_within_df$size!="Null",]
                  print(head(null_df))
                  quant = (quantile(non_null_within_df$value,probs=0.95))
                  limited = non_null_within_df[non_null_within_df$value>quant,]

                  permute_func = function(x, len_x, vector, n){
                  rand = sample(vector, n, replace=TRUE)
                  return(length(rand[rand>x])/n)}

                  non_null_within_df$p_value = sapply(non_null_within_df$value,
                  function(x) permute_func(x, length(non_null_within_df$value),
                  non_null_within_df$value, 10000))

                  non_null_between_df= m_between_df[m_between_df$size!="Null",]
                  non_null_between_df$p_value=sapply(non_null_between_df$value,
                  function(x) permute_func(x, length(non_null_within_df$value),
                  non_null_within_df$value, 10000))
                  non_null_between_df$p_value_adj =p.adjust(
                  non_null_between_df$p_value , method ="BH")

                  agg_between = aggregate(non_null_between_df$p_value_adj,
                  by=list(non_null_between_df$size,between_non_null_df$change),
                  FUN=function(x) length(x[x<0.05])/length(x))
                  colnames(agg_between) = c("size","change","power")

                  text_theme = element_text(size=20)
                  plot_theme = theme(axis.title.x = text_theme,
                  axis.title.y = text_theme, axis.text.x = text_theme,
                  axis.text.y = text_theme, legend.text = text_theme)

                  p = ggplot(agg_between,aes(x=as.numeric(change),
                  y=size,fill=power))+geom_tile()+xlim(-10,45)+ plot_theme
                  scale_fill_continuous(limits=c(0,1))
                  ggsave(paste0('%(base)s',"_change_vs_size_power_heatmap.png"),
                  plot = p, width = 5, height = 5)

                  p = ggplot(agg_between,aes(x=as.numeric(change),y=power,
                  col=size)) + geom_line() + xlim(-10,40) + plot_theme
                  ggsave(paste0('%(base)s',"_change_vs_size_power_lineplot.png"),
                  plot = p, width = 5, height = 5

                  agg_between_count= aggregate(between_non_null_df$p_value_adj,
                  by=list(between_non_null_df$size,between_non_null_df$change),length)
                  colnames(agg_between_count) = c("size","change","count")

                  p = ggplot(agg_between,aes(y=as.numeric(change),x=size,
                  col=count)) + geom_line() + ylim(-10,40) + plot_theme
                  ggsave(paste0('%(base)s',"_change_vs_size_count.png"),
                  plot = p, width = 5, height = 5))''' % locals())

    PowerPlot(r_df_between, r_df_within)
