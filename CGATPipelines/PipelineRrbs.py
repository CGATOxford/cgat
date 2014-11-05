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

import os
import shutil
import glob
import collections
import re
import gzip
import itertools
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
    scale_colour_discrete(guide = FALSE) +
    xlab("Read position") + ylab("Methylation (percentage)")
    ggsave('%(plot1_out)s', plot = p, width = 5, height = 5)

    p2 = ggplot(df, aes(x=position,
    y=count_methylated+count_unmethylated, col=context)) +
    geom_line() + geom_point(size=1) +
    facet_grid(context ~ ., scales="free") +
    scale_x_continuous(minor_breaks=seq(0,read_length),
    breaks=seq(0,read_length,5)) +
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
