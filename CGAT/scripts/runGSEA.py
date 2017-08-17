'''
runGSEA.py
=============================================

:Tags: Python

Usage
-----
This script will perform the enrichment analysis, by using gene set enrichment analysis
(GSEA) and leading edge analysis.
            "Leading edge are defined as genes that are common to multiple
             significantly enriched gene sets  and  coordinately  enriched
             in a phenotypic comparison of interest.They represent a rich
             source of biologically important  genes."
-----
To run this analysis with GSEA,you need to provide
two input files:
      1. Ranked list of genes (Expression data set file).
      2. Gene set
Notes:
      - It is important to make sure that the expression data set
        does not include duplicate ids (edit your rank list,so that
        all row ids are unique). Otherwise it will affect results.
        Expression data set file is a tab-delimited text file.
        First line of dataset will be considered as header.

      - A gene sets file defines one or more gene sets. For each gene
        set,the file contains the gene set name and the list of genes in
        that gene set. A gene sets file is a tab-delimited text file in
        gmx or gmt format. For descriptions and examples of each file
        format please refer to:
        http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

      - The Molecular Signatures Database (MSigDB)(http://software.broadinstitute.org/gsea/msigdb/index.jsp)
        is a collection of annotated gene sets, which can be used for gene set enrichment analysis.
        OR you can create your own gene set in gmt or gmx format.

This script will summarize the analysis in the following format:
       1. GSEA Statistics
       2. GSEA Report
       3. Leading edge analysis report

example
-------

The way the test is ran:
cgat runGSEA -f "Expression_data.tsv" -g "Gene_set.gmt" -n 10 -d 1 -l 4

Default run conditions:
cgat runGSEA -f "Expression_data.tsv" -g "Gene_set.gmt"

--------------
GSEA Statistics
---------------
It includes following statistics for GSEA(for each phenotype):
          - Enrichment Score (ES)
          - Normalized Enrichment Score (NES)
          - False Discovery Rate (FDR)
          - Nominal P Value

--------------
GSEA reports
--------------
	- Global Statistics and Plots include:
	    a) Enrichment plot,
            b) Three separate bar plots that provide a quick overview of top 20 (this number is user defined)
               enriched upregulated, downregulated and overall enriched genesets on the basis of their FDR values.
            c) Global distribution of normalized enrichment score
            d) Global distribution of normalized enrichment score with corresponding FDR q values and p values.
        - Reports:
		1 - Enrichment in Phenotype (of up and downregulated genes)
           	  This report provides summary of enrichment analysis of each phenotype.
            	  It includes details of which genesets are up and downregulated and a summary
	    	  of significant enriched gensets on the basis of FDR and p values.)
         	2 - Gene Set Details
                  This report provides summary of preprocessing steps of the genesets provided and
	   	  lists genes sets that were used in the anlysis and which one were discarded due to set thresholds
         	3 - Detailed Enrichment Results
         	  This report provides detail statistics of each geneset(for each phenotype). Three reports are
         	  generated. report for uoregulated genesets, downregulated genesets, and enriched genesets organised
         	  on the basis of their FDR values.

	    By default, enrichment plot for top 20 gene sets will be reported.

----------------------------
Leading edge analysis report
----------------------------
It will report graphs that help you visualize the overlap between the selected leading edge subsets. It also
summarises the analysis in the form of reports. By default top 10 enriched genesets will be used for leading edge analysis.
	- Leading edge plots include:
             a) Heat Map(unclustered)
          	This provides an overview of overlap between leading edge subsets
             b) Heat Map(clustered)
          	This heat map will be generated after hierarchical clustering of leading edge subset. It will
          	show you clustered genes among subsets
       	     c) Set-to-Set Heat Map
         	This plot help you to visualize intensity of overlap between subsets (i.e. the extent of overlap between two genesets)
	     d) Dendogram to illustrate the arrangement of the clusters produced by hierarchical clustering.

        - Reports:
           1- Leading_edge_summary_report: summary of genesets and corresponding enrichment statistics that were used for the leading edge analysis.
           2- Leading edge matrix (gmx) file provides detailed information on leading edge analysis genesets
              (i.e. participating genes in each gene set).
           3- Leading edge (gct,cluster format) files for unclustered and clustered gene set. It is a boolean matrix.
              that can be used as an input into other resources for additional analysis as this is ideal format for cluster representation
              (in GSEA)

For details on the algorithm please refer to
Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550)
                    and
Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273).

=============================================
coding: utf-8
--------------------
Command line options
--------------------
'''
import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import collections
import numpy as np
import matplotlib as mpl
mpl.use('Agg')  # noqa: E402
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as mlines
import csv
import random
import statsmodels.sandbox.stats.multicomp as sm
from matplotlib.colors import ListedColormap
from decimal import *
import pandas as pd
import random
import string
import itertools
import os
import scipy
from matplotlib import rc, font_manager
from scipy.cluster.hierarchy import dendrogram, linkage
##################################################
# PLOT CUSTOMIZATION
plt.ioff()
legend_properties = {'weight': 'bold'}
plt.rc('legend', fontsize=20)
plt.style.use('seaborn-white')
title_font = {
    'fontname': 'Ubuntu',
    'size': '20',
    'color': 'darkblue',
    'weight': 'bold',
    'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
axis_font = {'fontname': 'Ubuntu', 'size': '16', 'weight': 'bold'}
axis_font_h = {'fontname': 'Ubuntu', 'size': '12', 'weight': 'bold'}
title_font_h = {
    'fontname': 'Ubuntu',
    'size': '20',
    'color': 'darkblue',
    'weight': 'bold'}
title_font_ran = {
    'fontname': 'Ubuntu',
    'size': '16',
    'color': 'darkblue',
    'weight': 'bold'}
axis_font_col = {'fontname': 'Ubuntu', 'size': '10', 'weight': 'bold'}
# Set the font properties (for use in legend)
font_path = 'C:\Windows\Fonts\Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=12)

#######################################################


def intersect(a, b):
    return list(set(a) & set(b))


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title(
            'Hierarchical Clustering Dendrogram (truncated)',
            **title_font)
        plt.xlabel(
            '\nGene set index or (Cluster size)',
            fontsize=18,
            weight='bold')
        plt.ylabel('Distance', fontsize=18, weight='bold', labelpad=28)
        for i, d, c in zip(
                ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center', **axis_font_h)
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def preprocess_geneset(file, m, n, id_exp):
    c = 0
    t = 0
    t2 = 0
    f = open(file, "r")
    lines = list(csv.reader(f, delimiter="\t"))
    INDI = np.zeros((len(lines),), dtype=int)
    excluded_list = []
    included_list = []
    for item in lines:
        Gene_id = item[2:len(item)]
        F_id = intersect(Gene_id, id_exp)
        if len(F_id) in range(m, n + 1):
            included_list.append(
                [item[0], item[1], len(Gene_id), len(F_id), F_id])
            t = t + 1
            INDI[c] = t
        else:
            excluded_list.append([item[0], item[1], len(Gene_id), len(F_id)])
            t2 = t2 + 1
            INDI[c] = -(t2)
        c = c + 1
    f.close()
    del lines
    return excluded_list, included_list, INDI


def read_expression(file):
    f = open(file, "r")
    result = []
    lines = list(csv.reader(f, delimiter="\t"))
    lines.pop(0)
    entrez_id = [item[0] for item in lines]
    value = [item[1] for item in lines]
    value_arr = np.array(value, dtype=np.float)
    f.close()
    return entrez_id, value_arr


def calculate_enrichment_score(s_I, ma, value, s, STORE):
    ma[s_I] = 1
    rr1 = value[s_I]
    rr2 = np.sum(np.power(np.absolute(rr1), 0))
    A33 = np.divide(-1, (len(ma) - s))
    enrich_score = 0
    for ii in range(0, len(ma)):
        if(ma[ii] == 1):
            CC = np.power(np.absolute(value[ii]), 0)
            CC = np.divide(CC, rr2)
            CC = CC + enrich_score
            STORE[ii] = CC
            enrich_score = CC
        else:
            CC = enrich_score + A33
            STORE[ii] = CC
            enrich_score = CC
    del rr1
    del rr2
    return STORE


def generate_gen_set_report(EX, IN, m, n, geneset_indicator):
    # Create gene set details
    C = len(EX) + len(IN)
    with open("CGAT_Gene_set_details.tsv", 'w') as f:
        f.write("Gene Set details:" + "\n")
        f.write("Gene set size filters (min=" +
                str(m) +
                ",max=" +
                str(n) +
                ") resulted in filtering out " +
                str(len(EX)) +
                "/" +
                str(C) +
                " gene sets." +
                "\n")
        f.write("The remaining " + str(C - len(EX)) +
                " gene sets were used in the analysis." + "\n")
        f.write(
            "List of gene sets used and their sizes (restricted to features in the specified dataset) are mentioned here:" +
            "\n")
        f.write(
            "\n" +
            "NAME" +
            "\t" +
            "ORIGINAL SIZE" +
            "\t" +
            "AFTER RESTRICTING TO EXPRESSION DATASET" +
            "\t" +
            "STATUS" +
            "\n")
        for i in range(0, len(geneset_indicator)):
            if(geneset_indicator[i] < 0):
                f.write(EX[abs(geneset_indicator[i]) -
                           1][0] +
                        "\t" +
                        str(EX[abs(geneset_indicator[i]) -
                               1][2]) +
                        "\t" +
                        "\t" +
                        " Rejected!" +
                        "\n")
            else:
                f.write(IN[geneset_indicator[i] -
                           1][0] +
                        "\t" +
                        str(IN[geneset_indicator[i] -
                               1][2]) +
                        "\t" +
                        str(IN[geneset_indicator[i] -
                               1][3]) +
                        '  ' +
                        "\n")
    f.close()
    return


def generate_leading_edge_m(
        I_U,
        I_D,
        f_U,
        f_D,
        IN_PRO,
        hm,
        dict_new,
        dfp,
        ufp,
        matr,
        OE,
        OEI,
        id_for_leading,
        store_gene_leading_matrix,
        gene_set_a,
        store_gene_leading_info):
    '''
    #####################################################################
    I am writing this section for leading edge analysis. It is necessary to describe each and every parameter
    here (just for myself).
    I_U,I_D and f_U,f_D are sorted upregulated and downregulated indexes (on the basis of NES values) and
    coreesponding FDR values.
    IN_PRO is list of included geneset "in_list"
    dict_new is dictionary format of expression dataset "ind_dict"
    gene_set_a is stored included gene set "GG"
    id_for_leading is id of expression dataset. I am using array format "id_new"
    Apart from that, I am using store_gene_leading_info as SGLI, original_es_index as OEI, original_es as OE
    I am also using up_for_plot and down_for_plot as ufp and dfp, which stores information for retrieving
    original index in original_nes.
    ch is choice of user FDR and hm is selected number of genes for leading edge analysis
    matr is temp matrix, used by me earlier.
    This section is little bit complex.
    #######################################################################
    '''
    inn_u = np.append(f_U, f_D)
    in1 = np.argsort(np.append(f_U, f_D))
    in2 = np.sort(np.append(f_U, f_D))
    FINAL_L = in1[0:hm]
    for i in FINAL_L:
        CW = inn_u[i]
        if(i >= len(f_U)):
            aw = dfp[I_D[i - len(f_U)]]
        else:
            aw = ufp[I_U[i]]
        II = gene_set_a[aw]
        inter = intersect(dict_new, II[0])
        indices = sorted([dict_new[x] for x in inter])
        matr[indices] = 1
        if(OE[aw] < 0):
            v1 = np.where(matr[OEI[aw]:len(matr)] == 1)
        else:
            v1 = np.where(matr[0:OEI[aw]] == 1)
        v2 = id_for_leading[v1]
        store_gene_leading_matrix.append([IN_PRO[aw][0], len(II[0]), len(v1[0]), list(
            v2), store_gene_leading_info[0][aw], store_gene_leading_info[1][aw], CW])
        matr.fill(0)
        del indices
        del inter
        del II

    return store_gene_leading_matrix


def heatmap_plot(hh1, nl):
    fig, ax = plt.subplots()
    cMap = ListedColormap(['white',
                           'palegreen',
                           'yellowgreen',
                           'limegreen',
                           'lawngreen',
                           'mediumseagreen',
                           'springgreen',
                           'darkgreen',
                           'c',
                           'red'])
    heatmap = ax.pcolor(
        hh1,
        cmap=cMap,
        alpha=0.8,
        edgecolor='k',
        clip_on=False,
        facecolor='w',
        lw=0.4)
    cbar = plt.colorbar(heatmap, orientation='vertical', pad=0.06)
    cbar.ax.set_yticklabels(
        cbar.ax.get_yticklabels(),
        **axis_font_col)  # vertical colorbar
    cbar.ax.get_xaxis().set_tick_params(direction='out', width=1)
    cbar.ax.get_yaxis().set_tick_params(direction='out', width=1)

    # Format
    fig = plt.gcf()
    fig.set_size_inches(12, 10)
    # plt.figure(figsize=(8, 6), dpi=80)

    # turn on the frame
    ax.set_frame_on(True)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(hh1.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(hh1.shape[1]) + 0.5, minor=False)

    # want a more natural, table-like display
    # ax.invert_yaxis()
    ax.invert_xaxis()
    # ax.xaxis.tick_top()
    ax.xaxis.tick_bottom()
    # Set the labels

    ax.set_xticklabels(nl, minor=False, **axis_font_h)
    ax.set_yticklabels(nl, minor=False, **axis_font_h)
    ax.set_title(
        "Overlap intensity between Leading Edge subsets",
        **title_font)
    plt.ylabel('Gene sets', size=20, weight='bold', color='darkblue')
    # rotate the
    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(
        'Leading_Edge_Analysis/Leading_Edge_similarity_matrix_heatmap.jpeg',
        bbox_inches='tight')
    plt.close()
    return


def heatmap_leading_edge_subset(LM, subset_dict, jac_f):
    store_leading_h1 = np.zeros((len(LM), len(LM)), dtype=np.float)
    o_clus = []
    name_l = []
    for i in LM:
        name_l.append(i[0])
    idx = list(range(len(LM)))
    pairs = list(itertools.combinations(idx, 2))
    for i in range(0, len(pairs)):
        h = pairs[i]
        F_id = intersect(LM[h[0]][3], LM[h[1]][3])
        if (len(F_id) > 0):
            for vv in range(0, len(F_id)):
                subset_dict[F_id[vv]] = subset_dict[F_id[vv]] + 1
            store_leading_h1[h[0]][h[1]] = (
                len(F_id) / (LM[h[0]][2] + LM[h[1]][2]))
            jac_f.append((len(F_id) / (LM[h[0]][2] + LM[h[1]][2])))
        else:
            jac_f.append(0)
    heatmap_plot(store_leading_h1, name_l)
    return jac_f, subset_dict


def heatmap_plot_assign(ASS_M, name_x, name_y, file_to_save, oi, cMap):
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(
        ASS_M,
        cmap=cMap,
        alpha=0.8,
        edgecolor='k',
        lw=1,
        clip_on=False,
        facecolor="white")
    fig = plt.gcf()
    fig.set_size_inches(90, 10)
    ax.set_frame_on(True)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(ASS_M.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(ASS_M.shape[1]) + 0.5, minor=False)
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(name_x, minor=False, **axis_font_h)
    ax.set_yticklabels(name_y, minor=False, **axis_font_h)

    if(oi == 1):
        cbar = plt.colorbar(heatmap, orientation='vertical', pad=0.01)
        cbar.ax.set_yticklabels(
            cbar.ax.get_yticklabels(),
            **axis_font_col)  # vertical colorbar
        cbar.ax.get_xaxis().set_tick_params(direction='out', width=1)
        cbar.ax.get_yaxis().set_tick_params(direction='out', width=1)

    # rotate the
    plt.xticks(rotation=90)
    plt.ylabel('Gene sets', size=20, weight='bold', color='darkblue')
    plt.xlabel(
        '\nLeading Edge genes',
        size=20,
        weight='bold',
        color='darkblue')
    plt.tight_layout()
    plt.savefig(file_to_save, bbox_inches='tight')
    plt.close()
    return


def leading_edge_clustering(SCB, SCB_2, name_x, name_y):
    # For leading edge genes.
    te = scipy.spatial.distance.pdist(np.transpose(SCB), metric='correlation')
    Z = linkage(te, 'average')
    Z1 = dendrogram(Z, no_plot=1)
    idx1 = Z1['leaves']
    name_l_x_leaf = []
    for i in idx1:
        name_l_x_leaf.append(name_x[i])
    plot_dendrogram_for_cluster(Z, name_l_x_leaf, 'Leading edge genes')

    ind = scipy.cluster.hierarchy.fcluster(Z, 0.5 * te.max(), 'distance')
    name_l_x_clus = []
    yy = np.argsort(ind)
    for i in yy:
        name_l_x_clus.append(name_x[i])

    # For geneset.
    del te
    del Z
    del Z1
    del idx1
    te = scipy.spatial.distance.pdist(SCB, metric='correlation')
    Z = linkage(te, 'average')
    Z1 = dendrogram(Z, no_plot=1)
    idx1 = Z1['leaves']
    name_l_y_leaf = []
    for i in idx1:
        name_l_y_leaf.append(name_y[i])
    plot_dendrogram_for_cluster(Z, name_l_y_leaf, 'Leading edge genesets')

    ind_gene = scipy.cluster.hierarchy.fcluster(Z, 0.5 * te.max(), 'distance')
    yy_gene = np.argsort(ind_gene)
    name_l_y_clus = []
    for i in yy_gene:
        name_l_y_clus.append(name_y[i])

    FINAL_BOOLEAN_MA_copy = SCB_2[yy_gene, :]
    FINAL_BOOLEAN_MA = FINAL_BOOLEAN_MA_copy[:, yy]
    ##########################################################################
    # CREATE CLUSTERED AND UNCLUSTERED GCT FILES.
    f = open("Leading_Edge_Analysis/CGAT_leading_edge_matrix_unclustered.gct", "w")
    f2 = open("Leading_Edge_Analysis/CGAT_leading_edge_matrix_clustered.gct", "w")
    f.write("#1.2" + "\n")
    f2.write("#1.2" + "\n")
    f.write(str(len(name_y)) + "\t" + str(len(name_x)) + "\n")
    f2.write(str(len(name_y)) + "\t" + str(len(name_x)) + "\n")
    f.write("NAME" + "\t" + "DESCRIPTION" + "\t")
    f2.write("NAME" + "\t" + "DESCRIPTION" + "\t")

    for i in range(0, len(name_x)):
        f.write(name_x[i])
        f2.write(name_l_x_clus[i])
        if(i == len(name_x) - 1):
            f.write("\n")
            f2.write("\n")
        else:
            f.write("\t")
            f2.write("\t")

    for i in range(0, len(name_y)):
        f.write(name_y[i] + "\t" + "na" + "\t")
        f2.write(name_l_y_clus[i] + "\t" + "na" + "\t")
        for za in range(0, len(name_x)):
            f.write(str(SCB[i][za]))
            f2.write(str(SCB[yy_gene[i]][yy[za]]))
            if(za == len(name_x) - 1):
                f.write("\n")
                f2.write("\n")
            else:
                f.write("\t")
                f2.write("\t")

    f.close()
    f2.close()
    ##########################################################################
    cMap = plt.cm.Reds
    heatmap_plot_assign(
        FINAL_BOOLEAN_MA,
        name_l_x_clus,
        name_l_y_clus,
        'Leading_Edge_Analysis/Leading_Edge_heatmap_clustered.jpeg',
        1,
        cMap)
    return


def plot_enrichment_score(
        store_enrichment_score,
        original_es_index,
        original_es,
        in_list,
        A_W,
        l_id):
    plt.figure(figsize=(8, 6), dpi=80)
    plt.xlim(-200, l_id + 50)
    plt.plot(
        store_enrichment_score[A_W],
        color="firebrick",
        linewidth=3.5,
        linestyle="-")
    plt.plot(
        original_es_index[A_W],
        original_es[A_W],
        'go',
        markersize=15,
        mfc='none',
        mew=2.5,
        label="ES")
    plt.legend(
        bbox_to_anchor=(
            0.97,
            1),
        loc=2,
        prop=legend_properties,
        handletextpad=0)
    # plt.axvline(x= original_es_index[aw], ymin=0, ymax=original_es[aw],
    # linewidth=3, color="g",linestyle=":")
    plt.axhline(y=0, linewidth=2, color="k", linestyle="--")
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    q = in_list[A_W]
    file_to_report_pic = ".".join(["enplot", q[0], 'jpeg'])
    q2 = ''.join(['\nGene List Index', '\n\nNumber of genes: ',
                  str(l_id), '(in list), ', str(q[3]), '(in gene set)'])
    plt.xlabel(q2, **axis_font)
    plt.ylabel('\n\nRunning enrichment score(RES)', **axis_font, labelpad=28)
    # plt.title(''.join(["Gene Set", ":", q[0]]), **title_font)
    plt.title(q[0], **title_font_ran)
    plt.grid()
    plt.tight_layout()
    plt.savefig(file_to_report_pic, bbox_inches='tight')
    plt.close()
    return


def plot_random_ES(store_permute, A_W, in_list):
    plt.figure(figsize=(8, 6), dpi=80)
    Rans_plot = store_permute[:, A_W]
    n, bins, patches = plt.hist(
        Rans_plot, 80, normed=False, color='dodgerblue', alpha=0.9, histtype='bar')
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    plt.xlabel('\nEnrichment Score', **axis_font)
    plt.ylabel('Frequency', **axis_font, labelpad=28)
    q = in_list[A_W]
    file_to_report_pic = ".".join(["random_enplot", q[0], 'jpeg'])
    #plt.title(''.join([q[0], ":Random ES Distribution"]), **title_font_ran)
    plt.title("Random ES Distribution", **title_font_ran)
    plt.grid(b=True, which='both', linestyle='--')
    plt.tight_layout()
    plt.savefig(file_to_report_pic, bbox_inches='tight')
    plt.close()
    return


def plot_enrichment_score_subplot(store_enrichment_score, original_es_index,
                                  original_es, in_list, A_W,
                                  l_id, s_p, s_index, en_name, lk):
    if(s_index == 1):
        plt.figure(figsize=(22, 60), dpi=80)
    plt.subplot(s_p, 3, s_index)
    plt.xlim(-200, l_id + 50)
    plt.plot(
        store_enrichment_score[A_W],
        color="firebrick",
        linewidth=4,
        linestyle="-")
    plt.plot(
        original_es_index[A_W],
        original_es[A_W],
        'go',
        markersize=14,
        mfc='none',
        mew=2.5,
        label="ES")
    plt.axhline(y=0, linewidth=2, color="k", linestyle="--")
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    q = in_list[A_W]
    # file_to_report_pic = ".".join(["enplot",q[0],'jpeg'])
    if(s_index == 10):
        plt.ylabel(
            '\n\nR u n n i n g    e n r i c h m e n t    s c o r e (RES)',
            fontsize=20,
            weight='bold',
            labelpad=40)
    # plt.title(''.join(["Gene Set", ":", q[0]]),
        # fontsize=18, color='darkblue', weight='bold')
    plt.title(q[0],
              fontsize=16, color='darkblue', weight='bold')
    plt.grid()
    if(s_index == lk):
        q2 = ''.join(
            ['\n\nGene List Index', '\n\nNumber of genes: ', str(l_id), '(in list)'])
        plt.xlabel(q2, fontsize=20, weight='bold')
        plt.tight_layout(pad=0.1, w_pad=9, h_pad=1.0)
        plt.savefig(en_name, bbox_inches='tight')
        # plt.show()
        plt.close()
    return


def plot_dendrogram_for_cluster(ZZ, name_leaf, text_to_save):
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 10)
    fancy_dendrogram(
        ZZ,
        truncate_mode='lastp',
        p=20,                   # show only the last p merged clusters
        leaf_rotation=90.,
        leaf_font_size=14.,
        show_contracted=True,
        annotate_above=5,  # useful in small plots so annotations don't overlap
    )
    plt.yticks(fontsize=12, weight='bold')
    plt.xticks(weight='bold')
    plt.tight_layout()
    t_s = "Leading_Edge_Analysis/Last_20_merged_cluster_dendrogram_" + text_to_save + ".jpeg"
    plt.savefig(t_s, bbox_inches='tight')
    plt.close()
    fig, ax = plt.subplots()
    del t_s
    fig.set_size_inches(50, 10)
    dendrogram(
        ZZ,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=6.,  # font size for the x axis labels
        labels=name_leaf,
    )
    plt.yticks(fontsize=12, weight='bold')
    plt.xticks(weight='bold')
    plt.title('Hierarchical Clustering Dendrogram', **title_font)
    plt.xlabel(text_to_save, fontsize=18, weight='bold')
    plt.ylabel('Distance', fontsize=18, weight='bold', labelpad=28)
    plt.tight_layout()
    t_s = "Leading_Edge_Analysis/Hierarchical_cluster_dendrogram_" + text_to_save + ".jpeg"
    plt.savefig(t_s, bbox_inches='tight')
    plt.close()
    return


def plot_summary_report(nui, fg1, ufp, g_set, c):
    """
    nui = nes_up_index
    fg1 = original_nes
    ufp = up_for_plot
    g_set=in_list
    c="upregulated"
    """
    pl_x = []
    pl_y = []
    for i in range(0, 20):
        aw = ufp[nui[i]]
        itt = g_set[aw]
        pl_y.append(fg1[aw])
        pl_x.append(itt[0])
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 14)
    ax.set_frame_on(True)
    plt.barh(
        list(
            range(
                0,
                20)),
        pl_y,
        color='orange',
        align="center",
        height=0.6)
    cc = "Gene sets (" + c + ") by treatment"
    plt.ylabel(cc, labelpad=28, fontsize=25, weight='bold', color='darkblue')
    ax.xaxis.tick_top()
    plt.yticks(fontsize=16, weight='bold')
    plt.yticks(range(21), pl_x)
    plt.xticks(fontsize=16, weight='bold')
    ax.grid(which='major', linestyle='-', linewidth='0.3')
    plt.title("Normalized Enrichment Score(NES)", y=1.08,
              fontsize=25, weight='bold', color='darkblue')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    c2 = "Top_20_" + c + "_genesets_by_treatment.jpeg"
    plt.savefig(c2, bbox_inches='tight')
    plt.close()
    return


def plot_summary_report_for_fdr(pl_x, pl_y):
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 14)
    ax.set_frame_on(True)
    plt.barh(
        list(
            range(
                0,
                20)),
        pl_y,
        color='c',
        align="center",
        height=0.7,
        alpha=0.9)
    plt.axvline(x=0, linewidth=2, color="k")
    plt.yticks(range(21), pl_x)
    cc = "Top 20 enriched gene sets" + "\n(sorted by FDR-q values)"
    plt.ylabel(cc, labelpad=30, fontsize=25, weight='bold', color='darkblue')
    ax.xaxis.tick_top()
    plt.yticks(fontsize=16, weight='bold')
    ax.grid(which='major', linestyle='-', linewidth='0.4')
    plt.xticks(fontsize=16, weight='bold')
    plt.title("Normalized Enrichment Score(NES)", y=1.08,
              fontsize=25, weight='bold', color='darkblue')
    plt.gca().invert_yaxis()
    plt.grid(b=True, which='minor', linestyle='--')
    plt.tight_layout()
    plt.savefig(
        "Top_20_enriched_genesets(fdr_sorted).jpeg",
        bbox_inches='tight')
    plt.close()
    return


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-f",
        "--filename",
        dest="file_name",
        type="string",
        help="Expression data set(ranked list of genes) [default=%default].")

    parser.add_option(
        "-g",
        "--geneset",
        dest="geneset",
        type="string",
        help="Annotated gene sets database for enrichment analysis"
        "[default=%default].")

    parser.add_option(
        "-m",
        "--minimum",
        dest="min_gene",
        type="int",
        help="gene sets smaller than this are excluded from the analysis "
        "[default=%default].")

    parser.add_option(
        "-x",
        "--maximum",
        dest="max_gene",
        type="int",
        help="gene sets larger than this are excluded from the analysis "
        "[default=%default].")

    parser.add_option(
        "-s",
        "--randomseed",
        dest="seed",
        type="int",
        help="A number use to initialize a pseudorandom number generator"
        " [default=%default].")

    parser.add_option(
        "-n",
        "--permutation",
        dest="iteration",
        type="int",
        help="Number of permutations to perform in assessing the statistical significance of the enrichment score [default=%default].")

    parser.add_option(
        "-d",
        "--display",
        dest="plot_no",
        type="int",
        help="Displays enrichment plots for the specified no. of gene sets with the highest absolute normalized enrichment scores (each phenotype) [default=%default].")

    parser.add_option(
        "-l",
        "--num_leading",
        dest="fdr_num",
        type="int",
        help="Number of genesets for leading edge analysis, by default top 11 enriched genesets will be used for this analysis."
        "Minimum number of genesets should be 4. [default=%default].")

    parser.set_defaults(
        file_name=None,
        geneset=None,
        min_gene=25,
        max_gene=500,
        seed=42,
        iteration=1000,
        plot_no=20,
        fdr_num=10,
    )
    (options, args) = E.Start(parser, add_database_options=True)
    # Preprocess expression file.
    id, expression_value = read_expression(options.file_name)

    # Preprocess geneset
    ex_list, in_list, geneset_indicator = preprocess_geneset(
        options.geneset, options.min_gene, options.max_gene, id)
    generate_gen_set_report(
        ex_list,
        in_list,
        options.min_gene,
        options.max_gene,
        geneset_indicator)

    # Store filtered geneset. Because in_list contains annotation and
    # description also,for each geneset.
    GG = [item[4:len(item)] for item in in_list]

    # Store index of each id (from the ranked list of expression data as a
    # dicitionary)
    ind_dict = dict((k, i) for i, k in enumerate(id))

    # Calculate enrichment score for each geneset.
    # Create boolean array for total number of id in expression data and float
    # array for Enrichment Score of options.iteration permutation.
    temp = np.zeros((len(id),), dtype=np.bool)
    temp_2 = np.zeros((len(id),), dtype=np.float)
    store_permute = np.zeros((options.iteration, len(GG)), dtype=np.float)
    original_es = np.zeros((len(GG),), dtype=np.float)
    original_es_index = np.zeros((len(GG),), dtype=np.int)
    original_nes = np.zeros((len(GG),), dtype=np.float)
    permute_nes = np.zeros((options.iteration, len(GG)), dtype=np.float)
    nominal_p = np.zeros((len(GG),), dtype=np.float)
    nominal_p_nes = np.zeros((len(GG),), dtype=np.float)
    store_enrichment_score = []
    store_gene_leading_info = np.zeros((3, len(GG)), dtype=np.int)
    store_gene_leading_matrix = []

    count = 0
    for i in GG:
        inter = intersect(ind_dict, i[0])
        indices = sorted([ind_dict[x] for x in inter])
        S = len(i[0])
        enrich_score = calculate_enrichment_score(
            indices, temp, expression_value, S, temp_2)
        del indices
        del inter
        a1 = [
            np.absolute(
                np.max(enrich_score)), np.absolute(
                np.min(enrich_score))]
        if(np.absolute(np.max(enrich_score)) == np.absolute(np.min(enrich_score))):
            original_es[count] = np.max(enrich_score)
            original_es_index[count] = np.argmax(enrich_score)
        else:
            if(a1.index(max(a1)) == 0):
                original_es[count] = np.max(enrich_score)
                original_es_index[count] = np.argmax(enrich_score)
            else:
                original_es[count] = np.min(enrich_score)
                original_es_index[count] = np.argmin(enrich_score)
        store_enrichment_score.append(list(enrich_score))
        # This section has been added by me for "Leading Edge Analysis".
        t0 = original_es_index[count]
        if(original_es[count] < 0):
            se = temp[t0:len(id)]
        else:
            se = temp[0:t0 + 1]
        c0 = np.sum(se == 1)
        c1 = len(se)
        tag = (c0 * 100) / S
        tag2 = (c0 * 100) / len(id)
        gene_l = (c1 * 100) / len(id)
        signal_l = (tag / 100) * (1 - (gene_l / 100)) * \
            (len(id) / (len(id) - S))
        store_gene_leading_info[0][count] = tag
        store_gene_leading_info[1][count] = gene_l
        store_gene_leading_info[2][count] = signal_l * 100
        temp_2.fill(0)
        temp.fill(0)
        del a1
        del enrich_score
        count = count + 1

    print("Enrichment score calculation has been successfully completed")

    # Calculate Randon Background by random permutation of gene set.
    # calculate total element in gene set.
    s = 0
    for t in GG:
        s = s + len(t[0])
    size_info = [len(t[0]) for t in GG]

    # Random index generation
    np.random.seed(options.seed)
    id_new = np.array(id)

    # Calculate enrichment score for permuted genesets
    count = 0
    for per in range(0, options.iteration):
        if((per % 100) == 0):
            print(per)
        for i in range(0, len(size_info)):
            tt = size_info[i]
            tar = id_new[np.random.randint(len(id), size=(1, tt))]
            inter = intersect(ind_dict, tar[0])
            indices = sorted([ind_dict[x] for x in inter])
            S = len(tar[0])
            enrich_score = calculate_enrichment_score(
                indices, temp, expression_value, S, temp_2)
            del indices
            del inter
            a1 = [
                np.absolute(
                    np.max(enrich_score)), np.absolute(
                    np.min(enrich_score))]
            if(np.absolute(np.max(enrich_score)) == np.absolute(np.min(enrich_score))):
                store_permute[per, i] = np.max(enrich_score)
            else:
                if(a1.index(max(a1)) == 0):
                    store_permute[per, i] = np.max(enrich_score)
                else:
                    store_permute[per, i] = np.min(enrich_score)
            del a1
            temp.fill(0)
            temp_2.fill(0)
            del enrich_score

    print("Enrichment score calculation for permuted sets has been successfully completed")

    # Calculation of empirical p-value.
    for i in range(0, len(original_es)):
        if(original_es[i] >= 0):
            A1 = store_permute[:, i] >= original_es[i]
            A2 = store_permute[A1, i]
            nominal_p[i] = np.divide(len(A2), options.iteration)
        else:
            A1 = store_permute[:, i] <= original_es[i]
            A2 = store_permute[A1, i]
            nominal_p[i] = np.divide(len(A2), options.iteration)

    # Normalization of enrichment score

    for i in range(0, len(GG)):
        if(original_es[i] < 0):
            A1 = store_permute[:, i] < 0
            A2 = np.mean(np.absolute(store_permute[A1, i]))
            original_nes[i] = np.divide(original_es[i], A2)
        else:
            A1 = store_permute[:, i] >= 0
            A2 = np.mean(np.absolute(store_permute[A1, i]))
            original_nes[i] = np.divide(original_es[i], A2)

    # Normalization of enrichment score for each permutation.

    for per in range(0, options.iteration):
        for i in range(0, len(GG)):
            if(store_permute[per, i] < 0):
                A1 = store_permute[:, i] < 0
                A2 = np.mean(np.absolute(store_permute[A1, i]))
                permute_nes[per, i] = np.divide(store_permute[per, i], A2)
            else:
                A1 = store_permute[:, i] >= 0
                A2 = np.mean(np.absolute(store_permute[A1, i]))
                permute_nes[per, i] = np.divide(store_permute[per, i], A2)

    print("Normalization has been successfully completed")

    # Calculation of empirical p-value for normalized enrichment score.

    for i in range(0, len(original_nes)):
        if(original_nes[i] >= 0):
            A1 = permute_nes[:, i] >= original_nes[i]
            A2 = permute_nes[A1, i]
            nominal_p_nes[i] = np.divide(len(A2), options.iteration)
        else:
            A1 = permute_nes[:, i] <= original_nes[i]
            A2 = permute_nes[A1, i]
            nominal_p_nes[i] = np.divide(len(A2), options.iteration)

    nes_up_index = np.array([], dtype=np.int)
    nes_down_index = np.array([], dtype=np.int)
    # Extract two set NES>0 and NES<0
    if(sum(original_nes > 0) > 0):
        NN = original_nes[original_nes > 0]
        up_for_plot = np.where(original_nes >= 0)[0]
        nes_up = -np.sort(-NN)
        nes_up_index = np.argsort(-NN)
        del NN
    if(sum(original_nes < 0) > 0):
        NN = original_nes[original_nes < 0]
        down_for_plot = np.where(original_nes < 0)[0]
        nes_down = np.sort(NN)
        nes_down_index = np.argsort(NN)
        del NN

    # Calculate FDR
    fdr_upregulated = np.array([], dtype=np.float)
    fdr_downregulated = np.array([], dtype=np.float)
    if(len(nes_up_index) > 0):
        fdr_upregulated = np.zeros((len(nes_up),), dtype=np.float)
    if(len(nes_down_index) > 0):
        fdr_downregulated = np.zeros((len(nes_down),), dtype=np.float)

    # For upregulated:
    if(len(nes_up_index) > 0):
        p_value_up = nominal_p_nes[original_nes > 0]
        b, fdr_upregulated, w1, w2 = sm.multipletests(
            p_value_up[nes_up_index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    # For downregulated
    if(len(nes_down_index) > 0):
        p_value_down = nominal_p_nes[original_nes < 0]
        b, fdr_downregulated, w1, w2 = sm.multipletests(
            p_value_down[nes_down_index], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    print("FDR calculation has been successfully completed")

    # Generate graphical report and detail table report also
    kk = os.path.basename(options.file_name)
    part1, part2 = kk.split('.')
    kk = os.path.basename(options.geneset)
    part3, part4 = kk.split('.')
    nam1 = "CGAT_REPORT_FOR_upregulated"
    nam2 = "CGAT_REPORT_FOR_downregulated"
    nam3 = "CGAT_REPORT_FOR_fdr_sorted_up_down"
    nam4 = "CGAT_Summary_Report"
    file_to_report_1 = ".".join([nam1, part1, "xls"])
    file_to_report_2 = ".".join([nam2, part1, "xls"])
    file_to_report_3 = ".".join([nam3, part1, "xls"])
    file_to_report_4 = ".".join([nam4, part1, "txt"])

    # PREPARE FINAL SUMMARY REPORT.
    f = open(file_to_report_4, "w")
    f.write("Enrichment in phenotype (upregulated):" + "\n")
    f.write("@ " + str(len(nes_up_index)) + "/" +
            str(len(in_list)) + " gene sets are upregulated" + "\n")
    f.write("@ " + str(sum(fdr_upregulated < 0.05)) +
            " gene sets are significant at FDR < 5%" + "\n")
    f.write("@ " + str(sum(fdr_upregulated < 0.01)) +
            " gene sets are significant at FDR < 1%" + "\n")
    w1 = nominal_p[original_nes >= 0]
    w2 = nominal_p[original_nes < 0]
    f.write("@ " +
            str(sum(w1 < 0.01)) +
            " gene sets are significantly enriched at nominal pvalue < 1%" +
            "\n")
    f.write("@ " +
            str(sum(w1 < 0.05)) +
            " gene sets are significantly enriched at nominal pvalue < 5%" +
            "\n\n")
    f.write("Enrichment in phenotype (downregulated):" + "\n")
    f.write("@ " + str(len(nes_down_index)) + "/" +
            str(len(in_list)) + " gene sets are upregulated" + "\n")
    f.write("@ " + str(sum(fdr_downregulated < 0.05)) +
            " gene sets are significant at FDR < 5%" + "\n")
    f.write("@ " + str(sum(fdr_downregulated < 0.01)) +
            " gene sets are significant at FDR < 1%" + "\n")
    f.write("@ " +
            str(sum(w2 < 0.01)) +
            " gene sets are significantly enriched at nominal pvalue < 1%" +
            "\n")
    f.write("@ " +
            str(sum(w2 < 0.05)) +
            " gene sets are significantly enriched at nominal pvalue < 5%" +
            "\n\n")
    f.close()

    # Plot enrichment score of top gene set for each phenotype.

    # If specified number of top genesets for plotting enrichemnt score is higher than the total number of genesets.I will
    # plot enrichemnt score for all genesets.
    if(options.plot_no >= len(GG)):
        options.plot_no = len(GG)
    xcc = int(options.plot_no / 2) + 1
    for i in range(0, options.plot_no):
        if(len(nes_up_index) > 0):
            aw = up_for_plot[nes_up_index[i]]
            plot_enrichment_score(
                store_enrichment_score,
                original_es_index,
                original_es,
                in_list,
                aw,
                len(id))
            plot_random_ES(store_permute, aw, in_list)

        # Downregulated
        if(len(nes_down_index) > 0):
            aw = down_for_plot[nes_down_index[i]]
            plot_enrichment_score(
                store_enrichment_score,
                original_es_index,
                original_es,
                in_list,
                aw,
                len(id))
            plot_random_ES(store_permute, aw, in_list)
    for i in range(0, options.plot_no):
        if(len(nes_up_index) > 0):
            aw = up_for_plot[nes_up_index[i]]
            plot_enrichment_score_subplot(
                store_enrichment_score,
                original_es_index,
                original_es,
                in_list,
                aw,
                len(id),
                xcc,
                i + 1,
                "enplot_upregulated_summary.jpeg",
                options.plot_no)
    for i in range(0, options.plot_no):
        if(len(nes_down_index) > 0):
            aw = down_for_plot[nes_down_index[i]]
            plot_enrichment_score_subplot(
                store_enrichment_score,
                original_es_index,
                original_es,
                in_list,
                aw,
                len(id),
                xcc,
                i + 1,
                "enplot_downregulated_summary.jpeg",
                options.plot_no)
    # Plot summary of top 20 genesets of each phenotype
    if(len(nes_up_index) > 0):
        plot_summary_report(nes_up_index, original_nes,
                            up_for_plot, in_list, "upregulated")
    if(len(nes_down_index) > 0):
        plot_summary_report(nes_down_index, original_nes,
                            down_for_plot, in_list, "downregulated")

    # Generate pvalue vs ES graph.
    if(len(nes_down_index) > 0):
        aw = down_for_plot[nes_down_index]
    if(len(nes_up_index) > 0):
        aw2 = up_for_plot[nes_up_index]
    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(original_nes, nominal_p, 'ko')
    if(len(nes_up_index) > 0):
        plt.plot(nes_up, nominal_p[aw2], 'k-')
    if(len(nes_down_index) > 0):
        plt.plot(nes_down, nominal_p[aw], 'k-')
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    plt.xlabel("\nNormalized Enrichment Score(NES)", **axis_font)
    plt.ylabel('\np-value', **axis_font, labelpad=28)
    plt.title('NES vs Significance', **title_font)
    plt.grid()
    plt.tight_layout()
    ax2 = plt.twinx()
    if(len(nes_up_index) > 0):
        ax2.plot(
            nes_up,
            fdr_upregulated,
            marker='s',
            color='firebrick',
            linestyle='')
    if(len(nes_down_index) > 0):
        ax2.plot(
            nes_down,
            fdr_downregulated,
            marker='s',
            color='firebrick',
            linestyle='')
    ticks_font = font_manager.FontProperties(
        family='Helvetica',
        style='normal',
        size=12,
        weight='bold',
        stretch='normal')
    for label in ax2.get_yticklabels():
        label.set_fontproperties(ticks_font)
    ax2.tick_params('y', colors='firebrick')
    # ax2.tick_params('y', colors='firebrick',labelsize = 'x-large',width=20)
    ax2.yaxis.set_tick_params(
        'y',
        colors='firebrick',
        labelsize='x-large',
        width=20)
    ax2.set_ylabel('\nFDR q-value', **axis_font)
    plt.savefig('Genesets_null_distribution.jpeg', bbox_inches='tight')
    plt.close()

    # The histogram of the ES across all genesets
    plt.figure(figsize=(11, 6), dpi=80)
    n, bins, patches = plt.hist(original_nes, 50, normed=False, facecolor='darkgreen',
                                alpha=0.9, histtype='step', lw=4, color='darkgreen')
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    plt.xlabel('\nNormalized Enrichment scores across genesets', **axis_font)
    plt.ylabel('Frequency', **axis_font, labelpad=28)
    plt.title('Histogram of Normalized Enrichment Score', **title_font)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('Genesets_NES_Histogram.jpeg', bbox_inches='tight')
    plt.close()

    print("Graphical reports has been successfully completed")

    # Write detail report for upregulated gene.
    if(len(nes_up_index) > 0):
        with open(file_to_report_1, 'w') as f:
            f.write("NAME\tGENE SET(GS)\tGS DETAILS\tSIZE\tENRICHMENT SCORE(ES)\tNORMALIZED ENRICHMENT SCORE(NES)\tEMPERICAL p-value\tFDR-q value\tRANK AT MAX\tLEADING EDGE\n")
            for i in range(0, len(nes_up_index)):
                aw = up_for_plot[nes_up_index[i]]
                itt = in_list[aw]
                f.write(itt[0] +
                        "\t" +
                        itt[0] +
                        "\t" +
                        itt[1] +
                        "\t" +
                        str(itt[3]) +
                        "\t" +
                        str(original_es[aw]) +
                        "\t" +
                        str(original_nes[aw]) +
                        "\t" +
                        str(nominal_p[aw]) +
                        "\t" +
                        str(fdr_upregulated[i]) +
                        "\t" +
                        str(original_es_index[aw]) +
                        "\t" +
                        "tags=" +
                        str(store_gene_leading_info[0][aw]) +
                        "%, list=" +
                        str(store_gene_leading_info[1][aw]) +
                        "%" +
                        "\n")
            f.close()

    # Write detail report for downregulated gene.
    if(len(nes_down_index) > 0):
        with open(file_to_report_2, 'w') as f:
            f.write("NAME\tGENE SET(GS)\tGS DETAILS\tSIZE\tENRICHMENT SCORE(ES)\tNORMALIZED ENRICHMENT SCORE(NES)\tEMPERICAL p-value\tFDR-q value\tRANK AT MAX\tLEADING EDGE\n")
            for i in range(0, len(nes_down_index)):
                aw = down_for_plot[nes_down_index[i]]
                itt = in_list[aw]
                f.write(itt[0] +
                        "\t" +
                        itt[0] +
                        "\t" +
                        itt[1] +
                        "\t" +
                        str(itt[3]) +
                        "\t" +
                        str(original_es[aw]) +
                        "\t" +
                        str(original_nes[aw]) +
                        "\t" +
                        str(nominal_p[aw]) +
                        "\t" +
                        str(fdr_downregulated[i]) +
                        "\t" +
                        str(original_es_index[aw]) +
                        "\t" +
                        "tags=" +
                        str(store_gene_leading_info[0][aw]) +
                        "%, list=" +
                        str(store_gene_leading_info[1][aw]) +
                        "%" +
                        "\n")
            f.close()

    # Write detail report of up and down regulated genes in to a single file,
    # sorted by FDR values
    merge_mat = np.append(fdr_upregulated, fdr_downregulated)
    merge_index = np.argsort(np.append(fdr_upregulated, fdr_downregulated))
    merge_x = []
    merge_y = []
    count = 0
    with open(file_to_report_3, 'w') as f:
        f.write("NAME\tGENE SET(GS)\tGS DETAILS\tSIZE\tENRICHMENT SCORE(ES)\tNORMALIZED ENRICHMENT SCORE(NES)\tEMPERICAL p-value\tFDR-q value\tRANK AT MAX\tLEADING EDGE\n")
        for i in merge_index:
            CW = merge_mat[i]
            if(i >= len(fdr_upregulated)):
                aw = down_for_plot[nes_down_index[i - len(fdr_upregulated)]]
            else:
                aw = up_for_plot[nes_up_index[i]]
            itt = in_list[aw]
            ######################################
            # This section is for summary plot of enriched genesets.
            if(count < 20):
                merge_x.append(itt[0])
                merge_y.append(original_nes[aw])
            count = count + 1
            #####################################
            f.write(itt[0] +
                    "\t" +
                    itt[0] +
                    "\t" +
                    itt[1] +
                    "\t" +
                    str(itt[3]) +
                    "\t" +
                    str(original_es[aw]) +
                    "\t" +
                    str(original_nes[aw]) +
                    "\t" +
                    str(nominal_p[aw]) +
                    "\t" +
                    str(CW) +
                    "\t" +
                    str(original_es_index[aw]) +
                    "\t" +
                    "tags=" +
                    str(store_gene_leading_info[0][aw]) +
                    "%, list=" +
                    str(store_gene_leading_info[1][aw]) +
                    "%" +
                    "\n")
        f.close()
    # Plot top 20 enriched genset sorted by FDR values.
    plot_summary_report_for_fdr(merge_x, merge_y)
    del merge_mat
    del merge_index
    del merge_x
    del merge_y

    print("Reports have been successfully generated")

    '''
    ********************************************************************************************************************
                                            LEADING EDGE ANALYSIS
    *********************************************************************************************************************
    THIS SECTION IS FOR LEADING EDGE ANALYSIS.
    '''

    print("Leading edge analysis has been started")
    # Make a directory.
    newpath = r'Leading_Edge_Analysis'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    store_gene_leading_matrix = generate_leading_edge_m(
        nes_up_index,
        nes_down_index,
        fdr_upregulated,
        fdr_downregulated,
        in_list,
        options.fdr_num,
        ind_dict,
        down_for_plot,
        up_for_plot,
        temp,
        original_es,
        original_es_index,
        id_new,
        store_gene_leading_matrix,
        GG,
        store_gene_leading_info)

    # PREPARE LEADING EDGE SUMMARY FILE.
    f = open("Leading_Edge_Analysis/CGAT_LEADING_EDGE_ANALYSIS_SUMMARY.tsv", "w")
    f.write("Details of gene sets and signal used in this analysis" + "\n")
    f.write("There were " +
            str(options.fdr_num) +
            " gene sets used in the leading edge analysis (see below for details)" +
            "\n\n")
    f.write(
        "Gene set" +
        "\t" +
        "Members" +
        "\t" +
        "Members in signal" +
        "\t" +
        "Tag%" +
        "\t" +
        "List%" +
        "\t" +
        "FDR-q value" +
        "\n")
    for i in store_gene_leading_matrix:
        c = i
        f.write(c[0] + "\t" + str(c[1]) + "\t" + str(c[2]) + "\t" +
                str(c[4]) + "%\t" + str(c[5]) + "%\t" + str(c[6]) + "\n")
    f.close()

    # PREPARE LEADING EDGE MATRIX FILE.
    with open("Leading_Edge_Analysis/CGAT_leading_edge_matrix_for_results.gmx", 'w') as f:
        for_pd_in = []
        count = 0
        for i in store_gene_leading_matrix:
            c = i
            for_pd_in.append(c[3])
            f.write(c[0] + "_signal")
            if(count == (len(store_gene_leading_matrix) - 1)):
                f.write("\n")
            else:
                f.write("\t")
            count = count + 1

        count = 0
        for i in range(0, len(store_gene_leading_matrix)):
            f.write("na")
            if(count == (len(store_gene_leading_matrix) - 1)):
                f.write("\n")
            else:
                f.write("\t")
            count = count + 1
        v = pd.DataFrame(for_pd_in)
        v1 = v.transpose()
        v1.to_csv(f, encoding='utf-8', index=False, sep='\t', header=False)
        f.close()

    # PREPARE PLOT FOR GENE IN SUBSETS
    G1 = {}
    subset_dict = {}
    subset_dict_assign_matrix = {}
    name_l_x = []
    name_l_y = []
    # exp_for_clustering=[]

    for i in for_pd_in:
        G1 = i
        c1 = dict.fromkeys(G1, 1)
        subset_dict.update(c1)
        subset_dict_assign_matrix.update(c1)

    count = 0
    for j in subset_dict_assign_matrix.keys():
        subset_dict_assign_matrix[j] = count
        count = count + 1
        name_l_x.append(j)
    name_l_y = []
    for i in store_gene_leading_matrix:
        name_l_y.append(i[0] + "_signal")

    # CREATE ASSIGNMENT MATRIX FOR LEADING EDGE SUBSET
    store_unclustered_boolean = np.zeros(
        (len(store_gene_leading_matrix),
         len(subset_dict_assign_matrix)),
        dtype=np.int)
    store_clustered_boolean = np.zeros(
        (len(store_gene_leading_matrix),
         len(subset_dict_assign_matrix)),
        dtype=np.float)
    count = 0

    for i in store_gene_leading_matrix:
        c = i
        for jj in range(0, len(c[3])):
            xx = subset_dict_assign_matrix[c[3][jj]]
            store_unclustered_boolean[count][xx] = 1
            store_clustered_boolean[count][xx] = expression_value[ind_dict[c[3][jj]]]
        count = count + 1

    # PREPARE HEAT MAP FOR OVERLAPPING GENE SET UNCLUSTERED.
    jac_l = []
    # intensity_for_cluster = []
    jac_l, subset_dict = heatmap_leading_edge_subset(
        store_gene_leading_matrix, subset_dict, jac_l)

    # HEATMAP OF UNCLUSTERED AND CLUSTERED ASSIGNMENT MATRIX
    cMap = ListedColormap(['white', 'red'])
    heatmap_plot_assign(
        store_unclustered_boolean,
        name_l_x,
        name_l_y,
        'Leading_Edge_Analysis/Leading_Edge_heatmap_unclustered.jpeg',
        0,
        cMap)
    leading_edge_clustering(
        store_unclustered_boolean,
        store_clustered_boolean,
        name_l_x,
        name_l_y)

    # PLOT DISTRIBUTION oF GENES AMONG LEADING SUBSETS
    temp_dict_k = []
    temp_dict_v = []
    for i in subset_dict.keys():
        if(subset_dict[i] > 1):
            temp_dict_k.append(i)
            temp_dict_v.append(subset_dict[i])

    fig, ax = plt.subplots()
    ax.grid(zorder=3)
    ax.xaxis.grid()
    ax.yaxis.labelpad = 28
    plt.bar(range(len(temp_dict_k)), temp_dict_v, width=0.3, color='magenta')
    plt.xticks(range(len(temp_dict_k)), temp_dict_k)
    # Format
    fig = plt.gcf()
    fig.set_size_inches(12, 10)
    plt.yticks(fontsize=12, weight='bold')
    # plt.yticks(np.arange(min(temp_dict_v), max(temp_dict_v) + 1, 1.0))
    plt.xlabel('\nGenes', **axis_font)
    plt.ylabel('\nNumber Of Gene Sets', **axis_font)
    plt.title('Distribution of genes among leading edge subsets', **title_font)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=90, **axis_font_h)
    plt.tight_layout()
    plt.savefig(
        'Leading_Edge_Analysis/Genes_in_leading_edge_subset.jpeg',
        bbox_inches='tight')
    plt.close()
    del temp_dict_k
    del temp_dict_v

    # PLOT HISTOGRAM OF JACQUARD
    plt.figure(figsize=(11, 6), dpi=80)
    n, bins, patches = plt.hist(
        jac_l, 70, normed=1, color='darkgreen', alpha=0.9)
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')
    plt.xlabel('\nJacquard', **axis_font)
    plt.ylabel('Number of occurences', **axis_font, labelpad=28)
    plt.title('Overlapping among leading edge subset pairs', **title_font)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        'Leading_Edge_Analysis/Jacquard_distribution.jpeg',
        bbox_inches='tight')
    plt.close()

    print("Leading edge analysis has been finished")
    E.Stop()


if __name__ == "__main__":
    sys.exit(main())
