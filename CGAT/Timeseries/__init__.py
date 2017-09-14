'''
Timeseries.py - functions associated with analysis of time series data
======================================================================
:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Functions are split into::
 * clustering assessment
 * data transformation and normalisation
 * differential expression analysis
 * clustering and distance metric functions
'''
import sklearn.metrics.cluster.supervised as supervised
from math import log
import CGAT.Experiment as E
import numpy as np
import pandas as pd
import itertools
import os
import math
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
import rpy2.robjects as ro
import random
from CGAT.Timeseries import cmetrics as c2m


def get_r_path():
    """return path of R support functions.
    """
    return os.path.dirname(__file__)



def get_label_map(labels):
    '''
    return a dictionary with integer:string mapping
    '''
    label_set = set()
    map_dict = {}
    for val in labels:
        label_set.update(val)
    for lab, integer in enumerate(label_set):
        map_dict[integer] = lab

    return map_dict


def make_mapped_matrix(map_dict, input_frame):
    '''
    return a matrix with integer labels from mapping
    '''

    frame_index = input_frame.index.tolist()
    nindex = len(frame_index)
    ncols = len(input_frame.columns)
    integer_matrix = np.ndarray((nindex, ncols),
                                dtype=np.int32)

    E.info("mapping cluster labels")
    matrix_idx = [h for h, g in enumerate(frame_index)]
    for idx in matrix_idx:
        for col in range(ncols):
            mod = input_frame.iloc[idx][col+1]
            integer_matrix[idx][col] = map_dict[mod]

    return integer_matrix


def randIndexes(clustering_results):
    '''
    Calculate Rand index and adjusted Rand index over pairwise
    clustering comparisons.
    Use cythonised function to calculate indices
    '''

    # reassign module and gene labels with integer ids, integer comparison is
    # much faster than string comparison
    cluster_labels = clustering_results.values
    map_dict = get_label_map(cluster_labels)

    gene_map = {}
    for r, gene in enumerate(clustering_results.index):
        gene_map[gene] = r
    E.info("mapping gene ids")

    integer_matrix = make_mapped_matrix(map_dict, clustering_results)
    # take a small slice of the matrix for testing 5 genes, 3 clusterings

    E.info("counting clustering consensus")
    # use cythonized function to return rand index matrix
    cy_rand = c2m.consensus_metrics(integer_matrix)
    E.info("Rand Index calculated for all clusterings")

    return cy_rand


def unravel_arrays(metric_array):
    '''
    Unravel a numpy array such that only one half of the symmetrical
    matrix is output.  Do not output diagonal values.
    '''

    dim = metric_array.shape[0]
    flat_array = []

    for indx in itertools.combinations(range(0, dim), r=2):
        if indx[0] != indx[1]:
            flat_array.append(metric_array[indx[1], indx[0]])
        else:
            pass
    return flat_array


def mutualInformation(cluster1, cluster2):
    '''
    Calculate the mutual information for a given pair of clusterings.
    Assume clustering1 represents the reference ground truth.
    Code from scikit-learn.

    The mutual information of two clusterings, U and V is given by:

    MI(U, V) = sum(R)sum(C)P(i, j) [log(P(i, j)/(p(i)p'(j)))]

    This is the similarity of two clustering labels, where P(i) is
    the probability of a random sample in clustering Ui, and P'(j)
    the probability of a random sample in clustering Vj.
    '''

    cont = contingency(cluster1, cluster2)
    cont_sum = np.sum(cont)
    pi = np.sum(cont, axis=1)
    pj = np.sum(cont, axis=0)
    outer = np.outer(pi, pj)
    nnz = cont != 0

    cont_nm = cont[nnz]
    log_cont_nm = np.log(cont_nm)
    cont_nm /= cont_sum
    log_outer = -np.log(outer[nnz]) + log(pi.sum()) + log(pj.sum())
    mi = (cont_nm * (log_cont_nm - log(cont_sum)) + (cont_nm * log_outer))

    return mi.sum()


def contingency(cluster1, cluster2):
    '''
    Return a n x m matrix of clustering overlaps, where n
    is the number of clusters in clustering1 and m is the
    number of clusters in clustering2.  Return an np array.
    '''
    cont = pd.DataFrame(columns=cluster1.keys(), index=cluster2.keys())
    cont = cont.fillna(0.0)

    for x in itertools.product(cluster1.keys(), cluster2.keys()):
        set1 = cluster1[x[0]]
        set2 = cluster2[x[1]]
        intersect = len(set1.intersection(set2))
        cont[x[0]][x[1]] = intersect

    cont = cont.as_matrix()
    return cont


def entropy(cluster_labels):
    '''
    Calculate the entropy of a clustering.
    Entropy(X):
    H(X) = sum(p(i)log(1/(pi)))
    '''

    if len(cluster_labels) == 0:
        return 1.0
    else:
        pass

    cluster_prob = [len(cluster_labels[x]) for x in cluster_labels.keys()]
    pi = np.array(cluster_prob).astype(np.float)
    pi = pi[pi > 0]
    pi_sum = np.sum(pi)

    entropy = -np.sum((pi / pi_sum) * (np.log(pi) - log(pi_sum)))

    return entropy


def adjustedMutualInformation(cluster1, cluster2):
    '''
    Using the scikit-learn algorithms, calculate the adjusted mutual
    information for two clusterings.  Assume cluster1 is the
    reference/ground truth clustering.
    The adjusted MI accounts for higher scores by chance, particularly
    in the case where a larger number of clusters leads to a higher MI.

    AMI(U, V) = [MI(U, V) - E(MI(U, V))] / [max(H(U), H(V)) - E(MI(U, V))]
    where E(MI(U, V)) is the expected mutual information given the
    number of clusters and H(U), H(V) are the entropies of clusterings
    U and V.
    '''

    cont = contingency(cluster1, cluster2)
    mi = mutualInformation(cluster1, cluster2)
    sample_size = float(sum([len(cluster1[x]) for x in cluster1.keys()]))

    # Given the number of samples, what is the expected number
    # of overlaps that would occur by chance?
    emi = supervised.expected_mutual_information(cont, sample_size)

    # calculate the entropy for each clustering
    h_clust1, h_clust2 = entropy(cluster1), entropy(cluster2)

    if abs(h_clust1) == 0.0:
        h_clust1 = 0.0
    else:
        pass
    if abs(h_clust2) == 0.0:
        h_clust2 = 0.0
    else:
        pass

    ami = (mi - emi) / (max(h_clust1, h_clust2) - emi)

    # bug: entropy will return -0 in some instances
    # make sure this is actually 0 else ami will return None
    # instead of 0.0

    if np.isnan(ami):
        ami = np.nan_to_num(ami)
    else:
        pass

    return ami




def deseqNormalize(infile,
                   time_points,
                   reps,
                   conditions=None):
    '''
    Library size normalisation and variance stabilizing transformation of
    timeseries RNA-seq data

    :param infile: count table from NGS-seq experiment
    :type infile: str
    :param time_points: time point labels
    :type time_points: str list
    :param reps: replicates labels
    :type reps: str list
    :param conditions: if  multiple experimental conditions
    are to be normalised at the same time
    :type conditions: str list
    '''
    # MM: NB - this should be split into separate library size
    # normalisation and VST transformations
    # maybe add in different transformation options.

    pandas2ri.activate()
    reps = reps

    # load library
    R('''suppressMessages(library("DESeq"))''')

    # generates a lists for the design data frame
    # of the proper length
    # these need to be rpy2 objects to be parsed
    # properly in the string formatting

    E.info("converting to pandas dataframe object")

    if infile.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    data_frame = pd.read_table(infile,
                               index_col=0,
                               header=0,
                               sep="\t",
                               compression=comp)
    # py2ri requires activation
    pandas2ri.activate()
    rdf = pandas2ri.py2ri(data_frame)

    if not conditions:
        time_rep_comb = [x for x in itertools.product(time_points, reps)]
        time_cond = ro.StrVector([x[0] for x in time_rep_comb])
        rep_cond = ro.StrVector([x[1] for x in time_rep_comb])

        R.assign('countsTable', rdf)
        R('''design <- data.frame(row.names=colnames(countsTable),'''
          '''times=%s, replicates=%s)''' % (time_cond.r_repr(),
                                            rep_cond.r_repr()))
    elif conditions:
        design_dict = {}
        for x in data_frame.columns.values:
            sample_dict = {}
            sample_dict['condition'] = str(x).split(".")[0]
            sample_dict['times'] = int(str(x).split(".")[1])
            sample_dict['replicates'] = str(x).split(".")[2]
            design_dict[x] = sample_dict
            design_frame = pd.DataFrame(design_dict)
            design_frame = design_frame.T

        des_cond = design_frame['condition'].values.tolist()
        des_time = design_frame['times'].values.tolist()
        des_reps = design_frame['replicates'].values.tolist()

        cond_cond = ro.StrVector([x for x in des_cond])
        time_cond = ro.StrVector([x for x in des_time])
        rep_cond = ro.StrVector([x for x in des_reps])

        R.assign('countsTable', rdf)
        R.assign('design', design_frame)

    # create the count data set and normalize to library size
    # transform with variance stabilizing transformation
    # only select genes with an average of ten reads mapping

    E.info("calculating size factors and dispersion")
    R('''notZero <- (rowMeans(countsTable) > 1)''')
    R('''cds <- newCountDataSet(countsTable[notZero, ], design)''')
    R('''cds_size <- estimateSizeFactors(cds)''')
    R('''cds_disp <- estimateDispersions(cds_size, method="blind")''')

    E.info("applying variance stabilizing transformation")

    R('''vst <- varianceStabilizingTransformation(cds_disp)''')

    # format data set to long format with condition and replicate labels
    # convert to a numpy array

    R('''replicates <- c(%s)''' % rep_cond.r_repr())
    R('''times <- c(%s)''' % time_cond.r_repr())
    if conditions:
        R('''conditions <- c(%s)''' % cond_cond.r_repr())
        R('''trans_vst = data.frame(t(exprs(vst)), '''
          '''times, replicates, conditions)''')
    else:
        R('''trans_vst = data.frame(t(exprs(vst)), times, replicates)''')

    # load data and convert to pandas object
    data_file = pandas2ri.ri2py(R["trans_vst"])
    
    return data_file


def avTimeExpression(infile):
    '''
    Calculate average expression over replicates at each time point
    Requires genes as columns with 'replicates' and 'times' as additional
    columns
    '''

    # check file compression
    if infile.split(".")[-1] == "gz":
        comp = "gzip"
    else:
        comp = None

    df = pd.read_table(infile, sep="\t",
                       header=0, index_col=0,
                       compression=comp)

    # average over replicates at each time point
    # then recombine
    df_groups = df.groupby(by='times')
    data_frame = pd.DataFrame(index=df.columns,
                              columns=None)
    for names, groups in df_groups:
        _df = groups.drop(['times', 'replicates'], axis=1)
        _df = _df.apply(np.mean, axis=0)
        data_frame[names] = _df

    # check no extraneous columns remaining
    try:
        data_frame.drop(['replicates', 'times'],
                        inplace=True,
                        axis=0)
    except KeyError:
        pass

    return data_frame


def covarFilter(infile,
                time_points,
                replicates,
                quantile):
    '''
    Filter gene list based on the distribution of the
    sums of the covariance of each gene.  This is highly
    recommended to reduce the total number of genes used
    in the dynamic time warping clustering to reduce the
    computational time.  The threshold is placed at the
    intersection of the expected and observed value
    for the given quantile.
    '''

    time_points.sort()
    time_rep_comb = [x for x in itertools.product(time_points, replicates)]
    time_cond = ro.StrVector([x[0] for x in time_rep_comb])
    rep_cond = ro.StrVector([x[1] for x in time_rep_comb])
    df = pd.read_table(infile, sep="\t", header=0, index_col=0)

    df.drop(['replicates'], inplace=True, axis=1)
    df.drop(['times'], inplace=True, axis=1)
    df = df.fillna(0.0)

    # convert data frame and import into R namespace
    # py2ri requires activation
    pandas2ri.activate()
    R.assign('diff_data', pandas2ri.py2ri(df))

    E.info("loading data frame")

    # need to be careful about column headers and transposing data frames

    R('''trans_data <- data.frame(diff_data)''')
    R('''times <- c(%s)''' % time_cond.r_repr())
    R('''replicates <- c(%s)''' % rep_cond.r_repr())

    # calculate the covariance matrix for all genes
    # sum each gene's covariance vector

    E.info("calculating sum of covariance of expression")

    R('''covar.mat <- abs(cov(trans_data))''')
    R('''sum.covar <- rowSums(covar.mat)''')
    R('''exp.covar <- abs(qnorm(ppoints(sum.covar),'''
      '''mean=mean(sum.covar), sd=sd(sum.covar)))''')
    R('''sum.covar.quant <- quantile(sum.covar)''')
    R('''exp.covar.quant <- quantile(exp.covar)''')

    E.info("filter on quantile")

    R('''filtered_genes <- names(sum.covar[sum.covar > '''
      '''sum.covar.quant[%(quantile)i]'''
      ''' & sum.covar > exp.covar.quant[%(quantile)i]])''' % locals())
    R('''filtered_frame <- data.frame(diff_data[, filtered_genes],'''
      '''times, replicates)''')

    # load data and convert to pandas object
    filtered_frame = pandas2ri.ri2py(R["filtered_frame"]).T

    return filtered_frame


def clusterPCA(infile,
               cluster_file,
               image_dir):
    '''
    PCA for each module within an experimental condition across
    the time series.
    Take PC1 as the module eigengene and return the loadings and proportion
    of variance explained for the eigengene.
    The eigengene expression across the time series is taken to be the
    value for PC1 at each timepoint as a vector.
    This is basically what WGCNA moduleEigengenes does but it
    does not recover the PC loadings.

    Warning: this script will error if there is only one
    cluster. Make sure you have more than one cluster before
    trying to perform uninformative analyses.

    Parameters:

    rpath : string
        Path of R support libraries

    '''

    header = cluster_file.split("/")[-1].split("-")[0]

    # reshape data
    R('''sink(file='sink_file.txt')''')
    R('''suppressMessages(library("reshape2"))''')
    R('''suppressMessages(library("WGCNA"))''')
    # AH: these were hard-coded paths, parameterized them to point to the
    # directory of this module's location
    R('''source("%s")''' % os.path.join(get_r_path(), "summarySE.R"))
    R('''source("%s")''' % os.path.join(get_r_path(), "clusterEigengenes.R"))
    R('''cluster_match <- read.table('%(cluster_file)s', h=T, '''
      '''row.names=1)''' % locals())
    R('''express_data <- read.table('%(infile)s', '''
      '''h=T, row.names=1, stringsAsFactors=F)''' % locals())
    R('''sink(file=NULL)''')
    R('''colnames(cluster_match) <- c("genes", "cluster")''')
    R('''express_data <- data.frame(t(express_data))''')
    R('''express_data$times <- as.numeric(as.character(express_data$times))''')
    R('''data_melt <- melt(express_data, '''
      '''id.vars=c("times", "replicates"))''')

    # sometimes data is read in as a factor/string.
    # Explicitly convert to numeric

    R('''data_melt$value <- as.numeric(as.character(data_melt$value))''')
    R('''data_sum <- summarySE(data_melt, measurevar="value", '''
      '''groupvars=c("times", "variable"))''')
    R('''data_mod <- data.frame(data_sum$times,'''
      ''' data_sum$variable, data_sum$value)''')
    R('''colnames(data_mod) <- c("times", "gene", "value")''')
    R('''data_wide <- dcast(data_mod, gene ~ times, value.var="value")''')
    R('''rownames(data_wide) <- data_wide$gene''')
    R('''times <- as.numeric(as.character(unique(express_data$times)))''')
    R('''data_wide <- data.frame(data_wide[,-1])''')
    R('''colnames(data_wide) <- times''')

    # derive module eigengenes - return a dataframe of eigengene expression
    R('''eigen_clustered <- clusterPCA(cluster_frame=cluster_match, '''
      '''expression_frame=data_wide, n=times)''')
    R('''eigen_frame <- eigenExpress(eigen_clustered, n=times)''')

    # generate loadings plot for each eigengene
    R('''eigenLoad(clusterPCA(cluster_frame=cluster_match, '''
      '''expression_frame=data_wide, n=times), image.dir="%(image_dir)s", '''
      '''condition="%(header)s")''' % locals())

    # generate expression profile plots for all eigengenes
    R('''eigenPlot(eigen_frame, image.dir="%(image_dir)s", '''
      '''condition="%(header)s")''' % locals())

    eigen_frame = pandas2ri.ri2py(R["eigen_frame"])
    eigen_frame.index = eigen_frame['cluster']
    eigen_frame.drop(['cluster'], inplace=True, axis=1)

    return eigen_frame




def conditionDESeq2(data_frame, header, alpha, res_dir):
    '''
    Perform DESeq2-based analysis of condition:time interaction
    dependent differential expression
    '''

    E.info("Differential expression testing for %s" % header)
    cols = data_frame.columns

    # py2ri requires activation
    pandas2ri.activate()
    counts = pandas2ri.py2ri(data_frame)

    des_times = ro.IntVector([x.split(".")[1] for x in cols])
    des_reps = ro.StrVector([x.split(".")[2] for x in cols])
    des_cond = ro.StrVector([x.split(".")[0] for x in cols])
    genes = ro.StrVector([x for x in data_frame.index])

    # setup counts table and design frame

    R('''suppressPackageStartupMessages(library("DESeq2"))''')
    R('''sink(file="/dev/null")''')
    R('''times <- as.factor(%s)''' % des_times.r_repr())
    R('''reps <- c(%s)''' % des_reps.r_repr())
    R('''condition <- c(%s)''' % des_cond.r_repr())
    R('''design <- data.frame(times, reps, condition)''')
    R('''counts <- data.frame(%s)''' % counts.r_repr())
    R('''genes <- c(%s)''' % genes.r_repr())
    R('''rownames(counts) <- genes''')
    R('''rownames(design) <- colnames(counts)''')

    # use DESeq() with LRT and reduced formula.  Use effect
    # size moderation

    R('''dds <- DESeqDataSetFromMatrix(countData=counts, '''
      '''colData=design, '''
      '''design=~reps + times + condition + times:condition)''')
    R('''dds <- DESeq(dds, test="LRT", '''
      '''reduced=~reps + times + condition, betaPrior=T)''')
    R('''res <- results(dds)[order(results(dds)$padj, na.last=T), ]''')
    R('''res.df <- data.frame(res)''')

    # generate dispersion and MA plots
    R('''png("%s/%s-dispersions.png")''' % (res_dir,
                                            header))
    R('''plotDispEsts(dds)''')
    R('''dev.off()''')

    R('''png("%s/%s-MAplot.png")''' % (res_dir,
                                       header))
    R('''plotMA(res, alpha=%0.3f, ylim=c(-5,5))''' % alpha)
    R('''dev.off()''')
    R('''sink(file=NULL)''')

    df = pandas2ri.ri2py(R['res.df'])

    return df


def timepointDESeq2(data_frame, header, alpha, res_dir):
    '''
    Perform DESeq2-based differential expression analysis of condition:time
    interaction
    '''

    E.info("Differential expression testing for %s" % header)
    cols = data_frame.columns

    # py2ri requires activation
    pandas2ri.activate()
    counts = pandas2ri.py2ri(data_frame)

    des_times = ro.IntVector([x.split(".")[1] for x in cols])
    des_reps = ro.StrVector([x.split(".")[2] for x in cols])
    genes = ro.StrVector([x for x in data_frame.index])

    # setup counts table and design frame

    R('''suppressPackageStartupMessages(library("DESeq2"))''')
    R('''sink(file="/dev/null")''')
    R('''times <- as.factor(%s)''' % des_times.r_repr())
    R('''reps <- c(%s)''' % des_reps.r_repr())
    R('''design <- data.frame(times, reps)''')
    R('''counts <- data.frame(%s)''' % counts.r_repr())
    R('''genes <- c(%s)''' % genes.r_repr())
    R('''rownames(counts) <- genes''')
    R('''rownames(design) <- colnames(counts)''')

    # use DESeq() with LRT and reduced formula.  Use effect
    # size moderation

    R('''dds <- DESeqDataSetFromMatrix(countData=counts, '''
      '''colData=design, '''
      '''design=~reps + times )''')
    R('''dds <- DESeq(dds, betaPrior=T)''')
    R('''res <- results(dds)[order(results(dds)$padj, na.last=T), ]''')
    R('''res.df <- data.frame(res)''')

    # generate dispersion and MA plots
    R('''png("%s/%s-dispersions.png")''' % (res_dir,
                                            header))
    R('''plotDispEsts(dds)''')
    R('''dev.off()''')

    R('''png("%s/%s-MAplot.png")''' % (res_dir,
                                       header))
    R('''plotMA(res, alpha=%0.3f, ylim=c(-5,5))''' % alpha)
    R('''dev.off()''')
    R('''sink(file=NULL)''')

    df = pandas2ri.ri2py(R['res.df'])

    return df


def genSigGenes(file_list, alpha, out_dir):
    '''
    With a list of input files (results from DESeq2) generate a
    dictionary of genes that are statistically significantly differentially
    expressed between conditions and/or time points.
    '''

    alpha = float(alpha)
    deg_dict = {}
    for infle in file_list:
        if infle.split("/")[0].split(".")[0] == "diff_condition":
            header = infle.split("/")[1].split(".")[1]
            header = header.rstrip("-diff-cond.tsv")
            header = "%s" % header

        elif infle.split("/")[0].split(".")[0] == "diff_timepoints":
            header = infle.split("/")[-1].split("-")[0]
            header = "%s_%s" % (header.split("_")[0],
                                header.split("_")[-1])

        in_df = pd.read_table(infle, sep="\t", header=0, index_col=0)
        sig_genes = in_df[in_df['padj'] <= alpha]
        deg_dict[header] = sig_genes.index.tolist()

    if file_list[0].split("/")[0].split(".")[0] == "diff_condition":
        condition = file_list[0].split("/")[1].split(".")[0]
        condition = "%s-condition" % condition

    elif file_list[0].split("/")[0].split(".")[0] == "diff_timepoints":
        condition = file_list[0].split("/")[-1].split(".")[0]
        condition = "%s-time" % condition

    drawVennDiagram(deg_dict, condition, out_dir)


def drawVennDiagram(deg_dict, header, out_dir):
    '''
    Take a dictionary of gene IDs, with keys corresponding
    to timepoints/differential expression analyses and
    generate a Venn diagram.  Maximum of 5 overlapping sets
    possible using R package:: VennDiagram.
    '''

    keys = deg_dict.keys()
    try:
        keys = sorted(keys, key=lambda x: int(x.split("_")[1].rstrip("-time")))
    except IndexError:
        pass

    venn_size = len(keys)
    R('''suppressPackageStartupMessages(library("VennDiagram"))''')
    n1 = set(deg_dict[keys[0]])
    n2 = set(deg_dict[keys[1]])
    area1 = len(n1)
    area2 = len(n2)
    n12 = len(n1.intersection(n2))

    # for Venn > 2 sets
    if venn_size == 3:
        n3 = set(deg_dict[keys[2]])
        area3 = len(n3)
        n13 = len(n1.intersection(n3))
        n23 = len(n2.intersection(n3))
        n123 = len((n1.intersection(n2)).intersection(n3))
        cat1, cat2, cat3 = keys

        R('''png("%(out_dir)s/%(header)s-venn.png", '''
          '''width=480, height=480)''' % locals())
        R('''draw.triple.venn(%(area1)d, %(area2)d, %(area3)d, '''
          '''%(n12)d, %(n23)d, %(n13)d, %(n123)d, '''
          '''c('%(cat1)s', '%(cat2)s', '%(cat3)s'), '''
          '''col=c('red', 'yellow', 'skyblue'), '''
          '''fill=c('red', 'yellow', 'skyblue'), '''
          '''margin=0.05, alpha=0.5)''' % locals())
        R('''dev.off()''')

    elif venn_size == 4:
        n3 = set(deg_dict[keys[2]])
        area3 = len(n3)
        n13 = len(n1.intersection(n3))
        n23 = len(n2.intersection(n3))
        n123 = len((n1.intersection(n2)).intersection(n3))

        n4 = set(deg_dict[keys[3]])
        area4 = len(n4)
        n14 = len(n1.intersection(n4))
        n24 = len(n2.intersection(n4))
        n34 = len(n3.intersection(n4))
        n124 = len((n1.intersection(n2)).intersection(n4))
        n134 = len((n1.intersection(n3)).intersection(n4))
        n234 = len((n2.intersection(n3)).intersection(n4))
        n1234 = len(((n1.intersection(n2)).intersection(n3)).intersection(n4))
        cat1, cat2, cat3, cat4 = keys

        R('''png("%(out_dir)s/%(header)s-venn.png",'''
          '''width=480, height=480)''' % locals())
        R('''draw.quad.venn(%(area1)d, %(area2)d, %(area3)d, %(area4)d,'''
          '''%(n12)d, %(n13)d, %(n14)d, %(n23)d, %(n24)d, %(n34)d,'''
          '''%(n123)d, %(n124)d, %(n134)d, %(n234)d, %(n1234)d,'''
          '''c('%(cat1)s', '%(cat2)s', '%(cat3)s', '%(cat4)s'), '''
          '''col=c("red", "yellow", "skyblue", "orange"), '''
          '''fill=c("red", "yellow", "skyblue", "orange"), '''
          '''margin=0.05, alpha=0.5)''' % locals())
        R('''dev.off()''')

    elif venn_size == 5:
        n3 = set(deg_dict[keys[2]])
        area3 = len(n3)
        n13 = len(n1.intersection(n3))
        n23 = len(n2.intersection(n3))
        n123 = len((n1.intersection(n2)).intersection(n3))

        n4 = set(deg_dict[keys[3]])
        area4 = len(n4)
        n14 = len(n1.intersection(n4))
        n24 = len(n2.intersection(n4))
        n34 = len(n3.intersection(n4))
        n124 = len((n1.intersection(n2)).intersection(n4))
        n134 = len((n1.intersection(n3)).intersection(n4))
        n234 = len((n2.intersection(n3)).intersection(n4))
        n1234 = len(((n1.intersection(n2)).intersection(n3)).intersection(n4))

        n5 = set(deg_dict[keys[4]])
        area5 = len(n5)
        n15 = len(n1.intersection(n5))
        n25 = len(n2.intersection(n5))
        n35 = len(n3.intersection(n5))
        n45 = len(n4.intersection(n5))
        n125 = len((n1.intersection(n2)).intersection(n5))
        n135 = len((n1.intersection(n3)).intersection(n5))
        n145 = len((n1.intersection(n4)).intersection(n5))
        n235 = len((n2.intersection(n3)).intersection(n5))
        n245 = len((n2.intersection(n4)).intersection(n5))
        n345 = len((n3.intersection(n4)).intersection(n5))
        n1235 = len(((n1.intersection(n2)).intersection(n3)).intersection(n5))
        n1245 = len(((n1.intersection(n2)).intersection(n4)).intersection(n5))
        n1345 = len(((n1.intersection(n3)).intersection(n4)).intersection(n5))
        n2345 = len(((n2.intersection(n3)).intersection(n4)).intersection(n5))
        nstep = ((n1.intersection(n2)).intersection(n3))
        n12345 = len((nstep.intersection(n4)).intersection(n5))
        cat1, cat2, cat3, cat4, cat5 = keys

        R('''png("%(out_dir)s/%(header)s-venn.png", '''
          '''height=480, width=480)''' % locals())
        R('''draw.quintuple.venn(%(area1)d, %(area2)d, %(area3)d, '''
          '''%(area4)d, %(area5)d, %(n12)d, %(n13)d, %(n14)d,'''
          '''%(n15)d, %(n23)d, %(n24)d, %(n25)d, %(n34)d, %(n35)d,'''
          '''%(n45)d, %(n123)d, %(n124)d, %(n125)d, %(n134)d,'''
          '''%(n135)d, %(n145)d, %(n234)d, %(n235)d, %(n245)d,'''
          '''%(n345)d, %(n1234)d, %(n1235)d, %(n1245)d, %(n1345)d,'''
          '''%(n2345)d, %(n12345)d, '''
          '''c('%(cat1)s', '%(cat2)s', '%(cat3)s', '%(cat4)s', '%(cat5)s'),'''
          '''col=c("red", "yellow", "skyblue", "orange", "purple"),'''
          '''fill=c("red", "yellow", "skyblue", "orange", "purple"),'''
          '''alpha=0.05, margin=0.05, cex=rep(0.8, 31))''' % locals())
        R('''dev.off()''')

    elif venn_size > 5:
        raise ValueError("Illegal venn diagram size, must be <= 5")


def maSigPro(infile,
             order_terms=1,
             fdr=0.01,
             adjust="BH",
             stepwise="backward",
             include_p=0.01,
             rsq=0.2,
             var_group="all"):
    '''
    Generate differentially expressed genes for each experimental
    condition across a time series.  Uses the bioconductor
    package maSigPro to derive a set of genes of interest.
    '''

    ref_gtf = str(infile).split("-")[1]
    data_frame = pd.read_table(infile, sep="\t", index_col=0, header=0)
    design_dict = {}

    for x in data_frame.index.values:
        sample_dict = {}
        condition = str(x).split(".")[0]
        sample_dict[condition] = 1
        sample_dict['times'] = int(str(x).split(".")[1])
        sample_dict['replicates'] = str(x).split(".")[2]
        design_dict[x] = sample_dict

    design_frame = pd.DataFrame(design_dict)
    design_frame = design_frame.T
    cols = ['times', 'replicates', condition]
    design_frame = design_frame[cols]
    design_file = "deseq.dir/%s-%s-design.tsv" % (condition, ref_gtf)
    design_frame.to_csv(design_file, sep="\t")
    data_file = "deseq.dir/%s-%s-data.tsv" % (condition, ref_gtf)
    results_file = "deseq.dir/%s-%s-maSigPro.tsv" % (condition, ref_gtf)

    # data frame columns must be in the order time-replicate-condition
    # for maSigPro
    # define the numnber of higher-order terms included in the models

    masigpro_out = "deseq.dir/maSigPro.out"

    R('''suppressMessages(library("maSigPro"))''')
    R('''input_data <- read.table('%(infile)s', sep="\t", '''
      '''h=T, row.names=1)''' % locals())
    R('''input_data <- t(input_data[0:(length(input_data)-2)])''')

    E.info("constructing experimental design matrix")

    R('''input_design <- data.matrix(read.table('%(design_file)s', '''
      '''sep="\t", h=T, row.names=1))''' % locals())
    R('''%(condition)s_mat <- make.design.matrix(input_design, '''
      '''degree = %(order_terms)i )''' % locals())
    R('''sink(file = '%(masigpro_out)s')''' % locals())

    E.info("fitting linear model for each gene with "
           "%i polynomial terms" % order_terms)

    R('''%(condition)s_fit <- p.vector(input_data, %(condition)s_mat, '''
      '''Q = %(fdr)f, MT.adjust = '%(adjust)s')''' % locals())

    # fit a linear model to each of the genes called as
    # differentially expressed
    # report genes with model R-squared > threshold
    # maSigPro gives an un-suppressable output to stdout
    # therefore sink is used to shunt this to a temporary file 'maSigPro.out'

    R('''%(condition)s_step <- T.fit(%(condition)s_fit, '''
      '''step.method='%(stepwise)s', alfa=%(include_p)f)''' % locals())

    E.info("selecting significantly differentially "
           "expressed genes at FDR=%0.3f" % fdr)

    R('''sink(file=NULL)''')
    R('''%(condition)s_sigs <- get.siggenes(%(condition)s_step, '''
      '''rsq=%(rsq)f, vars='%(var_group)s')''' % locals())
    R('''write.table(%(condition)s_sigs$sig.genes$%(condition)s$group.coeffs'''
      ''',file="deseq.dir/%(condition)s-%(ref_gtf)s-coefficients.tsv", '''
      '''sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$sig.genes$%(condition)s$sig.pvalues,'''
      '''file="deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv",'''
      ''' sep="\t")''' % locals())
    R('''write.table(%(condition)s_sigs$summary, '''
      '''file='deseq.dir/%(condition)s-%(ref_gtf)s-geneids.tsv', '''
      '''sep="\t")''' % locals())
    # merge the p-value and coefficient results into a single file
    p_file = "deseq.dir/%(condition)s-%(ref_gtf)s-pvalues.tsv" % locals()
    coef_file = "deseq.dir/%s-%s-coefficients.tsv" % (condition,
                                                      ref_gtf)
    p_frame = pd.read_table(p_file, sep="\t")
    coef_frame = pd.read_table(coef_file, sep="\t")
    results_frame = pd.merge(coef_frame, p_frame,
                             how='right',
                             left_index=True,
                             right_index=True)

    results_frame.to_csv(results_file, sep="\t")

    R('''diff_genes <- data.frame(%(condition)s_fit$SELEC)''' % locals())
    diff_genes = pandas2ri.ri2py[R'diff_genes']

    return diff_genes




def splitReplicates(infile,
                    axis,
                    group_var,
                    outdir):
    '''
    Group data by replicates, new flat file for
    each.
    '''

    if axis == "row":
        axis = 1
    elif axis == "column":
        axis = 0

    inf_prefix = infile.split("/")[1].split("-")
    inf_prefix = inf_prefix[0] + "-" + inf_prefix[1]

    df = pd.read_table(infile,
                       sep="\t",
                       header=0,
                       index_col=0).T

    rep_groups = df.groupby(by=group_var,
                            axis=axis)
    for name, groups in rep_groups:
        outfile = outdir + "/" + inf_prefix + "-%s-expression.tsv" % name
        _df = groups.T
        _df.columns = _df.loc['times']
        _df.drop(['times'], axis=axis, inplace=True)
        _df.drop(['replicates'], axis=axis, inplace=True)
        _df.to_csv(outfile, sep="\t", index_label="gene_id")


def genResampleData(data_frame,
                    multiple_index,
                    replicates,
                    sample_reps,
                    times,
                    condition,
                    ref_gtf,
                    out_dir,
                    seed):
    '''
    Resample the data n-times with replacement - store
    in an sql database for ease of access
    '''

    vst_long = data_frame.T
    vst_long.index = multiple_index

    reps_dict = {}
    random.seed(seed)
    for it in range(1, replicates + 1):
        df = pd.DataFrame()
        df_list = []
        for i in times:
            k = str(random.randint(1,
                                   len(sample_reps)))
            series = vst_long.loc[str(i), 'R%s' % k]
            df_list.append((str(i), series))
        df = pd.DataFrame.from_items(df_list)
        cols = df.columns.tolist()
        cols = [int(x) for x in cols]
        cols.sort()
        cols = [str(x) for x in cols]
        df = df[cols]
        reps_dict[str(it)] = df
        table = "%s-%s-resample_%i-expression" % (condition,
                                                  ref_gtf,
                                                  it)
        seg_file = "%s/%s.tsv" % (out_dir, table)
        df.to_csv(seg_file, sep="\t")

    E.info("%i replicate datasets generated" % replicates)


def temporalCorrelate(series1, series2):
    '''
    Calculate the temporal correlation according to Chouakira & Nagabhushan
    Assumes both time series are of the same length
    '''

    series1 = list(series1)
    series2 = list(series2)

    sum_prod = []
    sum_usq = []
    sum_vsq = []
    for i in range(len(series1)-1):
        u = float(series1[i+1]) - float(series1[i])
        v = float(series2[i+1]) - float(series2[i])
        prod = u * v
        sum_prod.append(prod)
        sq_u = u**2
        sq_v = v**2
        sum_usq.append(sq_u)
        sum_vsq.append(sq_v)

    nume = sum(sum_prod)
    denom = math.sqrt(sum(sum_usq)) * math.sqrt(sum(sum_vsq))

    if denom != 0:
        return(nume/float(denom))
    else:
        return 0


def crossCorrelate(t, s, lag=0):
    '''
    Calculate the cross-correlation of two timeseries, s and t.
    Return the normalized correlation value at lag=n.
    Uses numpy.correlate; default is to return lag=0.
    TODO: return multiple lags?
    '''

    t_mean = np.mean(t)
    s_mean = np.mean(s)
    t_std = np.std(t)
    s_std = np.std(s)
    len_t = len(t)

    t_norm = [((x - t_mean)/(t_std * len_t)) for x in t]
    s_norm = [((y - s_mean)/s_std) for y in s]

    if lag == 0:
        xcorr = np.correlate(t_norm, s_norm)
    elif lag != 0:
        xcorr = np.correlate(t_norm, s_norm, mode=2)[len_t - 1 + lag]

    return xcorr


def adaptiveTune(value, k):
    '''
    Calculate the adaptive tuning function from Chouakira & Nagabhushan
    '''

    if k == 0:
        return 1.0
    else:
        return (2/(1 + math.exp(k*abs(value))))


def dtwWrapper(data, rows, columns, k):
    '''
    wrapper function for dynamic time warping.
    includes use of exponential adaptive tuning function
    with temporal correlation if k > 0
    '''

    # not explicitly called, but needs to be in R environment
    DTW = importr("dtw")

    # create a data frame of zeros of size number of ids x number of ids
    # fill it with the calculated distance metric for each pair wise comparison

    df_ = pd.DataFrame(index=rows,
                       columns=columns)
    df_ = df_.fillna(0.0).astype(np.float64)

    # fill the array with dtw-distance values
    pandas2ri.activate()

    for i in rows:
        E.info("DTW %s" % i)
        for j in columns:
            series1 = data.loc[i].values.tolist()
            series2 = data.loc[j].values.tolist()
            DTW_value = (R.dtw(series1,
                               series2)).rx('distance')[0][0]
            cort_value = temporalCorrelate(series1, series2)
            tuned_value = adaptiveTune(cort_value, k)
            time_dist = DTW_value * tuned_value
            df_.loc[i][j] = float(time_dist)
            df_[j][i] = float(time_dist)

    return df_


def correlateDistanceMetric(data, rows, columns, method, lag=0):
    '''
    wrapper for correlation coefficients as distance metrics
    for time-series clustering.
    Use either temporal correlation (analagous to template matching)
    or normalised cross correlation.
    '''

    # create blank (all 0's) dataframe to fill with correlation values

    df_ = pd.DataFrame(index=rows,
                       columns=columns)
    df_ = df_.fillna(0.0)

    if method == "cross-correlate":
        for i in rows:
            E.info("cross-correlation %s" % i)
            for j in columns:
                series1 = data.loc[i].values.tolist()
                series2 = data.loc[j].values.tolist()
                corr = crossCorrelate(series1, series2, lag=lag)
                df_.loc[i][j] = 1.0 - abs(corr)
                df_[j][i] = 1.0 - abs(corr)

    elif method == "temporal-correlate":
        for i in rows:
            E.info("temporal correlation %s" % i)
            for j in columns:
                series1 = data.loc[i].tolist()
                series2 = data.loc[j].tolist()
                corr = temporalCorrelate(series1, series2)
                df_.loc[i][j] = 1.0 - abs(corr)
                df_[j][i] = 1.0 - abs(corr)

    return df_


def splitFiles(infile, nchunks, out_dir):
    '''
    Give files names based on splitting into an arbitrary number of chunks
    '''

    df = pd.read_table(infile, sep="\t", header=0, index_col=0)
    total = len(df.index.tolist())

    # split into aribitrary number of chunks, or arbitrary chunk size?
    # small n bad for large input size, large n bad for small input size
    # set min/max chunk size, e.g. 100 genes minimum, 500 maximum?

    if total/nchunks < 100:
        step = 100
        E.warn("too few genes in each chunk, resetting to 100 genes per chunk")
    elif total/nchunks > 500:
        step = 500
        E.warn("too many genes per chunk, resetting to 500 genes per chunk")
    else:
        step = total/nchunks
        E.info("chunking input file into %i chunks" % step)

    file_pattern = infile.split("/")[1].rstrip("-expression.tsv")
    idx = 0
    for i in range(step, total, step):
        start = "%s" % idx
        end = "%s" % i
        file_name = "%s/%s-%s_%s-split.tsv" % (out_dir,
                                               file_pattern,
                                               start,
                                               end)
        with open(file_name, "w") as file_handle:
            file_handle.write(file_name + "\n")
        idx = i

    # final file
    start = "%s" % idx
    end = "%s" % total
    file_name = "%s/%s-%s_%s-split.tsv" % (out_dir,
                                           file_pattern,
                                           start,
                                           end)
    with open(file_name, "w") as file_handle:
        file_handle.write(file_name + "\n")


def mergeFiles(file_list, outfile):
    '''
    Merge files after split-transform operations
    '''

    # sort list by starting index of file
    res_list = sorted(file_list,
                      key=lambda x: int(x.split("/")[-1].split("-")[3].split("_")[0]))

    # merge files using pandas data frame merge method
    full_frame = pd.read_table(res_list[0], sep="\t", index_col=0, header=0)
    res_list.remove(res_list[0])
    for fle in res_list:
        df = pd.read_table(fle, sep="\t", index_col=0, header=0)
        full_frame = pd.merge(left=full_frame,
                              right=df,
                              how='inner',
                              left_index=True,
                              right_index=True)

    full_frame.to_csv(outfile, sep="\t")


def treeCutting(infile,
                expression_file,
                cluster_file,
                cluster_algorithm,
                deepsplit=False):
    '''
    Use dynamic tree cutting to derive clusters for each
    resampled distance matrix
    '''
    wgcna_out = "/dev/null"

    E.info("loading distance matrix")

    df = pd.read_table(infile, sep="\t",
                       header=0, index_col=0)
    df = df.fillna(0.0)
    genes = df.index
    genes_r = ro.StrVector([g for g in genes])

    # py2ri requires activation
    pandas2ri.activate()
    rdf = pandas2ri.py2ri(df)

    R.assign("distance_data", rdf)
    R.assign("gene_ids", genes_r)

    R('''sink(file='%(wgcna_out)s')''' % locals())
    R('''suppressPackageStartupMessages(library("WGCNA"))''')
    R('''suppressPackageStartupMessages(library("flashClust"))''')
    E.info("clustering data by %s linkage" % cluster_algorithm)
    R('''rownames(distance_data) <- gene_ids''')
    R('''clustering <- flashClust(as.dist(distance_data),'''
      ''' method='%(cluster_algorithm)s')''' % locals())
    if deepsplit:
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''minClusterSize=50, deepSplit=T)''')
    else:
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''minClusterSize=50, deepSplit=F)''')

    R('''color_cut <- labels2colors(cluster_cut)''')
    R('''write.table(color_cut, file = '%(cluster_file)s','''
      '''sep="\t")''' % locals())
    R('''cluster_matched <- data.frame(cbind(rownames(distance_data),'''
      '''color_cut))''')
    R('''colnames(cluster_matched) = c("gene_id", "cluster")''')
    R('''cluster_matched <- data.frame(cluster_matched$gene_id,'''
      '''cluster_matched$cluster)''')
    R('''sink(file=NULL)''')

    cluster_frame = pandas2ri.ri2py(R["cluster_matched"])
    cluster_frame.columns = ['gene_id', 'cluster']
    cluster_frame.index = cluster_frame['gene_id']
    cluster_frame.drop(['gene_id'], inplace=True, axis=1)

    return cluster_frame


def clusterAverage(file_list):
    '''
    Average distance measures across replicates
    '''

    # map replicate number on to dataframe as identifier
    # used to aggregate over replicates for each gene id
    # assumes filename is condition-reference-replicate

    df_dict = {}
    for fle in file_list:
        f = fle.split("/")[-1]
        rep = f.split("-")[2]
        _df = pd.read_table(fle, sep="\t",
                            header=0, index_col=0)
        df_dict[rep] = _df

    # group concatenated dataframes by gene_id
    concat_df = pd.concat(df_dict)
    group_df = concat_df.groupby(level=1)

    # aggregate/summarise distance metric over replicates for each gene_id
    agg_dict = {}
    for names, groups in group_df:
        agg_dict[names] = np.mean(groups, axis=0)

    agg_df = pd.DataFrame(agg_dict)
    return agg_df


def clusterAgreement(infile):
    '''
    calculate co-occurence of genes within resampled clusters
    '''
    # read in aggregated cluster assignment file, genes as rows,
    # iterations as columns

    df = pd.read_table(infile, sep="\t", header=0, index_col=0)
    genes = df.index.values

    # instantiate an empy matrix to count the number of times each
    # gene appears with others at each iteration

    dmat = pd.DataFrame(index=genes,
                        columns=genes)
    dmat = dmat.fillna(0)

    # generate all pair-wise combinations of genes, index into dmat using these
    reps = df.columns.values

    # alternative, faster code for consensus clustering - might be able to
    # improve it by using less nested for loops

    # generate a set for each cluster that contains the genes belonging to
    # each cluster count the number of times two genes occur in a cluster
    # and add to the dataframe repeat for each resampling iteration
    # cluster sets are generated at every iteration - time improvement over
    # directly accessing and storing as two dataframes is a factor of 7-8

    for i in reps:
        # from the input dataframe generate a list of sets, one set for each
        # cluster in that resampling iteration.
        # repeat for every resampling iteration.

        clusters = set(df[i].values.tolist())
        cluster_dict = {}
        for col in clusters:
            cluster_dict[col] = []
        for gene in genes:
            k_gene = df[i][gene]
            cluster_dict[k_gene].append(gene)
        rep_list = []

        # for each cluster add all the genes with that cluster ID to the set
        # add all of the cluster sets to a list container

        for col in clusters:
            col_set = set()
            clust_col = cluster_dict[col]
            gene_members = itertools.combinations_with_replacement(clust_col,
                                                                   2)
            col_set.add(gene_members)
            rep_list.append(col_set)

        # count if two genes occur in the same cluster, adding to both sides
        # of the symmetrical dataframe

        for cluster_set in rep_list:
            for combs in cluster_set:
                for x in combs:
                    if x[0] == x[1]:
                        dmat[x[0]][x[1]] += 1
                    else:
                        dmat[x[0]][x[1]] += 1
                        dmat[x[1]][x[0]] += 1

    # calculate the proportion of co-occurences

    prob = lambda x: x/float(len(reps))

    probs_df = dmat.applymap(prob)

    return probs_df


def consensusClustering(infile,
                        cutHeight,
                        cluster_algorithm,
                        min_size=30,
                        deepsplit=False):
    '''
    hierachichal clustering based on gene-cluster correlation across
    resampled datasets.  cut tree based with dynamic tree cut
    TODO: change this to cutHeight?  i.e. 0.2 = 80% clustering
    agreement OR use dynamic tree cut without deepsplit.
    '''
    condition = infile.split("/")[1].split("-")[0]
    wgcna_out = "tmp.dir/consensus-WGCNA.out"

    R('''sink(file='%(wgcna_out)s')''' % locals())
    R('''suppressMessages(library("WGCNA"))''')
    R('''suppressMessages(library("flashClust"))''')

    E.info("loading distance matrix")

    df = pd.read_table(infile, sep="\t", header=0, index_col=0)
    labels = df.index.tolist()
    labels_r = ro.StrVector([l for l in labels])

    # py2ri requires activation
    pandas2ri.activate()
    df_r = pandas2ri.py2ri(df)

    R.assign("distance.frame", df_r)
    R.assign("labels", labels_r)

    # large matricies/distance objects may need more
    # memory - allocate 1GB
    R('''memory.limit(10000)''')
    R('''rownames(distance.frame) <- labels''')
    R('''distance_data <- data.matrix(distance.frame)''')

    E.info("clustering data by %s linkage" % cluster_algorithm)

    R('''clustering <- flashClust(as.dist(1-distance_data),'''
      '''method='%(cluster_algorithm)s')''' % locals())

    if cutHeight > float(0.01):
        R('''cluster_cut <- cutreeStatic(dendro=clustering, '''
          '''minSize=%(min_size)i, cutHeight=%(cutHeight)s)''' % locals())

    elif deepsplit:
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''deepSplit=T, minClusterSize=%(min_size)i)''' % locals())
    else:
        R('''cluster_cut <- cutreeDynamic(dendro=clustering, '''
          '''deepSplit=F, minClusterSize=%(min_size)i)''' % locals())

    R('''color_cut <- labels2colors(cluster_cut)''')
    R('''cluster_matched <- data.frame(cbind(rownames(distance_data),'''
      '''color_cut))''')
    R('''colnames(cluster_matched) = c("gene_id", "cluster")''')
    R('''cluster_matched <- data.frame(cluster_matched$gene_id,'''
      '''cluster_matched$cluster)''')

    # plot and save dendrogram of clustering
    # AH: disabled, requires plots.dir to exist which might not be the case
    # AH: and thus causes this method to fail. Path names need to be parameterizable.
    # R('''png("plots.dir/%(condition)s-dendrogram-consensus_clustering.png")'''
    #   % locals())
    # R('''plotDendroAndColors(dendro=clustering, colors=color_cut,'''
    #   '''groupLabels="Dynamic tree cut",'''
    #   '''dendroLabels=F, addGuide=T, guideHang=0.05, '''
    #   '''hang=0.03, main="%(condition)s")''' % locals())
    # R('''dev.off()''')
    # R('''sink(file=NULL)''')
    cluster_frame = pandas2ri.ri2py(R["cluster_matched"])

    return cluster_frame
