'''
Counts.py - methods for manipulating counts data frames
==========================================================

:Tags: Python

Utility functions for dealing with counts data frames

Requirements:

Code
----

'''
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pandas as pd
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
import numpy as np
import numpy.ma as ma
import copy
import random
import sys
import sklearn.preprocessing as preprocessing

pandas2ri.activate()


def geometric_mean(array, axis=0):
    '''return the geometric mean of an array removing all zero-values but
    retaining total length
    '''
    non_zero = ma.masked_values(array, 0)
    log_a = np.log(non_zero)
    return log_a.mean(axis=axis)
#ma.exp(log_a.mean(axis=axis))


class Counts(object):
    """base class to store counts object

    Attributes
    -----------

    table : pandas DataFrame
       dataframe object with count data

    """

    def __init__(self, filename_or_table):
        # read in table in the constructor for Counts
        # e.g counts = Counts(pd.read_csv(...))

        if isinstance(filename_or_table, str):
            self.table = pd.read_csv(filename_or_table,
                                     sep="\t",
                                     index_col=0)
        else:
            self.table = filename_or_table

        assert self.table.shape, "counts table is empty"

    def clone(self):
        return copy.deepcopy(self)

    def restrict(self, design):
        ''' remove samples not in design '''

        self.table = self.table.ix[:, design.samples]

    def removeSamples(self, min_counts_per_sample=10):
        ''' remove samples without low counts '''

        max_counts_per_sample = self.table.max()

        low_samples = max_counts_per_sample < min_counts_per_sample
        high_samples = [not x for x in low_samples.tolist()]
        sample_names = self.table.columns
        nlow_samples = sum(low_samples)

        if nlow_samples:
            E.warn("%i empty samples are being removed: %s" %
                   (nlow_samples,
                    ",".join([sample_names[x] for x, y in

                              enumerate(low_samples) if y])))
            self.table = self.table.iloc[:, high_samples]

    def removeObservationsFreq(self, min_counts_per_row=1):
        '''remove Observations (e.g genes)

        * remove rows with less than x counts
        '''

        # Remove rows with low counts
        max_counts_per_row = self.table.max(1)
        self.table = self.table[
            max_counts_per_row >= min_counts_per_row]
        observations, samples = self.table.shape
        E.info("trimmed data: %i observations for %i samples" %
               (observations, samples))

    def removeObservationsPerc(self, percentile_rowsums=10):
        '''remove Observations (e.g genes)

        * remove the lowest percentile of rows in the table, sorted
           by total tags per row
        '''

        # percentile filtering
        percentile = float(percentile_rowsums) / 100.0
        sum_counts = self.table.sum(1)
        take = sum_counts >= sum_counts.quantile(percentile)
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (percentile_rowsums,
                sum(take),
                len(take) - sum(take)))
        self.table = self.table[take]

    def normalise(self, method="deseq-size-factors", row_title="total"):
        '''return a table with normalized count data.

        Implemented methods are:

        total-column

           Divide each value by the column total and multiply by 1,000,000

        deseq-size-factors

           Construct a "reference sample" by taking, for each row, the
           geometric mean of the counts in all samples.

           To get the sequencing depth of a column relative to the
           reference, calculate for each row the quotient of the counts in
           your column divided by the counts of the reference sample. Now
           you have, for each row and column, an estimate of the depth
           ratio.

           Simply take the median of all the quotients in a column to get
           the relative depth of the library.

           Divide all values in a column by the normalization factor. This
           normalization method removes all rows with a geometric mean of
           0.

        total-row

           Divide each value in a sample by the value in a particular row.
           The name of the row is given by `row_title`.

        total-count

           Normalise all values in a column by the ratio of the per-column
           sum of counts and the average column count across all rows.

        This method normalises the counts and returns the normalization
        factors that have been applied.

        '''

        if method == "deseq-size-factors":

            # compute row-wise geometric means
            geometric_means = geometric_mean(self.table, axis=1)

            # remove 0 geometric means
            take = geometric_means > 0
            geometric_means = geometric_means[take]
            self.table = self.table[take]

            normed = (self.table.T / geometric_means).T

            self.size_factors = normed.median(axis=0)

            normed = self.table / self.size_factors

        elif method == "edger":

            getCPM = R('''
            suppressMessages(library(edgeR))

            function(counts){
            countsTable = DGEList(counts)
            countsTable = calcNormFactors(countsTable)

            countsTable.cpm <- cpm(countsTable, normalized.lib.sizes=TRUE)

            return (as.data.frame(countsTable.cpm))
            }''' % locals())

            normed = pandas2ri.ri2py(getCPM(pandas2ri.py2ri(self.table)))
            print(normed)

        elif method == "total-count":

            # compute column-wise sum
            column_sums = self.table.sum(axis=0)
            column_sums_mean = np.mean(column_sums)

            self.size_factors = [(column_sums_mean / x) for x in column_sums]
            normed = self.table * self.size_factors

        elif method == "total-column":

            self.size_factors = self.table.sum(axis=0)
            normed = self.table * 1000000.0 / self.size_factors

        elif method == "total-row":

            self.size_factors = self.table.loc[row_title]
            self.table.drop(row_title, inplace=True)
            normed = self.table * 1000000.0 / self.size_factors
        else:
            raise NotImplementedError(
                "normalisation method '%s' not implemented" % method)

        # make sure we did not lose any rows or columns before replacing table
        assert normed.shape == self.table.shape

        self.table = normed

    def sort(self, sort_columns, inplace=True):
        ''' sort counts table by columns'''

        def sort_table(counts):
            counts.table = counts.table.sort(columns=sort_columns)

        if inplace:
            sort_table(self)
            return None

        else:
            tmp_counts = self.clone()
            sort_table(tmp_counts)
            return tmp_counts

    def log(self, base=2, pseudocount=1, inplace=True):
        ''' log transform the counts table '''

        if inplace:
            self.table = np.log(self.table + pseudocount)
            return None

        else:
            tmp_counts = self.clone()
            tmp_counts.table = np.log(tmp_counts.table + pseudocount)
            return tmp_counts

    def transform(self, method="vst", design=None, inplace=True, blind=True):
        '''
        perform transformation on counts table
        current methods are:
         - deseq2 variance stabalising transformation
        - deseq rlog transformation

        Need to supply a design table if not using "blind"
        '''

        assert method in ["vst", "rlog"], ("method must be one of"
                                           "[vst, rlog]")

        method2function = {"vst": "varianceStabilizingTransformation",
                           "rlog": "rlog"}

        t_function = method2function[method]

        r_counts = pandas2ri.py2ri(self.table)

        if not blind:
            assert design, ("if not using blind must supply a design table "
                            "(a CGAT.Expression.ExperimentalDesign object")

            # currently this only accepts "~group" design
            transform = R('''
            function(df, design){

            suppressMessages(library('DESeq2'))

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= df, colData = design, design = ~group))

            transformed <- suppressMessages(%(t_function)s(dds, blind=FALSE))
            transformed_df <- as.data.frame(assay(transformed))

            return(transformed_df)
            }''' % locals())

            r_design = pandas2ri.py2ri(design.table)
            df = pandas2ri.ri2py(transform(r_counts, r_design))

        else:

            transform = R('''
            function(df){

            suppressMessages(library('DESeq2'))

            design = data.frame(row.names = colnames(df),
                                group = seq(1, length(colnames(df))))

            dds <- suppressMessages(DESeqDataSetFromMatrix(
                     countData= df, colData = design, design = ~group))

            transformed <- suppressMessages(%(t_function)s(dds, blind=TRUE))
            transformed_df <- as.data.frame(assay(transformed))

            return(transformed_df)
            }''' % locals())

            df = pandas2ri.ri2py(transform(r_counts))

        # losing rownames for some reason during the conversion?!
        df.index = self.table.index

        if inplace:
            self.table = df
            # R replaces "-" in column names with ".". Revert back!
            self.table.columns = [x.replace(".", "-")
                                  for x in self.table.columns]
        else:
            tmp_counts = self.clone()
            tmp_counts.table = df
            tmp_counts.table.columns = [x.replace(".", "-")
                                        for x in tmp_counts.table.columns]
            return tmp_counts

    def plotTransformations(self, plot_filename="transformations.png"):
        '''perform transformations and plot to compare'''

        log = self.log(inplace=False)
        vst = self.transform(inplace=False, method="vst")
        rlog = self.transform(inplace=False, method="rlog")

        def counts2meanSdPlot(df):
            ''' takes a counts object and returns a dataframe with standard
            deviation vs. average expression'''

            mean_sd_df = pd.DataFrame({"mean": df.table.apply(np.mean, axis=1),
                                       "std": df.table.apply(np.std, axis=1)})

            mean_sd_df = mean_sd_df[mean_sd_df["mean"] > 0]

            mean_sd_df.sort("mean", inplace=True)
            return mean_sd_df

        plotTransformations = R('''
        function(log_df, rlog_df, vst_df, plot_outfile){

          library(ggplot2)
          library("gridExtra")

          plotMeanSd = function(df, ymax, title){
            df$index <- as.numeric(row.names(df))
            df$bin = .bincode(df$index, seq(1,max(df$index),1000))

            l_txt = element_text(size=20)
            m_txt = element_text(size=15)

            p = ggplot(df, aes(x=index, y=std)) +
                geom_point(alpha=0.1)  +
                stat_summary(mapping=aes(x=bin*1000, y=std),
                             col="red", geom = "point",
                             fun.y = "median", size=5) +
            theme_minimal() +
            ylim(0, ymax) +
            xlab("index") + ylab("sd") +
            theme(axis.text.x=element_text(size=15, angle=90,
                                           hjust=1, vjust=0.5),
                  axis.text.y=m_txt,
                  axis.title.y=m_txt,
                  title=m_txt) +
            ggtitle(title)

            return (p)
        }

        max_sd = max(log_df$std, rlog_df$std, vst_df$std)

        p1 = plotMeanSd(log_df, max_sd, "log")
        p2 = plotMeanSd(rlog_df, max_sd, "rlog")
        p3 = plotMeanSd(vst_df, max_sd, "vst")

        png("%(plot_filename)s",
            width=10, height=10, units="in", res=400, pointsize=1)
        grid.arrange(p1, p2, p3, ncol=3)
        dev.off()
        }''' % locals())

        mean_sd_log_df = counts2meanSdPlot(log)
        mean_sd_vst_df = counts2meanSdPlot(vst)
        mean_sd_rlog_df = counts2meanSdPlot(rlog)

        plotTransformations(mean_sd_log_df, mean_sd_rlog_df, mean_sd_vst_df)

    def zNormalise(self, inplace=True):
        ''' normalise each row to zero mean and unit variance
        (z-score) '''

        samples = self.table.columns
        genes = self.table.index

        z_df = pd.DataFrame(preprocessing.scale(
            self.table, axis=1, with_mean=True, with_std=True, copy=False))

        z_df.index = genes
        z_df.columns = samples

        if inplace:
            self.table = z_df
        else:
            return Counts(z_df)

    def plotDendogram(self, plot_filename=None,
                      distance_method="euclidean",
                      clustering_method="ward.D2"):

        r_counts = pandas2ri.py2ri(self.table)

        makeDendogram = R('''
        function(counts){
          png("%(plot_filename)s")
          par(mar = c(1,4,1,1))
          plot(hclust(dist(t(counts), method = "%(distance_method)s"),
                      method = "%(clustering_method)s"), main="")
          dev.off()
        }''' % locals())

        makeDendogram(r_counts)

    def plotPCA(self, design,
                variance_plot_filename=None, pca_plot_filename=None,
                x_axis="PC1", y_axis="PC2", colour="group", shape="group"):
        ''' use the prcomp function in base R to perform principal components
        analysis.

        Can specify colour and shape as either variables from design table
        or sample names (seperated into id_1, id_2, id_3 based on samples
        having names formated e.g Tissue-Treatment-Replicate)'''

        # TS: swap this for regexes
        assert (x_axis[0:2] == "PC" and y_axis[0:2] == "PC"),\
            "x_axis and y_axis names must start with 'PC'"

        r_counts = pandas2ri.py2ri(self.table)
        r_design = pandas2ri.py2ri(design.table)

        pc_number_1 = int(x_axis.replace("PC", ""))
        pc_number_2 = int(y_axis.replace("PC", ""))

        makePCA = R('''
        function(counts, design){

          suppressMessages(library(ggplot2))
          suppressMessages(library(grid))

          gene_pca <- prcomp(na.omit(t(counts), center = TRUE, scale=TRUE))

          m_text = element_text(size=12)
          s_text = element_text(size=8)


          variance = gene_pca$sdev^2
          variance_explained = round(variance/sum(variance), 5)

          variance_df = data.frame("Variance_explained" = variance_explained,
                                 "PC" = seq(1, length(variance)))
          p_variance = ggplot(variance_df, aes(x=PC, y=Variance_explained))+
          geom_point()+
          geom_line()+
          theme_bw()+
          ylab("Variance explained (%%)")+
          theme(axis.text.x = m_text,
                axis.title.y = m_text,
                axis.title.x = m_text,
                axis.text.y = m_text,
                aspect.ratio=1)

          ggsave("%(variance_plot_filename)s", width=10, height=10, unit="cm")

          PCs_df = data.frame(gene_pca$x)
          PCs_df['sample'] <- rownames(PCs_df)
          design['sample'] <- gsub("-", ".", rownames(design))

          PCs_df = merge(PCs_df, design)

          PCs_df$id_1 = sapply(strsplit(PCs_df$sample, "\\\."), "[", 1)
          PCs_df$id_2 = sapply(strsplit(PCs_df$sample, "\\\."), "[", 2)
          PCs_df$id_3 = sapply(strsplit(PCs_df$sample, "\\\."), "[", 3)

          p_pca = ggplot(PCs_df, aes(x=%(x_axis)s, y=%(y_axis)s)) +
          geom_point(size=3,
                     aes(shape=as.factor(%(shape)s),
                         colour=as.factor(%(colour)s))) +
          scale_colour_discrete(name=guide_legend(title='%(colour)s')) +
          scale_shape_discrete(name=guide_legend(title='%(shape)s')) +
          xlab(paste0('PC%(pc_number_1)i (Variance explained = ' ,
                       round(100 * variance_explained[%(pc_number_1)i], 1),
                       '%%)')) +
          ylab(paste0('PC%(pc_number_2)i (Variance explained = ' ,
                       round(100 * variance_explained[%(pc_number_2)i], 1),
                       '%%)')) +
          theme_bw() +
          theme(axis.text.x = s_text, axis.text.y = s_text,
                title = m_text, legend.text = m_text,
                legend.title = m_text, aspect.ratio=1)

          ggsave("%(pca_plot_filename)s", width=10, height=10, unit="cm")

        }''' % locals())

        makePCA(r_counts, r_design)

    def plotPairwise(self, outfile_scatter, outfile_heatmap, subset=False):
        ''' use the R base pairs function to plot all pairwise
        correlations between the samples

        subset will randomly subset n rows to speed up plotting'''

        plotGGpairs = R('''
        library(ggplot2)
        library(reshape2)

        function(df, df2){

        colnames(df) <- gsub("-", "_", colnames(df))

        width <- height <-  length(colnames(df)) * 100

        png("%(outfile_scatter)s", width=width, height=height, units = "px")

        panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          r <- abs(cor(x, y))
          txt <- format(c(r, 0.123456789), digits = digits)[1]
          txt <- paste0(prefix, txt)
          if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
          text(0.5, 0.5, txt, cex = cex.cor * r * 50)}

        panel.hist = function (x, ...) {
          par(new = TRUE)
          hist(x,
               breaks=30,
               col = "light blue",
               probability = TRUE,
               axes = FALSE,
               main = "")
          rug(x)}

        pairs(df, pch=20, cex=0.1,
              lower.panel = panel.smooth, upper.panel = panel.cor,
              diag.panel=panel.hist)

        dev.off()

        # use the full (not subsetted) df for heatmap correlation
        colnames(df2) <- gsub("-", "_", colnames(df2))
        df_cor <- as.matrix(cor(df2))
        df_cor_melt <- melt(df_cor)

        p = ggplot(df_cor_melt, aes(Var1, Var2, fill=value)) +
        geom_tile() +
        xlab("") + ylab("") +
        scale_fill_continuous(low="white", high="dodgerblue4", name="Correlation") +
        theme_bw() +
        theme(
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=0),
        aspect.ratio=1)

        ggsave("%(outfile_heatmap)s", width=10, height=10, unit="cm")
        }''' % locals())

        if subset:
            if len(self.table.index) > subset:
                rows = random.sample(self.table.index, subset)
                subset_counts = pandas2ri.py2ri(self.table.ix[rows])
            else:
                subset_counts = pandas2ri.py2ri(self.table)
        else:
            subset_counts = pandas2ri.py2ri(self.table)

        plotGGpairs(subset_counts, pandas2ri.py2ri(self.table))

    def heatmap(self, plotfile, zscore=False):
        ''' plots a heatmap
        set zscore=True to use a divergent colour scale
        '''
        # to do: add option to parse design file and add coloured row for
        # variable specified in design file.

        plotHeatmap = R('''
        function(df, zscore){

        library("Biobase")
        library("RColorBrewer")
        library("gplots")

        if(zscore[1]==TRUE){
        #hmcol <- colorRampPalette(colors = c("red", "white", "blue"))
        PuOr <- brewer.pal(11, "PuOr")
        hmcol <- c(colorRampPalette(c(PuOr[1], PuOr[6]))(100),
                   colorRampPalette(c(PuOr[6], PuOr[11]))(100)[-1])
        }
        else{
        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        }

        png("%(plotfile)s", width=1000, height=1000, units="px")

        heatmap.2(as.matrix(df),
                  col = hmcol,
                  scale="none", trace="none", margin=c(18, 10),
                  dendrogram="both", cexCol=2,
                  labRow = "",
                  hclustfun = function(x) hclust(x, method = 'average'),
                  distfun = dist)
        dev.off()
        }''' % locals())

        r_counts = pandas2ri.py2ri(self.table)

        plotHeatmap(r_counts, ro.BoolVector([zscore]))

    def shuffleRows(self,
                    min_cbin, max_cbin, width_cbin,
                    min_ibin, max_ibin, width_ibin,
                    tracks_map,  groups,
                    difference, s_max=100, i=1):
        '''take a dataframe and shuffle the rows to obtain spike in rows.
        return the indices to obtain the rows from the counts table
        and the counts per bin

        The shuffling of rows is parameterised as follows:
        * min_i_bins = minimum bin for initial values
        * width_i_bins = width of bins for initial values
        * max_i_bins = maximum bin for initial values
        * min_c_bins = minimum bin for change values
        * width_c_bins = width of bins for change values
        * max_c_bins = maximum bin for change values
        * tracks_map = dictionary mapping groups to tracks
        * difference = "relative", "logfold" or "abs_logfold"
        * s_max = maximum number of spikes per bin
        * i = number of iterations. More iterations = more filled bins
        '''

        # make bins with an extra bin at the end to capture spike-ins with
        # initial or changes values over the intended range.
        # these are ignored later
        c_bins = np.arange(min_cbin, max_cbin + width_cbin, width_cbin)
        i_bins = np.arange(min_ibin, max_ibin + width_ibin, width_ibin)

        bin_counts = np.zeros((len(i_bins) + 1, len(c_bins) + 1))

        indices = {(key1, key2): []
                   for key1 in np.digitize(i_bins, i_bins)
                   for key2 in np.digitize(c_bins, c_bins)}

        min_occup = min(bin_counts.flatten())
        for iteration in range(0,  i):
            E.info("performing shuffling iteration number %i.." % (
                iteration + 1))
            if min_occup < s_max:
                group1_rand = copy.copy(self.table.index.tolist())
                group2_rand = copy.copy(self.table.index.tolist())
                random.shuffle(group1_rand)
                random.shuffle(group2_rand)

                # subset the dataframe rows in the first random order
                # and by the column_ids in the first group
                # and return means across columns. Repeat for group2
                group1_mean = self.table.ix[group1_rand,
                                            tracks_map[groups[0]]].apply(
                                                np.mean, axis=1).tolist()
                group2_mean = self.table.ix[group2_rand,
                                            tracks_map[groups[1]]].apply(
                                                np.mean, axis=1).tolist()
                # retrieve the index for the bin in
                # which each index value falls
                change_idx, initial_idx = means2idxarrays(
                    group1_mean, group2_mean, i_bins, c_bins, difference)

                # for each initial and change value coordinate
                for idx, coord in enumerate(zip(initial_idx, change_idx)):
                    # ignore spike-in if change or initial fall into the final
                    # bin
                    if coord[0] < len(i_bins) and coord[1] < len(c_bins):
                        if coord in list(indices.keys()):
                            # if max fill of bin not reached
                            if bin_counts[coord] < s_max:
                                # ...append tuple of df indeces for groups
                                bin_counts[coord] += 1
                                indices[coord].append((group1_rand[idx],
                                                       group2_rand[idx]))

                    # find the minimum value in the counts array this should
                    # save time if all bins are filled before final iteration
                    min_occup = min(bin_counts.flatten())

        E.info("The largest bin has %i entries" % max(bin_counts.flatten()))

        return indices, bin_counts

    def outputSpikes(self, indices, tracks_map, groups,
                     output_method, spike_type,
                     min_cbin, width_cbin, max_cbin,
                     min_ibin, width_ibin, max_ibin,
                     min_sbin=1, width_sbin=1, max_sbin=1,
                     append=False):
        ''' method to output spike-ins generated by shuffling rows
        (counts.shuffleRows) or clusters of rows (counts.shuffleCluster)

        parameters:
        * indices - indices for spike-ins
        * min_i_bins = minimum bins for initial values
        * width_i_bins = width of bins for initial values
        * max_i_bins = maximum bin for initial values
        * min_c_bins = minimum bin for change values
        * width_c_bins = width of bins for change values
        * max_c_bins = maximum bin for change values
        * min_s_bins = minimum bin for size values
        * width_s_bins = width of bins for size values
        * max_s_bins = maximum bin for size values
        * tracks_map = dictionary mapping groups to tracks
        * spike_type = "relative" or "logfold"
        * output_method = "append" or "seperate"
        '''

        def makeHeader(tracks_map, groups, keep_columns=None):
            if keep_columns:
                header = keep_columns
            else:
                header = []
            header.extend(tracks_map[groups[0]])
            header.extend(tracks_map[groups[1]])

            return header

        if spike_type == "row":
            index = True
            # keep_columns = ["spike"]
            keep_columns = None

        elif spike_type == "cluster":
            index = False
            keep_columns = ["contig", "position"]

        header = makeHeader(tracks_map, groups, keep_columns=keep_columns)

        if output_method == "append":
            self.table = self.table.ix[:, header]
            self.table.to_csv(sys.stdout, index=index, header=True, sep="\t",
                              dtype={'position': int})
        else:
            sys.stdout.write("%s\t%s\n" % (
                "spike", "\t".join(map(str, header))))

        def getInitialChangeSize(key, width_ibin, min_ibin, width_cbin,
                                 min_cbin, width_s_bin, min_s_bin):
            initial_bin, change_bin, size_bin = key

            # initial and change values are the center of the bin
            initial = ((initial_bin * width_ibin) +
                       min_ibin - (width_ibin * 0.5))
            change = ((change_bin * width_cbin) +
                      min_cbin - (width_cbin * 0.5))
            size = ((size_bin * width_sbin) + min_sbin - 1)
            return initial, change, size

        def getInitialChange(key, width_ibin, min_ibin, width_cbin, min_cbin):
            initial_bin, change_bin = key

            # initial and change values are the center of the bin
            initial = ((initial_bin * width_ibin) +
                       min_ibin - (width_ibin * 0.5))
            change = ((change_bin * width_cbin) +
                      min_cbin - (width_cbin * 0.5))
            return initial, change

        n = 0

        if spike_type == "row":

            for key in indices:
                for pair in indices[key]:
                    initial, change = getInitialChange(
                        key, width_ibin, min_ibin, width_cbin, min_cbin)
                    row = ["_".join(map(str,
                                        ("spike-in", initial, change, n)))]
                    row.extend(self.table.ix[pair[0], tracks_map[groups[0]]])
                    row.extend(self.table.ix[pair[1], tracks_map[groups[1]]])
                    sys.stdout.write("%s\n" % "\t".join(map(str, row)))
                    n += 1

        elif spike_type == "cluster":
            for key in sorted(indices.keys()):
                initial, change, size = getInitialChangeSize(
                    key, width_ibin, min_ibin, width_cbin,
                    min_cbin, width_sbin, min_sbin)
                for values in sorted(indices[key]):
                    (c1s, c1e, c2s, c2e, c1rs, c1re,
                     c2rs, c2re) = values
                    cluster_id = "_".join(
                        map(str, ("spike-in", initial, change,
                                  size, c1rs - c1s, n)))

                    temp_cluster_df = self.table.ix[c1s:c1e, keep_cols]
                    temp_cluster_df['contig'] = cluster_id
                    temp_cluster_swap = self.table.ix[
                        c2rs:c2re, tracks_map[groups[1]]]
                    temp_cluster_swap.set_index(self.table.ix[c1rs:c1re].index,
                                                drop=True,  inplace=True)
                    temp_cluster_df.ix[c1rs:c1re, tracks_map[
                        groups[1]]] = temp_cluster_swap
                    temp_cluster_df.to_csv(sys.stdout, index=index,
                                           header=False, sep="\t",
                                           dtype={'position': int})
                    n += 1


def means2idxarrays(g1, g2, i_bins, c_bins, difference):
    '''take two arrays of values and return the initial values
    and differences as numpy digitised arrays'''

    if difference == "relative":
        # calculate difference between mean values for group1 and group2
        # g1 and g2 always the same length
        change = [g2[x] - g1[x] for x in range(0, len(g1))]
        initial = g1

    elif difference == "logfold":
        change = [np.log2((g2[x] + 1.0) / (g1[x] + 1.0))
                  for x in range(0, len(g1))]
        initial = [np.log2(g1[x] + 1.0) for x in range(0, len(g1))]

    elif difference == "abs_logfold":
        change = [abs(np.log2((g2[x] + 1.0) / (g1[x] + 1.0)))
                  for x in range(0, len(g1))]
        initial = [max(np.log2(g1[x] + 1.0), np.log2(g2[x] + 1.0))
                   for x in range(0, len(g1))]

    # return arrays of len(change) with the index position in c_bins
    # corresponding to the bin in which the value of change falls
    change_idx = np.digitize(change, c_bins, right=True)
    initial_idx = np.digitize(initial, i_bins, right=True)

    return(change_idx, initial_idx)


def findClusters(df, distance, size, tracks_map, groups):
    '''define clusters of genomic loci depending on thresholds for
    size and minimum members per cluster
    This was written with CpGs in mind but will work with any data frame
    containing "position" and "contig" columns'''

    positions = df['position'].tolist()
    contigs = df['contig'].tolist()
    current_pos = 0
    cluster_ix = []
    current_contig = ""
    cluster_dfs = {}
    n = 0
    for ix in range(0, len(positions)):
        next_pos = positions[ix]
        next_contig = contigs[ix]
        if (((next_pos < current_pos + distance) &
             (next_contig == current_contig))):
            cluster_ix.append(ix)
        else:
            if len(cluster_ix) >= size:
                start, end = (cluster_ix[0], cluster_ix[-1] + 1)
                cluster_dfs[n] = df.iloc[start:end]
                n += 1
            cluster_ix = []
            current_pos = next_pos
            current_contig = next_contig

    E.info("found %i clusters" % n)
    return (cluster_dfs)


def shuffleCluster(i_bins, c_bins, tracks_map, groups,
                   difference, s_max, i, clusters_dict,
                   s_bins_max, s_bins_min, s_bins_width):
    '''take a dictionary containing clusters (subdataframes) and shuffle
    subregions of clusters to obtain spike in clusters.
    return indeces from which the spike in clusters can be obtained from the
    original dataframe
    '''
    s_bins = list(range(s_bins_min, s_bins_max + 1, s_bins_width,))

    counts = np.zeros((len(i_bins) + 1, len(c_bins) + 1, len(s_bins) + 1))

    indices = {(key1, key2, key3): []
               for key1 in np.digitize(i_bins, i_bins)
               for key2 in np.digitize(c_bins, c_bins)
               for key3 in np.digitize(s_bins, s_bins)}

    for iteration in range(0,  i):
        E.info("performing shuffling iteration number %i.." % (iteration + 1))
        for size in s_bins:
            group1_mean = []
            group2_mean = []
            g1_rand_s = []
            g2_rand_s = []
            g1_rand = np.random.permutation(list(clusters_dict.keys()))
            g2_rand = np.random.permutation(list(clusters_dict.keys()))
            for perm in range(0, len(g1_rand)):
                cluster1 = clusters_dict[g1_rand[perm]].ix[
                    :, tracks_map[groups[0]]]
                cluster2 = clusters_dict[g2_rand[perm]].ix[
                    :, tracks_map[groups[1]]]
                c1_rand_s = random.randint(
                    min(cluster1.index), max(cluster1.index) - size)
                c1_e = int(c1_rand_s + size)
                c2_rand_s = random.randint(
                    min(cluster2.index), max(cluster2.index) - size)
                c2_e = int(c2_rand_s + size)

                c1_mean = np.mean(np.mean(cluster1.ix[c1_rand_s: c1_e]))
                c2_mean = np.mean(np.mean(cluster2.ix[c2_rand_s: c2_e]))
                group1_mean.append(c1_mean)
                group2_mean.append(c2_mean)
                g1_rand_s.append(c1_rand_s)
                g2_rand_s.append(c2_rand_s)

            change_idx, initial_idx,  = means2idxarrays(
                group1_mean, group2_mean, i_bins,
                c_bins,  difference)
            size_idx = np.digitize([size] * len(initial_idx), s_bins)
            for idx, coord in enumerate(zip(
                    initial_idx, change_idx, size_idx)):
                if counts[coord] < s_max:
                    counts[coord] += 1
                    cluster1_ix = clusters_dict[g1_rand[idx]].index.values
                    cluster1_start = cluster1_ix[0]
                    cluster1_end = cluster1_ix[-1]
                    cluster2_ix = clusters_dict[g2_rand[idx]].index.values
                    cluster2_start = cluster2_ix[0]
                    cluster2_end = cluster2_ix[-1]
                    c1_rand_s = g1_rand_s[idx]
                    c1_rand_e = int(c1_rand_s + size)
                    c2_rand_s = g2_rand_s[idx]
                    c2_rand_e = int(c2_rand_s + size)
                    indices[coord].append((
                        cluster1_start, cluster1_end, cluster2_start,
                        cluster2_end, c1_rand_s, c1_rand_e, c2_rand_s,
                        c2_rand_e))
    return indices, counts


def thresholdBins(indices, counts, s_min):
    '''use counts (np array) to remove bins from indices based on
    threshold (s_min)'''
    output_indices_keep = copy.copy(indices)
    for key in sorted(indices.keys()):
        if counts[key] < s_min:
            output_indices_keep.pop(key)

    E.info("%s/%s bins retained" % (len(list(output_indices_keep.keys())),
                                    len(list(indices.keys()))))
    return output_indices_keep


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
    counts_table = pd.read_table(tags_filename, sep="\t",
                                 index_col=0, comment="#")

    E.info("read data: %i observations for %i samples" % counts_table.shape)

    E.debug("sample names: %s" % list(counts_table.columns))

    inf = IOTools.openFile(design_filename)
    design_table = pd.read_csv(inf, sep="\t", index_col=0)
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

    low_samples = max_counts_per_sample < filter_min_counts_per_sample
    high_samples = [not x for x in low_samples.tolist()]
    sample_names = counts_table.columns
    nlow_samples = sum(low_samples)

    if nlow_samples:
        E.warn("%i empty samples are being removed: %s" %
               (nlow_samples,
                ",".join([sample_names[x] for x, y in

                          enumerate(low_samples) if y])))
        counts_table = counts_table.iloc[:, high_samples]

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

    # need to return altered design table based on filtering
    # in case samples are removed!

    design_table = design_table.ix[high_samples]

    nobservations, nsamples = counts_table.shape

    return nobservations, nsamples, counts_table, design_table


def normalizeTagData(counts, method="deseq-size-factors"):
    '''return a table with normalized count data.

    Implemented methods are:

    total-column

       Divide each value by the column total and multiply by 1,000,000

    deseq-size-factors

       Construct a "reference sample" by taking, for each row, the
       geometric mean of the counts in all samples.

       To get the sequencing depth of a column relative to the
       reference, calculate for each row the quotient of the counts in
       your column divided by the counts of the reference sample. Now
       you have, for each row and column, an estimate of the depth
       ratio.

       Simply take the median of all the quotients in a column to get
       the relative depth of the library.

       Divide all values in a column by the normalization factor.

    This method returns the normalized counts and any normalization
    factors that have been applied.

    '''

    if method == "deseq-size-factors":

        # compute row-wise geometric means
        geometric_means = geometric_mean(counts, axis=1)

        # remove 0 geometric means
        take = geometric_means > 0
        geometric_means = geometric_means[take]
        counts = counts[take]

        normed = (counts.T / geometric_means).T

        size_factors = normed.median(axis=0)

        normed = counts / size_factors

    elif method == "total-column":
        size_factors = counts.sum(axis=0)
        normed = counts * 1000000.0 / size_factors
    else:
        raise NotImplementedError(
            "normalization method '%s' not implemented" % method)

    # make sure we did not loose any rows or columns
    assert normed.shape == counts.shape

    return normed, size_factors
