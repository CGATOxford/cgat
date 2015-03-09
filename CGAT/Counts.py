##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
Counts.py - methods for manipulating counts data frames
==========================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

Utility functions for dealing with counts data frames

Requirements:
* fastq-dump >= 2.1.7

Code
----

'''
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pandas as pd
import numpy as np
import numpy.ma as ma
import copy
import random
import sys


def geometric_mean(array, axis=0):
    '''return the geometric mean of an array removing all zero-values but
    retaining total length
    '''
    non_zero = ma.masked_values(array, 0)
    log_a = ma.log(non_zero)
    return ma.exp(log_a.mean(axis=axis))


class Counts(object):
    '''base class to store counts object'''

    def __init__(self, table):
        # read in table in the constructor for Counts
        # e.g counts = Counts(pd.read_csv(...))
        self.table = table
        assert self.table.shape, "Counts table is empty"

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

        * remove rows with less than x number of counts
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
        take = sum_counts > sum_counts.quantile(percentile)
        E.info("percentile filtering at level %f: keep=%i, discard=%i" %
               (percentile_rowsums,
                sum(take),
                len(take) - sum(take)))
        self.table = self.table[take]

    def normalise(self, method="deseq-size-factors"):

        '''return a table with normalized count data.

        Implemented methods are:

        million-counts

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

        elif method == "million-counts":
            self.size_factors = self.table.sum(axis=0)
            normed = self.table * 1000000.0 / self.size_factors
        else:
            raise NotImplementedError(
                "normalization method '%s' not implemented" % method)

        self.table = normed
        # make sure we did not lose any rows or columns
        assert normed.shape == self.table.shape


########################################################################
# these functions for spike-in should be re-written to work with the ###
# counts class                                                       ###
########################################################################

def mapGroups(design_table):
    '''extract groups from the design table and return dictionaries
    mapping groups to tracks'''

    groups = list(set([design_table['group'][ix] for
                       ix in design_table.index]))
    groups_to_tracks = {}

    for group in groups:
        match_group = (design_table['group'] == group).tolist()
        tmp_design = design_table.iloc[match_group, ]
        groups_to_tracks[group] = tmp_design.index.tolist()

    return (groups,  groups_to_tracks)


def mapGroupsWithSuffix(design_table, shuffle_suffix, keep_suffix):
    '''use the regexes and suffixes supplied to extract groups from the
    design table and return dictionaries mapping each group to tracks
    which should be kept to tracks which should be shuffled
    '''

    design_table = design_table.ix[design_table['include'] == 1, ]
    design_table = design_table.ix[design_table['pair'] == 1, ]

    groups = list(set(design_table['group'].tolist()))

    groups_to_keep_tracks = {}
    groups_to_spike_tracks = {}

    keep_suffix = keep_suffix.split(",")
    for group in groups:
        match_group = (design_table['group'] == group).tolist()
        tmp_design = design_table.iloc[match_group, ]

        groups_to_spike_tracks[group] = [
            x + shuffle_suffix for x in tmp_design.index.tolist()]

        groups_to_keep_tracks[group] = copy.copy(groups_to_spike_tracks[group])
        groups_to_keep_tracks[group].extend(
            [x + y for x in tmp_design.index.tolist() for y in keep_suffix])

    return (groups,  groups_to_keep_tracks, groups_to_spike_tracks)


def means2idxarrays(g1, g2, i_bins, c_bins, difference):
    '''take two arrays of values and return the initial values
    and differences as numpy digitised arrays'''

    # note this currently returns a bin number for values
    # which fall outside of the bin. This is unwanted behaviour!
    # one solution would be to make one more bin than required and
    # throw it away at the end

    if difference == "relative":
        # calculate difference between mean values for group1 and group2
        # g1 and g2 always the same length
        change = [g2[x] - g1[x] for x in range(0, len(g1))]
        initial = g1
    elif difference == "logfold":
        change = [np.log2((g2[x]+1.0) / (g1[x]+1.0))
                  for x in range(0, len(g1))]
        initial = [np.log2(g1[x]+1.0) for x in range(0, len(g1))]

    # return arrays of len(change) with the index position in c_bins
    # corresponding to the bin in which the value of change falls
    change_idx = np.digitize(change, c_bins)
    initial_idx = np.digitize(initial, i_bins)

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


def shuffleRows(df, i_bins, c_bins, tracks_map,  groups,
                difference, s_max=100, i=1):
    '''take a dataframe and shuffle the rows to obtain spike in rows.
    return the indices to obtain the rows from the original dataframe'''
    counts = np.zeros((len(i_bins) + 1, len(c_bins) + 1))

    indices = {(key1, key2): []
               for key1 in np.digitize(i_bins, i_bins)
               for key2 in np.digitize(c_bins, c_bins)}

    min_occup = min(counts.flatten())

    for iteration in range(0,  i):
        E.info("performing shuffling iteration number %i.." % (iteration + 1))
        if min_occup < s_max:
            group1_rand = copy.copy(df.index.tolist())
            group2_rand = copy.copy(df.index.tolist())
            random.shuffle(group1_rand)
            random.shuffle(group2_rand)

            # subset the dataframe rows in the first random order
            # and by the column_ids in the first group
            # and return means across columns. Repeat for group2
            group1_mean = df.ix[group1_rand,
                                tracks_map[groups[0]]].apply(
                                    np.mean, axis=1).tolist()
            group2_mean = df.ix[group2_rand,
                                tracks_map[groups[1]]].apply(
                                    np.mean, axis=1).tolist()

            # retrieve the index for the bin in which each index value falls
            change_idx, initial_idx = means2idxarrays(
                group1_mean, group2_mean, i_bins, c_bins, difference)

            # for each initial and change value co-ordinate
            for idx, coord in enumerate(zip(initial_idx, change_idx)):
                if coord in indices.keys():
                    # if max fill of bin not reached
                    if counts[coord] < s_max:
                        counts[coord] += 1
                        # ...append tuple of df indeces for groups
                        indices[coord].append((group1_rand[idx],
                                               group2_rand[idx]))

            # find the minimum value in the counts array
            min_occup = min(counts.flatten())

    E.info("The largest bin has %i entries" % max(counts.flatten()))

    return indices, counts


def shuffleCluster(i_bins, c_bins, tracks_map, groups,
                   difference, s_max, i, clusters_dict,
                   s_bins_max, s_bins_min, s_bins_width):
    '''take a dictionary containing clusters (subdataframes) and shuffle
    subregions of clusters to obtain spike in clusters.
    return indeces from which the spike in clusters can be obtained from the
    original dataframe
    '''
    s_bins = range(s_bins_min, s_bins_max+1, s_bins_width, )

    counts = np.zeros((len(i_bins)+1, len(c_bins)+1, len(s_bins)+1))

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
            g1_rand = np.random.permutation(clusters_dict.keys())
            g2_rand = np.random.permutation(clusters_dict.keys())
            for perm in range(0, len(g1_rand)):
                cluster1 = clusters_dict[g1_rand[perm]].ix[
                    :, tracks_map[groups[0]]]
                cluster2 = clusters_dict[g2_rand[perm]].ix[
                    :, tracks_map[groups[1]]]
                c1_rand_s = random.randint(
                    min(cluster1.index), max(cluster1.index)-size)
                c1_e = int(c1_rand_s + size)
                c2_rand_s = random.randint(
                    min(cluster2.index), max(cluster2.index)-size)
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
            size_idx = np.digitize([size]*len(initial_idx), s_bins)
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
    s_min threshold'''
    output_indices_keep = copy.copy(indices)
    for key in indices.keys():
        if counts[key] < s_min:
            output_indices_keep.pop(key)

    E.info("%s/%s bins retained" % (len(output_indices_keep.keys()),
                                    len(indices.keys())))
    return output_indices_keep


def outputSpikes(df, id_column, indices, tracks_map, groups,
                 method, spike_type,
                 min_cbin, width_cbin, min_ibin, width_ibin,
                 min_sbin, width_sbin):

    def printHeader(keep_columns, groups):
        header = keep_columns
        header.extend(tracks_map[groups[0]])
        header.extend(tracks_map[groups[1]])
        sys.stdout.write("%s\n" % "\t".join(map(str, header)))
        return header

    if spike_type == "row":
        keep_cols = printHeader(id_column, groups)
    elif spike_type == "cluster":
        keep_cols = printHeader(["contig", "position"], groups)
    if method == "append":
        df = df.ix[:, keep_cols]
        df.to_csv(sys.stdout, index=False, header=False, sep="\t",
                  dtype={'position': int})

    def getInitialChangeSize(key, width_ibin, min_ibin, width_cbin, min_cbin,
                             width_s_bin, min_s_bin):
        initial_bin, change_bin, size_bin = key
        # initial and change values are the center of the bin
        initial = ((initial_bin*width_ibin) + min_ibin - (width_ibin*0.5))
        change = ((change_bin*width_cbin) + min_cbin - (width_cbin*0.5))
        size = ((size_bin*width_sbin) + min_sbin - 1)
        return initial, change, size

    def getInitialChange(key, width_ibin, min_ibin, width_cbin, min_cbin):
        initial_bin, change_bin = key
        # initial and change values are the center of the bin
        initial = ((initial_bin*width_ibin) + min_ibin - (width_ibin*0.5))
        change = ((change_bin*width_cbin) + min_cbin - (width_cbin*0.5))
        return initial, change

    n = 0
    if spike_type == "row":
        for key in indices:
            for pair in indices[key]:
                initial, change = getInitialChange(key, width_ibin, min_ibin,
                                                   width_cbin, min_cbin)
                row = ["_".join(map(str,
                                    ("spike-in", initial, change, n)))]
                row.extend(df.ix[pair[0], tracks_map[groups[0]]])
                row.extend(df.ix[pair[1], tracks_map[groups[1]]])
                sys.stdout.write("%s\n" % "\t".join(map(str, row)))
                n += 1

    elif spike_type == "cluster":
        for key in indices.keys():
            initial, change, size = getInitialChangeSize(
                key, width_ibin, min_ibin, width_cbin,
                min_cbin, width_sbin, min_sbin)
            for values in indices[key]:
                (c1s, c1e, c2s, c2e, c1rs, c1re,
                 c2rs, c2re) = values
                cluster_id = "_".join(map(str, ("spike-in", initial,
                                                change, size, c1rs-c1s, n)))
                temp_cluster_df = df.ix[c1s:c1e, keep_cols]
                temp_cluster_df['contig'] = cluster_id
                temp_cluster_swap = df.ix[c2rs:c2re, tracks_map[groups[1]]]
                temp_cluster_swap.set_index(df.ix[c1rs:c1re].index,
                                            drop=True,  inplace=True)
                temp_cluster_df.ix[c1rs:c1re, tracks_map[
                    groups[1]]] = temp_cluster_swap
                temp_cluster_df.to_csv(sys.stdout, index=False,
                                       header=False, sep="\t",
                                       dtype={'position': int})
                n += 1

##########################################################################
##########################################################################

##########################################################################
# remove the functions below (now methods for class Counts or          ###
# class ExpDesign (expression.py)                                      ###
##########################################################################


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

    million-counts

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

    elif method == "million-counts":
        size_factors = counts.sum(axis=0)
        normed = counts * 1000000.0 / size_factors
    else:
        raise NotImplementedError(
            "normalization method '%s' not implemented" % method)

    # make sure we did not loose any rows or columns
    assert normed.shape == counts.shape

    return normed, size_factors

##########################################################################
##########################################################################
