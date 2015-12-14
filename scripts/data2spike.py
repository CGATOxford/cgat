"""
data2spike.py - generate spike-in from a data frame
=====================================================

:Author: Tom smith
:Release: $Id: data2spike.py
:Date: |today|
:Tags: Genomics Statistics Power


Purpose
-------
This scripts can be used to generate spike-ins from a dataframe.
The output can be used to perform a downstream power analysis.

Usage
-----
Spike-ins are generated over a range of initial values and simulated
differences. Additionally, cluster spike-ins are generate over a range
of subcluster sizes. The range and bin width is parametised, as is the
the minimum and maximum number of spike ins required in each bin. The
number of randomisation iterations can be modified in order to
generate sufficient numbers of spike ins and/or optimise run time.

Input
+++++
The input to this script is a dataframe containing one item per row
(e.g gene, CpG).

The dataframe must contain rows with numerical data for shuffling. The
dataframe can contain additional columns which are shuffled alongside
the shuffling columns, e.g associated data columns or sample-specific
columns

If cluster spike-ins need to be generated, the dataframe must contain
columns called "contig" and "position" from which to derive the clusters.

The script further requires a design table describing the groups and columns.
The design table has four columns::

      track   include group   pair
      GS1        1       S       1
      GS2        1       S       1
      GD1        1       D       1
      GD2        1       D       1


track
     name of track - should correspond to column header in the counts
     table.
include
     flag to indicate whether or not to include this data in the output
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests)

Output
++++++

The original dataframe is outputted along with the spike-ins.
Optionally, the spike-ins can be outputted seperately.  If row
spike-ins need to be appended the user must specify an idendification
column for the dataframe.

"""

import copy
import random
import sys
import re
import pandas as pd
import numpy as np
import CGAT.Experiment as E


def groupMappers(design_table, spike_regex, shuffle_suffix, keep_suffix):
    '''from design table, extract groups and return dictionaries mapping
    groups to tracks'''

    design_table = design_table.ix[design_table['include'] == 1, ]
    design_table = design_table.ix[design_table['pair'] == 1, ]

    keep = [ix for ix in design_table.index if
            re.match(spike_regex, design_table['track'][ix])]
    design_table = design_table.ix[keep]

    groups = list(set([design_table['group'][ix] for
                       ix in design_table.index]))

    groups_to_keep_ix = {}
    groups_to_keep_tracks = {}
    groups_to_spike_tracks = {}

    keep_suffix = keep_suffix.split(",")
    for group in groups:
        groups_to_keep_ix[group] = [
            ix for ix in design_table.index
            if re.match(group, design_table['group'][ix])]
        groups_to_spike_tracks[group] = [
            x + shuffle_suffix
            for x in design_table.ix[
                groups_to_keep_ix[group], "track"].tolist()]
        groups_to_keep_tracks[group] = copy.copy(groups_to_spike_tracks[group])
        groups_to_keep_tracks[group].extend([
            x + y for
            x in design_table.ix[groups_to_keep_ix[group], "track"].tolist()
            for y in keep_suffix])
    return (groups,  groups_to_keep_tracks, groups_to_spike_tracks)


def means2idxarrays(g1, g2, i_bins, c_bins, difference):
    '''take two arrays of values and return the initial values
    and differences as numpy digitised arrays'''

    if difference == "relative":
        # calculate difference between mean values for group1 and group2
        # g1 and g2 will always be same length.
        change = [g2[x] - g1[x] for x in range(0, len(g1))]
        initial = g1
    elif difference == "logfold":
        # JJ bugfix list comprehension...
        change = [np.log2((g2[x]+1.0) /
                          (g1[x]+1.0)) for x in range(0, len(g1))]
        initial = [np.log2(g1[x]+1.0) for x in range(0, len(g1))]

    # Returns an array of len(change) with the index position in c_bins
    # corresponding to the bin in which that value of change falls.
    change_idx = np.digitize(change, c_bins)
    initial_idx = np.digitize(initial, i_bins)

    return(change_idx, initial_idx)


def findClusters(df, distance, size, tracks_map, groups):
    '''define clusters of genomic loci depending on thresholds for
    size and minimum members per cluster'''

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
                start, end = (cluster_ix[0], cluster_ix[-1])
                cluster_dfs[n] = df.loc[start:end]
                n += 1
            cluster_ix = []
            current_pos = next_pos
            current_contig = next_contig

    E.info("found %i clusters" % n)
    return (cluster_dfs)


def shuffleCluster(i_bins, c_bins, tracks_map, groups,
                   difference, s_max, i, clusters_dict,
                   s_bins_max, s_bins_min, s_bins_width):

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


def shuffleRows(df, i_bins, c_bins, tracks_map,  groups,
                difference, s_max=100, i=1):

    counts = np.zeros((len(i_bins) + 1, len(c_bins) + 1))
    # JJ
    E.info("Initial count table:")
    E.info(counts)

    indices = {(key1, key2): []
               for key1 in np.digitize(i_bins, i_bins)
               for key2 in np.digitize(c_bins, c_bins)}
    # JJ
    # E.info(indices)

    # find min value in the counts ndarray
    min_occup = min(counts.flatten())

    for iteration in range(0,  i):
        E.info("performing shuffling iteration number %i.." % (iteration + 1))
        # if min value is less than the specified maximum occupancy
        if min_occup < s_max:
            # permute the df index axes to get a random row order
            group1_rand = np.random.permutation(df.index)
            group2_rand = np.random.permutation(df.index)

            # subset the dataframe rows in the first random row order
            # subset the dataframe columns_ids in first group
            # return a list of the means across columns in the random order.
            group1_mean = df.ix[group1_rand,
                                tracks_map[groups[0]]].apply(
                                    np.mean, axis=1).tolist()
            # JJ
            # E.info( group1_mean )

            # same for second random order and column Ids in second group.
            group2_mean = df.ix[group2_rand,
                                tracks_map[groups[1]]].apply(
                                    np.mean, axis=1).tolist()

            # JJ
            # E.info( group2_mean )

            # retrieve the index for the bin in which each of the index
            # values falls.
            change_idx, initial_idx = means2idxarrays(
                group1_mean, group2_mean, i_bins, c_bins, difference)

            # JJ
            # E.info( change_idx )
            # E.info( initial_idx )
            # E.info( initial_idx.__class__ )
            E.info("Max %s" % str(initial_idx.max()))

            # for each of the initial and change values co-ordinates...
            for idx, coord in enumerate(zip(initial_idx, change_idx)):
                # if that combination is desired...
                if coord in indices.keys():
                    if counts[coord] < s_max:
                        counts[coord] += 1
                        # ...append tuple of the expression val for each group
                        indices[coord].append((group1_rand[idx],
                                               group2_rand[idx]))
            min_occup = min(counts.flatten())

    # JJ
    E.info("Final count table:")
    E.info(counts)
    # E.info( indices )

    return indices, counts


def thresholdBins(indices, counts, s_min):
    output_indices_keep = copy.copy(indices)
    for key in indices.keys():
        if counts[key] < s_min:
            output_indices_keep.pop(key)

    E.info("%s bins retained" % len(output_indices_keep.keys()))
    E.info("Just to check, the largest bin has %i subregions" %
           (max(counts.flatten())))

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
        keep_cols = printHeader([id_column], groups)
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

    parser.add_option("-t", "--design-file-tsv", dest="design_file",
                      type="string",
                      help="filename with design [default=%default].")

    parser.add_option("-i", "--change-bin-min", dest="min_cbin",
                      type="float",
                      help="minimum bin for change bins [default=%default].")

    parser.add_option("-x", "--change-bin-max", dest="max_cbin",
                      type="float",
                      help="maximum bin for change bins [default=%default].")

    parser.add_option("-w", "--change-bin-width", dest="width_cbin",
                      type="float",
                      help="bin width for change bins [default=%default].")

    parser.add_option("-s", "--initial-bin-min", dest="min_ibin",
                      type="float",
                      help="minimum bin for initial bins[default=%default].")

    parser.add_option("-e", "--initial-bin-max", dest="max_ibin",
                      type="float",
                      help="maximum bin for intitial bins[default=%default].")

    parser.add_option("-y", "--initial-bin-width", dest="width_ibin",
                      type="float",
                      help="bin width intitial bins[default=%default].")

    parser.add_option("-l", "--spike-minimum", dest="min_spike",
                      type="int",
                      help="minimum number of spike-ins required within each bin\
                      [default=%default].")

    parser.add_option("-u", "--spike-maximum", dest="max_spike",
                      type="int",
                      help="maximum number of spike-ins allowed within each bin\
                      [default=%default].")

    parser.add_option("-d", "--difference-method", dest="difference",
                      type="choice", choices=("relative", "logfold"),
                      help="method to use for calculating difference\
                      [default=%default].")

    parser.add_option("-r", "--iterations", dest="iterations", type="int",
                      help="number of iterations [default=%default].")

    parser.add_option("-a", "--id_columns", dest="id", action="append",  # JJ is actually id_column
                      help="name of identification column(s)\
                      [default=%default].")

    parser.add_option("-c", "--cluster-maximum-distance",
                      dest="cluster_max_distance", type="int",
                      help="maximum distance between adjacent loci in cluster\
                      [default=%default].")

    parser.add_option("-b", "--cluster-minimum-size",
                      dest="cluster_min_size", type="int",
                      help="minimum number of loci required per cluster\
                      [default=%default].")

    parser.add_option("-f", "--spike-type",
                      dest="spike_type", type="choice",
                      choices=("row", "cluster"),
                      help="spike in type [default=%default].")

    parser.add_option("-g", "--subcluster-min-size",
                      dest="min_sbin", type="int",
                      help="minimum size of subcluster\
                      [default=%default].")

    parser.add_option("-j", "--subcluster-max-size",
                      dest="max_sbin", type="int",
                      help="maximum size of subcluster\
                      [default=%default].")

    parser.add_option("-k", "--subcluster-bin-width",
                      dest="width_sbin", type="int",
                      help="bin width for subcluster size\
                      [default=%default].")

    parser.add_option("-o", "--output-method",
                      dest="output_method", type="choice",
                      choices=("append", "seperate"),
                      help="defines whether the spike-ins should be appended\
                      to the original table or seperately [default=%default].")

    parser.add_option("-p", "--shuffle-column-suffix",
                      dest="shuffle_suffix", type="string",
                      help="the suffix of the columns which are to be shuffled\
                      [default=%default].")

    parser.add_option("-q", "--keep-column-suffix",
                      dest="keep_suffix", type="string",
                      help="a list of suffixes for the columns which are to be\
                      keep along with the shuffled columns[default=%default].")

    parser.add_option("-z", "--spike-regex",
                      dest="spike_regex", type="string", action="append",
                      help="a regex to identify tracks for shuffling\
                      [default=%default].")

    parser.set_defaults(
        design_file=None,
        output_method="seperate",
        difference="relative",
        spike_type="row",
        genome_file=None,
        min_cbin=None,
        max_cbin=None,
        width_cbin=None,
        min_ibin=None,
        max_ibin=None,
        width_ibin=None,
        id_columns=None,
        max_spike=100,
        min_spike=None,
        iterations=1,
        cluster_max_distance=100,
        cluster_min_size=10,
        min_sbin=1,
        max_sbin=9,
        width_sbin=2,
        id_column=None,
        shuffle_suffix=None,
        keep_suffix=None,
        spike_regex=".*")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    df = pd.read_table(sys.stdin, sep="\t")

    if options.spike_type == "cluster":
        # JJ- THIS ASSUMES THESE COLUMNS EXIST, WHICH THE SUBSEQUENT SCRIPT DOES NOT
        # THEREFORE MOVED
        df_sort = df.sort(columns=["contig", "position"])
        df_sort.set_index(df.index, inplace=True)
        # print df_sort.head()
    else:
        df_sort = df

    if not options.design_file:
        raise ValueError("a design table is required.")

    if options.spike_type == "cluster":
        if options.max_sbin > options.cluster_min_size:
            raise ValueError("maximum size of subscluster bin: %s\
            is greater than minimum size of cluster: %s" % (
                options.max_sbin, options.cluster_min_size))
        df_columns = set(df.columns.values.tolist())
        msg = "cluster analysis requires columns named 'contig' and\
        'position' in the dataframe"
        assert "contig" in df_columns and "position" in df_columns, msg

        msg2 = '''the minimum cluster size (%s) must be larger than the maximum
        size of subcluster to be shuffled (%s)''' % (options.cluster_min_size,
                                                     options.max_sbin)
        assert options.max_sbin < options.cluster_min_size, msg2

        df_sort = df.sort(columns=["contig", "position"])
        df_sort.set_index(df.index, inplace=True)

    if not options.id_column:
        if options.output_method == "append" and options.spike_type == "row":
            raise ValueError("id column must be specified to append rows")
        elif options.output_method != "append":
            options.id_column = "id"

    if not options.min_spike:
        options.min_spike = options.max_spike

    design_table = pd.read_table(options.design_file, sep="\t")

    # get dictionaries to map group members to column names

    groups,  g_to_keep_tracks, g_to_spike_tracks = groupMappers(
        design_table, options.spike_regex,
        options.shuffle_suffix, options.keep_suffix)

    change_bins = np.arange(options.min_cbin, options.max_cbin,
                            options.width_cbin)
    # JJ
    E.info("Column boundaries are: %s" % str(change_bins))

    initial_bins = np.arange(options.min_ibin, options.max_ibin,
                             options.width_ibin)

    # JJ
    E.info("Row boundaries are: %s" % str(initial_bins))

    if options.spike_type == "cluster":
        E.info("looking for clusters...")
        clusters_dict = findClusters(df_sort, options.cluster_max_distance,
                                     options.cluster_min_size,
                                     g_to_spike_tracks, groups)
        if len(clusters_dict) == 0:
            raise Exception("no clusters were found, check parameters")

        E.info("repeatedly shuffling subcluster regions...")
        output_indices, counts = shuffleCluster(
            initial_bins, change_bins, g_to_spike_tracks, groups,
            options.difference, options.max_spike,
            options.iterations, clusters_dict,
            options.max_sbin, options.min_sbin, options.width_sbin)

    elif options.spike_type == "row":
        E.info("repeatedly shuffling rows...")
        output_indices, counts = shuffleRows(
            df, initial_bins, change_bins, g_to_spike_tracks, groups,
            options.difference, options.max_spike, options.iterations)

    filled_bins = thresholdBins(output_indices, counts, options.min_spike)
    if len(filled_bins) == 0:
        raise Exception("No bins contained a sufficient number of spike-ins")

    E.info("outputting spike-ins...")
    outputSpikes(df_sort, options.id_column, filled_bins, g_to_keep_tracks,
                 groups, method=options.output_method,
                 spike_type=options.spike_type,
                 min_cbin=options.min_cbin, width_cbin=options.width_cbin,
                 min_ibin=options.min_ibin, width_ibin=options.width_ibin,
                 min_sbin=options.min_sbin, width_sbin=options.width_sbin)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
