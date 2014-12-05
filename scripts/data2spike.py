"""
data2spike.py - generate spike-in from a data frame
=====================================================

:Author: Tom smith
:Release: $Id: data2spike.py 
:Date: |today|
:Tags: Genomics Statistics Power


Purpose
-------
This scripts can be used to generate spike-ins from a dataframe,
the primary purpose of which is to perform a power analysis.


Usage
-----

Input
+++++
The input to this script is a dataframe containing one item per row
(e.g gene, CpG).
The script further requires a design table describing the groups and columns.
The design table has for columns::

track    spike group   include
GS1_perc        1       S       1
GS2_perc        1       S       1
GD1_perc        1       D       1
GD2_perc        1       D       1
GS1_unmeth      0       S       1
GS2_unmeth      0       S       1
GD1_unmeth      0       D       1
GD2_unmeth      0       D       1
GS1_meth        0       S       1
GS1_meth        0       S       1
GS1_meth        0       D       1
GS1_meth        0       D       1
GS1_unwanted    0       S       0
GS2_unwanted    0       S       0
GD1_unwanted    0       D       0
GD2_unwanted    0       D       0


track
     name of track - should correspond to column header in the counts
     table.
spike
     flag to indicate which columns should be shuffled to generate the spik-ins

include
     flag to indicate whether or not to include this data in the output
group
     group indicator - experimental group
spike
     flag to indicate which columns should be shuffled to generate the spik-ins


MORE....


Output
++++++

The original dataframe is outputted along with the spike-ins



Command line options
--------------------

``--genome-file``
    required option. filename for genome fasta file


"""

import os
import copy
import random
import sys
import re
import optparse
import string
import collections
import array
import pandas as pd
import numpy as np

import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id: gtf2fasta.py 2861 2010-02-23 17:36:32Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-t", "--design-tsv-file", dest="design_file", type="string",
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

    parser.add_option("-a", "--id_columns", dest="id", action="append",
                      help="name of identification column(s)\
                      [default=%default].")

    parser.set_defaults(
        difference="relative",
        genome_file=None,
        min_cbin=None,
        max_cbin=None,
        bin_cwidth=None,
        min_ibin=None,
        max_ibin=None,
        bin_iwidth=None,
        id_columns=None,
        max_spike=100,
        min_spike=20,
        iterations=10)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    if not options.design_file:
        raise ValueError("a design table is required.")

    df = pd.read_table(sys.stdin, sep="\t")
    #print df.head()
    design_table = pd.read_table(options.design_file, sep="\t")
    groups = list(set([design_table['group'][ix]
                       for ix in design_table.index
                       if design_table['spike'][ix] == 1]))
    g_to_keep_ix = {}
    g_to_keep_tracks = {}
    g_to_spike_ix = {}
    g_to_spike_tracks = {}

    for group in groups:
        g_to_keep_ix[group] = [ix for ix in design_table.index
                               if re.match(group, design_table['group'][ix])
                               and design_table['include'][ix] == 1]
        g_to_keep_tracks[group] = design_table.ix[g_to_keep_ix[group], "track"]
        g_to_spike_ix[group] = [ix for ix in g_to_keep_ix[group]
                                if design_table['spike'][ix] == 1]
        g_to_spike_tracks[group] = design_table.ix[
            g_to_spike_ix[group], "track"]
    print g_to_keep_tracks[groups[0]]
    print g_to_keep_tracks[groups[1]]
    #print g_to_keep_ix
    #print g_to_keep_tracks
    #print g_to_spike_ix
    #print g_to_spike_tracks

    change_bins = [round(float(x), 2) for x in np.arange(
        options.min_cbin, options.max_cbin+options.width_cbin,
        options.width_cbin)]

    initial_bins = [round(float(x), 2) for x in np.arange(
        options.min_ibin, options.max_ibin+options.width_ibin,
        options.width_ibin)]

    output_counts = np.zeros(
        (len(initial_bins) + 1, len(change_bins) + 1))

    output_indices = {(key1, key2): []
                      for key1 in np.digitize(initial_bins, initial_bins)
                      for key2 in np.digitize(change_bins, change_bins)}
    output_indices_keys = output_indices.keys()
    #print "output_indices"
    #print output_indices

    min_occup = min(output_counts.flatten())

    for iteration in range(0, options.iterations):
        #print "min occupancy: %i" % min_occup
        if min_occup < options.max_spike:
            group1_rand = np.random.permutation(df.index)
            group2_rand = np.random.permutation(df.index)

            group1_mean = df.ix[group1_rand,
                                g_to_spike_tracks[groups[0]]].apply(
                                    np.mean, axis=1).tolist()
            group2_mean = df.ix[group2_rand,
                                g_to_spike_tracks[groups[1]]].apply(
                                    np.mean, axis=1).tolist()

            if options.difference == "relative":
                change = [group2_mean[x] - group1_mean[x]
                          for x in range(0, len(group1_mean))]
                initial = [group1_mean[x]
                           for x in range(0, len(group1_mean))]
            elif options.difference == "logfold":
                change = [np.log2((group2_mean[x]+1.0)/(group1_mean[x]+1.0))
                          for x in range(0, len(group1_mean))]
                initial = [np.log2(group1_mean[x]+1.0)
                           for x in range(0, len(group1_mean))]

            initial_idx = np.digitize(initial, initial_bins)
            change_idx = np.digitize(change, change_bins)

            for idx, coord in enumerate(zip(initial_idx, change_idx)):
                if coord in output_indices_keys:
                    if output_counts[coord] >= options.max_spike:
                        continue
                    output_counts[coord] += 1
                    output_indices[coord].append(
                        (group1_rand[idx], group2_rand[idx]))
        min_occup = min(output_counts.flatten())

    #print "min occupancy: %i" % min_occup

    output_indeces_keep = copy.copy(output_indices)
    for key in output_indices:
        if len(output_indices[key]) < options.min_spike:
            output_indeces_keep.pop(key)

    #print output_indices
    #print output_indeces_keep
    #print len(output_indices.keys())
    #print len(output_indeces_keep.keys())

    header = ["id"]
    #print g_to_keep_tracks[groups[0]]
    header.extend(g_to_keep_tracks[groups[0]])
    header.extend(g_to_keep_tracks[groups[1]])
    sys.stdout.write("%s\n" % "\t".join(map(str, header)))

    n = 0
    for key in output_indeces_keep:
        for pair in output_indeces_keep[key]:
            initial_bin, change_bin = key
            initial = ((initial_bin*options.width_ibin)+options.min_ibin)
            change = ((change_bin*options.width_cbin)+options.min_cbin)
            row = ["_".join(map(str, ("spike-in", initial, change, n)))]
            row.extend(df.ix[pair[0], g_to_keep_tracks[groups[0]]])
            row.extend(df.ix[pair[1], g_to_keep_tracks[groups[1]]])
            sys.stdout.write("%s\n" % "\t".join(map(str, row)))
            n += 1
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
