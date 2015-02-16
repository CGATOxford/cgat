'''
clustering.py - clustering metrics
====================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Provide clustering quality metrics.  Clustering metrics are:
* Rand index
* Adjusted Rand index
* clustering agreement rate
* cluster recall
* F-measure (with beta=1) for balance of recall/agreement
* Jaccard coefficient
* Adjusted mutual information

All metrics are calculated on a table of multiple clusterings (n > 1).

Options
-------

This script has no additional options.  All functions are executed on running.

Usage
-----

head clustering_results.tsv

gene_id              cluster   cluster   cluster
ENSMUSG00000001138   grey      grey      grey
ENSMUSG00000001305   grey      grey      grey
ENSMUSG00000003134   grey      grey      grey
ENSMUSG00000004110   grey      grey      grey
ENSMUSG00000004451   grey      grey      grey

python clustering.py clustering_results.tsv > metrics.tsv

head metrics.tsv

idx  AMI  AdjRand  F_measure  Jaccard  Positive  Rand  Recall  Total
0    0.0  0.019    0.5        0.33     0.5       0.34  0.5     3577.0
1    0.0  0.019    0.5        0.33     0.5       0.34  0.5     3577.0
2    0.0  0.019    0.5        0.33     0.5       0.34  0.5     3577.0

Type::

   python clustering.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import pandas as pd
import itertools
from math import factorial
import CGATPipelines.PipelineTimeseries as TS


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]
    E.info("loading input file: %s" % infile)
    assert infile

    df = pd.read_table(infile,
                       sep="\t",
                       header=None,
                       index_col=0)
    cluster_combs = (x for x in itertools.combinations(df.columns,
                                                       2))
    genes = df.index
    results_dict = {}
    all_clusters = {}
    for i in df.columns:
        clusters = set(df[i].values.tolist())
        cluster_dict = {}
        for clust in clusters:
            cluster_dict[clust] = []
        for gene in genes:
            cluster_dict[df[i][gene]].append(gene)

        for col in clusters:
            col_set = set()
            clust_col = cluster_dict[col]
            gene_members = itertools.combinations(clust_col,
                                                  2)
            col_set.update(gene_members)
            cluster_dict[col] = col_set
            all_clusters[i] = cluster_dict

    gene_len = len(genes)
    gene_factorial = factorial(gene_len)
    complete_pairs = int((gene_factorial/factorial(2))/factorial(gene_len - 2))

    E.info("generating all pair-wise cluster comparisons")
    for k in cluster_combs:
        E.info("calculating metrics for %s and %s" % (k[0],
                                                      k[1]))
        concord = TS.clusterConcordia(all_clusters[k[0]],
                                      all_clusters[k[1]],
                                      complete_pairs)
        metric_dict = TS.concordanceMetric(concord)
        metric_dict['AMI'] = TS.adjustedMutualInformation(all_clusters[k[0]],
                                                          all_clusters[k[1]])
        results_dict[k] = metric_dict

    E.info("aggregating results")
    res_frame = pd.DataFrame(results_dict).T
    res_frame = res_frame.reset_index()
    res_frame.drop(['level_0'], inplace=True, axis=1)
    res_frame.drop(['level_1'], inplace=True, axis=1)
    res_frame.to_csv(options.stdout,
                     sep="\t",
                     index_label='idx')

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
