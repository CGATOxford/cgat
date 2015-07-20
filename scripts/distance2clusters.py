'''
distance2clusters.py - RNAseq timeseries analysis
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
Script for clustering analysis of distance matrices.

Script functions are:
* combine multiple clusterings by averaging distance matrix or consensus of
repeated clusterings
* perform consensus clustering on average or aggregated distance matrix
* calculate PCA and eigengene for each cluster

Usage
-----

TODO: insert usage example

Options
-------

Each set of options is dependent on the task applied.  Not all options are
 appropriate for all tasks.  Please see details below for the exact intended
 usage of these options.

task::
  The task to perform on the input.  This will determine which options will
  be functional.

  Choices are::

    ``cluster`` - assign gene ids to clusters using the dynamic tree cutting
                  algorithm.  See WGCNA paper for details.
                  Active options::
                   --expression-file - file containing gene expression data
                   --cluster-file - file to output cluster labels to
                   --cluster-algorithm - hierarchical clustering algorithm
                     to apply
                   --split-clusters - apply deep splitting of dendrogram
                     after cutting.
                   --cluster-size - minimum number of items in a cluster
                     after tree cutting.
    ``clustagree`` - perform consensus clustering on either an averaged
                     distance matrix or across multiple clustering runs
                     Active options::
                     --method::
                       method applied to previous distance matrix calculation.
                       Either ``resample`` or ``replicate``.  ``resample`` will
                       perform consensus clustering over cluster labels and
                       cluster assignments. ``replicate`` will calculate the
                       average distance matrix across replicates.

    ``consensus-cluster`` - cut dendrogram from consensus clustering and
                            assign genes to discrete clusters.
                            Active options::
                            --cut-height - if not 0, position at which to cut
                              the dendrogram to form clusters.
                            --cluster-algorithm - agglomerative clustering
                              algorithm to use for tree building.  Choices
                              are single, average, complete or Wards.
                            --cluster-size - minimum number of genes in a
                              after tree cutting.  Clusters smaller than this
                              are merged with adjacent clusters.
                            --split-clusters - apply deep splitting of
                            dendrogram for large clusters.

    ``pca`` - perform principal components analysis on gene expression data
              within clusters.  Output PC1, representative expression profiles
              and loadings plot for each cluster.
              Active options::
              --cluster-file - input file of cluster assignments from tree
                tree cutting or other cluster assignments.
              --image-dir - directory to save plots to

Type::

   python distance2clusters.py --help

for command line help.

Command line options
--------------------
'''

import sys
import CGAT.Experiment as E
import CGAT.Timeseries as TS


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

    parser.add_option("--task", dest="task", type="string",
                      help="analysis task to be executed")

    parser.add_option("--infile", dest="infile", type="string",
                      help="input file path")

    parser.add_option("--method", dest="method", type="choice",
                      choices=("replicate", "resample"),
                      help="whether to use replicate or resample "
                      "for consensus clustering.")

    parser.add_option("--cluster-algorithm", dest="cluster", type="string",
                      help="hierarchical clustering algorithm")

    parser.add_option("--expression-file", dest="express", type="string",
                      help="matching expression data from input"
                      " distance matrix")

    parser.add_option("--cluster-file", dest="clustfile", type="string",
                      help="file to output cluster labels to")

    parser.add_option("--output-file", dest="outfile", type="string",
                      help="output file to write to")

    parser.add_option("--cut-height", dest="cutHeight", type="string",
                      help="threshold at which to define consensus clusters"
                      "as valid")

    parser.add_option("--split-clusters", dest="split", action="store_true",
                      help="switch for using deepSplit in tree cutting")

    parser.add_option("--cluster-size", dest="cluster_size", type="int",
                      help="minimum cluster size for tree cutting. Clusters "
                      "with fewer than this many objects will be merged with "
                      "nearest cluster. Default=30")

    parser.add_option("--image-dir", dest="images_dir", type="string",
                      help="directory to write plots/figures to")

# add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    parser.set_defaults(cutHeight=0,
                        conditions=None,
                        split=False,
                        cluster_size=30)

    if options.task == "cluster":

        data_frame = TS.treeCutting(infile=infile,
                                    expression_file=options.express,
                                    cluster_file=options.clustfile,
                                    cluster_algorithm=options.cluster,
                                    deepsplit=options.split)

    elif options.task == "clustagree":
        if options.method == "resample":
            data_frame = TS.clusterAgreement(infile)
        elif options.method == "replicate":
            file_list = infile.split(",")
            data_frame = TS.clusterAverage(file_list)

    elif options.task == "consensus-cluster":
        min_size = int(options.cluster_size)
        data_frame = TS.consensusClustering(infile=infile,
                                            cutHeight=float(options.cutHeight),
                                            cluster_algorithm=options.cluster,
                                            min_size=min_size,
                                            deepsplit=options.split)

    elif options.task == "pca":
        files = infile.split(",")
        infile = files[1]
        cluster_file = files[0]
        data_frame = TS.clusterPCA(infile=infile,
                                   cluster_file=cluster_file,
                                   image_dir=options.images_dir)

    else:
        pass

    data_frame.to_csv(options.stdout,
                      sep="\t",
                      header=True,
                      index_label="gene_id")

    # Write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
