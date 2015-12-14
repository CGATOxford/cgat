'''
expression2resample.py
=============================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Generates a distance matrix using one of three distance measures for time
series data.

The three distance measure options are dynamic time warping (not strictly a
metric), normalised cross-correlation or temporal correlation.

Dynamic time warping
====================

Dynamic time warping find the optimal alignment path between two time series
based on the sum of their Euclidean differences, allowing for constraints of
head and tail matching, but allowing for warping between any other time points.
It is not a distance metric in the sense that it does not always meet the
triangle inequality, however, it has been applied to hierarchical clustering
of time series in a number of applications, primarily those involving
signal processing.

Temporal correlation
====================

The temporal correlation measure and optional adaptive tuning function
are both derived from::
  Chouakria & Nagabhushan ADAC 2007, pp5-21

This correlates time series with respect to their behaviour, and is somewhat
analagous to template matching of time series expression profiles.

Cross correlation
=================

The normalised cross-correlation uses the numpy correlate function, normalised
for length of the leading time series.


Input is a single time-series expression data set with no replicates of
 normalised expression values.  If data are derived from RNAseq counts, it is
 recommended these are transformed on to an approximately normal scale using
 functions such as variance stabilising transformation, log transformation or
 other functions appropriate to the data.

The input file must be tab-delimited, and time points must be in order.

TODO:
* rename script to something more accurate: resample2distance.py


Options
-------
The main functionality of this script is determined by the distance metric
used.  Other options are used to adjust the behaviour of these measures
or are for use with a high performance computing cluster.

  --distance-metric - distance measure to use.  DTW is a dissimilarity
                      measure, whilst temporal correlation and cross
                      correlation are both similarity measures

  --expression-file - file containing gene expression data, with a single gene
                      along each row and the time points ordered along the
                      columns.  See example for details.

  --k - degree to which adaptive tuning will be applied.  This balances
        the dynamic time warping with the temporal correlation to give
        a dissimilarity measure that accounts for both the magnitude
        and behaviour of each time series expression profile.  See the paper
        by Chouakria & Nagabhushan for more details.

  --lag - lag to report for the cross correlation.  The default is 0.
          Negative numbers represent gene 'X' delaying with respect
          to 'Y', and positive values represent 'X' leading with respect
          to gene 'Y'.

  --parallel - this will split the distance expression data into chunks
               based on the co-ordinate provided within the file name.
               E.g.
               condition-R1-0_499-expression.tsv
               This will only calculate the distance matrix for all
               genes against genes 0-499 inclusive (0-based indexing).

  --out - output filename

Usage
-----

Example::

   head condition-expression.tsv
   gene_id             0      1      3      6      12      24      48  ...
   ENSMUSG00000001305  10.06  10.57  13.28  13.46  13.18   12.78   11.83
   ENSMUSG00000001674  11.97  11.99  14.37  14.32  14.011  13.57   12.76
   ENSMUSG00000003134  6.07   5.59   4.29   5.74   6.00    6.23    6.90
   ENSMUSG00000004110  13.15  13.54  12.57  8.6    5.56    5.16    7.50
   ...


   python expression2resample.py
   --distance-metric=dtw
   --k=0
   --out=distance_matrix.tsv
   condition-expression.tsv

   head distance_matrix.tsv

                      ENSMUSG00000001305    ENSMUSG00000001674   ...

  ENSMUSG00000001305  0.0                   27.35
  ENSMUSG00000001674  27.35                 0.0
  ...


Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import pandas as pd
import CGAT.Experiment as E
import CGAT.Timeseries as TS


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--k", dest="k", type="int", default=0,
                      help="value of k to adjust adaptive tuning function")

    parser.add_option("--out", dest="outfile", type="string",
                      help="output file name")

    parser.add_option("--expression-file", dest="expr", type="string",
                      help="file containing expression data")

    parser.add_option("--parallel", dest="parallel", action="store_true",
                      default=False, help="switches on parallel, will"
                      " split distance matrix into relevant number of"
                      " slices. Start-end positions are defined by"
                      " the file name.")

    parser.add_option("--distance-metric", dest="dist_metric", type="string",
                      help="distance metric to use for dissimilarity of time "
                      "series objects.  Choices: dtw, cross-correlate, "
                      "temporal-correlate. Default=dtw")

    parser.add_option("--lag", dest="lag", type="string",
                      help="cross correlation lag to report")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = args[-1]

    parser.set_defaults(lag=0,
                        k=0)
    if options.parallel:
        datfile = options.expr
    else:
        datfile = infile

    data = pd.read_table(datfile,
                         sep="\t",
                         index_col=0,
                         header=0)

    # data should already be sorted in time-series order
    # the time and replicate columns needs to be dropped to ensure only the
    # gene data is passed into the DTW function

    # drop header line(s) and non-numerical rows
    try:
        data.drop(['times'], inplace=True, axis=0)
        data.drop(['replicates'], inplace=True, axis=0)
    except ValueError:
        pass
    genes = data.index
    data = data.convert_objects(convert_numeric=True)

    # iterate over the genes list in nested loops to get
    # all pair-wise combinations.

    if options.dist_metric == "dtw":
        if options.parallel:
            start_idx = int(infile.split("-")[3].split("_")[0])
            end_idx = int(infile.split("-")[3].split("_")[1])
            slice_idx = genes[start_idx:end_idx]

            df_ = TS.dtwWrapper(data=data,
                                rows=genes,
                                columns=slice_idx,
                                k=options.k)
        else:
            df_ = TS.dtwWrapper(data=data,
                                rows=genes,
                                columns=genes,
                                k=options.k)

    elif options.dist_metric == "cross-correlate":

        if options.lag is None:
            options.lag = 0
        else:
            pass

        if options.parallel:
            start_idx = int(infile.split("/")[-1].split("-")[3].split("_")[0])
            end_idx = int(infile.split("/")[-1].split("-")[3].split("_")[1])
            slice_idx = genes[start_idx:end_idx]

            df_ = TS.correlateDistanceMetric(data=data,
                                             rows=genes,
                                             columns=slice_idx,
                                             method=options.dist_metric,
                                             lag=int(options.lag))
        else:
            df_ = TS.correlateDistanceMetric(data=data,
                                             rows=genes,
                                             columns=genes,
                                             method=options.dist_metric,
                                             lag=int(options.lag))

    elif options.dist_metric == "temporal-correlate":
        if options.parallel:
            start_idx = int(infile.split("/")[-1].split("-")[3].split("_")[0])
            end_idx = int(infile.split("/")[-1].split("-")[3].split("_")[1])
            slice_idx = genes[start_idx:end_idx]

            df_ = TS.correlateDistanceMetric(data=data,
                                             rows=genes,
                                             columns=slice_idx,
                                             method=options.dist_metric)
        else:
            df_ = TS.correlateDistanceMetric(data=data,
                                             rows=genes,
                                             columns=genes,
                                             method=options.dist_metric)

    if not options.outfile:
        df_.to_csv(options.stdout, sep="\t")
    else:
        df_.to_csv(options.outfile, sep="\t")

    # write footer and output benchmark information
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
