'''expression2expression.py - RNAseq filtering and normalisation
=============================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Script for normalisation and transformation of RNAseq timeseries data.

Script functions are:
* Normalise RNAseq libraries using DESeq size factors
* Transform counts tables by variance stablising transformation
* filter expression data on sum of absolute covariance
* fit polynomial model to data with significance testing
* calculate average expression over replicates for transformed time series
expression data

Usage
-----

::
   # applying library size normalisation and variance stabilising transformation
   head input_counts.tsv

   # zero-values and genes with mean counts <10 are filtered out
   gene_id            LPS.000.R1  LPS.000.R2  LPS.000.R3  LPS.001.R1  LPS.001.R2
   ENSMUSG00000000544 0           0           0           0           0
   ENSMUSG00000000817 13          8           12          5           20
   ENSMUSG00000001138 3878        3868        2484        1340        1266
   ENSMUSG00000001143 2891        2668        1573        931         1041

   python timeseries_analysis.py --task=deseq --replicates=R1,R2,R3
             --time=0,1,3,6,12,24,48 input_counts.tsv > output_file.tsv

   head output_file.tsv

   gene_id            LPS.000.R1  LPS.000.R2  LPS.000.R3  LPS.001.R1  LPS.001.R2
   ENSMUSG00000001138 11.20       11.06       11.16       10.73       10.35
   ENSMUSG00000001143 10.78       10.53       10.50       10.21       10.07


Input files is provided as the final argument to the script.  The exact form
 and nature of the input file is dependent on the task to be performed.

Options
-------

Each set of options is dependent on the task applied.  Not all options
are appropriate for all tasks.  Please see details below for the exact
intended usage of these options.

task::

  The task to perform on the input.  This will determine which options will
  be functional.

  Choices are::
    ``deseq`` - perform library size normalisation and variance
                stabilising transformation on expression count data.
                Active options:
                 --time - comma separated list of time points measured,
                   in numerical order
                 --replicates - comma separated list of labels associated
                   with each replicate

    ``masigpro`` - fit polynomial model to timeseries data to select
                   differentially expressed genes across the timeseries.
                   See maSigPro documentation for more details.
                   Active options::
                    --orders - order of polynomial model to fit
                    --fdr - false discovery rate for calling DEGs
                    --padjust - multiple testing correction to apply
                    --stepwise - stepwise regression to apply. forwards
                      or backwards
                    --pinclude - largest p-value required for inclusion
                      in the model by stepwise regression
                    --rsquared - variance explained cut-off for DEGs
                    --vargroup - variable group to report.  Either `each`,
                      all` or `group`.

    ``sumcovar`` - filter genes on the sum of their absolute covariance
                   Active options::
                    --reps - comma-separated list of replicate IDs
                    --timepoints - comma-separated list of time points
                    --quantile - quantile threshold to apply for filtering

    ``average_expression`` - calculate average expression over replicates
                             at each time point


Type::

   python expression2expression.py --help

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

    parser.add_option("--time", dest="timepoints", type="string",
                      help="a comma-separated list of time points measured")

    parser.add_option("--replicates", dest="reps", type="string",
                      help="a comma-separated list of replicate IDs")

    parser.add_option("--conditions", dest="conditions", type="string",
                      help="a comma-separated list of experimental conditions")

    parser.add_option("--orders", dest="orders", type="int",
                      help="order of polynomial terms to include in"
                      "maSigPro linear model")

    parser.add_option("--fdr", dest="fdr", type="string",
                      help="FDR for calling DEGs")

    parser.add_option("--padjust", dest="padjust", type="string",
                      help="multiple testing correction to apply to"
                      "control FDR")

    parser.add_option("--stepwise", dest="stepwise", type="string",
                      help="stepwise regression to use")

    parser.add_option("--pinclude", dest="pinclude", type="string",
                      help="p-value for inclusion in stepwise regression")

    parser.add_option("--rsquared", dest="rsquared", type="string",
                      help="rsquared cut-off for DEG reporting")

    parser.add_option("--var-group", dest="vargroup", type="string",
                      help="variable group reporting. each, all or"
                      "group")

    parser.add_option("--task", dest="task", type="string",
                      help="analysis task to be executed")

    parser.add_option("--infile", dest="infile", type="string",
                      help="input file path")

    parser.add_option("--quantile", dest="quantile", type="int",
                      help="see pipeline.ini for explanation")

# add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    parser.set_defaults(cutHeight=0,
                        conditions=None,
                        split=False,
                        cluster_size=30)

    if options.task == "deseq":
        timepoints = [int(x) for x in options.timepoints.split(",")]
        timepoints.sort()
        reps = [x for x in options.reps.split(",")]
        if not options.conditions:
            conditions = None
        else:
            conditions = [x for x in options.conditions.split(",")]

        data_frame = TS.deseqNormalize(infile=infile,
                                       time_points=timepoints,
                                       reps=reps,
                                       conditions=conditions)

    elif options.task == "masigpro":
        data_frame = TS.maSigPro(infile=infile,
                                 order_terms=int(options.orders),
                                 fdr=float(options.fdr),
                                 adjust=options.padjust,
                                 stepwise=options.stepwise,
                                 include_p=float(options.pinclude),
                                 rsq=float(options.rsquared),
                                 var_group=options.vargroup)

    elif options.task == "sumcovar":
        timepoints = [int(x) for x in options.timepoints.split(",")]
        reps = [x for x in options.reps.split(",")]
        data_frame = TS.covarFilter(infile=infile,
                                    time_points=timepoints,
                                    replicates=reps,
                                    quantile=int(options.quantile))

    elif options.task == "average_expression":
        data_frame = TS.avTimeExpression(infile)

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
