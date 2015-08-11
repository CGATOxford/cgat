'''
timeseries_deseq2.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Perform differential gene expression analysis using DESeq2 comparing
gene expression changes over two time points between 2 conditions or
between baseline and each time point (specified by the ``--method`` option).

Options
-------
Column names of input tables must be in the form condition.time.replicate

  --results-directory - directory to write all results/output files to

  --alpha - p-value threshold for statistical significance

  --method - analysis type.  Either ``timepoint`` or ``condition``.
    This will determine the design matrix and design formula.

    The standard form for the design formula for the condition analysis is:
      ~ replicates + times + condition + times:condition

    For the time point analysis the formula is:
      ~ replicates + condition + times

    The last term in the formula is the variable for which significance
    testing is performed.

  All options are mandatory.

Usage
-----

This script uses DESeq2 to call differentially expressed genes, either between
two time points (or between the baseline and each time point supplied), or
between two conditions for a given time point.  In the latter of these designs
the interaction between the condition and time point is the coefficient of
interest.

Output filenames are derived from the condition and the two time points, or
the test_reference conditions followed by time points.

For example::
  time point differential analysis will look like this::
    CD40L.0_3-time.tsv, CD40L.0_6-time.tsv, etc
  and condition differential analysis output will look like this::
    CD40L_IgM.0_3-diff-cond.tsv, CD40L_IgM.0_6-diff-cond.tsv, etc

  Output files are::
    [prefix]-MAplot.png - png image of MA plot differential expression analysis
    [prefix]-dispersions.png - mean-variance dependence plots with scaling
                               according to Bayes shrinkage (see DESeq2 manual)
    [prefix]-time/diff-cond.tsv - tab-delimited file of all results

Example::

   head counts_table.tsv

   gene_id       CD40L.000.R1  CD40L.000.R2  CD40L.000.R3  ...
   LNCGme00001   24            52            26
   LNCGme00002   22            25            10
   LNCGme00005   47            42            26
   LNCGme00009   150           204           120
   LNCGme00012   1419          1587          713

   cat counts_table.tsv |
   python timeseries_deseq2.py --method=timepoint
                               --alpha=0.01
                               --results-directory=timepoints.dir

   head CD40L.0_1-time.tsv

   gene_id              baseMean  log2FoldChange  lfcSE   stat    pvalue  padj
   LNCGse05190          20384.92   7.172          0.119   59.948  0.0     0.0
   ENSMUSG00000025981   12114.84   4.725          0.103   45.823  0.0     0.0
   ENSMUSG00000015312   59210.56   6.565          0.095   68.980  0.0     0.0
   ENSMUSG00000019850   165722.12  5.922          0.095   61.834  0.0     0.0


Type::

   python timeseries_deseq2.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import pandas as pd
import itertools
import re
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

    parser.add_option("--results-directory", dest="res_dir",
                      type="string", help="directory to write results"
                      "tables to")

    parser.add_option("--alpha", dest="alpha", type="string",
                      help="statistical significance p-value threshold")

    parser.add_option("--method", dest="method", type="string",
                      help="analysis design. "
                      "either timepoint or condition")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    try:
        infile = argv[-1]
        open(infile, "r")
        # check for compression state
        if infile.split(".")[-1] == "gz":
            comp = "gzip"
        else:
            comp = None

    except IOError:
        infile = options.stdin
        # check for compression state
        if infile.name.split(".")[-1] == "gz":
            comp = "gzip"
        else:
            comp = None

    alpha = float(options.alpha)
    res_dir = options.res_dir

    count_table = pd.read_table(infile,
                                sep="\t",
                                index_col=0,
                                header=0,
                                compression=comp)
    columns = count_table.columns
    conditions = set([x.split(".")[0] for x in columns])
    times = set([x.split(".")[1] for x in columns])

    data_dict = {}
    cond_times = [x for x in itertools.product(conditions, times)]
    base_col = {}
    time_dict = {}

    if options.method == "timepoint":

        # assumes all column names are in the form
        # `condition`:`time`:`replicate`
        # use `condition`.`time` as dictionary keys

        for x in cond_times:
            c_t = "%s.%s" % (x[0], x[1])
            cols = [k for k in count_table.columns if re.search(c_t, k)]
            if x[1] == '000':
                base_col[c_t] = count_table[cols]
            else:
                time_dict[c_t] = count_table[cols]

        for bt in itertools.product(base_col.keys(),
                                    time_dict.keys()):
            df = pd.merge(left=base_col[bt[0]],
                          right=time_dict[bt[1]],
                          how='outer',
                          left_index=True,
                          right_index=True)
            time = int(bt[1].split(".")[1])
            data_dict["%s_0_%i" % (bt[0].split(".")[0],
                                   time)] = df

        for each in data_dict.keys():
            df_ = data_dict[each]
            outfile = "%s/%s-time.tsv" % (res_dir,
                                          each)
            res_frame = TS.timepointDESeq2(df_,
                                           each,
                                           alpha,
                                           res_dir)
            res_frame.to_csv(outfile,
                             sep="\t",
                             index_label="gene_id")

    elif options.method == "condition":

        # assumes all column names are in the form
        # `condition`:`time`:`replicate`
        # use `condition`.`time` as dictionary keys

        for x in cond_times:
            c_t = "%s.%s" % (x[0], x[1])
            cols = [k for k in count_table.columns if re.search(c_t, k)]
            if int(x[1]) == 0:
                base_col[c_t] = count_table[cols]
            else:
                time_dict[c_t] = count_table[cols]

        # make a dataframe for each 0:time point combination
        # for all conditions, index on `condition:0_time`

        base_keys = base_col.keys()
        time_keys = time_dict.keys()
        for k in conditions:
            for x in itertools.product(base_keys, time_keys):
                if re.search(k, x[0]) and re.search(k, x[1]):
                    df = pd.merge(left=base_col[x[0]],
                                  right=time_dict[x[1]],
                                  how='outer',
                                  left_index=True,
                                  right_index=True)
                    time = int(x[1].split(".")[1])
                    data_dict["%s.0_%i" % (x[0].split(".")[0],
                                           time)] = df
                else:
                    pass

        time_span = set([x.split(".")[1] for x in data_dict.keys()])

        all_dict = {}
        for cond in itertools.combinations(conditions, 2):
            c1 = cond[0]
            c2 = cond[1]
            for x in time_span:
                key1 = "%s.%s" % (c1, x)
                key2 = "%s.%s" % (c2, x)
                df = pd.merge(left=data_dict[key1],
                              right=data_dict[key2],
                              how='outer',
                              left_index=True,
                              right_index=True)
                all_dict["%s_%s.%s-diff" % (c1, c2, x)] = df

        for each in all_dict.keys():

            df = all_dict[each]
            outfile = "%s/%s-cond.tsv" % (res_dir,
                                          each)
            res_frame = TS.conditionDESeq2(df,
                                           each,
                                           alpha,
                                           res_dir)
            res_frame.to_csv(outfile, sep="\t", index_label="gene_id")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
