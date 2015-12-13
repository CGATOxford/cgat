'''runExpression.py - wrap various differential expression tools
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script provides a convenience wrapper for differential expression
analysis for a variety of methods.

The aim of this script is to provide a common tabular output format
that is consistent between the different methods.

The script will call the selected method and output a variety of
diagnostic plots. Generally, the analysis aims to follow published
workflows for the individual method together with outputting diagnostic
plots to spot any problems. The script will also preprocess count
data to apply some common filtering methods.

The methods implemented are:

deseq2
   Application of DESeq2

deseq
   Application of DESeq

edger
   Application of EdgeR

ttest
   Application of Welch's ttest to FPKM values

mock
   A mock analysis. No differential analysis is performed,
   but fold changes are computed and output.

In addition, the script can process count data to output
modified couns or summary statistics.

summary
   Output summary statistics from a counts table.

dump
   Output the counts table after applying filtering.

spike
   Output a counts table with in-silico spike-ins. The spike-ins
   can be used to empirically check the power of any of the testing
   methods.

Usage
-----

Input
+++++

The input to this script is a table with of measurements reflecting
expression levels. For the tag counting methods such as DESeq or
EdgeR, these should be the raw counts, while for other methods such as
ttest, these can be normalized values such as FPKM values.

The script further requires a design table describing the tests to
be performed. The design table has for columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

These files should be tab separated as this is enforced in downstream
analyses and will cause the script to error.

track
     name of track - should correspond to column header in the counts
     table.
include
     flag to indicate whether or not to include this data
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.

Output
++++++

The script outputs a table with the following columns:

+------------------+------------------------------------------------------+
|*Column name*     |*Content*                                             |
+------------------+------------------------------------------------------+
|test_id           |Name of the test (gene name, ...                      |
+------------------+------------------------------------------------------+
|treatment_name    |Name of the treatment condition                       |
+------------------+------------------------------------------------------+
|treatment_mean    |Estimated expression value for treatment              |
+------------------+------------------------------------------------------+
|treatment_std     |Standard deviation                                    |
+------------------+------------------------------------------------------+
|control_name      |Name of the control condition                         |
+------------------+------------------------------------------------------+
|control_mean      |Estimated expression value for control                |
+------------------+------------------------------------------------------+
|control_std       |Standard deviation                                    |
+------------------+------------------------------------------------------+
|pvalue            |The p value for rejecting the null hypothesis         |
+------------------+------------------------------------------------------+
|qvalue            |Multiple testing correction                           |
+------------------+------------------------------------------------------+
|l2fold            |log2 foldchange of treatment/control                  |
+------------------+------------------------------------------------------+
|transformed_l2fold|a transformed log2 foldchange value.                  |
+------------------+------------------------------------------------------+
|fold              |foldchange of treatment/control                       |
+------------------+------------------------------------------------------+
|significant       |Flag, 1 if test called significant according to FDR   |
+------------------+------------------------------------------------------+
|status            |test status (OK|FAIL)                                 |
+------------------+------------------------------------------------------+

Additional plots and tables are generated and method specific.

Command line options
--------------------

'''

import sys
import os
import tempfile
from rpy2.robjects import r as R
import rpy2.robjects.numpy2ri

import CGAT.Experiment as E
import CGAT.Expression as Expression


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--tags-tsv-file", dest="input_filename_tags",
                      type="string",
                      help="input file with tag counts [default=%default].")

    parser.add_option(
        "--result-tsv-file", dest="input_filename_result",
        type="string",
        help="input file with results (for plotdetagstats) "
        "[default=%default].")

    parser.add_option("-d", "--design-tsv-file", dest="input_filename_design",
                      type="string",
                      help="input file with experimental design "
                      "[default=%default].")

    parser.add_option("-o", "--outfile", dest="output_filename", type="string",
                      help="output filename [default=%default].")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=(
                          "deseq", "edger", "deseq2",
                          "ttest",
                          "mock", "summary",
                          "dump", "spike",
                          "plottagstats",
                          "plotdetagstats"),
                      help="differential expression method to apply "
                      "[default=%default].")

    parser.add_option("--deseq-dispersion-method",
                      dest="deseq_dispersion_method",
                      type="choice",
                      choices=("pooled", "per-condition", "blind"),
                      help="dispersion method for deseq [default=%default].")

    parser.add_option("--deseq-fit-type", dest="deseq_fit_type", type="choice",
                      choices=("parametric", "local"),
                      help="fit type for deseq [default=%default].")

    parser.add_option("--deseq-sharing-mode",
                      dest="deseq_sharing_mode",
                      type="choice",
                      choices=("maximum", "fit-only", "gene-est-only"),
                      help="deseq sharing mode [default=%default].")

    parser.add_option(
        "--edger-dispersion",
        dest="edger_dispersion", type="float",
        help="dispersion value for edgeR if there are no replicates "
        "[default=%default].")

    parser.add_option("-f", "--fdr", dest="fdr", type="float",
                      help="fdr to apply [default=%default].")

    parser.add_option("-p", "--pseudocounts", dest="pseudo_counts",
                      type="float",
                      help="pseudocounts to add for mock analyis "
                      "[default=%default].")

    parser.add_option("-R", "--output-R-code", dest="save_r_environment",
                      type="string",
                      help="save R environment [default=%default].")

    parser.add_option("-r", "--reference-group", dest="ref_group",
                      type="string",
                      help="Group to use as reference to compute "
                      "fold changes against [default=$default]")

    parser.add_option("--filter-min-counts-per-row",
                      dest="filter_min_counts_per_row",
                      type="int",
                      help="remove rows with less than this "
                      "number of counts in total [default=%default].")

    parser.add_option("--filter-min-counts-per-sample",
                      dest="filter_min_counts_per_sample",
                      type="int",
                      help="remove samples with a maximum count per sample of "
                      "less than this number   [default=%default].")

    parser.add_option("--filter-percentile-rowsums",
                      dest="filter_percentile_rowsums",
                      type="int",
                      help="remove percent of rows with "
                      "lowest total counts [default=%default].")
    parser.add_option("--deseq2-design-formula",
                      dest="model",
                      type="string",
                      help="Design formula for DESeq2")
    parser.add_option("--deseq2-contrasts",
                      dest="contrasts",
                      type="string",
                      help=("contrasts for post-hoc testing writen"
                            " variable:control:treatment,..."))
    parser.add_option("--deseq2-plot",
                      dest="plot",
                      type="int",
                      help=("draw plots during deseq2 analysis"))

    parser.set_defaults(
        input_filename_tags=None,
        input_filename_result=None,
        input_filename_design=None,
        output_filename=sys.stdout,
        method="deseq",
        fdr=0.1,
        deseq_dispersion_method="pooled",
        deseq_fit_type="parametric",
        deseq_sharing_mode="maximum",
        edger_dispersion=0.4,
        ref_group=None,
        save_r_environment=None,
        filter_min_counts_per_row=1,
        filter_min_counts_per_sample=10,
        filter_percentile_rowsums=0,
        pseudo_counts=0,
        spike_foldchange_max=4.0,
        spike_expression_max=5.0,
        spike_expression_bin_width=0.5,
        spike_foldchange_bin_width=0.5,
        spike_max_counts_per_bin=50,
        model=None,
        plot=1
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    if options.input_filename_tags == "-":
        fh = tempfile.NamedTemporaryFile(delete=False)
        fh.write("".join([x for x in options.stdin]))
        fh.close()
        options.input_filename_tags = fh.name
    else:
        fh = None

    # load tag data and filter
    if options.method in ("deseq2", "deseq", "edger", "mock", "ttest"):
        assert options.input_filename_tags and os.path.exists(
            options.input_filename_tags)
        assert options.input_filename_design and os.path.exists(
            options.input_filename_design)

        Expression.loadTagData(options.input_filename_tags,
                               options.input_filename_design)

        nobservations, nsamples = Expression.filterTagData(
            filter_min_counts_per_row=options.filter_min_counts_per_row,
            filter_min_counts_per_sample=options.filter_min_counts_per_sample,
            filter_percentile_rowsums=options.filter_percentile_rowsums)

        if nobservations == 0:
            E.warn("no observations - no output")
            return

        if nsamples == 0:
            E.warn("no samples remain after filtering - no output")
            return

        sample_names = R('''colnames(countsTable)''')
        E.info("%i samples to test at %i observations: %s" %
               (nsamples, nobservations,
                ",".join(sample_names)))

    try:
        if options.method == "deseq2":
            Expression.runDESeq2(
                outfile=options.output_filename,
                outfile_prefix=options.output_filename_pattern,
                fdr=options.fdr,
                ref_group=options.ref_group,
                model=options.model,
                contrasts=options.contrasts,
                plot=options.plot
            )

        elif options.method == "deseq":
            Expression.runDESeq(
                outfile=options.output_filename,
                outfile_prefix=options.output_filename_pattern,
                fdr=options.fdr,
                dispersion_method=options.deseq_dispersion_method,
                fit_type=options.deseq_fit_type,
                sharing_mode=options.deseq_sharing_mode,
                ref_group=options.ref_group,
            )

        elif options.method == "edger":
            Expression.runEdgeR(
                outfile=options.output_filename,
                outfile_prefix=options.output_filename_pattern,
                fdr=options.fdr,
                ref_group=options.ref_group,
                dispersion=options.edger_dispersion)

        elif options.method == "mock":
            Expression.runMockAnalysis(
                outfile=options.output_filename,
                outfile_prefix=options.output_filename_pattern,
                ref_group=options.ref_group,
                pseudo_counts=options.pseudo_counts,
            )

        elif options.method == "summary":
            Expression.outputTagSummary(
                options.input_filename_tags,
                options.stdout,
                options.output_filename_pattern,
                filename_design=options.input_filename_design
            )

        elif options.method == "dump":
            assert options.input_filename_tags and os.path.exists(
                options.input_filename_tags)
            Expression.dumpTagData(options.input_filename_tags,
                                   options.input_filename_design,
                                   outfile=options.stdout)

        elif options.method == "plottagstats":
            assert options.input_filename_tags and os.path.exists(
                options.input_filename_tags)
            Expression.plotTagStats(
                options.input_filename_tags,
                options.input_filename_design,
                outfile_prefix=options.output_filename_pattern)

        elif options.method == "plotdetagstats":
            assert options.input_filename_result and os.path.exists(
                options.input_filename_result)
            Expression.plotDETagStats(
                options.input_filename_result,
                outfile_prefix=options.output_filename_pattern)

        elif options.method == "spike":
            Expression.outputSpikeIns(
                options.input_filename_tags,
                options.stdout,
                options.output_filename_pattern,
                filename_design=options.input_filename_design,
                foldchange_max=options.spike_foldchange_max,
                expression_max=options.spike_expression_max,
                max_counts_per_bin=options.spike_max_counts_per_bin,
                expression_bin_width=options.spike_expression_bin_width,
                foldchange_bin_width=options.spike_foldchange_bin_width,
            )

        elif options.method == "ttest":
            Expression.runTTest(
                outfile=options.output_filename,
                outfile_prefix=options.output_filename_pattern,
                fdr=options.fdr)

    except rpy2.rinterface.RRuntimeError:
        if options.save_r_environment:
            E.info("saving R image to %s" % options.save_r_environment)
            R['save.image'](options.save_r_environment)
        raise

    if fh and os.path.exists(fh.name):
        os.unlink(fh.name)

    if options.save_r_environment:
        R['save.image'](options.save_r_environment)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
