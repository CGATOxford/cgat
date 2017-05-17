'''counts2table.py - wrap various differential expression tools
=============================================================

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

sleuth
   Application of sleuth.

deseq2
   Application of DESeq2

edger
   Application of EdgeR

dexseq
   Application of DEXSeq

ttest
   Application of Welch's ttest to FPKM values

mock
   A mock analysis. No differential analysis is performed,
   but fold changes are computed and output.

Use --sleuth-genewise to test at gene rather than transcript level.
For genewise analysis, also require --gene-biomart option Use
following R code to identify the correct database, e.g
hsapiens_gene_ensembl) > library(biomaRt)
>listDatasets(useEnsembl(biomart="ensembl"))

Use the option --use-ihw to use the independent hypothesis weighting
method to calculate a weighted FDR. Note this will replace the
unweighted BH FDR in the final results table.

Usage
-----

Input
+++++

The input to this script is a table of measurements reflecting
expression levels. For the tag counting methods such as DESeq2 or
EdgeR, these should be the raw counts, while for other methods such as
ttest, these can be normalized values such as FPKM values. In
addition, sleuth does not use an expression table but rather the
directory of expression estimates from e.g kallisto.
See option --sleuth-counts-dir

The script further requires a design table describing the tests to
be performed. The design table has four columns::

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
     flag to indicate whether or not to include this data (0, 1)
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.

Note: additional columns included after pair can be used to specify
covariates (e.g replicate number etc)


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

To do:
 -- add some E.infos

Document!!!

'''

import sys
import pandas as pd

import CGAT.Experiment as E
import CGAT.Expression as Expression
import CGAT.IOTools as IOTools
import CGAT.Counts as Counts
import CGAT.R as R


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--tag-tsv-file", dest="input_filename_tags",
                      type="string",
                      help="input file with tag counts [default=%default].")

    parser.add_option("-d", "--design-tsv-file", dest="input_filename_design",
                      type="string",
                      help="input file with experimental design "
                      "[default=%default].")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("ttest", "sleuth", "edger", "deseq2", "mock",
                               "dexseq"),
                      help="differential expression method to apply "
                      "[default=%default].")

    parser.add_option("--deseq2-dispersion-method",
                      dest="deseq2_dispersion_method",
                      type="choice",
                      choices=("pooled", "per-condition", "blind"),
                      help="dispersion method for deseq2 [default=%default].")

    parser.add_option("--deseq2-fit-type", dest="deseq2_fit_type", type="choice",
                      choices=("parametric", "local"),
                      help="fit type for deseq2 [default=%default].")

    parser.add_option("--edger-dispersion",
                      dest="edger_dispersion", type="float",
                      help="dispersion value for edgeR if there are no "
                      "replicates [default=%default].")

    parser.add_option("-f", "--fdr", dest="fdr", type="float",
                      help="fdr to apply [default=%default].")

    # currently not implemented
    #parser.add_option("-R", "--output-R-code", dest="save_r_environment",
    #                  type="string",
    #                  help="save R environment to location [default=%default].")

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

    parser.add_option("--model",
                      dest="model",
                      type="string",
                      help=("model for GLM"))

    parser.add_option("--contrast",
                      dest="contrast",
                      type="string",
                      help=("contrast for differential expression testing"))

    parser.add_option("--sleuth-counts-dir",
                      dest="sleuth_counts_dir",
                      type="string",
                      help=("directory containing expression estimates"
                            "from sleuth. Sleuth expects counts"
                            "files to be called abundance.h5"))

    parser.add_option("--dexseq-counts-dir",
                      dest="dexseq_counts_dir",
                      type="string",
                      help=("directory containing counts for dexseq. DEXSeq "
                            "expects counts files to be called .txt and"
                            "to be generated by the DEXSeq_counts.py script"))

    parser.add_option("--dexseq-flattened-file",
                      dest="dexseq_flattened_file",
                      type="string",
                      help=("directory containing flat gtf for dexseq. DEXSeq "
                            "expects this to be generated by the"
                            "DEXSeq_prepare_annotations.py script"))

    parser.add_option("--outfile-sleuth-count",
                      dest="outfile_sleuth_count",
                      type="string",
                      help=("outfile for full count table generated by sleuth"))

    parser.add_option("--outfile-sleuth-tpm",
                      dest="outfile_sleuth_tpm",
                      type="string",
                      help=("outfile for full tpm table generated by sleuth"))

    parser.add_option("--use-ihw",
                      dest="use_ihw",
                      action="store_true",
                      help=("use the independent hypothesis weighting method "
                            "to obtain weighted FDR"))

    parser.add_option("--sleuth-genewise",
                      dest="sleuth_genewise",
                      action="store_true",
                      help=("run genewise, rather than transcript level testing"))

    parser.add_option("--gene-biomart",
                      dest="gene_biomart",
                      type="string",
                      help=("name of ensemble gene biomart"))

    parser.add_option("--de-test",
                      dest="DEtest",
                      type="choice",
                      choices=("wald", "lrt"),
                      help=("Differential expression test"))

    parser.add_option("--Rhistory",
                      dest="Rhistory",
                      type="string",
                      help=("Outfile for R history"))

    parser.add_option("--Rimage",
                      dest="Rimage",
                      type="string",
                      help=("Outfile for R image"))

    parser.set_defaults(
        input_filename_tags="-",
        input_filename_design=None,
        output_filename=sys.stdout,
        method="deseq2",
        fdr=0.1,
        deseq2_dispersion_method="pooled",
        deseq2_fit_type="parametric",
        edger_dispersion=0.4,
        ref_group=False,
        filter_min_counts_per_row=None,
        filter_min_counts_per_sample=None,
        filter_percentile_rowsums=None,
        spike_foldchange_max=4.0,
        spike_expression_max=5.0,
        spike_expression_bin_width=0.5,
        spike_foldchange_bin_width=0.5,
        spike_max_counts_per_bin=50,
        model=None,
        contrast=None,
        output_filename_pattern=None,
        sleuth_counts_dir=None,
        dexseq_counts_dir=None,
        dexseq_flattened_file=None,
        outfile_sleuth_count=None,
        outfile_sleuth_tpm=None,
        use_ihw=False,
        sleuth_genewise=False,
        gene_biomart=None,
        DEtest="wald",
        Rhistory=None,
        Rimage=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    RH = None
    if options.Rhistory or options.Rimage:
        RH = R.R_with_History()

    outfile_prefix = options.output_filename_pattern

    # Expression.py currently expects a refernce group for edgeR and
    # sleuth, regardless of which test is used
    if not options.ref_group and (
            options.method is "edger" or options.method is "sleuth"):
        raise ValueError("Must provide a reference group ('--reference-group')")

    # create Design object
    design = Expression.ExperimentalDesign(
        pd.read_csv(IOTools.openFile(options.input_filename_design, "r"),
                    sep="\t", index_col=0, comment="#"))

    if len(set(design.table[options.contrast])) > 2:

        if options.method == "deseq2" or options.method == "sleuth":
            if options.DEtest == "wald":
                raise ValueError(
                    "Factor must have exactly two levels for Wald Test. "
                    "If you have more than two levels in your factor, "
                    "consider LRT")
        else:
            E.info('''There are more than 2 levels for the contrast
            specified" "(%s:%s). The log2fold changes in the results table
            and MA plots will be for the first two levels in the
            contrast. The p-value will be the p-value for the overall
            significance of the contrast. Hence, some genes will have a
            signficant p-value but 0-fold change between the first two
            levels''' % (options.contrast, set(design[options.contrast])))

    # Sleuth reads in data itself so we don't need to create a counts object
    if options.method == "sleuth":
        assert options.sleuth_counts_dir, (
            "need to specify the location of the abundance.h5 counts files "
            " (--sleuth-counts-dir)")

        # validate design against counts and model
        design.validate(model=options.model)

        experiment = Expression.DEExperiment_Sleuth()
        results = experiment.run(design,
                                 base_dir=options.sleuth_counts_dir,
                                 model=options.model,
                                 contrast=options.contrast,
                                 outfile_prefix=outfile_prefix,
                                 counts=options.outfile_sleuth_count,
                                 tpm=options.outfile_sleuth_tpm,
                                 fdr=options.fdr,
                                 genewise=options.sleuth_genewise,
                                 gene_biomart=options.gene_biomart,
                                 DE_test=options.DEtest,
                                 ref_group=options.ref_group)

    # DEXSeq reads in data itself
    if options.method == "dexseq":
        assert options.dexseq_counts_dir, (
            "need to specify the location of the .txt counts files")

        # create Design object
        design = Expression.ExperimentalDesign(
            pd.read_csv(IOTools.openFile(options.input_filename_design, "r"),
                        sep="\t", index_col=0, comment="#"))

        # validate design against counts and model
        # design.validate(model=options.model)

        experiment = Expression.DEExperiment_DEXSeq()
        results = experiment.run(design,
                                 base_dir=options.dexseq_counts_dir,
                                 model=options.model,
                                 outfile_prefix=outfile_prefix,
                                 flattenedfile=options.dexseq_flattened_file,
                                 fdr=options.fdr)

    else:
        # create Counts object
        if options.input_filename_tags == "-":
            counts = Counts.Counts(pd.io.parsers.read_csv(
                sys.stdin, sep="\t", index_col=0, comment="#"))
        else:
            counts = Counts.Counts(pd.io.parsers.read_csv(
                IOTools.openFile(options.input_filename_tags, "r"),
                sep="\t", index_col=0, comment="#"))

        # validate design against counts and model
        design.validate(counts, options.model)

        # restrict counts to samples in design table
        counts.restrict(design)

        # remove sample with low counts
        if options.filter_min_counts_per_sample:
            counts.removeSamples(
                min_counts_per_sample=options.filter_min_counts_per_sample)

        # remove observations with low counts
        if options.filter_min_counts_per_row:
            counts.removeObservationsFreq(
                min_counts_per_row=options.filter_min_counts_per_row)

        # remove bottom percentile of observations
        if options.filter_percentile_rowsums:
            counts.removeObservationsPerc(
                percentile_rowsums=options.filter_percentile_rowsums)

        # check samples are the same in counts and design following counts
        # filtering and, if not, restrict design table and re-validate
        design.revalidate(counts, options.model)

        # set up experiment and run tests
        if options.method == "ttest":
            experiment = Expression.DEExperiment_TTest()
            results = experiment.run(counts, design)

        elif options.method == "edger":
            experiment = Expression.DEExperiment_edgeR()
            results = experiment.run(counts,
                                     design,
                                     model=options.model,
                                     contrast=options.contrast,
                                     outfile_prefix=outfile_prefix,
                                     ref_group=options.ref_group,
                                     fdr=options.fdr,
                                     dispersion=options.edger_dispersion)

        elif options.method == "deseq2":

            experiment = Expression.DEExperiment_DESeq2()
            results = experiment.run(counts,
                                     design,
                                     model=options.model,
                                     contrast=options.contrast,
                                     outfile_prefix=outfile_prefix,
                                     fdr=options.fdr,
                                     fit_type=options.deseq2_fit_type,
                                     ref_group=options.ref_group,
                                     DEtest=options.DEtest,
                                     R=RH)

    results.getResults(fdr=options.fdr)

    if options.use_ihw:
        results.calculateIHW(alpha=options.fdr)

    for contrast in set(results.table['contrast']):
        results.plotVolcano(contrast, outfile_prefix=outfile_prefix, R=RH)
        results.plotMA(contrast, outfile_prefix=outfile_prefix, R=RH)
        results.plotPvalueHist(contrast, outfile_prefix=outfile_prefix, R=RH)
        results.plotPvalueQQ(contrast, outfile_prefix=outfile_prefix, R=RH)

    results.table.to_csv(sys.stdout, sep="\t", na_rep="NA", index=False)

    results.summariseDEResults()

    # write out summary tables for each comparison/contrast
    for test_group in list(results.Summary.keys()):
        outf = IOTools.openFile("_".join(
            [outfile_prefix, test_group, "summary.tsv"]), "w")
        outf.write("category\tcounts\n%s\n"
                   % results.Summary[test_group].asTable())
        outf.close()

    if options.Rhistory:
        RH.saveHistory(options.Rhistory)
    if options.Rimage:
        RH.saveImage(options.Rimage)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
