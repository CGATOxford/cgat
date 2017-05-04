'''
assoc2coloc.py - Colocalization testing between association results for 2 traits
================================================================================

:Tags: Python

Purpose
-------

.. Genetic association results for different traits are often found to overlap in particular
regions.  This is particularly pertinent for overlapping trait-associations and eQTLs.  This
existence of an association and an eQTL at a locus does not however imply they are causally
linked by the same genetic variants (pleiotropy).  Testing for causal links can be performed
using only the summary statistics from association/eQTL analyses.

This script provides a Python wrapper and generalized input for the R `coloc` package
doi:10.1002/gepi.21765 - Wallace, arXiv (2013)
Statistical testing of shared genetic control for potentially related traits

Usage
-----

Input requires column headers:
SNP - SNP IDs, must be of the same type across traits, i.e. rsIDs
NMISS - number of non-missing individuals in analysis
P - p-value of association for the trait

e.g.

  SNP    NMISS   P
  snp1   1489    1E-10
  snp2   1478    0.579
  snp3   1490    3.1E-4
  .      .       .
  .      .       .

The MAF table must contain the headers:
SNP - this must be of the same ID type as the trait SNPs
MAF - minor allele frequency of the SNP in the population studied

e.g.

  SNP    MAF
  snp1   0.471
  snp2   0.044
  snp3   0.171
  .      .
  .      .

Options
-------
`--trait1/2-results` - summary statistics for trait1/2

`--maf-table` - Table of SNP IDs and MAFs

`--trait1/2-type` - trait type for trait1/2, either `cc` for case-control,
                  or `quant` for quantitative traits.

`--trait1/2-size` -  sample size for trait1/2.  Used to set NMISS column if it is
                  not present.

`--gene-list` - list of genes to restrict analysis to.  Only relevant when
                at least one of the traits is quantitative and relates to
                a QTL analysis, e.g. eQTL, pQTL.

`--trait1/2-snplist` - provide a file with a single SNP ID on each row, use
                     this to restrict the results from trait1/2 to those SNPs.
                     This is especially useful if there are multiple independent
                     association signals at the same locus.

`--chromosome` - chromosome to restrict analysis to.  This is particularly
                 useful when ised in conjunction with `--gene-list`.

`--restrict-from/to` - Restrict the colocalization testing to a specific
                     genomic region. Must also provide `--chromosome`.

`--trait1/2-p-column` - provide an alternative P-value column header

`--trait1/2-prevalence` - The population prevalence of the trait if binary

`--R-script` - Localtion of R script for running colocalisation analysis


Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri as py2ri
from rpy2.robjects import numpy2ri as np2ri
import rpy2.robjects as ro
import pandas as pd
import numpy as np


def testColoc(trait1, trait2, trait1_type, trait2_type,
              scriptsdir, maf_table, gene_list=None,
              trait1_prev=None, trait2_prev=None,
              chromosome=None, start=None, end=None):
    '''
    Perform colocalization testing between two traits.

    Arguments
    ------
    trait1: pandas.core.dataframe
      A data frame containing the summary statistics for
      trait 1

    trait2: pandas.core.dataframe
      A data frame containing the summary statistics for
      trait 2

    trait1_type: string
      Either `cc` or `quant`, denoting the type of trait 1

    trait2_type: string
      Either `cc` or `quant`, denoting the type of trait 2

    maf_table: pandas.core.dataframe
      Data frame containing SNP IDs and MAF

    gene_list: list
      A list of genes to restirct analysis to.  Either trait 1
      or trait 2 must be a quantitative trait

    snp_list: list
      A set of SNPs to restrict the primary trait results to.
      Usually this will be the SNPs that represent a single
      signal of association, particularly for regions where
      multiple independent signals exist.

    trait1_prev: float
      Prevalence of trait1 if binary

    trait2_prev: float
      Prevalence of trait2 if binary

    chromosome: int
      Chromosome to restrict the colocalisation analysis to

    start: int
      start co-ordinate to restrict analysis to.  Must also
    provide `chromosome`. 1-based index, closed [start, end]

    end: int
      end co-ordinate to restrict analysis to.  Must also
      provide `chromosome` and `start`. 1-based index, closed
      [start, end]

    Returns
    -------
    coloc_results: pandas.core.dataframe
      A data frame containing each region (e.g. genes) and
      the posterior probability in favour of each hypothesis:
      H0 - no association with trait1 or trait2, and no colocalisation
      H1 - association with trait 1, but no colocalisation
      H2 - association with trait2, but no colocalisation
      H3 - association with trait1 and 2, but no colocalisation
      H4 - association with trait1 and 2, and colocalised
    '''

    # push all elements into the R environment
    R('''sink(file="sink.text")''')
    R('''suppressPackageStartupMessages(library(coloc))''')
    R('''source("%(scriptsdir)s/coloQtl.R")''' % locals())

    E.info("Pushing results tables into R environment")
    py2ri.activate()
    np2ri.activate()
    r_trait1 = py2ri.py2ri(trait1)
    R.assign("r.trait1", r_trait1)

    r_trait2 = py2ri.py2ri(trait2)
    R.assign("r.trait2", r_trait2)

    r_maf = py2ri.py2ri(maf_table)
    R.assign("r.mafs", r_maf)

    if trait1_prev:
        R.assign("trait1.prev", trait1_prev)
    else:
        R('''trait1.prev <- NULL''')

    if trait2_prev:
        R.assign("trait2.prev", trait2_prev)
    else:
        R('''trait2.prev <- NULL''')

    E.info("Checking for gene list")
    if gene_list:
        E.info("Gene list contains {} genes".format(len(set(gene_list))))
        r_genes = ro.StrVector([rx for rx in set(gene_list)])
        R.assign("gene.list", r_genes)

        E.info("Iterating over gene list")
        R('''res.df <- geneListSnpColocQtl(gene_list=gene.list,'''
          '''results_table=r.trait1, MAF_table=r.mafs, '''
          '''eqtl_table=r.trait2, trait_type="%(trait1_type)s", '''
          '''prev=trait1.prev)''' % locals())
        R('''genes <- rownames(res.df)''')
        try:
            genes = [gx for gx in R["genes"]]
        except TypeError:
            E.warn("There are insufficient SNPs shared with the input "
                   "eQTL results to determine colocalisation")
            R('''genes <- dim(res.df)[1]''')
            genes = R["genes"]

    else:
        R('''res.df <- TwoTraitSnpColocQtl(trait1_table=r.trait1,'''
          '''trait2_table=r.trait2, MAF_table=r.mafs, '''
          '''trait1_type="%(trait1_type)s", trait2_type="%(trait2_type)s",'''
          '''prev1=trait1.prev, prev2=trait2.prev)''')

        R('''genes <- dim(res.df)[1]''')
        genes = R["genes"]
    try:
        try:
            coloc_results = pd.DataFrame(py2ri.ri2py(R["res.df"]))
        except NotImplementedError:
            coloc_results = pd.DataFrame(np.asarray(R["res.df"]))
        coloc_results.index = genes
    except TypeError:
        # if there are no results only a single vector of 0's
        # is returned.  Mush this into a dataframe to output
        E.warn("No colocalisation performed, outputting empty table")
        colo_vec = [fx for fx in R["res.df"]]
        coloc_results = pd.DataFrame(columns=["" for fx in colo_vec],
                                     index=[0, ])
        coloc_results.iloc[0, :] = colo_vec

    coloc_results.columns = ["nSNPs", "H0.PP", "H1.PP", "H2.PP", "H3.PP", "H4.PP"]

    R('''sink(file=NULL)''')

    return coloc_results


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--R-script", dest="scripts_r", type="string",
                      help="PATH to location of R scripts and functions")

    parser.add_option("--trait1-results", dest="trait1_res", type="string",
                      help="summary statistics for trait 1")

    parser.add_option("--trait2-results", dest="trait2_res", type="string",
                      help="summary statistics for trait 2")

    parser.add_option("--maf-table", dest="maf_table", type="string",
                      help="Table containing allele frequency info for "
                      "all SNPs")

    parser.add_option("--maf-snp-column", dest="maf_snpcol", type="string",
                      help="column header containing SNP IDs")

    parser.add_option("--trait1-snplist", dest="trait1_snplist", type="string",
                      help="restrict the analysis to this set of SNPs "
                      "for trait1")

    parser.add_option("--trait2-snplist", dest="trait2_snplist", type="string",
                      help="restrict the analysis to this set of SNPs "
                      "for trait2")

    parser.add_option("--gene-list", dest="gene_list", type="string",
                      help="list of genes to test eQTL-trait overlap with. "
                      "Either trait1 or trait2 must contain a GENE column.")

    parser.add_option("--trait1-type", dest="trait1_type", type="choice",
                      choices=["quant", "cc"], help="Trait 1 type, either "
                      "quantitative (quant) or binary (cc)")

    parser.add_option("--trait2-type", dest="trait2_type", type="choice",
                      choices=["quant", "cc"], help="Trait 2 type, either "
                      "quantitative (quant) or binary (cc)")

    parser.add_option("--trait1-size", dest="trait1_size", type="int",
                      help="sample size for trait1 analysis, only use this "
                      "if the NMISS column is missing")

    parser.add_option("--trait2-size", dest="trait2_size", type="int",
                      help="sample size for trait2 analysis, only use this "
                      "if the NMISS column is missing")

    parser.add_option("--trait1-p-column", dest="trait1_pcol", type="string",
                      help="Column header for P-value column in trait 1 "
                      "results file, if not `P`")

    parser.add_option("--trait2-p-column", dest="trait2_pcol", type="string",
                      help="Column header for P-value column in trait 2 "
                      "results file, if not `P`")

    parser.add_option("--trait1-prevalence", dest="trait1_prev", type="float",
                      help="Prevalence of trait 1 in the population.  Only "
                      "relevant for binary traits")

    parser.add_option("--trait2-prevalence", dest="trait2_prev", type="float",
                      help="Prevalence of trait 2 in the population. Only "
                      "relevant for binary traits")

    parser.add_option("--chromosome", dest="chrome", type="string",
                      help="Restrict analysis to this chromosome.")

    parser.add_option("--restrict-from", dest="restrict_from", type="int",
                      help="start co-ordinate to restrict analysis.  Must "
                      "provide `--chromosome` when restricting region")

    parser.add_option("--restrict-to", dest="restrict_to", type="int",
                      help="end co-ordinate to restrict analysis.  Must "
                      "provide `--chromosome` when restricting region")

    parser.set_defaults(chrome=None,
                        restrict_from=None,
                        restrict_to=None,
                        trait1_prev=None,
                        trait2_prev=None,
                        trait1_pcol="P",
                        trait2_pcol="P",
                        maf_snpcol="SNP")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # check all files contain necessary fields
    if options.trait1_res.endswith(".gz"):
        trait1_comp = "gzip"
    else:
        trait1_comp = None

    t1_nsize = False
    t2_nsize = False

    E.info("Parsing trait 1 file: {}".format(options.trait1_res))
    try:
        trait1_peek = pd.read_table(options.trait1_res, nrows=5,
                                    sep="\s*", header=0,
                                    index_col=None,
                                    compression=trait1_comp,
                                    engine='python')
        try:
            len_cols = len(set(trait1_peek.columns).intersection(["SNP",
                                                                  "NMISS",
                                                                  "{}".format(options.trait1_pcol)]))
            assert len_cols == 3
            trait1_sep = "\s*"
        except AssertionError:
            if options.trait1_size:
                t1_nsize = True
                trait1_sep = "\s*"
                E.warn("NMISS column is not present, "
                       "using input sample size n={}".format(options.trait1_size))
            else:
                raise IOError("Trait-1 input file does not contain "
                              "SNP, NMISS or P columns")
    except StopIteration:
        trait1_peek = pd.read_table(options.trait1_res, nrows=5,
                                    sep="\t", header=0,
                                    compression=trait1_comp,
                                    index_col=None)
        try:
            len_cols = len(set(trait1_peek.columns).intersection(["SNP",
                                                                  "NMISS",
                                                                  "{}".format(options.trait1_pcol)]))
            assert len_cols == 3
            trait1_sep = "\t"
        except AssertionError:
            if options.trait1_size:
                t1_nsize = True
                trait1_sep = "\t"
                E.warn("NMISS column is not present, "
                       "using input sample size n={}".format(options.trait1_size))
            else:
                raise IOError("Trait-1 input file does not contain "
                              "SNP, NMISS or P columns")

    if options.trait2_res.endswith(".gz"):
        trait2_comp = "gzip"
    else:
        trait2_comp = None

    E.info("Parsing trait 2 file: {}".format(options.trait2_res))
    try:
        trait2_peek = pd.read_table(options.trait2_res, nrows=5,
                                    sep="\s*", header=0,
                                    index_col=None,
                                    compression=trait2_comp,
                                    engine='python')
        try:
            len_cols = len(set(trait2_peek.columns).intersection(["SNP",
                                                                  "NMISS",
                                                                  "{}".format(options.trait2_pcol)]))
            assert len_cols == 3
            trait2_sep = "\s*"
        except AssertionError:
            if options.trait2_size:
                t2_nsize = True
                trait2_sep = "\s*"
                E.warn("NMISS column is not present, "
                       "using input sample size n={}".format(options.trait2_size))
            else:
                raise IOError("Trait-2 input file does not contain "
                              "SNP, NMISS or P columns")
    except StopIteration:
        trait2_peek = pd.read_table(options.trait2_res, nrows=5,
                                    sep="\t", header=0,
                                    compression=trait2_comp,
                                    index_col=None)
        try:
            len_cols = len(set(trait2_peek.columns).intersection(["SNP"
                                                                  "NMISS",
                                                                  "{}".format(options.trait2_pcol)]))
            assert len_cols == 3
            trait2_sep = "\t"
        except AssertionError:
            if options.trait2_size:
                t2_nsize = True
                trait2_sep = "\t"
                E.warn("NMISS column is not present, "
                       "using input sample size n={}".format(options.trait2_size))
            else:
                raise IOError("Trait-2 input file does not contain "
                              "SNP, NMISS or P columns")

    E.info("Parsing MAF table file: {}".format(options.maf_table))
    if options.maf_table.endswith(".gz"):
        maf_comp = "gzip"
    else:
        maf_comp = None
    try:
        maf_peek = pd.read_table(options.maf_table, nrows=5,
                                 sep="\s*", header=0,
                                 index_col=None,
                                 compression=maf_comp,
                                 engine='python')
    except StopIteration:
        maf_peek = pd.read_table(options.maf_table, nrows=5,
                                 sep="\t", header=0,
                                 compression=maf_comp,
                                 index_col=None)
    try:
        len_cols = len(set(maf_peek.columns).intersection(["{}".format(options.maf_snpcol),
                                                           "MAF"]))
        assert len_cols == 2
        maf_sep = "\s*"
    except AssertionError:
        raise IOError("Frequency table does not contain "
                      "SNP or MAF columns")

    trait1_results = pd.read_table(options.trait1_res,
                                   sep=trait1_sep,
                                   header=0,
                                   compression=trait1_comp,
                                   index_col=None)

    trait2_results = pd.read_table(options.trait2_res,
                                   sep=trait2_sep,
                                   header=0,
                                   compression=trait2_comp,
                                   index_col=None)
    if options.trait1_pcol != "P":
        trait1_results.loc[:, "P"] = trait1_results[:, options.trait1_pcol]
    else:
        pass

    if options.trait2_pcol != "P":
        trait2_results.loc[:, "P"] = trait2_results.loc[:, options.trait2_pcol]
    else:
        pass

    if t1_nsize:
        trait1_results.loc[:, "NMISS"] = options.trait1_size
    else:
        pass

    if t2_nsize:
        trait2_results.loc[:, "NMISS"] = options.trait2_size
    else:
        pass

    maf_table = pd.read_table(options.maf_table,
                              sep=maf_sep,
                              header=0,
                              compression=maf_comp,
                              index_col=None)

    if options.maf_snpcol != "SNP":
        maf_table.loc[:, "SNP"] = maf_table.loc[:, options.maf_snpcol]
    else:
        pass

    if options.gene_list:
        gene_list = set()
        with open(options.gene_list, "r") as gfile:
            for gene in gfile.readlines():
                gene_list.add(gene.rstrip("\n"))
    else:
        gene_list = None

    # restrict analysis to a specific set of SNP
    # good for picking just SNPs part of independent association
    # signals
    if options.trait1_snplist:
        t1_snplist = set()
        with open(options.trait1_snplist, "r") as t1_sfile:
            for t1snp in t1_sfile.readlines():
                t1_snplist.add(t1snp.rstrip("\n"))

        trait1_results = trait1_results.loc[trait1_results["SNP"].isin(t1_snplist)]
    else:
        pass

    if options.trait2_snplist:
        t2_snplist = set()
        with open(options.trait2_snplist, "r") as t2_sfile:
            for t2snp in t2_sfile.readlines():
                t2_snplist.add(t2snp.rstrip("\n"))

        trait2_results = trait2_results.loc[trait2_results["SNP"].isin(t2_snplist)]
    else:
        pass

    out_df = testColoc(trait1=trait1_results,
                       trait2=trait2_results,
                       trait1_type=options.trait1_type,
                       trait2_type=options.trait2_type,
                       scriptsdir=options.scripts_r,
                       gene_list=gene_list,
                       maf_table=maf_table,
                       trait1_prev=options.trait1_prev,
                       trait2_prev=options.trait2_prev,
                       chromosome=options.chrome,
                       start=options.restrict_from,
                       end=options.restrict_to)

    out_df.to_csv(options.stdout, index_label="Trait",
                  sep="\t")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
