'''
snps2architecture.py - Calculate & plot disease frequency by allele
===================================================================

:Author:
:Tags: Python

Purpose
-------
This script requires a Plink .ped file and .map file.  It calculates
two sets of summary statistics:
    1. frequency of cases carrying risk alleles
    2. cumulative frequency of cases carrying risk alleles

It also generates a plot of the desired cases frequency against the
number of risk alleles carried.

This allows one to observe whether there exists a threshold effect
at the genotype level for a trait.  For instance, if there is a
sudden jump in the number of cases explained at a given allele
carriage threshold.  This would provide evidence of a liability threshold
trait.

The steps are:
  * select SNPs that are also in gwas results
  * recode A1 alleles as risk alleles
  * set up case and control frequency arrays
  * count risk alleles for each individual
  * calculate P(Disease) for each risk allele bin
  * plot


Methods
-------
Options for `--method`.

`cases_explained` - the frequency of cases carrying the given
number of risk alleles.

`probability_phenotype` - the frequency of the (binary) phenotype
amongst all individuals, given the risk allele carriage

Reporting & Plotting
--------------------
This is defined by the `--plot-statistic` option.

`cumulative` - plot the cumulative frequency over the number of alleles

`frequency` - plot the observed frequency of cases for each number of
risk alleles

Example::

   python snps2architecture.py

Type::

   python snps2architecture.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.PipelineGWAS as gwas
import pandas as pd
import numpy as np


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

    parser.add_option("--method", dest="method", type="choice",
                      choices=["cases_explained", "probability_phenotype"],
                      help="Which results to report, either the proportion "
                      "of cases explained or the probability of the "
                      "phenotype given the number of alleles carried")

    parser.add_option("--map-file", dest="map_file", type="string",
                      help="plink .map file with SNP positions")

    parser.add_option("--ped-file", dest="ped_file", type="string",
                      help="plink ped file with phenotype and "
                      "genotype data - A2 major allele coded")

    parser.add_option("--gwas-file", dest="gwas", type="string",
                      help="gwas results file, assumes Plink "
                      "output format.  Must contain SNP, BP, "
                      "OR column headers.  Assumes results relate "
                      "to the A1 allele")

    parser.add_option("--flip-alleles", dest="flip", action="store_true",
                      help="force alleles to flip if OR < 1")

    parser.add_option("--plot-statistic", dest="plot_stat", type="choice",
                      choices=["frequency", "cumulative"],
                      help="plot either cases frequency or cumulative "
                      "frequency of cases")

    parser.add_option("--plot-path", dest="plot_path", type="string",
                      help="save path for plot")

    parser.add_option("--flag-explained-recessive", dest="explained",
                      action="store_true",
                      help="flag individuals explained by carriage of "
                      "2 risk alleles - NOT IMPLIMENTED")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # required files are .ped file, .map file and gwas results file
    E.info("reading GWAS results file: %s" % options.gwas)
    snp_df = pd.read_table(options.gwas, sep="\t", header=0,
                           index_col=None)
    snp_list = snp_df["SNP"].values

    # parse ped file
    E.info("Reading ped file: %s" % options.ped_file)
    ped_df = gwas.parsePed(options.ped_file,
                           compound_geno=True)

    # parse map file and get SNP indices that correspond to
    # ped file genotypes
    E.info("Fetching SNPs from map file: %s" % options.map_file)
    snp_index = gwas.getSNPs(options.map_file,
                             snp_list)

    E.info("SNPs found: %i" % len(snp_index))
    # extract SNPs and ORs as key, value pairs
    or_dict = snp_df.loc[:, ["SNP", "OR"]].to_dict(orient='list')
    snp_or = dict(zip(or_dict["SNP"], or_dict["OR"]))

    if options.flip:
        E.info("Flipping major alleles to risk alleles")
        flipped_genos = gwas.flipRiskAlleles(snp_index=snp_index,
                                             snp_results=snp_or,
                                             genos=ped_df["GENOS"].tolist())
        # merge flipped genotypes with pedigree frame to get phenotypes
        geno_df = pd.DataFrame(flipped_genos, index=ped_df["FID"])
    else:
        # split genos into a dataframe
        genos = np.array(ped_df["GENOS"].tolist())
        geno_df = pd.DataFrame(genos, index=ped_df["FID"])

    merged = pd.merge(geno_df, ped_df, left_index=True, right_on="FID")

    # need to discount missing genotypes > 1%

    # frequencies of number of risk alleles by trait frequency
    E.info("count #risk alleles per individual")
    risk_results = gwas.countRiskAlleles(ped_frame=merged,
                                         snp_index=snp_index.values(),
                                         report=options.method,
                                         flag=options.explained)
    risk_freqs = risk_results["freqs"]
    cumulative = risk_results["cumulative"]
    # select results upto and including cumulative freq = 1.0
    max_indx = [fx for fx, fy in enumerate(cumulative) if fy == 1.0][0]
    max_freqs = risk_freqs[:max_indx + 1]
    max_cum = cumulative[:max_indx + 1]
    bins = [ix for ix, iy in enumerate(cumulative)][:max_indx + 1]

    # plot!
    # need to add number of individuals into each bin as point size
    if options.plot_stat == "frequency":
        E.info("Generating plot of #risk alleles vs. P(Phenotype)")
        hist_df = gwas.plotRiskFrequency(bins=bins,
                                         frequencies=max_freqs,
                                         counts=risk_results["cases"][:max_indx + 1],
                                         savepath=options.plot_path,
                                         ytitle="P(Phenotype)")
    elif options.plot_stat == "cumulative":
        E.info("Generating plot of #risk alleles vs. cumulative frequency")
        hist_df = gwas.plotRiskFrequency(bins=bins,
                                         frequencies=max_cum,
                                         counts=risk_results["cases"][:max_indx + 1],
                                         savepath=options.plot_path,
                                         ytitle="Cumulative frequency cases")

    hist_df["freq"] = risk_results["freqs"][:max_indx + 1]
    hist_df["cumulative"] = risk_results["cumulative"][:max_indx + 1]
    hist_df["cases"] = risk_results["cases"][:max_indx + 1]
    hist_df["controls"] = risk_results["controls"][:max_indx + 1]
    hist_df["total"] = hist_df["cases"] + hist_df["controls"]
    hist_df.to_csv(options.stdout, sep="\t", index=None)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
