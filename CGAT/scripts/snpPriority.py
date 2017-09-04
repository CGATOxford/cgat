'''
snpPriority.py - score SNPs based on their LD score and SE weighted effect sizes
===============================================================================

:Tags: Python

Purpose
-------

.. Score SNPs based on their LD score and SE weighted effect sizes from
association analysis.

Usage
-----

.. Example use case

Example::

   python snpPriority.py

Type::

   python snpPriority.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.PipelineGWAS as gwas
import re
import pandas as pd
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--score-method", dest="method", type="choice",
                      choices=["PICS", "LDscore", "ABF", "R2_rank",
                               "get_eigen", "calc_prior", "credible_set",
                               "summarise"],
                      help="SNP scoring/prioritisation method to apply.")

    parser.add_option("--database", dest="database", type="string",
                      help="SQL database containing LD information "
                      "in table format. Expects columns SNP_A, "
                      "SNP_B, R2, BP_A and BP_B (Plink --r2 output)")

    parser.add_option("--ld-directory", dest="ld_dir", type="string",
                      help="directory containing tabix-index BGZIP "
                      "LD files.  Assumes Plink used to calculate LD")

    parser.add_option("--table-name", dest="table", type="string",
                      help="name of the SQL table containing the LD"
                      "values")

    parser.add_option("--chromosome", dest="chromosome", type="string",
                      help="chromosome to subset the association results "
                      "file on")

    parser.add_option("--ld-threshold", dest="ld_threshold", type="float",
                      help="the threshold of LD above which variants will "
                      "be taken forward.")

    parser.add_option("--rank-threshold", dest="rank_threshold", type="float",
                      help="the threshold in terms of the top n% SNPs to "
                      "output based on the ranking metric. e.g. "
                      "--rank-threshold=0.01 is the top 1% SNPs")

    parser.add_option("--credible-interval", dest="interval", type="float",
                      help="The credible set interval size to generate the "
                      "credible set of SNPs")

    parser.add_option("--prior-variance", dest="prior_var", type="float",
                      help="the prior variance used to weight the SNP "
                      "variance")

    parser.add_option("--fine-map-window", dest="map_window", type="int",
                      help="the region size to included around the index "
                      "SNP as the fine-mapping region.")

    parser.add_option("--eigen-score-directory", dest="eigen_dir", type="string",
                      help="PATH to directory containing tabix indexed "
                      "eigen score files")

    parser.add_option("--flat-prior", dest="flat_prior", action="store_true",
                      help="Ignore functional annotation information and "
                      "use an uninformative prior on each SNP")

    parser.add_option("--snp-set", dest="snp_set", type="string",
                      help="Pre-defined SNP set as a list of SNP IDs."
                      "If used to calculate priors contains column of scores.")

    parser.add_option("--distribution", dest="dist", type="choice",
                      choices=["normal", "t", "gamma", "lognormal",
                               "exponential"],
                      help="distribution from which to draw prior "
                      "probabilities")

    parser.add_option("--distribution-parameters", dest="dist_params", type="string",
                      help="distribution parameters as a comma-separated list")

    parser.add_option("--lead-snp-id", dest="lead_snp", type="int",
                      help="0-based item number in filename")

    parser.add_option("--filename-separator", dest="separator", type="string",
                      help="filename separator to extract information")

    parser.add_option("--snp-column", dest="snp_col", type="int",
                      help="0-based index of SNP ID column number")

    parser.add_option("--probability-column", dest="prob_col", type="int",
                      help="0-based index of posterior probabilities column"
                      " number")

    parser.set_defaults(ld_dir=None,
                        dist="normal",
                        dist_params=None,
                        snp_set=None,
                        prior_var=0.04,
                        interval=0.99,
                        eigen_dir=None,
                        map_window=100000,
                        ld_threshold=0.5,
                        database=None,
                        table=None,
                        flat_prior=False,
                        lead_snp=2,
                        separator="_",
                        snp_col=0,
                        prob_col=1,
                        )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    if len(infile.split(",")) > 1:
        pass
    else:
        peek = pd.read_table(infile, nrows=5, sep="\s*", header=0)
        try:
            if len(peek["TEST"] != "ADD"):
                clean = False
            else:
                clean = True
        except KeyError:
            clean = True

    if options.method == "LDscore":
        snpscores = gwas.snpPriorityScore(gwas_results=infile,
                                          database=options.database,
                                          table_name=options.table,
                                          chromosome=options.chromosome,
                                          ld_dir=options.ld_dir,
                                          clean=clean)
        # take top 1%, all SNPs doesn't achieve anything useful
        ranks = int(len(snpscores.index) * 0.01)
        snpscores = snpscores.iloc[:ranks]

    elif options.method == "PICS":
        snp_list = {}
        if options.snp_set and not options.flat_prior:
            with IOTools.openFile(options.snp_set, "r") as sfile:
                for line in sfile.readlines():
                    snp = line.split("\t")[0]
                    try:
                        score = float(line.split("\t")[-1].rstrip("\n"))
                    except ValueError:
                        score = 0
                    snp_list[snp] = float(score)

            # get the parameter estimates for the distribution
            # if they have not been provided
            if not options.dist_params:
                dist_params = gwas.estimateDistributionParameters(data=snp_list.values(),
                                                                  distribution=options.dist)
            else:
                dist_params = tuple([float(fx) for fx in options.dist_params.split(",")])

            E.info("Calculating priors on SNPs")
            priors = gwas.calcPriorsOnSnps(snp_list=snp_list,
                                           distribution=options.dist,
                                           params=dist_params)

        elif options.snp_set and options.flat_prior:
            with IOTools.openFile(options.snp_set, "r") as sfile:
                for line in sfile.readlines():
                    snp = line.split("\t")[0]
                    snp_list[snp] = 1.0

            priors = snp_list

        else:
            # allow for no priors or scores to be set,
            # use of priors will be ignored,
            # i.e. when prior and likelihood are not from
            # conjugate distributions
            priors = None

        # PICS scores expects the gwas results file to
        # only contain the region of interest, which
        # represents an independent association signal
        # if a SNP has not been genotyped,
        # but it is in strong LD, it will cause problems
        # downstream <- only allow SNPs that
        # are present in the analysis
        snpscores = gwas.PICSscore(gwas_results=infile,
                                   database=options.database,
                                   table_name=options.table,
                                   chromosome=options.chromosome,
                                   priors=priors,
                                   clean=clean,
                                   ld_dir=options.ld_dir,
                                   ld_threshold=options.ld_threshold)

        snpscores.columns = ["SNP", "PICS"]
        posterior_sum = 0
        snpscores.sort_values(ascending=False,
                              inplace=True)
        post_snps = []
        for snp in snpscores.index:
            if posterior_sum < 99.0:
                posterior_sum += snpscores.loc[snp]
                post_snps.append(snp)
            else:
                break

        snpscores = snpscores.loc[post_snps]

        snpscores.drop_duplicates(inplace=True)

    elif options.method == "R2_rank":
        # rank SNPs based on their LD with the lead
        # SNP, take the top n% SNPs
        snpscores = gwas.LdRank(gwas_results=infile,
                                database=options.database,
                                table_name=options.table,
                                ld_dir=options.ld_dir,
                                chromosome=options.chromosome,
                                ld_threshold=options.ld_threshold,
                                top_snps=options.rank_threshold,
                                clean=clean)

    elif options.method == "ABF":
        snpscores = gwas.ABFScore(gwas_results=infile,
                                  region_size=options.map_window,
                                  chromosome=options.chromosome,
                                  prior_variance=options.prior_var,
                                  clean=clean)
    elif options.method == "get_eigen":
        E.info("Fetching Eigen scores")
        snpscores = gwas.getEigenScores(eigen_dir=options.eigen_dir,
                                        bim_file=infile,
                                        snp_file=options.snp_set)
        snpscores = pd.DataFrame(snpscores).T

    elif options.method == "credible_set":
        E.info("Creating credible set")

        snpscores = gwas.makeCredibleSet(probs_file=infile,
                                         credible_set=options.interval,
                                         lead_snp_indx=options.lead_snp,
                                         filename_sep=options.separator,
                                         snp_column=options.snp_col,
                                         probs_column=options.prob_col)

    elif options.method == "summarise":
        E.info("Collating SNP prioritisation resuslts")
        file_list = infile.split(",")
        snpscores = gwas.summariseResults(file_list=file_list)

    snpscores.to_csv(options.stdout, index_label="SNP",
                     sep="\t")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
