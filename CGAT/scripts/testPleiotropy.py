'''
testPleiotropy - implementation of the Pleiotropy Estimation and Testing method
===============================================================================

:Tags: Python

Purpose
-------

Test genetic variants for pleiotropy on 2 traits using the PET-B method
described in Zhang et al Genetic Epidemiology 2014.  This calculates
a pleiotropy correlation coefficient (PCC), and calculates p-values
by bootstrapping.

The method tests a composite null of no pleiotropy versus the
alternative hypothesis of pleiotropy.  The composite null contains
both the complete null (b1 = b2 = 0) for genetic coefficients on
traits 1 and 2 (b1 and b2), respectively, and the incomplete
null (b1 = 0, b2 =/=0 | b1 =/= 0, b2 = 0).  Bootstraping is
used to generate empirical p-values as the distribution of the
PCC under the composite null distribution is unknown.

Usage
-----

.. Example use case

Example::

   python testPleiotropy.py

Type::

   python testPleiotropy.py --help

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
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri as py2ri
from rpy2.robjects.packages import importr


def pythonWrapper4Pet(dataframe, snps, covars,
                      trait1, trait2, model1,
                      scriptsdir,
                      model2, resamples=999):
    '''
    This is just Python wrapper around the R code
    for the PET calculations

    '''
    py2ri.activate()

    E.info("Checking regression models")
    if model1 == "logistic":
        R('''trait1.mod <- binomial''')
        R('''trait1.link <- "logit" ''')
    elif model1 == "linear":
        R('''trait1.mod <- gaussian''')
        R('''trait1.link <- "identity" ''')

    if model2 == "logistic":
        R('''trait2.mod <- binomial''')
        R('''trait2.link <- "logit" ''')
    elif model2 == "linear":
        R('''trait2.mod <- gaussian''')
        R('''trait2.link <- "identity" ''')
    E.info("Running {} regression for trait 1: {}".format(model1,
                                                          trait1))
    E.info("Running {} regression for trait 2: {}".format(model2,
                                                          trait2))

    R('''source("%(scriptsdir)s/PET_functions.R")''' % locals())
    E.info("Pushing data objects into the R environment")
    # push everything into the R environment
    r_df = py2ri.py2ri_pandasdataframe(dataframe)
    R.assign("data.df", r_df)

    r_snps = ro.StrVector([sp for sp in snps])
    R.assign("snp.list", r_snps)

    E.info("Parsing covariates")
    covars = covars.split(",")
    r_covar = ro.StrVector([cv for cv in covars])
    R.assign("covar.list", r_covar)
    E.info("{} covariates found to adjust "
           "in regression  models".format(len(covars)))

    # clean up, replacing "missing values" with NAs for R
    R('''data.df[data.df == -9] <- NA''')
    R('''pet_results <- list()''')

    # loop over all SNP, calculate PCC and p-value
    # this takes a long time <- need to think of speed ups
    # possible Python-pure implementation, i.e. with LIMIX?
    E.info("Iteratively calculating PCC for all SNPs")
    R('''results <- loopPET(data.df=data.df, trait1="%(trait1)s", trait2="%(trait2)s", '''
      '''trait1.link=trait1.link, trait2.link=trait2.link, '''
      '''trait1.mod=trait1.mod, trait2.mod=trait2.mod, covars=covar.list,'''
      '''resamples=%(resamples)i, snp.list=snp.list)''' % locals())

    R('''out.res <- data.frame(do.call(rbind, results))''')
    R('''colnames(out.res) <- c("PCC", "pvalue")''')
    py_out = py2ri.ri2py_dataframe(R["out.res"])

    return py_out


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--R-scripts", dest="scripts_r", type="string",
                      help="PATH to R scripts and functions")

    parser.add_option("--trait1", dest="trait1", type="string",
                      help="name/column header of trait 1 in the "
                      "input data table")

    parser.add_option("--trait2", dest="trait2", type="string",
                      help="name/column header of trait 2 in the "
                      "input data table")

    parser.add_option("--snp-list", dest="snp_list", type="string",
                      help="optional list of snps on which to "
                      "restrict analysis.")

    parser.add_option("--covariates", dest="covars", type="string",
                      help="column headers that refer to covariates "
                      "to adjust primary traits for")

    parser.add_option("--resamples", dest="resample", type="int",
                      help="number of resamples with replacement "
                      "to use for bootstrapping")

    parser.add_option("--trait1-model", dest="trait1_mod", type="choice",
                      choices=["logistic", "linear"],
                      help="model to use to fit covariates and trait")

    parser.add_option("--trait2-model", dest="trait2_mod", type="choice",
                      choices=["logistic", "linear"],
                      help="model to use to fit covariates and trait")

    parser.set_defaults(resample=999,
                        )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    # snp headers are assumed to start with 'rs'
    # read the dataframe with pandas then push
    # it into R

    df = pd.read_table(infile, sep="\t", header=0,
                       index_col=None)

    E.info("Parsing SNP IDs")
    if options.snp_list:
        snp_set = set()
        with open(options.snp_list, "r") as sfile:
            snp_list = [s.rstrip("\n") for s in sfile.readlines()]
            snp_list = set(snp_list)

        for snp in snp_list:
            snp_re = re.compile(snp)
            snp_set.update([sx for sx in df.columns if re.search(snp_re, sx)])
        snps = [st for st in snp_set]
    else:
        snp_re = re.compile("^rs")
        snps = [sx for sx in df.columns if re.search(snp_re, sx)]

    E.info("{} SNPs found in data table".format(len(snps)))

    out_df = pythonWrapper4Pet(dataframe=df,
                               snps=snps,
                               covars=options.covars,
                               scriptsdir=options.scripts_r,
                               trait1=options.trait1,
                               trait2=options.trait2,
                               model1=options.trait1_mod,
                               model2=options.trait2_mod,
                               resamples=options.resample)

    out_df.to_csv(options.stdout, sep="\t", index_label="SNP")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
