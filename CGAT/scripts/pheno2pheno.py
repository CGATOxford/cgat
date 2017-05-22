'''
pheno2pheno.py - format and manipulate phenotype files
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Project specific operations on phenotype files

Usage
-----

.. Example use case

Example::

   python pheno2pheno.py

Type::

   python pheno2pheno.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
import pandas as pd
import re
import numpy as np
import CGAT.PipelineGWAS as gwas


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

    parser.add_option("--task", dest="task", type="choice",
                      choices=["set_factors", "dichotimise_phenotype",
                               "plink_format", "select_ethnicity",
                               "merge_covariates", "subset_phenotypes"],
                      help="task to execute on phenotype file(s)")

    parser.add_option("--R-script", dest="r_script", type="string",
                      help="R script for table reformatting")

    parser.add_option("--adjustment", dest="adjust", type="choice",
                      choices=["snp"],
                      help="adjustements to make pre- or post mergeing")

    parser.add_option("--pheno-id", dest="dichot_var", type="string",
                      help="column header of variable to be dichotimised")

    parser.add_option("--reference-variable", dest="ref_level", type="string",
                      help="level of variable to be dichotimised toi set to 1")

    parser.add_option("--missing-var-label", dest="missing_label", type="string",
                      help="missing/unobserved value labels")

    parser.add_option("--id-variable", dest="id_var", type="string",
                      help="ID variable column header")

    parser.add_option("--ethnicity-id", dest="ethnic_var", type="string",
                      help="column header for variable containing "
                      "ethnicity data")

    parser.add_option("--ethnicity-label", dest="ethnic", type="string",
                      help="ethnicity label to select samples on")

    parser.add_option("--covariate-file", dest="covar_file", type="string",
                      help="a comma-separated list of files to be merged, or "
                      "a single file")

    parser.add_option("--fam-file", dest="fam_file", type="string",
                      help="Plink .fam file that specifies which samples "
                      "to subset from the phenotypes file")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]
    if options.task == "set_factors":
        pandas2ri.activate()
        R('''source("%s")''' % options.r_script)
        R('''format <- format_phenotypes("%s")''' % infile)
        pheno_df = pandas2ri.ri2py_dataframe(R["format"])
        pheno_df["IID"] = pheno_df["f.eid"]
        cols = pheno_df.columns.tolist()
        cols = [xc if not re.search("f.eid", xc) else "FID" for xc in cols]
        # columns need to be FID, IID, ...
        cols.remove("FID")
        cols.remove("IID")
        new_cols = cols
        new_cols.insert(0, "FID")
        new_cols.insert(1, "IID")
        pheno_df.columns = new_cols
        pheno_df.to_csv(options.stdout, sep="\t", index_col=None)

    elif options.task == "dichotimise_phenotype":
        # catch situation where delimiter does is not tab
        try:
            df = pd.read_table(infile, sep="\t", header=0, index_col=None)
            assert len(df.columns) > 1
        except AssertionError:
            df = pd.read_table(infile, sep="\s+", index_col=None)

        var = pd.Series(df[options.dichot_var].copy(), dtype=np.int64)
        ref = np.int64(options.ref_level)

        # there maybe multiple missing/unobserved data categories to deal with
        # output in plink format
        missing = options.missing_label.split(",")
        miss_mask = var.isin([int(x) for x in missing])
        var[miss_mask] = np.nan

        mask = var.isin([ref])
        var[mask] = 2
        nas = np.isnan(var)
        var[~(mask) & ~(nas)] = 1

        # set NA or unobserved to missing, assume missing value is -9 (Plink standard)
        var[nas] = -9
        p_df = df.loc[:, ("FID", "IID", options.dichot_var)]
        p_df[options.dichot_var] = pd.Series(var, dtype=np.int64)
        p_df.index = p_df["FID"]
        p_df.drop(labels="FID", axis=1, inplace=True)
        p_df.to_csv(options.stdout, sep="\t", index_col=None)

    elif options.task == "plink_format":
        # add IID and FID columns based on individual IDs
        pheno_df = pd.read_table(infile, sep="\t", header=0, index_col=None)
        pheno_df["IID"] = pheno_df[options.id_var]
        pheno_df["FID"] = pheno_df[options.id_var]
        cols = pheno_df.columns.tolist()
        cols = [xc for xc in cols if not re.search(options.id_var, xc)]
        # columns need to be FID, IID, ...
        cols.remove("FID")
        cols.remove("IID")
        new_cols = cols
        new_cols.insert(0, "FID")
        new_cols.insert(1, "IID")
        resort_df = pheno_df[new_cols]
        resort_df.index = resort_df["FID"]
        resort_df.drop(labels="FID", axis=1, inplace=True)
        resort_df.to_csv(options.stdout, sep="\t")

    elif options.task == "select_ethnicity":
        # select ethnicity
        pheno_df = pd.read_table(infile, sep="\t", header=0, index_col=None)
        ethnic_var = pheno_df.loc[:, options.ethnic_var].copy()
        ethnic_mask = ethnic_var == int(options.ethnic)
        select_indv = ethnic_var[ethnic_mask].index
        filter_df = pheno_df.loc[select_indv, :]
        filter_df.index = filter_df["FID"]
        filter_df.drop(labels="FID", axis=1, inplace=True)
        filter_df.to_csv(options.stdout, sep="\t", index_col=None)

    elif options.task == "merge_covariates":
        if len(options.covar_file.split(",")) > 1:
            filelist = options.covar_file.split(",")
            df = pd.read_table(filelist.pop(0), sep="\t",
                               index_col=None, header=0)
            if options.adjust == "snp":
                re_snp = re.compile(".raw")
                snp_file = [fil for fil in filelist if re.search(re_snp,
                                                                 fil)][0]
                _df = pd.read_table(snp_file, sep="\t", header=0,
                                    index_col=None)

                cols = _df.columns[6:]
                real_cols = list(_df.columns[:6])
                snp_cols = [sc.split("_")[:-1][0] for sc in cols]

                # list methods work in place, don't assign as a new variable
                real_cols.extend(snp_cols)
                _df.columns = real_cols
                df = pd.merge(left=df, right=_df,
                              on=["FID", "IID"],
                              how='inner')
                try:
                    filelist.remove(snp_file)
                except:
                    pass

            for fle in filelist:
                _df = pd.read_table(fle, sep="\t", header=0,
                                    index_col=None)
                df = pd.merge(left=df, right=_df,
                              on=["FID", "IID"],
                              how='inner')
            # python outputs NA as blank when writing to stdout,
            # plink expects values, use string NAs
            df = df.fillna("NA")
            df.index = df["FID"]
            df.drop(["FID"], inplace=True, axis=1)

            # IID also need to be one of the first 2 columns
            cols = df.columns.tolist()

            # columns need to be FID, IID, ...
            cols.remove("IID")
            new_cols = cols
            new_cols.insert(0, "IID")
            df = df[new_cols]
            df.to_csv(options.stdout, index_col=0,
                      index_label="FID", sep="\t")
        else:
            E.warn("only a single covariates file provided."
                   "No merging possible, exiting")

    elif options.task == "subset_phenotypes":
        fam_df = pd.read_table(options.fam_file, sep=None,
                               index_col=None, header=None)
        fam_df.columns = ["FID", "IID", "PAT", "MAT", "SEX",
                          "PHENO"]

        pheno_df = pd.read_table(infile, sep=None,
                                 index_col=0, header=0)
        fam_ids = fam_df["FID"]
        sub_pheno = pheno_df.loc[fam_ids]

        sub_pheno.to_csv(options.stdout, index_col=0,
                         index_label="FID", sep="\t")
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
