'''
geno2geno.py - format, manipulate and process genotying data files
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Format, manipulation and processing of genotyping files

Usage
-----

.. Example use case

Example::

   python geno2geno.py

Type::

   python geno2geno.py --help

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
import os


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
                      choices=["mafs", "penetrance",
                               "detect_duplicates", "allele_diff"],
                      help="task to perform")

    parser.add_option("--ped-file", dest="ped_file", type="string",
                      help="plink format .ped file")

    parser.add_option("--map-file", dest="map_file", type="string",
                      help="plink format .map file")

    parser.add_option("--freq-file", dest="mafs", type="string",
                      help="text file containing populations minor "
                      "allele frequencies of variants.  One row per "
                      "variant with ID MAF")

    parser.add_option("--groups-file", dest="group_file", type="string",
                      help="file containing group labels for individuals "
                      "in the provided ped file")

    parser.add_option("--ref-label", dest="ref_label", type="string",
                      help="group label to be used as the reference case")

    parser.add_option("--test-label", dest="test_label", type="string",
                      help="group label to be used as the test case")

    parser.add_option("--subset", dest="subset", type="choice",
                      choices=["cases", "gender"], help="subset the "
                      "data by either case/control or gender")

    parser.add_option("--take-last", dest="take", action="store_true",
                      help="if use duplicates will take the last variant, "
                      "default behaviour is to take the first")

    parser.add_option("--outfile-pattern", dest="out_pattern", type="string",
                      help="outfile pattern to use for finding duplicates "
                      "and triallelic variants")

    parser.add_option("--snp-set", dest="snp_subset", type="string",
                      help="list of SNPs to include")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(mafs=None,
                        subset=None,
                        take_last=False)

    if options.task == "mafs":
        mafs = gwas.countByVariantAllele(options.ped_file,
                                         options.map_file)

        mafs.to_csv(options.stdout,
                    index_col=None,
                    sep="\t")

    elif options.task == "penetrance":
        summary, pens = gwas.calcPenetrance(options.ped_file,
                                            options.map_file,
                                            subset=options.subset,
                                            mafs=options.mafs,
                                            snpset=options.snp_subset)

        pens.to_csv(options.stdout, sep="\t", index_label="SNP")
        summary.to_csv("/".join([os.getcwd(), "penetrance_summary.txt"]),
                       sep="\t", index_label="SNP")

    elif options.task == "allele_diff":
        allele_diffs = gwas.calcMaxAlleleFreqDiff(ped_file=options.ped_file,
                                                  map_file=options.map_file,
                                                  group_file=options.group_file,
                                                  test=options.test_label,
                                                  ref=options.ref_label)

        allele_diffs.to_csv(options.stdout, sep="\t")

    elif options.task == "detect_duplicates":
        # find variants with duplicated position and shared reference
        # allele indicative of triallelic variants - also same ID
        # ouput to a filter list
        infile = argv[-1]
        dups, tris, oves = gwas.findDuplicateVariants(bim_file=infile,
                                                      take_last=options.take)

        if os.path.isabs(options.out_pattern):
            with open(options.out_pattern + ".triallelic",
                      "w") as otfile:
                for tvar in tris:
                    otfile.write("%s\n" % tvar)

            with open(options.out_pattern + ".duplicates",
                      "w") as odfile:
                for dvar in dups:
                    odfile.write("%s\n" % dvar)

            with open(options.out_pattern + ".overlapping",
                      "w") as ovfile:
                for ovar in oves:
                    ovfile.write("%s\n" % ovar)
        else:
            outpattern = os.path.abspath(options.out_pattern)
            with open(outpattern + ".triallelic",
                      "w") as otfile:
                for tvar in tris:
                    otfile.write("%s\n" % tvar)

            with open(outpattern + ".duplicates",
                      "w") as odfile:
                for dvar in dups:
                    odfile.write("%s\n" % dvar)

            with open(outpattern + ".overlapping",
                      "w") as ovfile:
                for ovar in oves:
                    ovfile.write("%s\n" % ovar)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
