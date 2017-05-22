'''
qc2qc.py - format and manipulate output files from GWA QC steps
===============================================================

:Author:
:Tags: Python

Purpose
-------

.. Filtering, merging and transforming QC output files from GWA

Usage
-----

.. Example use case

Example::

   python qcs2qc.py

Type::

   python qc2qcs.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
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
                      choices=["merge_exclusions", "flag_hets",
                               "find_inbreds", "flag_relations",
                               "discordant_gender"],
                      help="task to execute on phenotype file(s)")

    parser.add_option("--gender-check-file", dest="gender_check",
                      type="string", help="output from gender checking "
                      "by Plink, suffix should be .sexcheck")

    parser.add_option("--relationship-file", dest="relations",
                      type="string", help="output file from IBS "
                      "calculation.  Should contain all pairwise "
                      "relationships.")

    parser.add_option("--inbreeding-coef-file", dest="inbreed_file",
                      type="string", help="file containing either Plink "
                      "or GCTA estimates of F, inbreeding coefficient")

    parser.add_option("--inbreeding-coefficient", dest="inbred_coeff", type="choice",
                      choices=["Fhat1", "Fhat2", "Fhat3", "F",
                               "ibc"], help="inbreeding coefficient "
                      "to use to identify highly inbred individuals")

    parser.add_option("--inbred-cutoff", dest="inbred_cutoff", type="float",
                      help="threshold above which individuals are classed "
                      "as inbred.")

    parser.add_option("--ibs-cutoff", dest="ibs_cutoff", type="float",
                      help="IBS threshold to flag individuals as being "
                      "closely related")

    parser.add_option("--trimmed-relationships", dest="rel_cutoff",
                      type="string", help="output file from Plink "
                      "--rel-cutoff with trimmed data set of unrelated "
                      "individuals.")

    parser.add_option("--heterozygotes-file", dest="hets_file", type="string",
                      help="file from heterozygote analysis containing observed "
                      "homozygosity and F coefficients")

    parser.add_option("--auxillary-file", dest="aux_file", type="string",
                      help="a file of IIDs and FIDs for individuals that are "
                      "to be removed from analysis, unrelated to QC")

    parser.add_option("--plotting-path", dest="plot_path", type="string",
                      help="PATH to save any plots to")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.task == "flag_hets":
        # calculate heterozygosity rates, find and flag
        # individuals > 3 s.d. away from mean value
        # rate = (nonissing - homs) / nonmissing
        # i.e. non-homozygote rate
        flags = gwas.flagExcessHets(options.hets_file,
                                    plot=True,
                                    plot_path=options.plot_path)
        flags.to_csv(options.stdout, index=None, sep="\t")

    elif options.task == "merge_exclusions":
        exclusions = gwas.mergeQcExclusions(hets_file=options.hets_file,
                                            inbred_file=options.inbreed_file,
                                            related_file=options.relations,
                                            gender_file=options.gender_check,
                                            mask_file=options.aux_file)
        exclusions.to_csv(options.stdout, index=None, sep="\t")
    elif options.task == "find_inbreds":
        inbreds = gwas.flagInbred(inbred_file=options.inbreed_file,
                                  inbreeding_coefficient=options.inbred_coeff,
                                  ibc_threshold=options.inbred_cutoff,
                                  plot=True,
                                  plot_path=options.plot_path)
        inbreds.to_csv(options.stdout, sep="\t", index=None)
    elif options.task == "flag_relations":
        # the input file is likely to be huge! Ergo, read the file in chunks
        # calculate any related individuals and store them, store
        # an array of IBD values for plotting, drop the rest
        relate = gwas.flagRelated(ibd_file=options.relations,
                                  chunk_size=500000,
                                  threshold=options.ibs_cutoff,
                                  plot=True,
                                  plotting_path=options.plot_path)
    elif options.task == "discordant_gender":
        sex_discord = gwas.flagGender(gender_file=options.gender_check,
                                      plot=True,
                                      plot_path=options.plot_path)
        sex_discord.to_csv(options.stdout, index=None, sep="\t")
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
