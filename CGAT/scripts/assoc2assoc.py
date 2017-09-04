'''
assoc2assoc.py filter, transform and process results from genome-wide analyses
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Format, manipulation and processing of genome-wide association results

Usage
-----

* `extract_results` - extract GWAS results for a specific SNP set
* `merge_freq` - merge results with a .bim file to get allele information
Example::

   python assoc2assoc.py

Type::

   python assoc2assoc.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.PipelineGWAS as gwas
import CGAT.IOTools as IOTools
import re


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--task", dest="task", type="choice",
                      choices=["get_hits", "extract_results",
                               "merge_freq"],
                      help="task to perform")

    parser.add_option("--p-threshold", dest="p_threshold", type="float",
                      help="threshold for association p-value, below "
                      "which results will be output")

    parser.add_option("--output-directory", dest="outdir", type="string",
                      help="output file directory")

    parser.add_option("--snp-set", dest="snpset", type="string",
                      help="file containing list of SNP per row to "
                      "extract from GWAS results")

    parser.add_option("--frequency-directory", dest="freq_dir", type="string",
                      help="Directory containing plink .frq files corresponding"
                      " to all chromosomes")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # if the input is a list of files, split them
    infile = argv[-1]
    infiles = infile.split(",")
    if len(infiles) > 1:
        results = gwas.GWASResults(assoc_file=infiles)
    elif len(infiles) == 1:
        results = gwas.GWASResults(assoc_file=infile)
    else:
        raise IOError("no input files detected, please specifiy association "
                      "results files as the last command line argument")

    if options.task == "get_hits":
        hits = results.getHits(float(options.p_threshold))
        for name, region in hits:
            try:
                try:
                    top_reg = region.sort_values(by="CHISQ",
                                                 ascending=False)
                    top_bp = top_reg.iloc[0]["BP"]
                    top_snp = top_reg.iloc[0]["SNP"]
                except KeyError:
                    top_reg = region
                    top_reg.loc[:, "STAT"] = abs(top_reg["STAT"])
                    top_reg = top_reg.sort_values(by="STAT",
                                                  ascending=False)
                    top_bp = top_reg.iloc[0]["BP"]
                    top_snp = top_reg.iloc[0]["SNP"]
            except KeyError:
                top_reg = region
                top_reg.loc[:, "STAT"] = abs(top_reg["T"])
                top_reg = top_reg.sort_values(by="T",
                                              ascending=False)
                top_bp = top_reg.iloc[0]["BP"]
                top_snp = top_reg.iloc[0]["SNP"]

            outname = "_".join(["chr%s" % str(name),
                                str(top_bp),
                                top_snp,
                                "significant"])

            outfile = outname + ".tsv"
            out_file = "/".join([options.outdir, outfile])
            E.info("output association results from Chr%s to %s" %
                   (str(name), out_file))
            # this keeps outputing the first column as unamed: 0,
            # need to remove this
            try:
                if region.columns[0] != "A1":
                    region.drop([region.columns[0]], inplace=True, axis=1)
            except:
                pass

            region.to_csv(out_file, sep="\t", index=None)

    elif options.task == "extract_results":
        with IOTools.openFile(options.snpset, "r") as sfile:
            snpset = sfile.readlines()
            snpset = [snp.rstrip("\n") for snp in snpset]

        snp_df = results.extractSNPs(snpset)
        snp_df.dropna(axis=0, how='all', inplace=True)
        snp_df.drop_duplicates(subset=["SNP"], inplace=True)
        snp_df.to_csv(options.stdout, sep="\t", index=None)

    elif options.task == "merge_freq":
        # sequentially merge GWAS result with frequency data
        # to make file for GCTA joint analysis
        regex = re.compile("(\S+).frq$")
        cojo_df = results.mergeFrequencyResults(options.freq_dir,
                                                file_regex=regex)
        cojo_df.to_csv(options.stdout, sep="\t", index=None)
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
