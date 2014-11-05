'''
codeml2tsv.py - analyze results from codeml kaks run
=============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Example::

   python codeml2tsv.py --help

Type::

   python codeml2tsv.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import tempfile
import subprocess
import optparse

from types import *

import CGAT.Genomics as Genomics
import CGAT.Experiment as E
import CGAT.WrapperCodeML as WrapperCodeML
import scipy
import scipy.stats
import CGAT.TreeTools as TreeTools

import CGAT.Stats as Stats


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: codeml2tsv.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--methods", dest="methods", type="string",
                      help="""methods for analysis.
write-ks-tree: write out ks tree(s).
write-ka-tree: write out ka tree(s).
                      """  )

    parser.add_option("--column-prefix", dest="prefix", type="string",
                      help="prefix for rows.")

    parser.add_option("--input-pattern", dest="pattern_input_filenames", type="string",
                      help="input pattern.")

    parser.add_option("--filter-probability", dest="filter_probability", type="float",
                      help="threshold for probability above which to include positive sites.")

    parser.add_option("--filter-omega", dest="filter_omega", type="float",
                      help="threshold for omega above which to include positive sites.")

    parser.add_option("--models", dest="models", type="string",
                      help="restrict output to set of site specific models.")

    parser.add_option("--significance-threshold", dest="significance_threshold", type="float",
                      help="significance threshold for log-likelihood test.")

    parser.add_option("--mode", dest="mode", type="choice",
                      choices=("pairs", "1xn"),
                      help="analysis mode.")

    parser.set_defaults(
        methods="",
        prefix=None,
        filter_probability=0,
        filter_omega=0,
        models="",
        significance_threshold=0.05,
        mode="pairs",
    )

    (options, args) = E.Start(parser)

    options.methods = options.methods.split(",")
    options.models = options.models.split(",")

    codeml = WrapperCodeML.CodeML()

    results = []
    if len(args) == 0:
        # read from stdin, if no arguments are given
        results.append(codeml.parseOutput(sys.stdin.readlines()))
    else:
        # read multiple results
        for f in args:
            try:
                results.append(codeml.parseOutput(open(f, "r").readlines()))
            except WrapperCodeML.ParsingError, msg:
                options.stdlog.write(
                    "# parsing error in file %s: %s.\n" % (f, msg))
                continue

    if options.prefix:
        prefix_tree = ">%s\n" % options.prefix
        prefix_header = "prefix\t"
        prefix_row = "%s\t" % options.prefix
    else:
        prefix_tree = ""
        prefix_header = ""
        prefix_row = ""

    for method in options.methods:

        if method == "write-ks-tree":
            for result in results:
                options.stdout.write(
                    prefix_tree + TreeTools.Tree2Newick(result.mTreeKs) + "\n")

        elif method == "write-ka-tree":
            for result in results:
                options.stdout.write(
                    prefix_tree + TreeTools.Tree2Newick(result.mTreeKa) + "\n")

        elif method == "write-kaks-tree":
            for result in results:
                options.stdout.write(
                    prefix_tree + TreeTools.Tree2Newick(result.mTreeKaks) + "\n")

        elif method == "lrt":
            # perform log-likelihood ratio test between successive models
            # Assumption is that the models are nested with the previous model
            # being the less complex model.
            first_result = results[0]
            last_result = results[0]
            x = 1

            options.stdout.write(
                "%sm1\tm2\tstatus\tlnL1\tnp1\tlnl2\tnp2\tP-value\n" % prefix_header)

            for result in results[1:]:

                if options.mode == "pairs":
                    reference_result = last_result
                    reference_id = x - 1
                elif options.mode == "1xn":
                    reference_result = first_result
                    reference_id = 0

                if reference_result.mNumParameters >= result.mNumParameters:
                    if options.loglevel >= 1:
                        options.stdlog.write("number of parameters of full model not increased (null=%i, full=%i).\n" % (
                            reference_result.mNumParameters, result.mNumParameters))
                    continue

                lrt = Stats.doLogLikelihoodTest(
                    result.mLogLikelihood, result.mNumParameters,
                    reference_result.mLogLikelihood, reference_result.mNumParameters,
                    options.significance_threshold)

                if lrt.mPassed:
                    c = "passed"
                else:
                    c = "failed"

                options.stdout.write("%s%i\t%i\t%s\t%f\t%i\t%f\t%i\t%5.2e\n" % (prefix_row, reference_id, x, c,
                                                                                lrt.mFullLogLikelihood, lrt.mFullNumParameters,
                                                                                lrt.mNullLogLikelihood, lrt.mNullNumParameters,
                                                                                lrt.mProbability,
                                                                                ))

                last_result = result
                x += 1

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
