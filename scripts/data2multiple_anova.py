'''
data2multiple_anova.py - compute multiple ANOVA between samples
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: data2multiple_anova.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms.")
    parser.add_option("-t", "--tree-nh-file", dest="filename_tree", type="string",
                      help="filename with tree(s).")
    parser.add_option("--skip-header", dest="add_header", action="store_false",
                      help="do not add header to flat format.")
    parser.add_option("--output-with-header", dest="write_header", action="store_true",
                      help="write header and exit.")
    parser.add_option("--debug", dest="debug", action="store_true",
                      help="debug mode")
    parser.add_option("--display-tree", dest="display_tree", action="store_true",
                      help="display the tree")

    parser.add_option("-m", "--method", dest="methods", type="choice", action="append",
                      choices=("contrasts", "spearman", "pearson", "compute"),
                      help="methods to perform on contrasts.")

    parser.set_defaults(
        columns="all",
        filename_tree=None,
        add_header=True,
        write_header=False,
        debug=False,
        methods=[],
        value_format="%6.4f",
        pvalue_format="%e",
        display_tree=False,
    )

    (options, args) = E.Start(parser, quiet=True)

    if options.columns not in ("all", "all-but-first"):
        options.columns = map(lambda x: int(x) - 1, options.columns.split(","))

    data = []

    options.filenames = args

    for filename in options.filenames:

        infile = open(filename, "r")
        table, headers = IOTools.readTable(
            infile, take=options.columns, headers=False)
        infile.close()

        data.append(table)

    fields = ["Df", "Sum Sq", "F value", "Pr(>F)", "Mean Sq"]

    options.stdout.write("set1\tset2")
    for field in fields:
        options.stdout.write("\t%s" % field)
    options.stdout.write("\n")

    # CODE needs to be refactored for rpy2 usage

    for x in range(len(data)):

        for y in range(x + 1, len(data)):

            rpy.set_default_mode(rpy.NO_CONVERSION)

            factors = ["x"] * len(data[x][:, 0]) + ["y"] * len(data[y][:, 0])
            values = list(data[x][:, 0]) + list(data[y][:, 0])

            linear_model = R.lm(
                R("y ~ x"), data=R.data_frame(x=factors, y=values))
            rpy.set_default_mode(rpy.BASIC_CONVERSION)
            result = R.anova(linear_model)

            options.stdout.write(
                "%s\t%s" % (options.filenames[x], options.filenames[y]))
            for field in fields:
                options.stdout.write("\t%s" % str(result[field]))
            options.stdout.write("\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
