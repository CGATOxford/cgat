'''
histogram2histogram.py - operations on histograms
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a histogram from stdin and
outputs a modified histogram on stdout. Modifications
are normalization, the cumulative, etc.

Usage
-----

Example::

   python histogram2histogram.py --help

Type::

   python histogram2histogram.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import numpy
from functools import reduce

USAGE = """python %s < stdin > stdout

read in data and append columns to a density histogram

-> relative frequencies
-> cumulative counts and frequencies in both directions

'#' at start of line is a comment

""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: histogram2histogram.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--is-int", dest="is_ints", action="store_true",
                      help="categories are integers.")
    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("append", "cumul", "rcumul", "normalize"),
                      help="method(s) to apply.")
    parser.add_option("--no-headers", dest="headers", action="store_false",
                      help="histogram has no headers.")
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for plotting.")
    parser.add_option("", "--truncate", dest="truncate", type="string",
                      help="truncate at range.")
    parser.add_option("", "--no-out-of-range", dest="cumulate_out_of_range", action="store_false",
                      help="add up bins out of range.")
    parser.add_option("--bin-format", dest="format_bin", type="string",
                      help="format for bins.")
    parser.add_option("--value-format", dest="format_val", type="string",
                      help="format for vals.")

    parser.set_defaults(
        is_ints=False,
        method="append",
        columns="all",
        headers=True,
        truncate=None,
        cumulate_out_of_range=True,
        format_bin="%6.4f",
        format_val="%6.4f",
    )

    (options, args) = E.Start(parser)

    # old histogram2histogram.py semantics - need to merged with newer
    # code below.
    if options.method == "append":

        vals = []

        # retrieve histogram
        lines = [x for x in sys.stdin.readlines() if x[0] != "#"]

        # check if first line contains a header
        d = lines[0][:-1].split("\t")[0]
        try:
            if options.is_ints:
                value = int(d)
            else:
                value = float(d)
        except ValueError:
            print("\t".join(
                (d, "counts", "frequency",
                 "cumulative counts", "increasing cumulative frequency",
                 "cumulative counts", "decreasing cumulative frequency")))
            del lines[0]

        data = [list(map(float, x[:-1].split("\t"))) for x in lines]

        if len(data) == 0:
            raise ValueError("no data found")

        total = float(reduce(lambda x, y: x + y, [x[1] for x in data]))

        cumul_down = int(total)
        cumul_up = 0

        if options.is_ints:
            form = "%i\t%i\t%6.4f\t%i\t%6.4f\t%i\t%6.4f"
        else:
            form = "%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f"

        for bin, val in data:
            percent = float(val) / total
            cumul_up += val
            percent_cumul_up = float(cumul_up) / total
            percent_cumul_down = float(cumul_down) / total

            print(form %
                  (bin, val, percent, cumul_up, percent_cumul_up,
                   cumul_down, percent_cumul_down))

            cumul_down -= val

    else:

        if options.truncate:
            options.truncate = list(map(float, options.truncate.split(",")))

        options.method = options.method.split(",")
        data, legend = IOTools.readTable(sys.stdin,
                                         numeric_type=numpy.float32,
                                         take=options.columns,
                                         headers=options.headers,
                                         truncate=options.truncate,
                                         cumulate_out_of_range=options.cumulate_out_of_range)

        nfields = len(legend)

        # note: because of MA, iteration makes copy of slices
        # Solution: inplace edits.
        nrows, ncols = data.shape

        for method in options.method:
            if method == "cumul":
                l = [0] * ncols
                for x in range(nrows):
                    for y in range(1, ncols):
                        data[x, y] += l[y]
                        l[y] = data[x, y]

            elif method == "rcumul":
                l = [0] * ncols
                for x in range(nrows - 1, 0, -1):
                    for y in range(1, ncols):
                        data[x, y] += l[y]
                        l[y] = data[x, y]

            elif method == "normalize":
                m = [0] * ncols
                for x in range(nrows):
                    for y in range(1, ncols):
                        # the conversion to float is necessary
                        m[y] = max(m[y], float(data[x, y]))

                for y in range(1, ncols):
                    if m[y] == 0:
                        m[y] = 1.0

                for x in range(nrows):
                    for y in range(1, ncols):
                        data[x, y] = data[x, y] / m[y]
            else:
                raise "unknown method %s" % method

        print("\t".join(legend))

        format = options.format_bin + "\t" + \
            "\t".join([options.format_val] * (nfields - 1))

        for d in data:
            print(format % tuple(d))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
