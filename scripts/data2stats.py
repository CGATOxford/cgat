'''
data2stats.py - summary statistics on table columns/rows
========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script computes summary statistics (min, mean, median)
on rows/columns of a table.

Usage
-----

Example::

   python data2stats.py < table.in > stats.out

Type::

   python data2stats.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import scipy
import scipy.stats
import CGAT.Experiment as E


def PrintValues(outfile, values,  options, prefix="", titles=None):

    if options.flat or options.aggregate_column:

        if options.add_header:
            if prefix:
                outfile.write("prefix\t")

            if titles:
                outfile.write("column\t")

            print "\t".join(("nval", "min", "max", "mean", "median", "stddev", "sum", "q1", "q3"))

        for x in range(len(values)):

            vals = values[x]

            if len(vals) == 0:

                if options.output_empty:
                    if titles:
                        outfile.write(titles[x] + "\t")
                    if prefix:
                        outfile.write(prefix + "\t")

                    outfile.write("0" + "\tna" * 8 + "\n")

                continue

            if titles:
                outfile.write(titles[x] + "\t")
            if prefix:
                outfile.write(prefix + "\t")

            vals.sort()
            if len(vals) > 4:
                q1 = options.value_format % vals[len(vals) // 4]
                q3 = options.value_format % vals[len(vals) * 3 // 4]
            else:
                q1 = options.value_format % vals[0]
                q3 = options.value_format % vals[-1]

            outfile.write("\t".join(("%i" % len(vals),
                                     options.value_format % float(
                                         min(vals)),
                                     options.value_format % float(
                                         max(vals)),
                                     options.value_format % scipy.mean(
                                         vals),
                                     options.value_format % scipy.median(
                                         vals),
                                     options.value_format % scipy.std(vals),
                                     options.value_format % reduce(
                                         lambda x, y: x + y, vals),
                                     q1, q3,
                                     )) + "\n")

    else:

        if titles:
            outfile.write("category\t%s" % string.join(titles, "\t") + "\n")

        outfile.write(
            "count\t%s" % (string.join(
                map(lambda v: "%i" % len(v), values), "\t")) + "\n")
        outfile.write(
            "min\t%s" % (string.join(
                map(lambda v: options.value_format % min(v),
                    values), "\t")) + "\n")
        outfile.write(
            "max\t%s" % (string.join(
                map(lambda v: options.value_format % max(v),
                    values), "\t")) + "\n")
        outfile.write(
            "mean\t%s" % (string.join(
                map(lambda v: options.value_format % scipy.mean(v),
                    values), "\t")) + "\n")
        outfile.write(
            "median\t%s" % (string.join(
                map(lambda v: options.value_format % scipy.median(v),
                    values), "\t")) + "\n")
        outfile.write(
            "stddev\t%s" % (string.join(
                map(lambda v: options.value_format % scipy.std(v),
                    values), "\t")) + "\n")
        outfile.write(
            "sum\t%s" % (string.join(
                map(lambda v: options.value_format %
                    reduce(lambda x, y: x + y, v), values), "\t")) + "\n")
        outfile.write(
            "q1\t%s" % (string.join(
                map(lambda v: options.value_format %
                    scipy.stats.scoreatpercentile(v, per=25),
                    values), "\t")) + "\n")
        outfile.write(
            "q3\t%s" % (string.join(
                map(lambda v: options.value_format %
                    scipy.stats.scoreatpercentile(v, per=75),
                    values), "\t")) + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take for calculating histograms.")
    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum value for histogram.")
    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum value for histogram.")
    parser.add_option("--scale", dest="scale", type="float",
                      help="scale values.")
    parser.add_option("-a", "--aggregate-column", dest="aggregate_column",
                      type="int",
                      help="use column to aggregate.")
    parser.add_option("-i", "--no-title", dest="titles", action="store_false",
                      help="do not use supplied column titles.")
    parser.add_option("-e", "--header-names", dest="headers", type="string",
                      help="headers.")
    parser.add_option("-r", "--rows", dest="rows", action="store_true",
                      help="data is in rows.")
    parser.add_option("--ignore-zeros", dest="ignore_zeros",
                      action="store_true",
                      help="ignore zero values.")
    parser.add_option("-f", "--format", dest="value_format", type="string",
                      help="number format.")
    parser.add_option("-x", "--flat-output", dest="flat", action="store_true",
                      help="flat format.")
    parser.add_option("--skip-header", dest="add_header", action="store_false",
                      help="do not add header to flat format.")
    parser.add_option("--output-with-header", dest="write_header",
                      action="store_true",
                      help="write header and exit.")
    parser.add_option("--skip-empty", dest="output_empty",
                      action="store_false",
                      help="do not output empty columns.")
    parser.add_option("--output-empty", dest="output_empty",
                      action="store_true",
                      help="output empty columns.")

    parser.set_defaults(
        columns="all",
        min_value=None,
        max_value=None,
        scale=None,
        aggregate_column=None,
        titles=True,
        headers=None,
        rows=False,
        value_format="%6.4f",
        flat=False,
        test="%5.2f",
        add_header=True,
        write_header=False,
        ignore_zeros=False,
        output_empty=False,
        separator="\t"
    )

    (options, args) = E.Start(parser, quiet=True)

    if options.columns not in ("all", "all-but-first", "variable"):
        options.columns = map(lambda x: int(x) - 1, options.columns.split(","))

    if options.headers:
        options.headers = options.headers.split(",")

    # write header for flat output
    if options.write_header:
        options.stdout.write("\t".join(("nval", "min", "max",
                                        "mean", "median", "stddev",
                                        "sum", "q1", "q3")) + "\n")
        return

    # retrieve histogram
    lines = filter(lambda x: x[0] != "#", sys.stdin.readlines())

    outfile = options.stdout

    if len(lines) > 0:

        ncols = len(string.split(lines[0][:-1], "\t"))

        if options.columns == "all":
            options.columns = range(0, ncols)
        elif options.columns == "all-but-first":
            options.columns = range(1, ncols)
        elif options.columns == "variable":
            pass

        if options.rows:

            # ignore first value: is row title
            if options.columns != "variable":
                del options.columns[0]

            if options.titles:
                del lines[0]

            # write header for flat output
            if options.flat:
                if options.headers:
                    head = options.headers[0]
                else:
                    head = "row"
                options.stdout.write("\t".join(
                    (head, "nval", "min", "max", "mean", "median",
                     "stddev", "sum", "q1", "q3")) + "\n")
                options.add_header = False

            for l in lines:
                data = l[:-1].split(options.separator)

                if options.columns == "variable":
                    vals = data[1:]
                else:
                    vals = map(lambda x: data[x], options.columns)

                # remove unknown values
                vals = [float(x)
                        for x in vals if x and x.lower() not in ("na", "nan")]

                if options.ignore_zeros:
                    vals = [x for x in vals if x != 0.0]

                # now convert to float
                vals = map(float, vals)

                PrintValues(outfile, [vals], options, data[0])

        else:

            last_aggregate = None

            if options.titles:
                data = lines[0][:-1].split("\t")

                if not options.headers:
                    options.headers = map(lambda x: data[x], options.columns)
                    del lines[0]

                if options.aggregate_column is not None:
                    outfile.write(
                        "category\t%s" % "\t".join(options.headers) + "\n")

            vals = [[] for x in range(len(options.columns))]

            for l in lines:

                data = l[:-1].split("\t")

                for c in range(len(options.columns)):

                    try:
                        val = string.atof(data[options.columns[c]])

                    except IndexError:
                        E.warn("IndexError in line: %s" % l[:-1])
                        continue
                    except ValueError:
                        continue

                    if options.aggregate_column is not None:

                        if last_aggregate != data[options.aggregate_column]:

                            if last_aggregate:
                                PrintValues(
                                    outfile, vals, options, last_aggregate)

                            vals = [[] for x in range(len(options.columns))]
                            last_aggregate = data[options.aggregate_column]

                    if options.scale:
                        val *= options.scale

                    if options.max_value is not None \
                       and val > options.max_value:
                        val = options.max_value

                    if options.min_value is not None \
                       and val < options.min_value:
                        val = options.min_value

                    vals[c].append(val)

            lines = None

            # remove empty columns
            nvals = []
            titles = []
            for c in range(len(options.columns)):
                if vals[c] or options.output_empty:
                    nvals.append(vals[c])

                    if options.headers:
                        titles.append(options.headers[c])

            vals = nvals

            PrintValues(outfile, vals, options, last_aggregate, titles)

    else:
        if options.titles:
            titles = ["missing", ]
        else:
            titles = []

        if options.output_empty:
            PrintValues(outfile, [[], ], options, None, titles)

    if options.loglevel >= 1:
        E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
