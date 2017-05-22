'''
modify_table.py - 
======================================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python modify_table.py --help

Type::

   python modify_table.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import CGAT.Experiment as E
from functools import reduce

"""read in data and append columns to a density histogram

-> relative frequencies
-> cumulative counts and frequencies in both directions

'#' at start of line is a comment
"""

parser = E.OptionParser(
    version="%prog version: $Id$")



def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take from table.")

    parser.add_option("-m", "--method", dest="methods", type="string",
                      help="methods to apply to columns.",
                      action="append")

    parser.add_option("-e", "--echo", dest="echo",
                      action="store_true",
                      help="echo columns not taken.")

    parser.add_option("-r", "--replace", dest="replace",
                      action="store_true",
                      help="replace orginial values.")

    parser.set_defaults(
        columns="1",
        echo=False,
        replace=False,
        format="%5.2f",
        methods=[])

    (options, args) = E.Start(parser)

    options.columns = [int(x) - 1 for x in options.columns.split(",")]

    print(E.GetHeader())
    print(E.GetParams())

    vals = []

    # retrieve histogram
    lines = [x for x in sys.stdin.readlines() if x[0] != "#"]

    headers = lines[0][:-1].split("\t")
    del lines[0]

    notcolumns = [x for x in range(len(headers)) if x not in options.columns]

    data = [[] for x in range(len(headers))]

    for l in lines:
        d = l[:-1].split("\t")
        for c in options.columns:
            data[c].append(float(d[c]))
        for c in notcolumns:
            data[c].append(d[c])

    if len(data) == 0:
        raise ValueError("no data found")

    totals = [0] * len(headers)

    for c in options.columns:
        totals[c] = reduce(lambda x, y: x + y, data[c])

    new_columns = []
    new_headers = []

    if options.echo:
        for c in notcolumns:
            new_headers.append(headers[c])
            new_columns.append(data[c])

    for c in options.columns:
        if not options.replace:
            new_columns.append(data[c])
            new_headers.append(headers[c])

        for method in options.methods:
            if method == "normalize":
                new_columns.append([d / totals[c] for d in data[c]])
                new_headers.append("normalized")

    print(string.join(new_headers, "\t"))

    for d in zip(*new_columns):
        print(string.join(list(map(str, d)), "\t"))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
