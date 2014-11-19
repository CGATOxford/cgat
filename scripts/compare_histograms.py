'''
histograms2distance.py - compute Kullback-Leibler distance
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Output the Kullback-Leibler distance between two or
more distributions given as histograms.

Usage
-----

Example::

   python histograms2kl.py --help

Type::

   python histograms2kl.py --help

for command line help.

This script was formerly called :file:`compare_histograms.py`.

Command line options
--------------------

'''
import sys
import math
import numpy
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      help="method to use [kl=kullback-leibler]",
                      choices=("kl",))
    parser.add_option("-n", "--no-normalize", dest="normalize",
                      action="store_false",
                      help="do not normalize data")
    parser.add_option("-p", "--pseudocounts", dest="pseudocounts",
                      type="int",
                      help="pseudocounts to add.")
    parser.add_option("-f", "--number-format", dest="number_format",
                      type="string",
                      help="number format.")

    parser.set_defaults(
        method="kl",
        columns="all",
        headers=True,
        xrange=None,
        pseudocounts=1,
        normalize=True,
        number_format="%6.4f"
    )

    (options, args) = E.Start(parser,
                              add_pipe_options=True)

    if options.xrange:
        options.xrange = map(float, options.xrange.split(","))

    data, legend = IOTools.readTable(sys.stdin,
                                     numeric_type=numpy.float32,
                                     take=options.columns,
                                     headers=options.headers,
                                     truncate=options.xrange)

    nrows, ncols = data.shape

    # first: normalize rows
    for y in range(1, ncols):
        for x in range(nrows):
            data[x, y] = data[x, y] + float(options.pseudocounts)
        if options.normalize:
            t = numpy.sum(data[:, y])
            for x in range(nrows):
                data[x, y] = data[x, y] / t

    for x in range(1, len(legend) - 1):
        for y in range(x + 1, len(legend)):

            if options.method == "kl":
                d1 = 0.0
                d2 = 0.0
                for bin in range(nrows):
                    p = data[bin, x]
                    q = data[bin, y]
                    d1 += p * math.log(p / q)
                    d2 += q * math.log(q / p)

                options.stdout.write("%s\t%s\t%s\n" %
                                     (legend[x], legend[y],
                                      options.number_format % d1))
                options.stdout.write("%s\t%s\t%s\n" %
                                     (legend[y], legend[x],
                                      options.number_format % d2))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
