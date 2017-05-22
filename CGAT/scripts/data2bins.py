'''
data2bins.py - split table according to values in a column
==========================================================

:Tags: Python

Purpose
-------

Separate data in table according to a column.

Missing values (na) are ignored.

Reads data from stdin unless infile is given.

Usage
-----

Example::

   python data2bins.py --help

Type::

   python data2bins.py --help

for command line help.

Command line options
--------------------

'''
import sys
import os
import math
import bisect
import CGAT.Experiment as E
import CGAT.CSV as CSV


class Outputter:

    def __init__(self, filename, headers=None):
        self.mFilename = filename
        self.mOutfile = IOTools.openFile(filename, "w")
        self.mCounts = 0
        if headers:
            self.mOutfile.write("\t".join(headers) + "\n")

    def write(self, data):
        self.mOutfile.write("\t".join(map(str, data)) + "\n")
        self.mCounts += 1

    def __del__(self):
        self.mOutfile.close()
        if self.mCounts == 0:
            os.remove(self.mFilename)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: data2bins.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("--column", dest="column", type="int",
                      help="column to split on.")

    parser.add_option("--num-bins", dest="num_bins", type="int",
                      help="number of bins to create.")

    parser.add_option("--method", dest="method", type="choice",
                      choices=("equal-sized-bins",),
                      help="method to use to bin data.")

    parser.add_option("--no-headers", dest="has_headers", action="store_false",
                      help="matrix has no row/column headers.")

    parser.add_option("-p", "--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="OUTPUT filename with histogram information on aggregate coverages [%default].")

    parser.set_defaults(
        has_headers=True,
        method="equal-sized-bins",
        column=1,
        num_bins=4,
        output_filename_pattern="bin%i",
    )

    (options, args) = E.Start(parser)
    options.column -= 1

    if args:
        if args[0] == "-":
            infile = sys.stdin
        else:
            infile = IOTools.openFile(args[0], "r")
    else:
        infile = sys.stdin

    fields, data = CSV.readTable(infile)

    c = options.column
    values = [float(x[c]) for x in data]

    bins = []

    if options.method == "equal-sized-bins":
        increment = int(math.floor(float(len(values)) / options.num_bins))
        indices = list(range(0, len(values)))
        indices.sort(key=lambda x: values[x])
        for x in range(len(values)):
            values[indices[x]] = x
        bins = list(range(0, len(values) - increment, increment))

    elif options.method == "pass":
        pass

    E.debug("bins=%s" % str(bins))

    outputters = []
    for x in range(0, len(bins)):
        outputters.append(
            Outputter(options.output_filename_pattern % x, fields))

    # output tables
    for x in range(0, len(data)):
        bin = bisect.bisect(bins, values[x]) - 1
        outputters[bin].write(data[x])

    # stats
    if options.loglevel >= 1:
        options.stdlog.write("# bin\tstart\tcounts\tfilename\n")
        for x in range(0, len(bins)):
            options.stdlog.write("# %i\t%f\t%i\t%s\n" % (
                x, bins[x], outputters[x].mCounts, outputters[x].mFilename))

    E.info("ninput=%i, noutput=%i" %
           (len(data), sum((x.mCounts for x in outputters))))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
