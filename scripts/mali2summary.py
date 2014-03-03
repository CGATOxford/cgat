'''
mali2summary.py - compute summary stats on a mali
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

output summary information on a multiple alignment.

* column/row occupancy

Usage
-----

Example::

   python mali2summary.py --help

Type::

   python mali2summary.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import optparse
import math
import time

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali
import scipy
import scipy.stats


def analyzeMali(mali, options, prefix_row=""):

    if len(mali) == 0:
        raise "not analyzing empty multiple alignment"

    # count empty sequences
    row_data = map(lambda x: Mali.MaliData(
        x.mString, options.gap_chars, options.mask_chars), mali.values())
    col_data = map(lambda x: Mali.MaliData(
        x, options.gap_chars, options.mask_chars), mali.getColumns())

    if len(row_data) == 0 or len(col_data) == 0:
        return False

    if options.loglevel >= 2:
        for row in row_data:
            options.stdlog.write("# row: %s\n" % str(row))
        for col in col_data:
            options.stdlog.write("# col: %s\n" % str(col))

    options.stdout.write(prefix_row)

    # calculate average column occupancy
    col_mean = scipy.mean(map(lambda x: x.mNChars, col_data))
    col_median = scipy.median(map(lambda x: x.mNChars, col_data))
    length = mali.getLength()

    if float(int(col_median)) == col_median:
        options.stdout.write("%5.2f\t%5.2f\t%i\t%5.2f" % (col_mean, 100.0 * col_mean / length,
                                                          col_median, 100.0 * col_median / length))
    else:
        options.stdout.write("%5.2f\t%5.2f\t%5.1f\t%5.2f" % (col_mean, 100.0 * col_mean / length,
                                                             col_median, 100.0 * col_median / length))

    row_mean = scipy.mean(map(lambda x: x.mNChars, row_data))
    row_median = scipy.median(map(lambda x: x.mNChars, row_data))
    width = mali.getWidth()

    if float(int(row_median)) == row_median:
        options.stdout.write("\t%5.2f\t%5.2f\t%i\t%5.2f" % (row_mean, 100.0 * row_mean / width,
                                                            row_median, 100.0 * row_median / width))
    else:
        options.stdout.write("\t%5.2f\t%5.2f\t%5.1f\t%5.2f" % (row_mean, 100.0 * row_mean / width,
                                                               row_median, 100.0 * row_median / width))

    options.stdout.write("\n")

    return True

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: mali2summary.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm"),
                      help="input format of multiple alignment")

    parser.add_option("-a", "--alphabet", dest="alphabet", type="choice",
                      choices=("aa", "na"),
                      help="alphabet to use [default=%default].", )

    parser.add_option("-p", "--pattern-mali", dest="pattern_mali", type="string",
                      help="filename pattern for input multiple alignment files.")

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        mask_chars="nN",
        gap_chars="-.",
        alphabet="na",
        pattern_mali=None,
    )

    (options, args) = E.Start(parser)

    if options.pattern_mali:
        prefix_header = "prefix\t"
        prefix_row = "\t"
    else:
        prefix_header = ""
        prefix_row = ""

    options.stdout.write(
        "%sncol_mean\tpcol_mean\tncol_median\tpcol_median\tnrow_mean\tprow_mean\tnrow_median\tprow_median\n" % (prefix_header, ))

    ninput, nskipped, noutput, nempty = 0, 0, 0, 0

    if options.pattern_mali:

        ids, errors = IOTools.ReadList(sys.stdin)

        E.debug("read %i identifiers.\n" % len(ids))

        nsubstitutions = len(re.findall("%s", options.pattern_mali))

        for id in ids:

            filename = options.pattern_mali % tuple([id] * nsubstitutions)
            ninput += 1

            if not os.path.exists(filename):
                nskipped += 1
                continue

            # read multiple alignment in various formats
            mali = Mali.Mali()
            mali.readFromFile(open(filename, "r"), format=options.input_format)

            if mali.isEmpty():
                nempty += 1
                continue

            E.debug("read mali with %i entries from %s.\n" %
                    (len(mali), filename))

            if analyzeMali(mali, options, prefix_row="%s\t" % id):
                noutput += 1

    else:

        # read multiple alignment in various formats
        mali = Mali.Mali()
        mali.readFromFile(sys.stdin, format=options.input_format)
        ninput += 1

        if mali.isEmpty():
            nempty += 1
        else:
            E.debug("read mali with %i entries." % (len(mali)))

            if analyzeMali(mali, options, prefix_row=""):
                noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i, nempty=%i." %
           (ninput, noutput, nskipped, nempty))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
