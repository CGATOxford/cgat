'''
coverage2stats.py
=============================================

:Tags: Python

Purpose
-------

The output from bedtools genomecov is converted to a set of statistics
for each contig that is output. Statistics output are the mean and sd
of coverage for each contig. There is also the option to output a matrix
that represents a histogram of coverage over each contig.

Usage
-----

Example::

   python coverage2stats.py --help

Type::

   python coverage2stats.py --help

for command line help.

Documentation
-------------

The script assumes the input is from stdin and outputs the results to stdout.
The input is expected to be in tab delimited text format (no header line)

          <contig><base-position><coverage>

Command line options
--------------------

'''

import sys
import optparse
import numpy as np
import collections
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser(version="%prog version: $Id",
                                   usage=globals()["__doc__"])

    parser.add_option("--bin", dest="bin", action="store_true",
                      help="output average in bins across the interval")
    parser.add_option("-n", "--num-bins", dest="bin_number", type=int,
                      help="number of bins for coverage profile")
    parser.add_option("-o", "--output-filename-prefix",
                      dest="output_filename_prefix",
                      help="pattern to write coverage bins to")

    parser.set_defaults(
        bin=False, bin_number=10)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    inf = options.stdin

    coverage_result = collections.defaultdict(list)
    E.info("reading in coverage data")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contig, coverage = data[0], data[2]
        coverage_result[contig].append(coverage)
    E.info("read %i contigs" % len(list(coverage_result.keys())))

    options.stdout.write("contig\tcov_mean\tcov_sd\n")
    if options.bin:
        outf = IOTools.openFile(options.output_filename_prefix + ".binned",
                                "w")
        outf.write("%s" % "\t".join(
            [str(i) for i in range(1, options.bin_number + 1, 1)]) + "\n")
    for contig, coverage in coverage_result.items():
        coverage = list(map(float, coverage))
        options.stdout.write(
            "%s\t%s\t%s\n" % (contig,
                              str(np.mean(coverage)),
                              str(np.std(coverage))))
        if options.bin:
            bin_means = []
            bins = np.linspace(0, len(coverage), options.bin_number + 1)
            if len(coverage) < len(bins) - 1:
                E.warn("will not calculate coverage means for %s: too short" %
                       contig)
                continue
            for i in range(len(bins)):
                try:
                    bin_mean = np.mean(coverage[int(bins[i]):int(bins[i + 1])])
                except IndexError:
                    continue
                bin_means.append(bin_mean)
            outf.write(contig + "\t" + "\t".join(map(str, bin_means)) + "\n")
    outf.close()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
