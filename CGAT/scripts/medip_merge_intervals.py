'''
medip_merge_intervals.py - merge differentially methylated regions
==================================================================

:Tags: Python

Purpose
-------

This script takes the output of DESeq or EdgeR and merges
adjacent intervals that show a similar expression change.

Input is data like this::

   contig start end treatment_name  treatment_mean  treatment_std   control_name    control_mean    control_std     pvalue  qvalue  l2fold  fold    significant     status                                                                                 
   chr1 10000 11000        CD14    32.9785173324   0       CD4     41.7117152603   0       0.199805206526  1.0     0.338926100945  1.26481475319   0       OK                                                                                   
   chr1 14000 15000        CD14    9.32978709019   0       CD4     9.31489982941   0       1.0     1.0     -0.00230390372974       0.998404330063  0       OK                                                                                   
   chr1 15000 16000        CD14    9.04603350905   0       CD4     9.01484414416   0       1.0     1.0     -0.00498279072069       0.996552150193  0       OK                                                                                   
   chr1 16000 17000        CD14    0.457565479197  0       CD4     0.14910378845   0       0.677265200643  1.0     -1.61766129852  0.325863281276  0       OK                                                                                   

The second and third window would be merged, as

1. Their methylation levels are within 10% of each other.
2. They are both not differentially methylated.

It aggregates the following:

* mean values: average
* std values: max
* pvalue: max
* qvalue: max
* fold: min/max (depending on enrichment/depletion)
* l2fold: min/max (depending on enrichment/depletion)

The analysis outputs bed files with intervals that are
potentially activated in one of the conditions. Windows
with a positive fold change are collected in the ``treatment``,
while windows with a negative fold change are collected in the
``control``.

For methylation analysis, it might be more interesting
to report windows that are depleted (instead of enriched)
of signal. Thus, if the option ``--invert`` is given,
windows with a negative l2fold change are labeled ``treatment``.
Less methylation means that this region is "active" in the
``treatment`` condition.

Note that the input is assumed to be sorted by coordinate.

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
import re
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

DATA = collections.namedtuple(
    "DATA",
    "test_id contig start end treatment_name  treatment_mean  treatment_std "
    "control_name    control_mean    control_std     pvalue  qvalue  "
    "l2fold  fold    significant     status nintervals")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-o", "--min-overlap", dest="min_overlap", type="int",
                      help="minimum overlap")

    parser.add_option(
        "-w", "--pattern-window",
        dest="pattern_window", type="string",
        help="regular expression to extract window coordinates from "
        "test id [%default]")

    parser.add_option(
        "-i", "--invert", dest="invert", action="store_true",
        help="invert direction of fold change [%default]")

    parser.set_defaults(min_overlap=10,
                        invert=False,
                        pattern_window="(\S+):(\d+)-(\d+)"),

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    outfiles = IOTools.FilePool(options.output_filename_pattern)

    if options.invert:
        test_f = lambda l2fold: l2fold < 0
    else:
        test_f = lambda l2fold: l2fold > 0

    def read():

        rx_window = re.compile(options.pattern_window)
        # filter any of the DESeq/EdgeR message that end up at the top of the
        # output file

        for data in IOTools.iterate(options.stdin):

            contig, start, end = rx_window.match(data.test_id).groups()
            start, end = list(map(int, (start, end)))

            yield DATA._make((data.test_id,
                              contig, start, end,
                              data.treatment_name,
                              float(data.treatment_mean),
                              float(data.treatment_std),
                              data.control_name,
                              float(data.control_mean),
                              float(data.control_std),
                              float(data.pvalue),
                              float(data.qvalue),
                              float(data.l2fold),
                              float(data.fold),
                              int(data.significant),
                              data.status,
                              0))

    def grouper(data, distance=10):

        last = next(data)
        entries = [last]

        while 1:
            d = next(data)
            if d is None:
                break
            if d.contig == last.contig and d.start < last.start:
                raise ValueError("error not sorted by start")

            if ((d.contig != last.contig) or
                    (d.start - last.end > distance) or
                    (d.status != last.status) or
                    (d.significant != last.significant) or
                    (d.l2fold * last.l2fold < 0)):
                yield entries
                entries = []

            entries.append(d)
            last = d

        yield entries

    counter = E.Counter()

    options.stdout.write("\t".join(DATA._fields) + "\n")

    # set of all sample names - used to create empty files
    samples = set()

    # need to sort by coordinate
    all_data = list(read())
    all_data.sort(key=lambda x: (x.contig, x.start))

    group_id = 0

    for group in grouper(iter(all_data), distance=options.min_overlap):
        group_id += 1

        start, end = group[0].start, group[-1].end
        assert start < end, 'start > end: %s' % str(group)
        n = float(len(group))
        counter.input += n

        g = group[0]

        if g.l2fold < 0:
            l2fold = max([x.l2fold for x in group])
            fold = max([x.fold for x in group])
        else:
            l2fold = min([x.l2fold for x in group])
            fold = min([x.fold for x in group])

        outdata = DATA._make((
            str(group_id),
            g.contig, start, end,
            g.treatment_name,
            sum([x.treatment_mean for x in group]) / n,
            max([x.treatment_std for x in group]),
            g.control_name,
            sum([x.control_mean for x in group]) / n,
            max([x.control_std for x in group]),
            max([x.pvalue for x in group]),
            max([x.qvalue for x in group]),
            l2fold,
            fold,
            g.significant,
            g.status,
            int(n)))

        samples.add(g.treatment_name)
        samples.add(g.control_name)
        if g.significant:
            if test_f(g.l2fold):
                # treatment lower methylation than control
                outfiles.write(
                    g.treatment_name, "%s\t%i\t%i\t%i\t%f\n" % (
                        g.contig, g.start, g.end,
                        group_id,
                        sum([x.treatment_mean for x in group]) / n))

            else:
                outfiles.write(
                    g.control_name, "%s\t%i\t%i\t%i\t%f\n" % (
                        g.contig, g.start, g.end,
                        group_id,
                        sum([x.control_mean for x in group]) / n))

        options.stdout.write("\t".join(map(str, outdata)) + "\n")

        counter.output += 1

    # create empty files
    for sample in samples:
        outfiles.write(sample, "")

    outfiles.close()
    E.info("%s" % counter)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
