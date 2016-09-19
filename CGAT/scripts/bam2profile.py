'''
bam2profile.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a :term:`bam` file as input and counts the number of reads
in features.  Output is a :term:`bed` formatted file
for each interval with the relevant counts of reads mapping to each interval.
Reads can be selected on the basis of :term:`tags`, such as NH or NM and a
threshold for the proportion of reads mapping to an interval before
it is output.  Read counting is done using mapping position within an interval.

The :term:`bam` file must have an :term:`index` file (.bai) in the same
directory.

Options
-------

+-------------------+--------------------------------------------------------+
|--condition        |used to specify whether a range or categorical filter   |
|                   |is used to select reads for counting.  The default      |
|                   |behaviour with `cat` is to output by contig.            |
+-------------------+--------------------------------------------------------+
|--interval         |size of intervals over which to count reads             |
+-------------------+--------------------------------------------------------+
|--tag              |SAM format tag used to count reads, e.g. to only ouput  |
|                   |reads that are multi-mapping use the NH tag and supply  |
|                   |a float value to the '--threshold' option               |
+-------------------+--------------------------------------------------------+
|--threshold        |The cut-off value above which regions with this number  |
|                   |of reads mapping to it will be ouput.  e.g. for regions |
|                   |where > 1% of reads are multimapping use --tag NH and   |
|                   |--threshold 0.1                                         |
+-------------------+--------------------------------------------------------+
|--clip             |Counts reads that are clipped based on the CIGAR        |
|                   |string. Combine with `--threshold` to select regions    |
|                   |with proportion reads clipped > threshold.  Choices are |
|                   |soft, hard, or both                                     |
+-------------------+--------------------------------------------------------+


Usage
-----

For example to output only 1000bp intervals that have >10% of reads that
are multi-mapping::

  python bam2feature.py example.bam
         --condition range
         --interval 1000
         --tag NH
         --threshold 0.1
         > example.bed
  
To just count reads within 100bp intervals you can supply the command thus::

  python bam2feature.py example.bam
  --condition range
  --interval 100
  > example.bed

The output files are a :term:bedGraph format file with each interval that has
>0 reads aligned to it::

  track type=bedGraph name="NH intervals" description="Intervals with NH > 1.00%" visibility=full priority=20
  chr1   24690001  24700001   7.00
  chr11  58110001  58120001   12.00
  chr12  63790001  6938000    12.00
  chr17  24860001  2487000    13.00
  chr18  67400001  6741000    13.00

The default behaviour is to output the sum of reads across regions in
10kbp windows

Type::

   python cgat_script_template.py --help

for command line help.


Command line options
--------------------

'''

import sys
import re
import pysam
import numpy
from functools import partial
from collections import deque
import CGAT.Experiment as E


class ProfileGenerator(object):

    '''basic profile generator from samfile
    using a list of Indicator functions to
    check if a read shoud be recorded or not.

    This generator currently assumes that all
    alignments have the same length.
    '''

    def __init__(self, indicators):
        self.indicators = indicators

    def __call__(self, samfile):

        last_pos = None
        last_tid = -1
        indicators = self.indicators

        running_totals = numpy.zeros(len(indicators),
                                     dtype=numpy.int)

        read_ends = [deque() for x in range(len(indicators))]

        def wrapUp(read_ends, running_totals, last_tid, last_pos, pos):
            '''once a new position is encountered, output
            profiles up to this position.'''

            results = []
            try:
                mi = min([x[0] for x in read_ends if x])
            except ValueError:
                # all queues are empty
                mi = pos

            if pos <= mi:
                results.append(
                    ((last_tid,
                      last_pos,
                      pos),
                     running_totals))

            while mi < pos:
                if last_pos != mi:
                    results.append(
                        ((last_tid,
                          last_pos,
                          pos),
                         running_totals))
                # remove elements
                for idx, ends in enumerate(read_ends):
                    while ends and ends[0] == mi:
                        ends.popleft()
                        running_totals[idx] -= 1
                last_pos = mi

                try:
                    mi = min([x[0] for x in read_ends if x])
                except ValueError:
                    # all queues are empty
                    break
            return results

        # iterate over reads in profiles
        for read in samfile:
            pos = read.pos

            if last_tid >= 0 and (last_tid != read.tid or pos != last_pos):
                # at end of chromosome
                if last_tid != read.tid:
                    # set pos to maxmimum end position
                    pos = max([x[-1] for x in read_ends if x]) + 1

                for result in wrapUp(read_ends, running_totals,
                                     last_tid, last_pos, pos):
                    yield result

                # reset at end of chromosome
                if last_tid != read.tid:
                    running_totals = numpy.zeros(len(indicators),
                                                 dtype=numpy.int)
                    read_ends = [deque() for x in range(len(indicators))]

            # count reads with indicators and record
            # end of alignments
            for idx, indicator in enumerate(indicators):
                if indicator(read):
                    running_totals[idx] += 1
                    read_ends[idx].append(read.aend)

            last_tid = read.tid
            last_pos = read.pos

        # should never fail if first indicator is the
        # number of reads.
        pos = max([x[-1] for x in read_ends if x]) + 1
        for result in wrapUp(read_ends, running_totals,
                             last_tid, last_pos, pos):
            yield result


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--indicator", dest="indicators", type="string",
                      action='append',
                      help="indicators [default=%default]")

    parser.set_defaults(indicators=[],
                        )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        args.append("-")

    # set up indicator functions
    indicators = []
    for indicator in options.indicators:
        if re.search("[=<>]", indicator):
            tag, operator, value = re.match(
                "([^=><]+)\s*([=><]+)\s*([^=><]+)", indicator).groups()
            # add more tags here
            # or do type conversion based on type of value
            if tag in ("NM", "NH"):
                value = int(value)

            if operator == "=":
                indicators.append(partial(
                    lambda tag, value, read: read.opt(tag) == value,
                    tag, value))
            elif operator == "<":
                indicators.append(partial(
                    lambda tag, value, read: read.opt(tag) < value,
                    tag, value))
            elif operator == ">":
                indicators.append(partial(
                    lambda tag, value, read: read.opt(tag) > value,
                    tag, value))
            elif operator == "<=":
                indicators.append(partial(
                    lambda tag, value, read: read.opt(tag) <= value,
                    tag, value))
            elif operator == ">=":
                indicators.append(partial(
                    lambda tag, value, read: read.opt(tag) >= value,
                    tag, value))

    if len(indicators) == 0:
        raise ValueError('please specify at least one indicator')

    # return total counts as first column
    indicators.insert(0, lambda read: True)
    E.info("reading input file %s" % args[0])
    samfile = pysam.Samfile(args[0], "r")

    profile_generator = ProfileGenerator(indicators)

    outfile = options.stdout

    outfile.write("contig\tstart\tend\tnreads\t%s\n" %
                  "\t".join(options.indicators))

    # simply output profile, alternatively, profiles
    # could be aggregated here as well.
    for coords, counts in profile_generator(samfile):
        tid, start, end = coords
        outfile.write("%i\t%i\t%i\t%s\n" %
                      (tid, start, end,
                       "\t".join(map(str, counts))))

    # write footer and output benchmarkinformation.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
