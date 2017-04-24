'''fastqs2fastqs.py - manipulate (merge/reconcile) fastq files
=============================================================

:Tags: Genomics NGS FASTQ FASTQ Manipulation

Purpose
-------

This script manipulates multiple fastq files and outputs
new fastq files. Currently only the method ``reconcile``
is implemented.

reconcile
+++++++++

Reconcile reads from a pair of fastq files.

This method takes two fastq files and outputs two fastq files such
that all reads in the output are present in both output files.

The typical use case is that two fastq files containing the first and
second part of a read pair have been independently filtered, for
example by quality scores, truncation, etc. As a consequence some
reads might be missing from one file but not the other. The reconcile
method will output two files containing only reads that are common to
both files.

The two files must be sorted by read identifier.

Example input, read2 and read3 are only present in either of the
files:

   # File1        # File 2

   @read1         @read1
   AAA            AAA
   +              +
   !!!            !!!
   @read2         @read3
   CCC            TTT
   +              +
   !!!            !!!
   @read4         @read4
   GGG            GGG
   +              +
   !!!            !!!

Example output, only the reads common to both files are output::

   # File1        # File 2

   @read1         @read1
   AAA            AAA
   +              +
   !!!            !!!
   @read4         @read4
   GGG            GGG
   +              +
   !!!            !!!

Usage
-----

Example::

   python fastqs2fastqs.py \
            --method=reconcile \
            --output-filename-pattern=myReads_reconciled.%s.fastq \
            myReads.1.fastq.gz myReads.2.fastq.gz

In this example we take a pair of fastq files, reconcile by read
identifier and output 2 new fastq files named
``myReads_reconciled.1.fastq.gz`` and
``myReads_reconciled.2.fastq.gz``.

Type::

   python fastqs2fastqs.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


class PatternGetter:

    def __init__(self, pattern):
        self.pattern = re.compile(pattern)

    def __call__(self, id):
        return self.pattern.search(id).groups()[0]


def plain_getter(id):
    return id


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=('reconcile',),
                      help="method to apply [default=%default].")

    parser.add_option(
        "-c", "--chop-identifier", dest="chop", action="store_true",
        help="whether or not to trim last character of the  "
        "sequence name. For example sometimes ids in the first "
        "file in the pair will end with \1 and the second "
        "with \2. If --chop-identifier is not specified "
        "then the results will be wrong [default=%default].")

    parser.add_option(
        "-u", "--unpaired", dest="unpaired", action="store_true",
        help="whether or not to write out unpaired reads "
        "to a separate file")

    parser.add_option(
        "--id-pattern-1", dest="id_pattern_1",
        help="If specified will use the first group from the"
        "pattern to determine the ID for the first read",
        default=None)

    parser.add_option(
        "--id-pattern-2", dest="id_pattern_2",
        help="As above but for read 2",
        default=None)

    parser.add_option(
        "-o", "--output-filename-pattern",
        dest="output_pattern", type="string",
        help="pattern for output files [default=%default].")

    parser.set_defaults(
        method="reconcile",
        chop=False,
        unpaired=False,
        output_pattern="%s.fastq.gz",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError(
            "please supply at least two fastq files on the commandline")

    fn1, fn2 = args
    c = E.Counter()

    if options.id_pattern_1:
        id1_getter = PatternGetter(options.id_pattern_1)
    else:
        id1_getter = plain_getter

    if options.id_pattern_2:
        id2_getter = PatternGetter(options.id_pattern_2)
    else:
        id2_getter = plain_getter

    if options.method == "reconcile":

        # IMS: switching to no store second set of read names and only use
        # lazily. Since generators don't have a size must keep track
        id_lengths = {fn1: 0, fn2: 0}

        def getIds(infile, id_getter=plain_getter):
            '''return ids in infile.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]:
                    break
                r = id_getter(l[0].split()[0])
                # decide if to chop read number off
                id_lengths[infile.name] += 1
                if options.chop:
                    yield r[:-1]
                else:
                    yield r

        def write(outfile, infile, take, unpaired_file=None,
                  id_getter=plain_getter):
            '''filter fastq files with ids in take.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]:
                    break
                r = id_getter(l[0].split()[0])
                if options.chop:
                    r = r[:-1]
                if r not in take:
                    if unpaired_file is None:
                        continue
                    else:
                        unpaired_file.write("\n".join(l) + "\n")
                else:
                    outfile.write("\n".join(l) + "\n")

        E.info("reading first in pair")
        inf1 = IOTools.openFile(fn1)
        ids1 = set(getIds(inf1, id1_getter))

        E.info("reading second in pair")
        inf2 = IOTools.openFile(fn2)
        # IMS: No longer keep as a set, but lazily evaluate into intersection
        # leads to large memory saving for large inf2, particularly if
        # inf1 is small.
        ids2 = getIds(inf2, id2_getter)
        take = ids1.intersection(ids2)

        E.info("first pair: %i reads, second pair: %i reads, "
               "shared: %i reads" %
               (id_lengths[fn1],
                id_lengths[fn2],
                len(take)))

        if options.unpaired:
            unpaired_filename = IOTools.openFile(
                options.output_pattern % "unpaired", "w")
        else:
            unpaired_filename = None

        with IOTools.openFile(options.output_pattern % "1", "w") as outf:
            inf = IOTools.openFile(fn1)
            E.info("writing first in pair")
            write(outf, inf, take, unpaired_filename, id1_getter)

        with IOTools.openFile(options.output_pattern % "2", "w") as outf:
            inf = IOTools.openFile(fn2)
            E.info("writing second in pair")
            write(outf, inf, take, unpaired_filename, id2_getter)

        if options.unpaired:
            unpaired_filename.close()

    # write footer and output benchmark information.
    E.info("%s" % str(c))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
