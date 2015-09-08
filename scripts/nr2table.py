'''
nr2table.py - convert description in nr fasta file to table
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts converts headers in ncbi `nr` fasta files to
a tabular format.

Usage
-----

Example::

   python nr2table.py --help

Type::

   python nr2table.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E
import CGAT.FastaIterator as FastaIterator

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: nr2table.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.set_defaults()

    (options, args) = E.Start(parser)

    iterator = FastaIterator.FastaIterator(sys.stdin)

    sequences = []

    ninput, noutput, nentries = 0, 0, 0

    options.stdout.write("gid\tsrc\tacc\tprotid\tannotation\tspecies\n")
    while 1:

        cur_record = iterator.next()

        if cur_record is None:
            break

        ninput += 1

        records = cur_record.title.split(chr(1))
        for record in records:

            a, anno = re.search("(\S+)\s+(.+)", record).groups()

            vals = a.split("|")

            try:
                gi, gid, src, acc, protid = vals
            except ValueError:
                raise "parsing error for record: %s" % record

            try:
                annotation, species = re.search(
                    "(.+)\s+\[(.*)\]", anno).groups()
            except AttributeError:
                annotation = anno.strip()
                species = ""

            annotation = re.sub("\t", " ", annotation)

            options.stdout.write("\t".join((gid, src, acc, protid,
                                            annotation, species)) + "\n")

            nentries += 1

        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write(
            "# ninput=%i, noutput=%i, nentries=%i\n" % (ninput, noutput, nentries))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
