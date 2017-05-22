"""
bed2graph.py - compute the overlap graph between two bed files
==============================================================

:Tags: Python

Purpose
-------

This script ouputs a list of the names of all overlapping intervals 
between two bed files.

Usage
-----

Type::

   python bed2graph.py A.bed.gz B.bed.gz > graph.out

for command line help.

Command line options
--------------------

"""

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-o", "--output-section", dest="output", type="choice",
                      choices=("full", "name"),
                      help="output either ``full`` overlapping entries, only the ``name``s. [default=%default].")

    parser.set_defaults(
        output="full",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError("two arguments required")

    if args[0] == "-":
        infile1 = options.stdin
    else:
        infile1 = IOTools.openFile(args[0], "r")

    infile2 = IOTools.openFile(args[1], "r")

    idx = Bed.readAndIndex(infile2, with_values=True)

    output = options.output
    outfile = options.stdout

    if output == "name":
        outfile.write("name1\tname2\n")
        outf = lambda x: x.fields[0]
    else:
        outf = str

    for bed in Bed.iterator(infile1):
        try:
            overlaps = idx[bed.contig].find(bed.start, bed.end)
        except (KeyError, IndexError):
            # ignore missing contig and zero length intervals
            continue

        for o in overlaps:
            outfile.write("\t".join((outf(bed), outf(o[2]))) + "\n")

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
