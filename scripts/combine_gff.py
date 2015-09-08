'''
combine_gff.py - merge overlapping intervals
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

combine gff features as overlapping regions.

Usage
-----

Example::

   python combine_gff.py --help

Type::

   python combine_gff.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: combine_gff.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--format", dest="format",
                      help="output format.", type="choice", choices=("flat", "full", "first"))

    parser.set_defaults(
        format="flat"
    )

    (options, args) = E.Start(parser)

    last_e = None
    for line in sys.stdin:
        if line[0] == "#":
            continue
        if options.format in ("full", "first"):
            last_e = GTF.Entry()
        else:
            last_e = GTF.Entry()
        last_e.Read(line)
        break

    for line in sys.stdin:

        if line[0] == "#":
            continue

        if options.format in ("full", "first"):
            e = GTF.Entry()
        else:
            e = GTF.Entry()
        e.Read(line)

        if not GTF.Overlap(last_e, e):
            print str(last_e)
            last_e = e
        else:
            last_e.start = min(last_e.start, e.start)
            last_e.end = max(last_e.end, e.end)
            if options.format == "full":
                last_e.mInfo += " ; " + e.mInfo
            continue

    print str(last_e)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
