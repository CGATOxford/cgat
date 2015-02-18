'''
merge_tables.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

merge two tables with the same number of lines

Usage
-----

Example::

   python merge_tables.py --help

Type::

   python merge_tables.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import CGAT.Experiment as E

parser = E.OptionParser(
    version="%prog version: $Id$")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser.add_option("-t", "--table", dest="tables", type="string",
                      help="tables to merge.",
                      action="append")

    parser.set_defaults(
        tables=[])

    (options, args) = E.Start(parser)

    if len(options.tables) < 1:
        raise "please specify at least one table."

    files = []
    for t in options.tables:
        files.append(open(t, "r"))

    while 1:
        frags = []

        stop = False

        for f in files:
            l = f.readline()
            if not l:
                stop = True
                break
            frags.append(l[:-1])

        if stop:
            break

        print string.join(frags, "\t")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
