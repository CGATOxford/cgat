'''
jalview.py - build annotation for viewing in jalview
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will build annotations for a multiple alignment
such that both can be visualized with :file:`jalview`.

Usage
-----

Example::

   python jalview.py --help

Type::

   python jalview.py --help

for command line help.

Command line options
--------------------

'''
import sys
import time

import CGAT.Experiment as E
import CGAT.Mali as Mali


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: jalview.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("list2annotation", ),
                      help="methods.")

    parser.add_option("--filename-mali", dest="filename_mali", type="string",
                      help="filename with multiple alignment used for calculating sites - used for filtering")

    parser.add_option("--jalview-title", dest="jalview_title", type="string",
                      help="title for jalview annotation.")

    parser.set_defaults(
        method=None,
        jalview_symbol="*",
        jalview_title="anno",
        filename_mali=None,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if not options.filename_mali:
        raise "please specify a multiple alignment."

    mali = Mali.Mali()
    mali.readFromFile(open(options.filename_mali, "r"))

    if options.method == "list2annotation":

        options.stdout.write("JALVIEW_ANNOTATION\n")
        options.stdout.write("# Created: %s\n\n" %
                             (time.asctime(time.localtime(time.time()))))

        codes = [""] * mali.getWidth()

        first = True
        for line in sys.stdin:
            if line[0] == "#":
                continue
            if first:
                first = False
                continue

            position = int(line[:-1].split("\t")[0])
            codes[position - 1] = options.jalview_symbol

        options.stdout.write("NO_GRAPH\t%s\t%s\n" %
                             (options.jalview_title, "|".join(codes)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
