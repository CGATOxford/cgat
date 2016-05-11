'''
csv2csv.py - operate on tables
==============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

operate on tables.

Usage
-----

Example::

   python csv2csv.py --help

Type::

   python csv2csv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import csv
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$")

    parser.add_option(
        "-s", "--method=sort --sort-order", dest="sort", type="string",
        help="fields to take (in sorted order).")

    (options, args) = E.Start(parser, add_csv_options=True)

    reader = csv.DictReader(E.stdin, dialect=options.csv_dialect)

    if options.sort:
        fields = options.sort.split(",")
    else:
        fields = None

    writer = csv.DictWriter(E.stdout,
                            fields,
                            dialect=options.csv_dialect,
                            lineterminator=options.csv_lineterminator,
                            extrasaction='ignore')

    E.stdout.write("\t".join(fields) + "\n")

    for row in reader:
        row = IOTools.convertDictionary(row)
        writer.writerow(row)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
