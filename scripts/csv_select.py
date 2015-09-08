'''
csv_select.py - select rows from a table
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

extract rows from a csv-formatted table.

The select statement is a one-line, for example::

   csv_select.py "int(r['mC-foetal-sal-R4']) > 0" < in > out

Note the required variable name r for denoting field names. Please
also be aware than numeric values need to be converted first before
testing.

Usage
-----

Type::

   python csv_select.py --help

for command line help.

Command line options
--------------------

'''
import sys
import csv
import _csv
import CGAT.Experiment as E
import CGAT.CSV as CSV


class CommentStripper:

    """iterator class for stripping comments from file.
    """

    def __init__(self, file):
        self.mFile = file

    def __iter__(self):
        return self

    def next(self):
        while 1:
            line = self.mFile.readline()
            if not line:
                raise StopIteration
            if line[0] != "#":
                return line


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: csv_cut.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-r", "--remove", dest="remove", action="store_true",
                      help="remove specified columns, keep all others.")

    parser.add_option("-u", "--unique", dest="unique", action="store_true",
                      help="output rows are uniq.")

    parser.add_option("-l", "--large", dest="large", action="store_true",
                      help="large columns. Do not use native python CSV module [default=%default].")

    parser.add_option("-f", "--filename-fields", dest="filename_fields", type="string",
                      help="filename with field information.")

    parser.set_defaults(
        remove=False,
        unique=False,
        filename_fields=None,
    )

    (options, args) = E.Start(parser,
                              add_csv_options=True,
                              quiet=True)

    statement = " ".join(args)

    if options.large:
        reader = CSV.DictReaderLarge(CommentStripper(sys.stdin),
                                     dialect=options.csv_dialect)
    else:
        reader = csv.DictReader(CommentStripper(sys.stdin),
                                dialect=options.csv_dialect)

    exec "f = lambda r: %s" % statement in locals()

    counter = E.Counter()
    writer = csv.DictWriter(options.stdout,
                            reader.fieldnames,
                            dialect=options.csv_dialect,
                            lineterminator=options.csv_lineterminator)

    writer.writerow(dict((fn, fn) for fn in reader.fieldnames))

    while 1:
        counter.input += 1
        try:
            row = reader.next()
        except _csv.Error, msg:
            options.stderr.write("# error while parsing: %s\n" % (msg))
            counter.errors += 1
            continue
        except StopIteration:
            break

        if not row:
            break

        if f(row):
            writer.writerow(row)
            counter.output += 1
        else:
            counter.filtered += 1

    E.info("%s" % counter)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
