'''csv_cut.py - select columns from a table
========================================

:Tags: Python

Purpose
-------

extract named columns from a csv formatted table


.. todo::

   describe purpose of the script.

Usage
-----

Extract the two columns gene and length from a table in standard input::

   python csv_cut.py gene length < stdin

The script permits the use of patterns. For example, the command will
select the column gene and all columns that contain the part 'len'::

   python csv_cut.py gene %len% < stdin

Type::

   python csv_cut.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E
import csv
import six
import _csv
import hashlib
import CGAT.CSV as CSV


class UniqueBuffer:
    mKeys = {}

    def __init__(self, outfile):
        self.mOutfile = outfile

    def write(self, out):
        key = hashlib.md5(out).digest()
        if key not in self.mKeys:
            self.mKeys[key] = True
            self.mOutfile.write(out)


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
        large=False,
        filename_fields=None,
    )

    (options, args) = E.Start(parser,
                              add_csv_options=True,
                              quiet=True)

    input_fields = args

    if options.filename_fields:
        input_fields = [x[:-1].split("\t")[0] for x in [x for x in IOTools.openFile(options.filename_fields, "r").readlines() if x[0] != "#"]]

    if options.unique:
        outfile = UniqueBuffer(options.stdout)
    else:
        outfile = options.stdout

    while 1:
        line = options.stdin.readline()

        if not line:
            E.Stop()
            sys.exit(0)

        if line[0] == "#":
            continue

        first_line = line
        break

    old_fields = first_line[:-1].split("\t")

    fields = []
    for f in input_fields:
        # do pattern search
        if f[0] == "%" and f[-1] == "%":
            pattern = re.compile(f[1:-1])
            for o in old_fields:
                if pattern.search(o) and o not in fields:
                    fields.append(o)
        else:
            if f in old_fields:
                fields.append(f)

    if options.remove:
        fields = set(fields)
        fields = [x for x in old_fields if x not in fields]

    if options.large:
        reader = CSV.DictReaderLarge(CSV.CommentStripper(options.stdin),
                                     fieldnames=old_fields,
                                     dialect=options.csv_dialect)
    else:
        reader = csv.DictReader(CSV.CommentStripper(options.stdin),
                                fieldnames=old_fields,
                                dialect=options.csv_dialect)

    writer = csv.DictWriter(outfile,
                            fields,
                            dialect=options.csv_dialect,
                            lineterminator=options.csv_lineterminator,
                            extrasaction='ignore')

    print("\t".join(fields))

    first_row = True
    ninput, noutput, nerrors = 0, 0, 0

    while 1:
        ninput += 1
        try:
            row = six.next(reader)
        except _csv.Error as msg:
            options.stderr.write("# error while parsing: %s\n" % (msg))
            nerrors += 1
            continue
        except StopIteration:
            break
        if not row:
            break
        writer.writerow(row)
        noutput += 1

    E.info("ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput, nerrors))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
