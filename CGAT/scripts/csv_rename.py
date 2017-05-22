'''
csv_rename.py - rename columns in a table
=========================================

:Tags: Python

Purpose
-------

rename columns in a csv file

Usage
-----

Example::

   csv_rename.py gene=id < stdin

Type::

   python csv_rename.py --help

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
        version="%prog version: $Id: csv_rename.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-r", "--remove", dest="remove", action="store_true",
                      help="remove specified columns, keep all others.")

    parser.add_option("-u", "--unique", dest="unique", action="store_true",
                      help="output rows are uniq.")

    parser.add_option("-f", "--filename-fields", dest="filename_fields", type="string",
                      help="filename with field information.")

    parser.set_defaults(
        filename_fields=None,
    )

    (options, args) = E.Start(parser,
                              add_csv_options=True)
    mapper = {}
    for x in args:
        a, b = x.split("=")
        mapper[a.strip()] = b.strip()

    while 1:
        line = options.stdin.readline()

        if not line:
            E.Stop()
            sys.exit(0)

        if line[0] == "#":
            options.stdout.write(line)
            continue

        break

    header = []
    nreplaced = 0
    for x in line[:-1].split():
        if x in mapper:
            nreplaced += 1
            header.append(mapper[x])
        else:
            header.append(x)

    options.stdout.write("\t".join(header) + "\n")
    nlines = 0
    for line in options.stdin:
        nlines += 1
        options.stdout.write(line)

    if options.loglevel >= 1:
        ninput = len(header)
        noutput = ninput
        options.stdout.write("# ninput=%i, noutput=%i, nreplaced=%i, nlines=%i\n" % (
            ninput, noutput, nreplaced, nlines))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
