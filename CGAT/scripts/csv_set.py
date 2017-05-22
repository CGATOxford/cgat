'''
csv_set.py - set operations on a table
======================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python csv_set.py --help

Type::

   python csv_set.py --help

for command line help.

Command line options
--------------------

'''
import sys

import CGAT.Experiment as E
import CGAT.CSV as CSV
import hashlib


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

    parser = E.OptionParser(
        version="%prog version: $Id: csv_set.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-u", "--unique", dest="unique", action="store_true",
                      help="output rows are uniq.")

    parser.add_option("-1", "--join-fields1", dest="join_fields1", type="string",
                      help="join fields in first table.")
    parser.add_option("-2", "--join-fields2", dest="join_fields2", type="string",
                      help="join fields in second table.")
    parser.add_option("-m", "--method", dest="method", type="choice",
                      help="set operation to perform.", choices=("intersection", "rest", "union"))

    parser.set_defaults(
        remove=False,
        unique=False,
        join_fields1=None,
        join_fields2=None,
        method="intersection",
    )

    (options, args) = E.Start(parser, add_csv_options=True)

    if len(args) != 2:
        raise ValueError("please specify two files to join")

    if not options.join_fields1 or not options.join_fields2:
        raise ValueError("please specify at least one join field per table")

    options.join_fields1 = options.join_fields1.split(",")
    options.join_fields2 = options.join_fields2.split(",")

    options.filename1, options.filename2 = args

    fields1, table1 = CSV.readTable(open(options.filename1, "r"))
    fields2, table2 = CSV.readTable(open(options.filename2, "r"))

    if options.unique:
        outfile = UniqueBuffer(sys.stdout)
    else:
        outfile = options.stdout

    nfields1 = []
    for x in range(len(fields1)):
        if fields1[x] in options.join_fields1:
            nfields1.append(x)
    nfields2 = []
    for x in range(len(fields2)):
        if fields2[x] in options.join_fields2:
            nfields2.append(x)

    # calculate row indices: double keys are not taken care of here
    keys = {}
    for row1 in table1:
        v = [row1[x] for x in nfields1]
        key = hashlib.md5("".join(v)).digest()
        keys[key] = row1

    if options.method == "intersection":
        # build new field list
        take = list(range(len(fields1)))
        c = len(take)
        for x in fields2:
            if x not in options.join_fields2:
                take.append(c)
            c += 1

        t = fields1 + fields2

        new_fields = [t[x] for x in take]

        print("\t".join(new_fields))

        for row2 in table2:
            v = [row2[x] for x in nfields2]
            key = hashlib.md5("".join(v)).digest()
            if key in keys:
                new_row = keys[key] + row2
                outfile.write(
                    "\t".join([new_row[x] for x in take]) + "\n")

    elif options.method == "rest":

        new_fields = fields2
        print("\t".join(new_fields))

        for row2 in table2:
            v = [row2[x] for x in nfields2]
            key = hashlib.md5("".join(v)).digest()
            if key not in keys:
                outfile.write("\t".join(row2) + "\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
