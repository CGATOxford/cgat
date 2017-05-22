'''
join_tables.py - join tables
======================================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python join_tables.py --help

Type::

   python join_tables.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import CGAT.Experiment as E

USAGE = """python %s < stdin > stdout

OPTIONS:

-t, --titles=           column titles
-m, --missing-value=          missing value
-h, --header-names=          add headers for files
-s, --method=sort --sort-order=             sort by column titles (given by sort order)
'#' at start of line is a comment
""" % sys.argv[0]


def GetData(f):

    while 1:
        line = f.readline()
        if not line:
            return None
        if line[0] == "#":
            continue
        break

    return line[:-1].split("\t")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: join_tables.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to do join on. Files have to sorted alphabeticlly and incrementally on these columns.")
    parser.add_option("-m", "--missing-value", dest="missing_value", type="string",
                      help="value for missing entries.")

    parser.set_defaults(
        headers=False,
        pattern_filename=None,
        title="",
        footer="",
        columns="1",
        missing_value="",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) < 2:
        raise ValueError("there have to be at least two tables")

    if options.columns:
        options.columns = [int(x) - 1 for x in options.columns.split(",")]

    # open all files
    files = []
    for filename in args:

        if os.path.exists(filename):
            files.append(open(filename, "r"))

    if len(files) <= 1:
        raise ValueError("less than two files opened")

    nfiles = len(files)

    # go through all files and iteratively find smallest entry
    data = []
    entries = []
    lengths = []
    takes = []

    for f in files:
        d = GetData(f)
        l = len(d)
        t = []
        for x in range(l):
            if x not in options.columns:
                t.append(x)

        takes.append(t)
        data.append([d[x] for x in t])
        entries.append([d[x] for x in options.columns])
        lengths.append(len(t))

    activa = list(range(nfiles))
    # take only first entry for duplicate entries
    # columns need to be sorted incrementally on all columns
    while len(activa) > 0:

        line = []

        min_field = min([entries[x] for x in activa])

        for f in range(nfiles):

            if files[f] and entries[f] == min_field:

                line += data[f]

                while entries[f] == min_field:

                    d = GetData(files[f])

                    if not d:
                        files[f].close()
                        files[f] = None
                        activa.remove(f)
                        break

                    data[f] = [d[x] for x in takes[f]]
                    entries[f] = [d[x] for x in options.columns]

            else:
                line += [options.missing_value for x in range(lengths[f])]

        options.stdout.write(
            "\t".join(min_field) + "\t" + "\t".join(line) + "\n")

    # close all files
    for f in files:
        if f:
            f.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
