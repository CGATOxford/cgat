'''
csvs2csv.py - join tables
=========================

:Tags: Python

Purpose
-------

This script reads several tab-separated tables and joins them.

.. note:: 
   working with multiple columns per table and sorting is
   not implemented correctly and likely to fail.

Usage
-----

Example::

   python combine_tables.py --help

Type::

   python combine_tables.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import os
import glob

import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-t", "--no-titles", dest="titles", action="store_false",
        help="no titles in input.")

    parser.add_option(
        "-i", "--skip-titles", dest="skip_titles", action="store_true",
        help="skip output of titles.")

    parser.add_option(
        "-m", "--missing-value", dest="missing_value", type="string",
        help="entry to use for missing values.")

    parser.add_option("--header-names", dest="headers", type="string",
                      help="add headers for files.")

    parser.add_option(
        "-c", "--columns", dest="columns", type="string",
        help="columns to use for joining. Multiple columns can be specified "
        "as a comma-separated list [default=%default].")

    parser.add_option(
        "-g", "--glob", dest="glob", type="string",
        help="wildcard expression for table names.")

    parser.add_option(
        "-s", "--sort-order", dest="sort", type="string",
        help="sort by column titles alphabetical|numeric|list of columns.")

    parser.add_option(
        "-e", "--merge-overlapping", dest="merge", action="store_true",
        help="simply merge tables without matching up rows. "
        "[default=%default].")

    parser.add_option(
        "--sort-keys", dest="sort_keys", type="choice",
        choices=("numeric", "alphabetic"),
        help="sort key columns by value.")

    parser.add_option(
        "--keep-empty", dest="ignore_empty", action="store_false",
        help="keep empty tables. The default is to ignore them.")

    parser.add_option(
        "--add-file-prefix", dest="add_file_prefix", action="store_true",
        help="add file prefix to columns headers in multi-column tables "
        "[default=%default]")

    parser.add_option(
        "--regex-filename", dest="regex_filename", type="string",
        help="pattern to apply to filename to build prefix "
        "[default=%default]")

    parser.set_defaults(
        titles=True,
        skip_titles=False,
        missing_value="na",
        headers=None,
        sort=None,
        glob=None,
        columns="1",
        sort_keys=False,
        merge=False,
        ignore_empty=True,
        add_file_prefix=False,
        regex_filename="(.*)"
    )

    (options, args) = E.Start(parser)

    if options.headers:
        if "," in options.headers:
            options.headers = options.headers.split(",")
        else:
            options.headers = re.split("\s+", options.headers.strip())

    if options.sort and options.sort not in ("numeric", "alphabetic"):
        if "," in options.sort:
            options.sort = options.sort.split(",")
        else:
            options.sort = re.split("\s+", options.sort)

    if options.merge:
        options.columns = []
    else:
        options.columns = [int(x) - 1 for x in options.columns.split(",")]

    options.filenames = []

    if options.glob:
        options.filenames += glob.glob(options.glob)

    options.filenames += args

    if len(options.filenames) < 1:
        print(USAGE, "no tables specified/found.")
        sys.exit(1)

    if options.loglevel >= 1:
        options.stdlog.write("# combining %i tables.\n" %
                             len(options.filenames))
        sys.stdout.flush()
        if len(options.filenames) == 1:
            for line in IOTools.openFile(options.filenames[0]):
                options.stdout.write(line)
            E.Stop()
            sys.exit(0)

    if options.headers and options.headers[0] != "auto" and \
       len(options.headers) != len(options.filenames):
        raise "number of provided headers (%i) is not equal to number filenames (%i)." %\
              (len(options.headers), len(options.filenames))

    tables = []
    keys = {}
    sorted_keys = []
    sizes = {}
    if options.merge:
        titles = ["count"]
    else:
        titles = []

    for filename in options.filenames:

        prefix = os.path.basename(filename)

        if os.path.exists(filename):
            file = IOTools.openFile(filename, "r")
            lines = [x for x in file if x[0] != "#"]

        else:
            lines = []

        if len(lines) == 0 and options.ignore_empty:
            continue

        table = {}
        sizes = {}
        max_size = 0
        ncolumns = 0

        if options.titles:
            data = lines[0][:-1].split("\t")
            if not titles:
                key = "-".join([data[x] for x in options.columns])
                titles = [key]
            for x in range(len(data)):
                if x in options.columns:
                    continue
                ncolumns += 1
                if options.add_file_prefix:
                    p = re.search(options.regex_filename, prefix).groups()[0]
                    titles.append("%s_%s" % (p, data[x]))
                else:
                    titles.append(data[x])

            del lines[0]
        else:
            ncolumns = 1

        n = 0
        for line in lines:
            data = line[:-1].split("\t")
            row_keys = [data[x] for x in options.columns]
            if options.sort_keys:
                if options.sort_keys == "numeric":
                    row_keys.sort(lambda x, y: cmp(float(x), float(y)))
                else:
                    row_keys.sort()
            if options.merge:
                key = n
            else:
                key = "-".join(row_keys)

            if key not in keys:
                sorted_keys.append(key)
                keys[key] = 1
                sizes[key] = 0

            max_size = max(len(data) - len(options.columns), max_size)
            table[key] = [data[x]
                          for x in [x for x in range(0, len(data)) if x not in options.columns]]
            n += 1

        # enter columns of "na" for empty tables.
        if max_size == 0:
            max_size = ncolumns

        tables.append((max_size, table))

    if len(tables) == len(titles) - 1:

        if options.headers:
            headers = ["bin"]
            if options.headers[0] == 'auto':
                for t in range(len(tables)):
                    headers.append(os.path.basename(options.filenames[t]))
                    headers += [""] * (tables[t][0] - 1)

            else:
                for t in range(len(tables)):
                    headers.append(options.headers[t])
                    headers += [""] * (tables[t][0] - 1)

            # use headers as titles, if headers is given and skip-titles is
            # turned on
            if options.titles and options.skip_titles:
                titles = headers
            else:
                # otherwise: print the headers out right away
                sys.stdout.write("\t".join(headers) + "\n")

        order = list(range(0, len(tables) + 1))

        if options.titles:

            if options.sort:
                sort_order = []

                if options.sort == "numeric":
                    t = list(zip(list(map(int, titles[1:])), list(range(1, len(titles) + 1))))
                    t.sort()

                    for tt in t:
                        sort_order.append(titles[tt[1]])

                elif options.sort == "alphabetical":
                    t = list(zip(titles[1:], list(range(1, len(titles) + 1))))
                    t.sort()

                    for tt in t:
                        sort_order.append(titles[tt[1]])
                else:
                    sort_order = options.sort

                map_title2pos = {}
                for x in range(1, len(titles)):
                    map_title2pos[titles[x]] = x

                order = [0, ]
                for x in sort_order:
                    if x in map_title2pos:
                        order.append(map_title2pos[x])

            else:
                order = list(range(0, len(titles)))

            sys.stdout.write(
                "\t".join([titles[order[x]] for x in range(len(titles))]))
            sys.stdout.write("\n")

        if options.sort_keys:
            if options.sort_keys:
                if options.sort_keys == "numeric":
                    sorted_keys.sort(lambda x, y: cmp(float(x), float(y)))
                else:
                    sorted_keys.sort()

        for key in sorted_keys:

            sys.stdout.write("%s" % key)

            for x in order[1:]:
                max_size, table = tables[x - 1]
                c = 0
                if key in table:
                    sys.stdout.write("\t")
                    sys.stdout.write("\t".join(table[key]))
                    c = len(table[key])

                assert(max_size == 1)

                sys.stdout.write(
                    "\t%s" % options.missing_value * (max_size - c))

            sys.stdout.write("\n")

    else:

        # for multi-column table, just write
        if options.titles:
            sys.stdout.write(
                "\t".join([titles[x] for x in range(len(titles))]))
            sys.stdout.write("\n")

        for key in sorted_keys:

            sys.stdout.write("%s" % key)

            for x in range(len(tables)):

                max_size, table = tables[x]
                c = 0
                if key in table:
                    sys.stdout.write("\t")
                    sys.stdout.write("\t".join(table[key]))
                    c = len(table[key])

                sys.stdout.write(
                    "\t%s" % options.missing_value * (max_size - c))

            sys.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
