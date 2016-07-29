'''combine_tables.py - join tables
===============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads several tab-separated tables and joins them into a
single one.

.. todo::

   * Rename to tables2table.py
   * Use pandas dataframes for fast IO and merging/joining

Usage
-----

The option ``--header-names`` sets the column titles explicitely. Add
``--skip-titles`` if you want to avoid echoing the original title in
the input files.


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
import string
import os
import glob
import collections

import CGAT.IOTools as IOTools
import CGAT.Experiment as E


def readTable(filename, options):
    '''read table and filter.
    '''

    if os.path.exists(filename):
        lines = IOTools.openFile(filename, "r").readlines()
    else:
        lines = []

    # extract table by regular expression
    if options.regex_start:
        rx = re.compile(options.regex_start)
        for n, line in enumerate(lines):
            if rx.search(line):
                E.info("reading table from line %i/%i" % (n, len(lines)))
                lines = lines[n:]
                break
        else:
            E.info("start regex not found - no table")
            lines = []

    if options.regex_end:
        rx = re.compile(options.regex_end)
        for n, line in enumerate(lines):
            if rx.search(line):
                break
        lines = lines[:n]

    # remove comments and empty lines
    lines = [x for x in lines if not x.startswith("#") and x.strip()]

    return lines


def concatenateTables(outfile, options, args):
    '''concatenate tables.'''

    missing_value = options.missing_value

    rx = re.compile(options.regex_filename)

    if options.headers is None or options.headers == "auto":
        row_headers = [
            [y for y in rx.search(x).groups()] for x in options.filenames]
    else:
        row_headers = [options.headers]

    tables = []
    # read all tables
    for filename in options.filenames:
        table = readTable(filename, options)
        if len(table) == 0:
            E.warn("table '%s' is empty" % filename)
            continue
        tables.append(table)

    if options.cat is None:
        if len(row_headers) == 1:
            row_head_titles = ["filename"]
        else:
            row_head_titles = ["pattern" + str(x)
                               for x in range(len(row_headers))]
    else:
        row_head_titles = [x.strip() for x in options.cat.split(",")]
        if len(row_headers[0]) != len(row_head_titles):
            raise ValueError(
                "row header (%i) has different number of fields in "
                "regular expression than supplied by the --cat option (%i)" %
                (len(row_headers[0]), len(row_head_titles)))

    # collect titles
    if options.input_has_titles:
        titles = collections.OrderedDict()
        for table in tables:
            for key in table[0][:-1].split("\t"):
                # skip any titles that conflict with
                # the newly added titles
                if key in row_head_titles:
                    continue
                titles[key] = 1

        outfile.write("%s\t%s\n" %
                      ("\t".join([x for x in row_head_titles]),
                       "\t".join(titles.keys())))

        map_title2column = collections.defaultdict(lambda: None)
        for x, title in enumerate(titles.keys()):
            map_title2column[title] = x
    else:
        ncolumns = [len(table[0].split('\t')) for table in tables]
        if min(ncolumns) != max(ncolumns):
            raise ValueError('tables have unequal number of columns '
                             '(min=%i, max=%i)' %
                             (min(ncolumns),  max(ncolumns)))
        # create a pseudo dictionary of columns
        titles = collections.OrderedDict(
            [(x, x) for x in range(min(ncolumns))])

    all_titles = set(titles.keys())
    for nindex, table in enumerate(tables):
        if options.input_has_titles:
            titles = table[0][:-1].split("\t")
            map_old2new = [map_title2column[t] for t in titles]
            del table[0]
        else:
            map_old2new = list(range(len(all_titles)))

        for l in table:
            data = [missing_value] * len(all_titles)
            for x, value in enumerate(l[:-1].split("\t")):
                if map_old2new[x] is None:
                    continue

                data[map_old2new[x]] = value

            row = "\t".join([str(x) for x in row_headers[nindex]] +
                            data) + "\n"
            outfile.write(row)


def joinTables(outfile, options, args):
    '''join tables.'''

    if options.headers and options.headers[0] != "auto" and \
            len(options.headers) != len(options.filenames):
        raise ValueError("number of provided headers (%i) "
                         "is not equal to number filenames (%i)." %
                         (len(options.headers), len(options.filenames)))

    tables = []
    keys = {}
    sorted_keys = []
    sizes = {}

    if options.merge:
        titles = ["count"]
    else:
        titles = []

    headers_to_delete = []

    if options.prefixes:
        prefixes = [x.strip() for x in options.prefixes.split(",")]
        if len(prefixes) != len(options.filenames):
            raise ValueError(("number of prefixes (%i) and tables (%i) "
                              "do not match") %
                             (len(prefixes),
                              len(options.filenames)))
    else:
        prefixes = None

    E.debug("joining on columns %s and taking columns %s" %
            (options.columns, options.take))

    for nindex, filename in enumerate(options.filenames):

        E.info("processing %s (%i/%i)" %
               (filename, nindex + 1, len(options.filenames)))

        prefix = os.path.basename(filename)

        lines = readTable(filename, options)

        # skip (or not skip) empty tables
        if len(lines) == 0 and options.ignore_empty:
            E.warn("%s is empty - skipped" % filename)
            headers_to_delete.append(nindex)
            continue

        table = {}
        sizes = {}
        max_size = 0
        ncolumns = 0

        if options.input_has_titles:
            data = string.split(lines[0][:-1], "\t")
            # no titles have been defined so far
            if not titles:
                key = "-".join([data[x] for x in options.columns])
                titles = [key]

            # set take based on column titles or numerically
            if options.take:
                take = []
                # convert numeric columns for filtering
                for x in options.take:
                    try:
                        take.append(int(x) - 1)
                    except ValueError:
                        # will raise error if x is not present
                        take.append(data.index(x))
            else:
                # tables with max 100 columns
                take = None

            for x in range(len(data)):
                if x in options.columns or (take and x not in take):
                    continue
                ncolumns += 1
                if options.add_file_prefix:
                    try:
                        p = re.search(
                            options.regex_filename, prefix).groups()[0]
                    except AttributeError:
                        E.warn("can't extract title from filename %s" % prefix)
                        p = "unknown"
                    titles.append("%s_%s" % (p, data[x]))
                elif options.use_file_prefix:
                    try:
                        p = re.search(
                            options.regex_filename, prefix).groups()[0]
                    except:
                        E.warn("can't extract title from filename %s" % prefix)
                        p = "unknown"
                    titles.append("%s" % p)
                elif prefixes:
                    titles.append("%s_%s" % (prefixes[nindex], data[x]))
                else:
                    titles.append(data[x])

            del lines[0]
        else:

            # set take based on numeric columns if no titles are present
            if options.take:
                take = []
                # convert numeric columns for filtering
                for x in options.take:
                    take.append(int(x) - 1)
            else:
                # tables with max 100 columns
                take = None

            # IMS: We might still want filename titles even if the input
            # columns don't have titles.
            if options.add_file_prefix:
                if not titles:
                    titles = ["ID"]
                try:
                    p = re.search(options.regex_filename, prefix).groups()[0]
                except AttributeError:
                    E.warn("can't extract title from filename %s" % prefix)
                    p = "unknown"
                titles.append("%s_%s" % (p, data[x]))
            elif options.use_file_prefix:
                if not titles:
                    titles = ["ID"]
                try:
                    p = re.search(options.regex_filename, prefix).groups()[0]
                except:
                    E.warn("can't extract title from filename %s" % prefix)
                    p = "unknown"
                titles.append("%s" % p)
            ncolumns = 1

        n = 0
        for line in lines:
            data = string.split(line[:-1], "\t")
            try:
                row_keys = [data[x] for x in options.columns]
            except IndexError, msg:
                raise IndexError(
                    "error while parsing %s: %s" % (filename, msg))
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

            if take:
                max_size = len(take)
                table[key] = [data[x] for x in take]
            else:
                max_size = max(len(data) - len(options.columns), max_size)
                table[key] = [data[x]
                              for x in range(0, len(data))
                              if x not in options.columns]
            n += 1

        # enter columns of "na" for empty tables.
        if max_size == 0:
            max_size = ncolumns

        tables.append((max_size, table))

    # delete in reverse order
    if options.headers:
        for nindex in headers_to_delete[::-1]:
            del options.headers[nindex]

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
            if options.input_has_titles and options.skip_titles:
                titles = headers
            else:
                # otherwise: print the headers out right away
                outfile.write(string.join(headers, "\t") + "\n")

        order = range(0, len(tables) + 1)

        if options.input_has_titles or \
           (options.use_file_prefix or options.add_file_prefix):

            if options.sort:
                sort_order = []

                if options.sort == "numeric":
                    t = zip(map(int, titles[1:]), range(1, len(titles) + 1))
                    t.sort()

                    for tt in t:
                        sort_order.append(titles[tt[1]])

                elif options.sort == "alphabetical":
                    t = zip(titles[1:], range(1, len(titles) + 1))
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
                order = range(0, len(titles))

            outfile.write(
                "\t".join(map(lambda x: titles[order[x]], range(len(titles)))))
            outfile.write("\n")

        if options.sort_keys:
            if options.sort_keys:
                if options.sort_keys == "numeric":
                    sorted_keys.sort(lambda x, y: cmp(float(x), float(y)))
                else:
                    sorted_keys.sort()

        for key in sorted_keys:

            outfile.write("%s" % key)

            for x in order[1:]:

                max_size, table = tables[x - 1]
                c = 0
                if key in table:
                    outfile.write("\t")
                    outfile.write(string.join(table[key], "\t"))
                    c = len(table[key])

                assert(max_size == 1)

                outfile.write("\t%s" % options.missing_value * (max_size - c))

            outfile.write("\n")

    else:

        # for multi-column table, just write
        if options.input_has_titles:
            outfile.write(
                "\t".join(map(lambda x: titles[x], range(len(titles)))))
            outfile.write("\n")

        for key in sorted_keys:

            outfile.write("%s" % key)

            for x in range(len(tables)):

                max_size, table = tables[x]
                c = 0
                if key in table:
                    outfile.write("\t")
                    outfile.write("\t".join(table[key]))
                    c = len(table[key])

                outfile.write("\t%s" % options.missing_value * (max_size - c))

            outfile.write("\n")

# ------------------------------------------------------------------------


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--no-titles",
                      dest="input_has_titles",
                      action="store_false",
                      help="no titles in input [%default].")

    parser.add_option("--ignore-titles",
                      dest="ignore_titles",
                      action="store_true",
                      help="ignore titles in input [%default]")

    parser.add_option("-i", "--skip-titles",
                      dest="skip_titles",
                      action="store_true",
                      help="skip output of titles.")

    parser.add_option("-m", "--missing-value",
                      dest="missing_value",
                      type="string",
                      help="entry to use for missing values.")

    parser.add_option("--header-names", dest="headers", type="string",
                      help="add headers for files as a ,-separated "
                      "list [%default].")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for joining. Multiple columns "
                      "can be specified as a comma-separated list "
                      "[default=%default].")

    parser.add_option("-k", "--take",
                      dest="take",
                      type="string",
                      action="append",
                      help="columns to take. If not set, all columns "
                      "except for "
                      "the join columns are taken [%default]")

    parser.add_option("-g", "--glob", dest="glob", type="string",
                      help="wildcard expression for table names.")

    parser.add_option(
        "-s", "--sort-order", dest="sort", type="string",
        help="sort by column titles in particular given order: "
        "alphabetical|numeric|list of columns.")

    parser.add_option(
        "-e", "--merge-overlapping", dest="merge", action="store_true",
        help="simply merge tables without matching up "
        "rows. [default=%default].")

    parser.add_option("-a", "--cat", dest="cat", type="string",
                      help="simply concatenate tables. Adds an "
                      "additional column called X with the filename "
                      " [default=%default].")

    parser.add_option("--sort-keys", dest="sort_keys", type="choice",
                      choices=("numeric", "alphabetic"),
                      help="sort key columns by value.")

    parser.add_option("--keep-empty", dest="ignore_empty",
                      action="store_false",
                      help="keep empty tables. The default is "
                      "to ignore them.")

    parser.add_option("--ignore-empty",
                      dest="ignore_empty",
                      action="store_true",
                      help="ignore empty tables - this is "
                      "the default [%default].")

    parser.add_option("--add-file-prefix",
                      dest="add_file_prefix",
                      action="store_true",
                      help="add file prefix to "
                      "columns headers. Suitable for multi-column"
                      "tables [default=%default]")

    parser.add_option("--use-file-prefix",
                      dest="use_file_prefix",
                      action="store_true",
                      help="use file prefix as column headers. "
                      "Suitable for two-column tables "
                      "[default=%default]")

    parser.add_option("--prefixes", dest="prefixes", type="string",
                      help="list of prefixes to use. "
                      ", separated list of prefixes. "
                      "The number of prefixes need to correspond to the "
                      "number of input files [default=%default]")

    parser.add_option("--regex-filename", dest="regex_filename",
                      type="string",
                      help="pattern to apply to filename to "
                      "build prefix [default=%default]")

    parser.add_option("--regex-start",
                      dest="regex_start",
                      type="string",
                      help="regular expression to start "
                      "collecting table in a file [default=%default]")

    parser.add_option("--regex-end",
                      dest="regex_end",
                      type="string",
                      help="regular expression to end collecting "
                      "table in a file [default=%default]")

    parser.add_option("--test", dest="test",
                      type="int",
                      help="test combining tables with "
                      "first X rows [default=%default]")

    parser.set_defaults(
        input_has_titles=True,
        skip_titles=False,
        missing_value="na",
        headers=None,
        sort=None,
        glob=None,
        columns="1",
        sort_keys=False,
        merge=False,
        ignore_empty=True,
        regex_start=None,
        regex_end=None,
        add_file_prefix=False,
        use_file_prefix=False,
        cat=None,
        take=[],
        regex_filename="(.*)",
        prefixes=None,
        test=0,
    )

    (options, args) = E.Start(parser, argv=argv)

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
        options.columns = map(lambda x: int(x) - 1, options.columns.split(","))

    options.filenames = []

    if options.glob:
        options.filenames += glob.glob(options.glob)

    options.filenames += args

    if len(options.filenames) < 1:
        raise ValueError("no tables found.")

    E.info("combining %i tables" % len(options.filenames))

    if options.cat:
        concatenateTables(options.stdout, options, args)
    else:
        joinTables(options.stdout, options, args)

    E.Stop()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
