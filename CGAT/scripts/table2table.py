'''table2table.py - operate on tables
==================================

:Tags: Python

Purpose
-------

This script implements a few methods for manipulating tables.

Methods working on all tables:
++++++++++++++++++++++++++++++

transpose
   transpose a table

split-fields
   Split muliple-value fields in each row at ``--separator``. Output
   multiple rows with all combinations.

group
   Group values by column

join-column
   Join rows in a table by columns

expand-table
   If a field in a row contains multiple values,
   the row is expanded into multiple rows such
   that all values have space.

flatten-table
   Output a table as row/column/value tuples.

as-column
   Output table as a single column. Colums in the original table are
   appended and output.

collapse-table
   Collapse a table of two columns with row names in the first
   column. Outputs a table with multiple columns for each row name.

Methods for numerical columns
+++++++++++++++++++++++++++++

Some methods make only sense for columns containing numerical values.
If a table contains both numerical and non-numerical data, the
numerical columns can be specified by the ``--columns`` option.

normalize-by-value
   divide all cells in a table by a value

multiply-by-value
   multiply all cells in a table by a value

lower-bound
   replace all cells with a value of less than lower bound with the lower
   bound.

upper-bound
   replace all cells with a value of more than upper bound with the upper
   bound.

normalize-by-table
   divide each cell in a table with the corresponding entry in a secondary
   table.

normalize-by-max
   divide table columns by maximum per column

kullback-leibler
   compute kullback-leibler divergence between two columns. Compute
   both D(a||b), D(b||a) and (D(a||b) + D(b||a)) / 2

rank
   substitute cells with their ranks in a column

fdr
   compute an FDR over selected columns. Replaces the columns
   with the qvalues.

Usage
-----

Example::

   python table2table.py --help

Type::

   python table2table.py --help

for command line help.

Command line options
--------------------

'''
import sys
import math
import itertools
import collections

import CGAT.Experiment as E
import CGAT.CSV as CSV
import CGAT.Stats as Stats
import CGAT.IOTools as IOTools
import scipy
from functools import reduce


def getColumns(fields, columns="all"):
    '''return columns to take.'''
    if columns == "all":
        return list(range(len(fields)))
    elif columns == "all-but-first":
        return list(range(1, len(fields)))
    else:
        map_field2column = dict([(y, x) for x, y in enumerate(fields)])
        c = []
        for x in columns.split(","):
            if x in map_field2column:
                c.append(map_field2column[x])
            else:
                c.append(int(x) - 1)
        return c


def readAndTransposeTable(infile, options):
    """read table from infile and transpose
    """
    rows = []
    if options.transpose_format == "default":

        for line in infile:
            if line[0] == "#":
                continue
            rows.append(line[:-1].split("\t"))

    elif options.transpose_format == "separated":
        for line in infile:
            if line[0] == "#":
                continue
            key, vals = line[:-1].split("\t")
            row = [key] + vals.split(options.separator)
            rows.append(row)

    ncols = max([len(x) for x in rows])
    nrows = len(rows)

    new_rows = [[""] * nrows for x in range(ncols)]

    for r in range(0, len(rows)):
        for c in range(0, len(rows[r])):
            new_rows[c][r] = rows[r][c]

    if options.set_transpose_field:
        new_rows[0][0] = options.set_transpose_field

    for row in new_rows:
        options.stdout.write("\t".join(row) + "\n")


def readAndGroupTable(infile, options):
    """read table from infile and group.
    """
    fields, table = CSV.readTable(
        infile, with_header=options.has_headers, as_rows=True)
    options.columns = getColumns(fields, options.columns)
    assert options.group_column not in options.columns

    converter = float
    new_fields = [fields[options.group_column]] + [fields[x]
                                                   for x in options.columns]

    if options.group_function == "min":
        f = min
    elif options.group_function == "max":
        f = max
    elif options.group_function == "sum":
        f = lambda z: reduce(lambda x, y: x + y, z)
    elif options.group_function == "mean":
        f = scipy.mean
    elif options.group_function == "cat":
        f = lambda x: ";".join([y for y in x if y != ""])
        converter = str
    elif options.group_function == "uniq":
        f = lambda x: ";".join([y for y in set(x) if y != ""])
        converter = str
    elif options.group_function == "stats":
        f = lambda x: str(Stats.DistributionalParameters(x))
        # update headers
        new_fields = [fields[options.group_column]]
        for c in options.columns:
            new_fields += list(["%s_%s" %
                                (fields[c], x) for x in Stats.DistributionalParameters().getHeaders()])

    # convert values to floats (except for group_column)
    # Delete rows with unconvertable values and not in options.columns
    new_table = []
    for row in table:
        skip = False
        new_row = [row[options.group_column]]

        for c in options.columns:
            if row[c] == options.missing_value:
                new_row.append(row[c])
            else:
                try:
                    new_row.append(converter(row[c]))
                except ValueError:
                    skip = True
                    break
        if not skip:
            new_table.append(new_row)
    table = new_table

    new_rows = CSV.groupTable(table,
                              group_column=0,
                              group_function=f)

    options.stdout.write("\t".join(new_fields) + "\n")
    for row in new_rows:
        options.stdout.write("\t".join(map(str, row)) + "\n")


def readAndExpandTable(infile, options):
    '''splits fields in table at separator.

    If a field in a row contains multiple values,
    the row is expanded into multiple rows such
    that all values have space.
    '''

    fields, table = CSV.readTable(
        infile, with_header=options.has_headers, as_rows=True)

    options.stdout.write("\t".join(fields) + "\n")

    for row in table:

        data = []
        for x in range(len(fields)):
            data.append(row[x].split(options.separator))

        nrows = max([len(d) for d in data])

        for d in data:
            d += [""] * (nrows - len(d))

        for n in range(nrows):
            options.stdout.write("\t".join([d[n] for d in data]) + "\n")


def readAndCollapseTable(infile, options, missing_value=""):
    '''collapse a table.

    Collapse a table of two columns with row names in the first
    column. Outputs a table with multiple columns for each row name.
    '''

    fields, table = CSV.readTable(
        infile, with_header=options.has_headers, as_rows=True)

    if len(fields) != 2:
        raise NotImplementedError("can only work on tables with two columns")

    values = collections.defaultdict(list)

    # column header after which to add
    separator = table[0][0]
    row_names = set([x[0] for x in table])

    row_name, value = table[0]

    values[row_name].append(value)
    added = set([row_name])
    for row_name, value in table[1:]:
        if row_name == separator:
            for r in row_names:
                if r not in added:
                    values[r].append(missing_value)
            added = set()

        values[row_name].append(value)
        added.add(row_name)

    for r in row_names:
        if r not in added:
            values[r].append(missing_value)

    sizes = set([len(x) for x in list(values.values())])
    assert len(sizes) == 1, "unequal number of row_names"
    size = list(sizes)[0]

    options.stdout.write(
        "row\t%s\n" % ("\t".join(["column_%i" % x for x in range(size)])))

    for key, row in list(values.items()):
        options.stdout.write("%s\t%s\n" % (key, "\t".join(row)))


def computeFDR(infile, options):
    '''compute FDR on a table.
    '''

    fields, table = CSV.readTable(
        infile, with_header=options.has_headers, as_rows=True)

    options.stdout.write("\t".join(fields) + "\n")

    for row in table:

        data = []
        for x in range(len(fields)):
            data.append(row[x].split(options.separator))

        nrows = max([len(d) for d in data])

        for d in data:
            d += [""] * (nrows - len(d))

        for n in range(nrows):
            options.stdout.write("\t".join([d[n] for d in data]) + "\n")


def readAndJoinTable(infile, options):

    fields, table = CSV.readTable(
        infile, with_header=options.has_headers, as_rows=True)

    join_column = options.join_column - 1
    join_name = options.join_column_name - 1

    join_rows = list(set([x[join_column] for x in table]))
    join_rows.sort()

    join_names = list(set([x[join_name] for x in table]))
    join_names.sort()

    join_columns = list(
        set(range(len(fields))).difference(set((join_column, join_name))))
    join_columns.sort()

    new_table = []
    map_old2new = {}

    map_name2start = {}
    x = 1
    for name in join_names:
        map_name2start[name] = x
        x += len(join_columns)

    row_width = len(join_columns) * len(join_names)
    for x in join_rows:
        map_old2new[x] = len(map_old2new)
        new_row = [x, ] + ["na"] * row_width
        new_table.append(new_row)

    for row in table:
        row_index = map_old2new[row[join_column]]
        start = map_name2start[row[join_name]]
        for x in join_columns:
            new_table[row_index][start] = row[x]
            start += 1

    # print new table
    options.stdout.write(fields[join_column])
    for name in join_names:
        for column in join_columns:
            options.stdout.write(
                "\t%s%s%s" % (name, options.separator, fields[column]))
    options.stdout.write("\n")

    for row in new_table:
        options.stdout.write("\t".join(row) + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-m", "--method", dest="methods", type="choice", action="append",
        choices=("transpose", "normalize-by-max", "normalize-by-value",
                 "multiply-by-value",
                 "percentile", "remove-header", "normalize-by-table",
                 "upper-bound", "lower-bound", "kullback-leibler",
                 "expand", "compress", "fdr", "grep"),
        help="""actions to perform on table.""")

    parser.add_option("-s", "--scale", dest="scale", type="float",
                      help="factor to scale matrix by.")

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="output number format [default]")

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="Parameters for various functions.")

    parser.add_option(
        "-t", "--header-names", dest="has_headers", action="store_true",
        help="matrix has row/column headers.")

    parser.add_option("--transpose", dest="transpose", action="store_true",
                      help="transpose table.")

    parser.add_option(
        "--set-transpose-field", dest="set_transpose_field", type="string",
        help="set first field (row 1 and col 1) to this value [%default].")

    parser.add_option(
        "--transpose-format", dest="transpose_format", type="choice",
        choices=("default", "separated", ),
        help="input format of un-transposed table")

    parser.add_option(
        "--expand", dest="expand_table", action="store_true",
        help="expand table - multi-value cells with be expanded over "
        "several rows.")

    parser.add_option("--no-headers", dest="has_headers", action="store_false",
                      help="matrix has no row/column headers.")

    parser.add_option("--columns", dest="columns", type="string",
                      help="columns to use.")

    parser.add_option("--file", dest="file", type="string",
                      help="columns to test from table.",
                      metavar="FILE")

    parser.add_option("-d", "--delimiter", dest="delimiter", type="string",
                      help="delimiter of columns.",
                      metavar="DELIM")

    parser.add_option(
        "-V", "--invert-match", dest="invert_match",
        action="store_true",
        help="invert match.")

    parser.add_option("--sort-by-rows", dest="sort_rows", type="string",
                      help="output order for rows.")

    parser.add_option("-a", "--value", dest="value", type="float",
                      help="value to use for various algorithms.")

    parser.add_option(
        "--group", dest="group_column", type="int",
        help="group values by column. Supply an integer column "
        "[default=%default]")

    parser.add_option("--group-function", dest="group_function", type="choice",
                      choices=(
                          "min", "max", "sum", "mean", "stats", "cat", "uniq"),
                      help="function to group values by.")

    parser.add_option("--join-table", dest="join_column", type="int",
                      help="join rows in a table by columns.")

    parser.add_option(
        "--collapse-table", dest="collapse_table", type="string",
        help="collapse a table. Value determines the missing variable "
        "[%default].")

    parser.add_option(
        "--join-column-name", dest="join_column_name", type="int",
        help="use this column as a prefix.")

    parser.add_option(
        "--flatten-table", dest="flatten_table", action="store_true",
        help="flatten a table [%default].")

    parser.add_option("--as-column", dest="as_column", action="store_true",
                      help="output table as a single column.")

    parser.add_option(
        "--split-fields", dest="split_fields", action="store_true",
        help="split fields.")

    parser.add_option(
        "--separator", dest="separator", type="string",
        help="separator for multi-valued fields [default=%default].")

    parser.add_option(
        "--fdr-method", dest="fdr_method", type="choice",
        choices=(
            "BH", "bonferroni", "holm", "hommel", "hochberg", "BY"),
        help="method to perform multiple testing correction by controlling "
        "the fdr [default=%default].")

    parser.add_option(
        "--fdr-add-column", dest="fdr_add_column", type="string",
        help="add new column instead of replacing existing columns. "
        "The value of the option will be used as prefix if there are "
        "multiple columns [%default]")

    # IMS: add option to use a column as the row id in flatten
    parser.add_option(
        "--id-column", dest="id_column", type="string",
        help="list of column(s) to use as the row id when flattening "
        "the table. If None, then row number is used. [default=%default].")

    parser.add_option(
        "--variable-name", dest="variable_name", type="string",
        help="the column header for the 'variable' column when flattening "
        "[default=%default].")

    parser.add_option(
        "--value-name", dest="value_name", type="string",
        help="the column header for the 'value' column when flattening "
        "[default=%default].")

    parser.set_defaults(
        methods=[],
        scale=1.0,
        has_headers=True,
        format=None,
        value=0.0,
        parameters="",
        columns="all",
        transpose=False,
        set_transpose_field=None,
        transpose_format="default",
        group=False,
        group_column=0,
        group_function="mean",
        missing_value="na",
        sort_rows=None,
        flatten_table=False,
        collapse_table=None,
        separator=";",
        expand=False,
        join_column=None,
        join_column_name=None,
        compute_fdr=None,
        as_column=False,
        fdr_method="BH",
        fdr_add_column=None,
        id_column=None,
        variable_name="column",
        value_name="value",
        file=None,
        delimiter="\t",
        invert_match=False,
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    options.parameters = options.parameters.split(",")

    if options.group_column:
        options.group = True
        options.group_column -= 1

    ######################################################################
    ######################################################################
    ######################################################################
    # if only to remove header, do this quickly
    if options.methods == ["remove-header"]:

        first = True
        for line in options.stdin:
            if line[0] == "#":
                continue
            if first:
                first = False
                continue
            options.stdout.write(line)

    elif options.transpose or "transpose" in options.methods:

        readAndTransposeTable(options.stdin, options)

    elif options.flatten_table:
        # IMS: bug fixed to make work. Also added options for keying
        # on a particular and adding custom column headings

        fields, table = CSV.readTable(
            options.stdin, with_header=options.has_headers, as_rows=True)

        options.columns = getColumns(fields, options.columns)

        if options.id_column:
            id_columns = [int(x) - 1 for x in options.id_column.split(",")]
            id_header = "\t".join([fields[id_column]
                                   for id_column in id_columns])
            options.columns = [
                x for x in options.columns if x not in id_columns]
        else:
            id_header = "row"

        options.stdout.write(
            "%s\t%s\t%s\n" % (id_header, options.variable_name,
                              options.value_name))

        for x, row in enumerate(table):

            if options.id_column:
                row_id = "\t".join([row[int(x) - 1]
                                    for x in options.id_column.split(",")])
            else:
                row_id = str(x)

            for y in options.columns:
                options.stdout.write(
                    "%s\t%s\t%s\n" % (row_id, fields[y], row[y]))

    elif options.as_column:

        fields, table = CSV.readTable(
            options.stdin, with_header=options.has_headers, as_rows=True)
        options.columns = getColumns(fields, options.columns)
        table = list(zip(*table))

        options.stdout.write("value\n")

        for column in options.columns:
            options.stdout.write("\n".join(table[column]) + "\n")

    elif options.split_fields:

        # split comma separated fields
        fields, table = CSV.readTable(options.stdin,
                                      with_header=options.has_headers,
                                      as_rows=True)

        options.stdout.write("%s\n" % ("\t".join(fields)))

        for row in table:
            row = [x.split(options.separator) for x in row]
            for d in itertools.product(*row):
                options.stdout.write("%s\n" % "\t".join(d))

    elif options.group:
        readAndGroupTable(options.stdin, options)

    elif options.join_column:
        readAndJoinTable(options.stdin, options)

    elif options.expand_table:
        readAndExpandTable(options.stdin, options)

    elif options.collapse_table is not None:
        readAndCollapseTable(options.stdin, options, options.collapse_table)

    elif "grep" in options.methods:

        options.columns = [int(x) - 1 for x in options.columns.split(",")]

        patterns = []

        if options.file:
            infile = IOTools.openFile(options.file, "r")
            for line in infile:
                if line[0] == "#":
                    continue
                patterns.append(line[:-1].split(options.delimiter)[0])
        else:
            patterns = args

        for line in options.stdin:

            data = line[:-1].split(options.delimiter)
            found = False

            for c in options.columns:

                if data[c] in patterns:
                    found = True
                    break

            if (not found and options.invert_match) or (found and not options.invert_match):
                print(line[:-1])
    else:

        ######################################################################
        ######################################################################
        ######################################################################
        # Apply remainder of transformations
        fields, table = CSV.readTable(
            options.stdin, with_header=options.has_headers, as_rows=False)
        # convert columns to list
        table = [list(x) for x in table]

        ncols = len(fields)
        if len(table) == 0:
            raise ValueError("table is empty")

        nrows = len(table[0])

        E.info("processing table with %i rows and %i columns" % (nrows, ncols))

        options.columns = getColumns(fields, options.columns)

        # convert all values to float
        for c in options.columns:
            for r in range(nrows):
                try:
                    table[c][r] = float(table[c][r])
                except ValueError:
                    continue

        for method in options.methods:

            if method == "normalize-by-value":

                value = float(options.parameters[0])
                del options.parameters[0]

                for c in options.columns:
                    table[c] = [x / value for x in table[c]]

            elif method == "multiply-by-value":

                value = float(options.parameters[0])
                del options.parameters[0]

                for c in options.columns:
                    table[c] = [x * value for x in table[c]]

            elif method == "normalize-by-max":

                for c in options.columns:
                    m = max(table[c])
                    table[c] = [x / m for x in table[c]]

            elif method == "kullback-leibler":
                options.stdout.write("category1\tcategory2\tkl1\tkl2\tmean\n")
                format = options.format
                if format is None:
                    format = "%f"

                for x in range(0, len(options.columns) - 1):
                    for y in range(x + 1, len(options.columns)):
                        c1 = options.columns[x]
                        c2 = options.columns[y]
                        e1 = 0
                        e2 = 0
                        for z in range(nrows):
                            p = table[c1][z]
                            q = table[c2][z]
                            e1 += p * math.log(p / q)
                            e2 += q * math.log(q / p)

                        options.stdout.write("%s\t%s\t%s\t%s\t%s\n" % (
                            fields[c1], fields[c2],
                            format % e1,
                            format % e2,
                            format % ((e1 + e2) / 2)))
                E.Stop()
                sys.exit(0)

            elif method == "rank":

                for c in options.columns:
                    tt = table[c]
                    t = list(zip(tt, list(range(nrows))))
                    t.sort()
                    for i, n in zip([x[1] for x in t], list(range(nrows))):
                        tt[i] = n

            elif method in ("lower-bound", "upper-bound"):

                boundary = float(options.parameters[0])
                del options.parameters[0]
                new_value = float(options.parameters[0])
                del options.parameters[0]

                if method == "upper-bound":
                    for c in options.columns:
                        for r in range(nrows):
                            if isinstance(table[c][r], float) and \
                                    table[c][r] > boundary:
                                table[c][r] = new_value
                else:
                    for c in options.columns:
                        for r in range(nrows):
                            if isinstance(table[c][r], float) and \
                                    table[c][r] < boundary:
                                table[c][r] = new_value

            elif method == "fdr":
                pvalues = []
                for c in options.columns:
                    pvalues.extend(table[c])

                assert max(pvalues) <= 1.0, "pvalues > 1 in table: max=%s" % \
                    str(max(pvalues))
                assert min(pvalues) >= 0, "pvalue < 0 in table: min=%s" % \
                    str(min(pvalues))

                # convert to str to avoid test for float downstream
                qvalues = list(map(
                    str, Stats.adjustPValues(pvalues,
                                             method=options.fdr_method)))

                if options.fdr_add_column is None:
                    x = 0
                    for c in options.columns:
                        table[c] = qvalues[x:x + nrows]
                        x += nrows
                else:
                    # add new column headers
                    if len(options.columns) == 1:
                        fields.append(options.fdr_add_column)
                    else:
                        for co in options.columns:
                            fields.append(options.fdr_add_column + fields[c])

                    x = 0
                    for c in options.columns:
                        # add a new column
                        table.append(qvalues[x:x + nrows])
                        x += nrows
                    ncols += len(options.columns)

            elif method == "normalize-by-table":

                other_table_name = options.parameters[0]
                del options.parameters[0]
                other_fields, other_table = CSV.readTable(
                    IOTools.openFile(other_table_name, "r"),
                    with_header=options.has_headers,
                    as_rows=False)

                # convert all values to float
                for c in options.columns:
                    for r in range(nrows):
                        try:
                            other_table[c][r] = float(other_table[c][r])
                        except ValueError:
                            continue

                # set 0s to 1 in the other matrix
                for c in options.columns:
                    for r in range(nrows):
                        if isinstance(table[c][r], float) and \
                                isinstance(other_table[c][r], float) and \
                                other_table[c][r] != 0:
                            table[c][r] /= other_table[c][r]
                        else:
                            table[c][r] = options.missing_value

        # convert back
        if options.format is not None:
            for c in options.columns:
                for r in range(nrows):
                    if isinstance(table[c][r], float):
                        table[c][r] = format % table[c][r]

        options.stdout.write("\t".join(fields) + "\n")
        if options.sort_rows:
            old2new = {}
            for r in range(nrows):
                old2new[table[0][r]] = r
            for x in options.sort_rows.split(","):
                if x not in old2new:
                    continue
                r = old2new[x]
                options.stdout.write(
                    "\t".join(map(str,
                                  [table[c][r] for c in range(ncols)])) + "\n")
        else:
            for r in range(nrows):
                options.stdout.write(
                    "\t".join(map(str,
                                  [table[c][r] for c in range(ncols)])) + "\n")

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
