"""CSV.py - Tools for parsing CSV files
========================================

The methods in this module provide utility functions for
working with :term:`CSV` or :term:`TSV` formatted files.

With pandas providing fast and flexible access to :term:`CSV`
formatted files, most of the functionaly here is now superfluous.

:class:`DictReader` is derived from :py:class:`csv.DictReader`
and adds the capability to skip comment characters.

"""

import six
import csv
import sys

PY3 = sys.version_info > (3,)


def getMapColumn2Type(rows, ignore_empty=False, get_max_values=False):
    """map fields to types based on rows.

    Preference is Int to Float to String.

    If get_max_values is set to true, the maximum values for integer
    columns are returned in a dictionary.
    """

    headers = list(rows[0].keys())
    map_column2type = {}
    is_full = {}

    max_values = {}

    for row in rows:
        for h in headers:

            if row[h] == "":
                continue

            is_full[h] = True

            if isinstance(row[h], int):
                t = int
                if h not in max_values:
                    max_values[h] = int(row[h])
                else:
                    max_values[h] = max(max_values[h], int(row[h]))

            elif isinstance(row[h], float):
                t = float
            else:
                continue

            map_column2type[h] = t

    ignored = []
    for h in headers:
        if h not in map_column2type:
            if h in is_full or not ignore_empty:
                map_column2type[h] = bytes
            else:
                ignored.append(h)
    if get_max_values:
        return map_column2type, ignored, max_values
    else:
        return map_column2type, ignored


class CommentStripper(six.Iterator):
    """Iterator for stripping comments from file.

    This iterator will skip any lines beginning with ``#``
    or any empty lines at the beginning of the output.
    """

    def __init__(self, infile):
        self.infile = infile

    def __iter__(self):
        return self

    def __next__(self):
        while 1:
            line = next(self.infile)
            if line is None:
                raise StopIteration
            if line.strip() != "" and not line.startswith("#"):
                return line


class UnicodeCsvReader(object):

    def __init__(self, f, encoding="utf-8", **kwargs):
        self.csv_reader = csv.reader(f, **kwargs)
        self.encoding = encoding

    def __iter__(self):
        return self

    def __next__(self):
        # read and split the csv row into fields
        row = next(self.csv_reader)
        # now decode
        if PY3:
            return [str(cell, self.encoding) for cell in row]
        else:
            return [str(cell) for cell in row]

    def next(self):
        return self.__next__()

    @property
    def line_num(self):
        return self.csv_reader.line_num


class DictReader(csv.DictReader):
    """Like csv.DictReader, but skip lines starting with ``#``.
    """

    def __init__(self, infile, *args, **kwargs):
        csv.DictReader.__init__(self,
                                CommentStripper(infile),
                                *args, **kwargs)


class UnicodeDictReader(csv.DictReader):

    def __init__(self, f, encoding="utf-8", fieldnames=None, **kwds):
        csv.DictReader.__init__(self, f, fieldnames=fieldnames, **kwds)
        self.reader = UnicodeCsvReader(f, encoding=encoding, **kwds)


class DictReaderLarge:
    """Substitute for :py:class:`csv.DictReader` that handles very large
    fields.

    :py:mod:`csv` is implemented in C and limits the number of columns
    per table. This class has no such limit, but will not be as fast.

    This class is only a minimal implementation. For example, it does
    not handle dialects.
    """

    def __init__(self, infile, fieldnames, *args, **kwargs):
        self.mFile = infile
        self.mFieldNames = fieldnames
        self.mNFields = len(fieldnames)

    def __iter__(self):
        return self

    def __next__(self):

        line = next(self.mFile)
        if not line:
            raise StopIteration
        data = line[:-1].split("\t")
        assert len(data) == self.mNFields
        return dict(list(zip(self.mFieldNames, data)))


def readTable(infile,
              as_rows=True,
              with_header=True,
              ignore_incomplete=False,
              dialect="excel-tab"):
    """read a table from infile

    Arguments
    ---------
    infile : File
       File or list of lines
    as_rows : bool
       If true, return table as a list of rows. Otherwise,
       return as a list of columns.
    ignore_incomplete : bool
       If true, incomplete rows are ignored.
    dialect : string
       CSV dialect.

    Returns
    -------
    fields : list
        List of field names
    rows : list
        List of rows
    """

    lines = [x for x in infile if x[0] != "#"]

    if len(lines) == 0:
        return [], []

    if with_header:
        fields = lines[0][:-1].split("\t")
        del lines[0]
    else:
        fields = lines[0][:-1].split("\t")
        fields = list(map(str, list(range(len(fields)))))

    nfields = len(fields)

    try:
        reader = csv.reader(lines.__iter__(),
                            dialect=dialect)
    except TypeError:
        reader = csv.reader(lines.__iter__())

    table = list(reader)

    if ignore_incomplete:
        table = [x for x in table if len(x) == nfields]
    else:
        for r, row in enumerate(table):
            if len(row) != nfields:
                if not ignore_incomplete:
                    raise ValueError(
                        "missing elements in line %s, received=%s, "
                        "expected=%s" %
                        (r, str(row),  str(fields)))

                raise ValueError

    if not as_rows:
        table = list(zip(*table))

    return fields, table


def readTables(infile, *args, **kwargs):
    """read a set of csv tables from a single file.

    Tables within the file should be separated by `//`.

    See :func:`readTable` for additional arguments.

    Arguments
    ---------
    infile : File
       File or list of lines

    Returns
    -------
    tables : list
        A list of tuples (fields, data), one for each table.
    """

    lines = [x for x in infile if x[0] != "#"]
    chunks = [x for x in range(len(lines)) if lines[x][:2] == "//"]
    if not lines[-1].startswith("//"):
        chunks.append(len(lines))

    result = []

    start = 0
    for end in chunks:
        fields, table = readTable(lines[start:end], *args, **kwargs)
        result.append((fields, table))
        start = end + 1

    return result


def __DoGroup(rows, group_column, group_function, missing_value="na"):

    values = []
    for x in range(len(rows[0])):
        if x == group_column:
            values.append(rows[0][x])
        else:
            v = [x for x in [y[x] for y in rows] if x != missing_value]
            if len(v) == 0:
                values.append(missing_value)
            else:
                values.append(group_function([y[x] for y in rows]))

    return values


def groupTable(table,
               group_column=0,
               group_function=min,
               missing_value="na"):
    '''group table by *group_column*.

    The table need not be sorted.
    Arguments
    ---------
    table : list
        List of rows
    group_column : int
        Column to group on
    group_function : function
        Function to apply on grouped values
    missing_value : string
        String to use for missing values.
    '''

    table.sort(lambda x, y: cmp(x[group_column], y[group_column]))

    rows = []
    last_value = None
    new_table = []

    for row in table:
        if row[group_column] != last_value:

            if last_value is not None:
                new_table.append(
                    __DoGroup(rows, group_column, group_function,
                              missing_value))

            rows = []
            last_value = row[group_column]

        rows.append(row)

    if last_value is not None:
        new_table.append(
            __DoGroup(rows, group_column, group_function, missing_value))

    return new_table


def convertTable(table, columns, function=float,
                 skip_errors=False):
    """convert values in particular columns of a table to a new type.

    Arguments
    ---------
    table : list
        Rows in the table. Each row is a list or tuple.
    columns : list
        Indices of columns to convert
    function : function
        Function to apply for conversion
    skip_errors : bool
        If True, errors are ignored and rows with unconvertible
        values are ignored. The default is to raise a ValueError.

    Returns
    -------
    table : list
        A new table with the converted values.

    """
    new_table = []
    for row in table:
        skip = False
        for c in columns:
            try:
                row[c] = float(row[c])
            except ValueError as msg:
                if skip_errors:
                    skip = True
                    break
                else:
                    raise ValueError(msg)

        if not skip:
            new_table.append(row)

    return new_table
