'''
MatrixTools.py - 
======================================================

:Tags: Python

Code
----

'''
import string
import collections
import numpy
import pandas
from CGAT import IOTools as IOTools


class Matrix:

    def __init__(self):
        pass


def addOptions(parser):
    """add matrices to option parser."""

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="format.")

    parser.add_option("--row-names", dest="row_names", type="string",
                      help="list of row names separated by ','.")

    parser.add_option("--col-names", dest="col_names", type="string",
                      help="list of col names separated by ','.")

    parser.add_option("--file-row-names", dest="file_row_names", type="string",
                      help="filename with row names.")

    parser.add_option("--file-col-names", dest="file_col_names", type="string",
                      help="filename with col names.")

    parser.add_option("--numeric", dest="numeric", action="store_true",
                      help="row and column titles are numeric.")

    parser.add_option("--is-asymmetric", dest="asymmetric", action="store_true",
                      help="matrix is asymmetric.")

    parser.add_option("--default-value", dest="default", type="string",
                      help="default value for missing values.")

    parser.add_option("--default-diagonal", dest="default_diagonal", type="string",
                      help="default value for missing values on diagonal.")

    parser.add_option("--input-format", dest="input_format", type="choice",
                      choices=(
                          "row-col-weight", "", "row-col-weight-replicates"),
                      help="input format.")

    parser.set_defaults(
        default="0",
        default_diagonal="0",
        format="string",
        asymmetric=False,
        is_numeric=False,
        row_names=None,
        col_names=None,
        file_row_names=None,
        file_col_names=None,
        input_format="row-col-weight",
    )


def getMapTokens(options):

    map_token2row, map_token2col = {}, {}

    if options.file_row_names:
        row_tokens = [string.split(
            x[:-1], "\t")[0] for x in open(options.file_row_names, "r").readlines()]
        for row_token in row_tokens:
            map_token2row[row_token] = len(map_token2row)

    if options.row_names:
        for x in options.row_names.split(","):
            map_token2row[x] = len(map_token2row)

    if options.file_col_names:
        col_tokens = [string.split(
            x[:-1], "\t")[0] for x in open(options.file_col_names, "r").readlines()]
        for col_token in col_tokens:
            map_token2col[col_token] = len(map_token2col)

    if options.col_names:
        for x in options.col_names.split(","):
            map_token2col[x] = len(map_token2col)

    if not options.asymmetric:
        if map_token2row and not map_token2col:
            map_token2col = map_token2row
        elif map_token2col and not map_token2row:
            map_token2row = map_token2col

    return map_token2row, map_token2col


def getMatrixFromEdges(lines, options,
                       in_map_token2row={},
                       in_map_token2col={}):
    """read matrix from lines
    """

    # remove comments
    lines = [x for x in lines if x[0] != "#" and len(x[:-1]) > 0]

    if in_map_token2row:
        map_token2row = in_map_token2row
    else:
        map_token2row = {}

    if in_map_token2col:
        map_token2col = in_map_token2col
    else:
        map_token2col = {}

    if options.format == "string":

        has_row_names = len(map_token2row) > 0
        has_col_names = len(map_token2col) > 0

        # if either row/column names are not given:
        if not map_token2row or not map_token2col:

            row_tokens = [string.split(x[:-1], "\t")[0] for x in lines]
            col_tokens = [string.split(x[:-1], "\t")[1] for x in lines]

            if options.is_numeric:
                row_tokens = list(map(float, row_tokens))
                col_tokens = list(map(float, col_tokens))
                row_tokens.sort()
                col_tokens.sort()
                row_tokens = list(map(str, row_tokens))
                col_tokens = list(map(str, col_tokens))
            else:
                row_tokens.sort()
                col_tokens.sort()

            if not has_row_names:
                for row_token in row_tokens:
                    if row_token not in map_token2row:
                        map_token2row[row_token] = len(map_token2row)
            if not has_col_names:
                for col_token in col_tokens:
                    if col_token not in map_token2col:
                        map_token2col[col_token] = len(map_token2col)

        if not options.asymmetric:
            for col_token in list(map_token2col.keys()):
                if col_token not in map_token2row:
                    map_token2row[col_token] = len(map_token2row)
            map_token2col = map_token2row

        matrix = [[options.default for j in range(len(map_token2col))] for i in range(
            len(map_token2row))]

        if len(map_token2col) == len(map_token2row):
            for j in range(len(map_token2col)):
                matrix[j][j] = options.default_diagonal

        # return matrix
        m = Matrix()

        if options.input_format == "row-col-weight":
            for line in lines:
                row_token, col_token, weight = string.split(
                    line[:-1], "\t")[:3]
                matrix[map_token2row[row_token]][
                    map_token2col[col_token]] = weight
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][
                        map_token2row[row_token]] = weight

        elif options.input_format == "row-col-weight-replicates":
            replicates = [[0 for j in range(len(map_token2col))]
                          for i in range(len(map_token2row))]
            for line in lines:
                row_token, col_token, weight, n = string.split(
                    line[:-1], "\t")[:4]
                matrix[map_token2row[row_token]][
                    map_token2col[col_token]] = weight
                replicates[map_token2row[row_token]][
                    map_token2col[col_token]] = int(n)
                if not options.asymmetric:
                    matrix[map_token2col[col_token]][
                        map_token2row[row_token]] = weight
                    replicates[map_token2col[col_token]][
                        map_token2row[row_token]] = int(n)
            m.mReplicates = replicates

        col_tokens = list(map_token2col.items())
        col_tokens.sort(lambda x, y: cmp(x[1], y[1]))
        row_tokens = list(map_token2row.items())
        row_tokens.sort(lambda x, y: cmp(x[1], y[1]))

        m.mMatrix = matrix
        m.mMapRow2Token = row_tokens
        m.mMapCol2Token = col_tokens
        m.mMapToken2Row = map_token2row
        m.mMapToken2Col = map_token2col

    return matrix


def buildMatrixFromLists(lists, dtype=numpy.float, default=None):
    '''build a matrix from a list of lists.

    Each list is a list of tuples (row, value).
    The columns are given by order of the lists.

    Returns matrix, row_headers
    '''

    all_rows = collections.defaultdict()
    for l in lists:
        all_rows.update([(x[0], 0) for x in l])
    for x, v in enumerate(all_rows.items()):
        all_rows[v[0]] = x

    matrix = numpy.zeros((len(all_rows), len(lists)), dtype=dtype)
    if default:
        matrix.fill(default)

    for col, l in enumerate(lists):
        for row, value in l:
            matrix[all_rows[row], col] = value
    return matrix, list(all_rows.keys())


def buildMatrixFromTables(infiles, column, column_header=0, dtype=numpy.float,
                          default=None):
    '''build a matrix from a column called *column* in a series of input files.

    If column_value is None, the first column is taken as the name of the row.

    The columns are given by order of the input files.

    returns matrix, row_headers
    '''

    lists = []
    for infile in infiles:
        data = pandas.read_table(IOTools.openFile(infile))
        lists.append(list(zip(list(data[column_header]), list(data[column]))))

    return buildMatrixFromLists(lists, dtype=dtype, default=default)


def buildMatrixFromEdges(edges,
                         in_map_token2row={},
                         in_map_token2col={},
                         is_symmetric=False,
                         missing_value=0,
                         diagonal_value=0,
                         dtype=numpy.int):
    """build a matrix from an edge-list representation.

    For example, the following list of tuples::

       [('A', 'B', 1),
        ('A', 'C', 2),
        ('B', 'C', 3)]

    will be converted to the following matrix::

         A B C
       A   1 2
       B     3
       C 

    If *is_symmetric* is set to True, the matrix is assumed to be
    symmetric and missing values will automatically be filled::

         A B C
       A   1 2
       B 1   3
       C 2 3

    If edge list may contain four elements, in which case the
    fourth element is expected to be the value of the lower
    diagonal in a symmetric matrix::

       [('A', 'B', 1, 4),
        ('A', 'C', 2, 5),
        ('B', 'C', 3, 6)]

    will yield::

         A B C
       A   1 2
       B 4   3
       C 5 6


    returns a numpy matrix and lists of row and column names.
    """

    if in_map_token2row:
        map_token2row = in_map_token2row
    else:
        map_token2row = {}

    if in_map_token2col:
        map_token2col = in_map_token2col
    else:
        map_token2col = {}

    has_row_names = len(map_token2row) > 0
    has_col_names = len(map_token2col) > 0

    # if either row/column names are not given:
    if not map_token2row or not map_token2col:

        row_tokens = sorted(list(set([x[0] for x in edges])))
        col_tokens = sorted(list(set([x[1] for x in edges])))

        if not has_row_names:
            for row_token in row_tokens:
                if row_token not in map_token2row:
                    map_token2row[row_token] = len(map_token2row)
        if not has_col_names:
            for col_token in col_tokens:
                if col_token not in map_token2col:
                    map_token2col[col_token] = len(map_token2col)

        if is_symmetric:
            for col_token in list(map_token2col.keys()):
                if col_token not in map_token2row:
                    map_token2row[col_token] = len(map_token2row)
            map_token2col = map_token2row

        matrix = numpy.matrix([missing_value] *
                              len(map_token2row) * len(map_token2col),
                              dtype=numpy.dtype)
        matrix.shape = (len(map_token2row), len(map_token2col))

        if len(map_token2col) == len(map_token2row):
            for j in range(len(map_token2col)):
                matrix[j][j] = diagonal_value

        if len(edges[0] == 3):
            if is_symmetric:
                for row, col, value in edges:
                    matrix[map_token2row[row], map_token2col[col]] = \
                        matrix[map_token2row[col], map_token2col[row]] = value
            else:
                for row, col, value in edges:
                    matrix[map_token2row[row], map_token2col[col]] = value
        elif len(edges[0] == 4):
            for row, col, value1, value2 in edges:
                matrix[map_token2row[row], map_token2col[col]] = value1
                matrix[map_token2row[col], map_token2col[row]] = value2
        else:
            raise ValueError(
                "unexpected number of elements in list, expected 3 or 4, "
                "got %i" % (len(edges[0])))

    col_tokens = list(map_token2col.items())
    col_tokens.sort(lambda x, y: cmp(x[1], y[1]))
    row_tokens = list(map_token2row.items())
    row_tokens.sort(lambda x, y: cmp(x[1], y[1]))

    return matrix, row_tokens, col_tokens
