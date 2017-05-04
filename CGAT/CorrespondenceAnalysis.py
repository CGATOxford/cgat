'''
CorrespondenceAnalysis.py - 
======================================================

:Tags: Python

Code
----

'''
import numpy
import numpy.linalg
import numpy.linalg.linalg


def GetIndices(matrix):
    """return order (1st eigenvector) of row and column indicies.

    """

    original_nrows, original_ncols = matrix.shape
    # calculate row and column sums
    row_sums = numpy.array(numpy.sum(matrix, 1)).flatten()
    col_sums = numpy.array(numpy.sum(matrix, 0)).flatten()

    # check for empty rows/columns
    # remove rows/columns that are empty
    if 0 in row_sums or 0 in col_sums:
        # map rows with data to original rows
        row_map_nonempty = numpy.arange(original_nrows)[row_sums != 0]
        row_map_nonempty = dict(enumerate(row_map_nonempty))
        col_map_nonempty = numpy.arange(original_ncols)[col_sums != 0]
        col_map_nonempty = dict(enumerate(col_map_nonempty))

        # truncate matrix
        matrix = matrix[row_sums != 0, :][:, col_sums != 0]

        # filter row and column sums
        row_sums = row_sums[row_sums != 0]
        col_sums = col_sums[col_sums != 0]

    nrows, ncols = matrix.shape

    a = numpy.zeros((nrows, nrows), numpy.float)
    for x in range(0, nrows):
        a[x, x] = 1.0 / float(row_sums[x])

    b = numpy.zeros((ncols, ncols), numpy.float)
    for x in range(0, ncols):
        b[x, x] = 1.0 / float(col_sums[x])

    M = numpy.dot(
        a, numpy.dot(
            matrix, numpy.dot(
                b, numpy.transpose(matrix))))

    try:
        row_eigenvector = numpy.linalg.eig(M)[1][:, 1]
    except numpy.linalg.linalg.LinAlgError as msg:
        raise ValueError(msg)

    M = numpy.dot(
        b, numpy.dot(
            numpy.transpose(matrix), numpy.dot(
                a, matrix)))

    try:
        col_eigenvector = numpy.linalg.eig(M)[1][:, 1]
    except numpy.linalg.linalg.LinAlgError as msg:
        raise ValueError(msg)

    row_eigenvector = row_eigenvector.astype(numpy.float)
    col_eigenvector = col_eigenvector.astype(numpy.float)

    # insert columns ignored at the computation and give them the lowest
    # eigenvalue
    row_values = [row_eigenvector.min() - 1.0] * original_nrows
    for x, y in enumerate(row_eigenvector):
        row_values[row_map_nonempty[x]] = y
    col_values = [col_eigenvector.min() - 1.0] * original_ncols
    for x, y in enumerate(col_eigenvector):
        col_values[col_map_nonempty[x]] = y

    return row_values, col_values


def GetPermutatedMatrix(matrix,
                        map_row_new2old, map_col_new2old,
                        row_headers=None, col_headers=None):
    """return a permuted matrix. Note, that currently this is very
    inefficient, as I do not know how to do this in numpy.
    """

    nrows, ncols = matrix.shape

    result = numpy.zeros((nrows, ncols), matrix.dtype)
    for r in range(0, nrows):
        for c in range(0, ncols):
            result[r, c] = matrix[map_row_new2old[r], map_col_new2old[c]]

    if not row_headers or not col_headers:
        return result

    rows = []
    for x in map_row_new2old:
        rows.append(row_headers[x])

    cols = []
    for x in map_col_new2old:
        cols.append(col_headers[x])

    return result, rows, cols


def PermuteRows(matrix):
    pass


if __name__ == "__main__":

    num_rows = 6
    num_cols = 5
    matrix = numpy.zeros((num_rows, num_cols), numpy.int)

    matrix[0, 2] = 1
    matrix[0, 3] = 1
    matrix[1, 0] = 1
    matrix[1, 1] = 1
    matrix[1, 4] = 1
    matrix[2, 1:5] = 1
    matrix[3, 2:4] = 1
    matrix[4, 0] = 1
    matrix[4, 4] = 1
    matrix[5, 0] = 1
    matrix[5, 2:5] = 1

    print("matrix=", matrix)

    row_indices, col_indices = GetIndices(matrix)

    map_row_new2old = numpy.argsort(row_indices)
    map_col_new2old = numpy.argsort(col_indices)

    print(map_row_new2old)
    print(map_col_new2old)

    print(GetPermutatedMatrix(matrix, map_row_new2old, map_col_new2old))
