'''
MatlabTools.py - 
======================================================

:Tags: Python

Code
----

'''
import sys
import string
import re
import numpy


def WriteMatrixOld(matrix, separator="\t"):
    # this is clumsy

    lines = []
    for x in range(0, matrix.shape[0]):
        lines.append(string.join(list(map(str, matrix[x, ])), separator))

    return string.join(lines, "\n")


def WriteMatrix(matrix, outfile=sys.stdout, separator="\t", format="%f",
                row_headers=None, col_headers=None):
    """write matrix to stream.
    """
    if col_headers:
        outfile.write(separator + separator.join(col_headers) + "\n")

    for x in range(0, matrix.shape[0]):
        if row_headers:
            outfile.write(row_headers[x] + separator)
        outfile.write(
            string.join([format % x for x in matrix[x, ]], separator) + "\n")


def ReadMatrix(file,
               separator="\t",
               numeric_type=numpy.float,
               take="all",
               headers=False
               ):
    """read a matrix. There probably is a routine for this in Numpy, which
    I haven't found yet.
    """

    lines = [x for x in file.readlines() if x[0] != "#"]

    row_headers, col_headers = [], []

    if headers:
        col_headers = lines[0][:-1].split("\t")[1:]
        del lines[0]

    num_rows = len(lines)

    if take != "all":
        num_cols = len(take)
    else:
        l = len(string.split(lines[0][:-1], "\t"))
        if headers:
            take = list(range(1, l))
        else:
            take = list(range(0, l))

    num_cols = len(take)

    matrix = numpy.zeros((num_rows, num_cols), numeric_type)

    nrow = 0
    for l in lines:
        data = l[:-1].split("\t")
        if headers:
            row_headers.append(data[0])

        try:
            data = [float(data[x]) for x in take]
        except ValueError:
            print("error parsing data", data)
            raise

        matrix[nrow] = data
        nrow += 1

    return matrix, row_headers, col_headers


def ReadSparseMatrix(filename,
                     separator="\t",
                     numeric_type=numpy.float,
                     is_symmetric=None):
    """read sparse matrix."""

    lines = [string.split(x[:-1], separator)[:3]
             for x in [x for x in open(filename, "r").readlines() if x[0] != "#"]]

    data = [(string.atoi(x[0]), string.atoi(x[1]), string.atof(x[2]))
            for x in lines]

    num_rows = max([x[0] for x in data])
    num_cols = max([x[1] for x in data])

    if (is_symmetric):
        num_rows = max(num_rows, num_cols)
        num_cols = max(num_rows, num_cols)

    matrix = numpy.zeros((num_rows, num_cols), numeric_type)

    for row, col, weight in data:
        matrix[row - 1, col - 1] = weight

    if is_symmetric:
        for col, row, weight in data:
            matrix[row - 1, col - 1] = weight

    return matrix


def ReadBinarySparseMatrix(filename,
                           separator="\t",
                           numeric_type=numpy.float,
                           is_symmetric=None):
    """read sparse matrix."""

    lines = [string.split(x[:-1], separator)[:2]
             for x in [x for x in open(filename, "r").readlines() if x[0] != "#"]]

    data = [(string.atoi(x[0]), string.atoi(x[1])) for x in lines]

    num_rows = max([x[0] for x in data])
    num_cols = max([x[1] for x in data])

    if (is_symmetric):
        num_rows = max(num_rows, num_cols)
        num_cols = max(num_rows, num_cols)
    matrix = numpy.zeros((num_rows, num_cols), numeric_type)

    for row, col in data:
        matrix[row - 1, col - 1] = 1

    if is_symmetric:
        for col, row in data:
            matrix[row - 1, col - 1] = 1

    return matrix


def readMatrix(infile,
               format="full",
               separator="\t",
               numeric_type=numpy.float,
               take="all",
               headers=True,
               missing=None,
               ):
    """read a matrix from file ane return a numpy matrix.

    formats accepted are:
    * full
    * sparse
    * phylip
    """

    row_headers, col_headers = [], []

    lines = [x for x in infile.readlines() if x[0] != "#"]

    if len(lines) == 0:
        raise IOError("no input")

    if format == "full":

        if headers:
            col_headers = lines[0][:-1].split("\t")[1:]
            del lines[0]

        num_rows = len(lines)

        if take != "all":
            num_cols = len(take)
        else:
            l = len(string.split(lines[0][:-1], "\t"))
            if headers:
                take = list(range(1, l))
            else:
                take = list(range(0, l))

        num_cols = len(take)

        matrix = numpy.zeros((num_rows, num_cols), numeric_type)

        nrow = 0
        for l in lines:
            data = l[:-1].split("\t")
            if headers:
                row_headers.append(data[0])

            if missing is None:
                try:
                    data = [float(data[x]) for x in take]
                except ValueError as msg:
                    raise ValueError("error %s: data=%s" % (msg, str(data)))
                except IndexError as msg:
                    raise IndexError("error %s: data=%s" % (msg, str(data)))

            else:
                d = []
                for x in take:
                    try:
                        d.append(float(data[x]))
                    except ValueError:
                        d.append(missing)
                    except IndexError as msg:
                        raise IndexError(
                            "error %s: data=%s" % (msg, str(data)))

                data = d

            matrix[nrow] = data

            nrow += 1

    elif format == "phylip":
        # read in symmetric phylip matrices
        # note: they can wrap around
        if take != "all":
            raise ValueError(
                "phylip matrix does not support take - "
                "only full matrices are processed.")

        if not headers:
            raise ValueError("phylip matrix always has headers.")

        num_rows = int(lines[0].strip())

        num_cols = num_rows

        matrix = numpy.zeros((num_rows, num_cols), numeric_type)
        take = list(range(1, num_rows))

        nrow = 0
        ncol = 0
        for l in lines[1:]:

            data = re.split("\s+", l[:-1])

            if ncol == 0:
                row_headers.append(data[0])

            try:
                data = list(map(float, data[1:len(data)]))
            except ValueError:
                raise ValueError(
                    "parsing error in conversion to "
                    "float in line %s" % l)

            for x in range(len(data)):
                matrix[nrow][ncol] = data[x]
                ncol += 1

            # deal with wrapping around
            if ncol == num_cols:
                ncol = 0
                nrow += 1

        col_headers = row_headers

    return matrix, row_headers, col_headers


def writeMatrix(outfile, matrix,
                format="full",
                separator="\t",
                value_format="%f",
                row_headers=None,
                col_headers=None):
    """write matrix to stream.
    """

    if format == "full":
        if col_headers:
            outfile.write(separator + separator.join(col_headers) + "\n")

        for x in range(0, matrix.shape[0]):
            if row_headers:
                outfile.write(row_headers[x] + separator)
            outfile.write(
                string.join([value_format % x for x in matrix[x, ]], separator) + "\n")

    elif format == "phylip":
        if not row_headers:
            raise ValueError("phylip output requires row headers")

        nrows = len(row_headers)
        outfile.write("%i\n" % nrows)

        for x in range(0, nrows):
            outfile.write(row_headers[x] + separator)
            outfile.write(
                separator.join([value_format % y for y in matrix[x, ]]) + "\n")

    else:
        raise ValueError("unknown output format %s" % format)
