"""
Histogram2D.py - functions for handling two-dimensional histograms.
===================================================================

:Tags: Python

"""

import string


def Calculate(values, mode=0, bin_function=None):
    """Return a list of (value, count) pairs, summarizing the input values.
    Sorted by increasing value, or if mode=1, by decreasing count.

    If bin_function is given, map it over values first.

    """

    if bin_function:
        values = list(map(bin_function, values))

    bins = {}
    for val in values:
        v = "%f-%f" % tuple(val)
        bins[v] = bins.get(v, 0) + 1

    bb = list(bins.items())

    if mode:
        bb.sort(lambda x, y: cmp(y[1], x[1]))
    else:
        bb.sort()

    r = []
    for v, n in bb:
        x, y = list(map(string.atof, string.split(v, "-")))
        r.append((x, y, n))

    return r


def Print(h, bin_function=None):
    """print a histogram.

    A histogram can either be a list/tuple of values or
    a list/tuple of lists/tuples where the first value contains
    the bin and second contains the values (which can again be
    a list/tuple).

    :param format: output format. 
        0 = print histogram in several lines,
        1 = print histogram on single line
    """

    if bin_function:
        h = list(map(bin_function, h))

    for hh in h:
        print(string.join(list(map(str, hh)), "\t"))
