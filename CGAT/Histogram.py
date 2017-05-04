"""
Histogram.py - Various functions to deal with histograms
===========================================================

:Author:
:Tags: Python

Histograms can be calculated from a list/tuple/array of
values. The histogram returned is then a list of tuples
of the format [(bin1,value1), (bin2,value2), ...].

"""

import sys
import re
import math
import scipy
import scipy.stats
import bisect
import numpy
from functools import reduce


def CalculateFromTable(dbhandle,
                       field_name,
                       from_statement,
                       num_bins=None,
                       min_value=None,
                       max_value=None,
                       intervals=None,
                       increment=None):
    """get a histogram using an SQL-statement.
    Intervals can be either supplied directly or are build
    from the data by providing the number of bins and optionally
    a minimum or maximum value.

    If no number of bins are provided, the bin-size is 1.

    This command uses the INTERVAL command from MYSQL, i.e. a bin value
    determines the upper boundary of a bin.
    """

    if not min_value:
        min_value = int(math.floor(dbhandle.Execute(
            "SELECT MIN(%s) %s" % (field_name, from_statement)).fetchone()[0]))

    if not max_value:
        max_value = int(math.ceil(dbhandle.Execute(
            "SELECT MAX(%s) %s" % (field_name, from_statement)).fetchone()[0]))

    if increment:
        step_size = increment
    elif num_bins:
        step_size = int(float(max_value - min_value) / float(num_bins))
    else:
        step_size = 1

    if not intervals:
        intervals = list(range(min_value, max_value, step_size))

    i_string = ",".join(list(map(str, intervals)))

    statement = "SELECT INTERVAL( %s, %s )-1 AS i, COUNT(*) %s GROUP BY i" % (
        field_name, i_string, from_statement)

    return Convert(dbhandle.Execute(statement).fetchall(), intervals)


def CalculateConst(values,
                   num_bins=None,
                   min_value=None,
                   max_value=None,
                   intervals=None,
                   increment=None,
                   combine=None):
    """calculate a histogram based on a list or tuple of values.
    """

    if not min_value:
        min_value = int(math.floor(min(values)))

    if not max_value:
        max_value = int(math.ceil(max(values))) + 1

    if increment:
        step_size = increment
    elif num_bins:
        step_size = int(float(max_value - min_value) / float(num_bins))
    else:
        step_size = 1

    if not intervals:
        intervals = list(range(min_value, max_value, step_size))

    histogram = [0] * len(intervals)
    for v in values:
        i = 0
        while i < len(intervals) and v > intervals[i]:
            i += 1
        if i < len(intervals):
            histogram[i] += 1

    return intervals, histogram


def Calculate(values,
              num_bins=None,
              min_value=None,
              max_value=None,
              intervals=None,
              increment=None,
              combine=None,
              no_empty_bins=0,
              dynamic_bins=False,
              ignore_out_of_range=True):
    """calculate a histogram based on a list or tuple of values.

    use scipy for calculation.
    """

    if len(values) == 0:
        return []

    if not intervals:

        if min_value is None:
            min_value = min(values)

        if max_value is None:
            max_value = max(values)

        if dynamic_bins:
            intervals = list(
                set([x for x in values if min_value <= x <= max_value]))
            intervals.sort()
        else:
            if increment:
                step_size = increment
            elif num_bins and max_value:
                step_size = float(max_value - min_value) / float(num_bins)
            else:
                step_size = 1.0

            num_bins = int(
                math.ceil((float(max_value) - float(min_value)) / float(step_size)))
            intervals = [float(min_value) + float(x) * float(step_size)
                         for x in range(num_bins + 1)]

    if not ignore_out_of_range:
        new_values = []
        for v in values:
            if v < min_value:
                v = min_value
            elif v > max_value:
                v = max_value
            new_values.append(v)

        values = new_values

    return Convert(scipy.stats.histogram2(values, intervals), intervals, no_empty_bins)


def Scale(h, scale=1.0):
    """rescale bins in histogram.
    """
    n = []
    for b, v in h:
        n.append((b * scale, v))
    return n


def Convert(h, i, no_empty_bins=0):
    """add bins to histogram.
    """
    n = []
    for x in range(0, len(h)):
        if no_empty_bins and h[x] == 0:
            continue
        n.append((i[x], h[x]))
    return n


def Combine(source_histograms, missing_value=0):
    """combine a list of histograms
    Each histogram is a sorted list of bins and counts.
    The counts can be tuples.
    """

    new_bins = {}

    # get all bins
    for h in source_histograms:
        for data in h:
            new_bins[data[0]] = []

    # add values
    length = 1
    l = 0
    for h in source_histograms:

        # add data for used bins
        for data in h:
            bin = data[0]
            if len(data) != 2:
                v = data[1:]
            else:
                v = data[1]

            if isinstance(v, list) or isinstance(v, tuple):
                l = len(v)
                for x in v:
                    new_bins[bin].append(x)
            else:
                l = 1
                new_bins[bin].append(v)

        # add missing value for unused bins
        for b in list(new_bins.keys()):
            if len(new_bins[b]) < length:
                for x in range(0, l):
                    new_bins[b].append(missing_value)

        length += l

    return __ConvertToList(new_bins)


def Print(h, intervalls=None, format=0, nonull=None, format_value=None, format_bin=None):
    """print a histogram.

    A histogram can either be a list/tuple of values or
    a list/tuple of lists/tuples where the first value contains
    the bin and second contains the values (which can again be
    a list/tuple).

    format
       0 = print histogram in several lines
       1 = print histogram on single line

    """

    Write(sys.stdout, h, intervalls, format, nonull, format_value, format_bin)


def Write(outfile, h, intervalls=None, format=0, nonull=None,
          format_value=None, format_bin=None):
    """print a histogram.

    A histogram can either be a list/tuple of values or
    a list/tuple of lists/tuples where the first value contains
    the bin and second contains the values (which can again be
    a list/tuple).

    :param format: output format.
        0 = print histogram in several lines,
        1 = print histogram on single line

    """

    lines = []

    if len(h) == 0:
        return

    if format_value:
        def fv(x):
            if x == "na":
                return x
            else:
                return format_value % x
    else:
        fv = str

    if format_bin:
        fb = lambda x: format_bin % x
    else:
        fb = str

    if not isinstance(h[0], list) and not isinstance(h[0], tuple):
        for x in range(0, len(h)):

            if intervalls:
                lines.append(fb(intervalls[x]) + "\t" + fv(h[x]))
            else:
                lines.append(fb(x) + "\t" + fv(h[x]))
    else:
        for x, v in h:
            if isinstance(v, list) or isinstance(v, tuple):
                val = "\t".join(list(map(fv, v)))
            else:
                val = fv(v)

            if intervalls:
                lines.append(fb(intervalls[x]) + "\t" + val)
            else:
                lines.append(fb(x) + "\t" + val)

    # print values
    if nonull is not None:
        for l in range(0, len(lines)):
            lines[l] = re.sub("\t0", "\t%s" % nonull, lines[l])

    if format == 0:
        outfile.write("\n".join(lines) + "\n")
    elif format == 1:
        outfile.write(" | ".join(lines) + "\n")
    else:
        outfile.write("\n".join(lines) + "\n")


def Fill(h):
    """fill every empty value in histogram with
    previous value.
    """

    new_h = []

    x, v = h[0]

    if isinstance(v, list) or isinstance(v, tuple):
        l = len(v)
        previous_v = [0] * l
    else:
        previous_v = 0

    for x, v in h:
        if isinstance(v, list) or isinstance(v, tuple):
            for i in range(0, l):
                if v[i] == 0:
                    v[i] = previous_v[i]
                else:
                    previous_v[i] = v[i]
        else:
            if v == 0:
                v = previous_v
            else:
                previous_v = v

        new_h.append((x, v))

    return new_h


def Normalize(h):

    # first count totals
    if isinstance(h[0][1], list) or isinstance(h[0][1], tuple):
        l = len(h[0][1])
        is_list = 1
    else:
        is_list = 0
        l = 1

    totals = [0.0] * l

    for bin, v in h:
        if is_list:
            for x in range(0, l):
                try:
                    totals[x] += v[x]
                except TypeError:
                    pass
        else:
            totals[0] += v

    for x in range(0, l):
        if totals[x] == 0:
            totals[x] = 1

    # first count totals
    new_histogram = []
    for bin, v in h:
        if is_list:
            vv = []
            for x in range(0, l):
                vv.append(float(v[x]) / totals[x])
        else:
            vv = float(v) / totals[0]

        new_histogram.append((bin, vv))

    return new_histogram


def Add(h1, h2):
    """adds values of histogram h1 and h2 and
    returns a new histogram
    """

    new_bins = {}
    # get all bins
    for h in (h1, h2):
        if h:
            for bin, v in h:
                new_bins[bin] = 0

    for h in (h1, h2):
        if h:
            for bin, v in h:
                new_bins[bin] += v

    return __ConvertToList(new_bins)


def __ConvertToList(new_bins):
    """converts a hash to a histogram.
    """

    # convert to list
    keys = list(new_bins.keys())
    keys.sort()
    new_histogram = []

    for k in keys:
        new_histogram.append((k, new_bins[k]))

    return new_histogram


def SmoothWrap(histogram, window_size):
    """smooth histogram by sliding window-method, where
    the window is wrapped around the borders. The sum of
    all values is entered at center of window.
    """

    new_histogram = [0] * len(histogram)

    half_window_size = window_size / 2
    length = len(histogram)

    # 1. start with window
    cumul = 0
    for i in range(length - half_window_size, length):
        cumul = cumul + histogram[i]
    for i in range(0, half_window_size + 1):
        cumul = cumul + histogram[i]

    # 2. iterate over histogram and add values over windows_size
    y = length - half_window_size
    z = half_window_size
    for i in range(0, length):

        new_histogram[i] = cumul

        y = y + 1
        z = z + 1
        if y >= length:
            y = 0
        if z >= length:
            z = 0
        cumul = cumul - histogram[y] + histogram[z]

    return new_histogram


def GetMaximumIndex(histogram):
    return histogram.index(max(histogram))


def GetMinimumIndex(histogram):
    return histogram.index(min(histogram))


def PrintAscii(histogram, step_size=1):
    """print histogram ascii-style.
    """

    l = len(histogram)
    m = max(histogram)

    if m == 0:
        print("----> histogram is empty")
        return

    f = 100.0 / m

    print("----> histogram: len=%i, max=%i" % (l, m))
    for x in range(1, l, step_size):
        s = "|"
        s += " " * (int(histogram[x] * f) - 1) + "*"

        print("%5i" % x, s)


def Count(data):
    """count categorized data. Returns a list
    of tuples with (count, token).
    """

    counts = []
    data.sort()

    last_c = None
    for c in data:
        if last_c != c:
            if last_c:
                counts.append((t, last_c))
            last_c = c
            t = 0
        t += 1

    counts.append((t, last_c))

    return counts


def Accumulate(h, num_bins=2, direction=1):
    """add successive counts in histogram.
    Bins are labelled by group average.
    """

    if len(h) == []:
        return []

    new_histogram = []

    if isinstance(h[0][1], list) or isinstance(h[0][1], tuple):
        l = len(h[0][1])
        is_list = 1
    else:
        is_list = 0
        l = 1

    if direction != 1:
        h.reverse()

    i = 0
    if is_list:
        vv = [0] * l
        bb = 0
        for b, v in h:

            bb += b
            for x in range(0, l):
                vv[x] += v[x]

            i += 1
            if i % num_bins == 0:
                new_histogram.append((float(bb) / float(num_bins), vv))
                vv = [0] * l
                bb = 0

        if (i % num_bins):
            new_histogram.append((float(bb) / float(i % num_bins), vv))
    else:
        vv = 0
        bb = 0
        for b, v in h:
            bb += b
            vv += v

            i += 1
            if i % num_bins == 0:
                new_histogram.append(float(bb) / float(num_bins), vv)
                vv = 0
                bb = 0

        if (i % num_bins):
            new_histogram.append(float(bb) / float(i % num_bins), vv)

    # reorder h and new_histogram
    if direction != 1:
        h.reverse()
        new_histogram.reverse()

    return new_histogram


def Cumulate(h, direction=1):
    """calculate cumulative distribution.
    """

    if len(h) == []:
        return []

    new_histogram = []

    if isinstance(h[0][1], list) or isinstance(h[0][1], tuple):
        l = len(h[0][1])
        is_list = 1
    else:
        is_list = 0
        l = 1

    if direction != 1:
        h.reverse()

    if is_list:
        vv = [0] * l
        for b, v in h:
            for x in range(0, l):
                vv[x] += v[x]
            new_histogram.append((b, [x for x in vv]))
    else:
        vv = 0
        for b, v in h:
            vv += v
            new_histogram.append((b, vv))

    # reorder h and new_histogram
    if direction != 1:
        h.reverse()
        new_histogram.reverse()

    return new_histogram


def AddRelativeAndCumulativeDistributions(h):
    """adds relative and cumulative percents to a histogram.
    """

    if len(h) == []:
        return []

    new_histogram = []
    total = float(reduce(lambda x, y: x + y, [x[1] for x in h]))

    cumul_down = int(total)
    cumul_up = 0

    for bin, val in h:
        percent = float(val) / total
        cumul_up += val
        percent_cumul_up = float(cumul_up) / total
        percent_cumul_down = float(cumul_down) / total

        new_histogram.append(
            (bin, (val, percent, cumul_up, percent_cumul_up, cumul_down, percent_cumul_down)))
        cumul_down -= val

    return new_histogram


def histogram(values, mode=0, bin_function=None):
    """Return a list of (value, count) pairs, summarizing the input values.
    Sorted by increasing value, or if mode=1, by decreasing count.
    If bin_function is given, map it over values first.
    Ex: vals = [100, 110, 160, 200, 160, 110, 200, 200, 220]
    histogram(vals) ==> [(100, 1), (110, 2), (160, 2), (200, 3), (220, 1)]
    histogram(vals, 1) ==> [(200, 3), (160, 2), (110, 2), (100, 1), (220, 1)]
    histogram(vals, 1, lambda v: round(v, -2)) ==> [(200.0, 6), (100.0, 3)]"""

    if bin_function:
        values = list(map(bin_function, values))
    bins = {}
    for val in values:
        bins[val] = bins.get(val, 0) + 1
        if mode:
            return sort(list(bins.items()), lambda x, y: cmp(y[1], x[1]))
        else:
            return sort(list(bins.items()))


def cumulate(histogram):
    """cumulate histogram in place.

    histogram is list of (bin, value) or (bin, (values,) )
    """

    if len(histogram) == 0:
        return

    if isinstance(h[0][1], list) or isinstance(h[0][1], tuple):
        n = len(histogram[0][1])
        l = [0] * n
        for bin, vv in histogram:
            for x in range(n):
                vv[x] += l[x]
                l[x] = vv[x]
    else:
        l = 0
        for x in range(len(histogram)):
            histogram[x] = (histogram[x][0], histogram[x][1] + l)
            l = histogram[x][1]


def normalize(histogram):
    """normalize histogram in place.

    histogram is list of (bin, value) or (bin, (values,) )
    """

    if len(histogram) == 0:
        return

    if isinstance(histogram[0][1], list) or isinstance(histogram[0][1], tuple):
        n = len(histogram[0][1])
        m = [0] * n
        for bin, d in histogram:
            for x in range(n):
                m[x] = max(m[x], d[x])
        for bin, d in histogram:
            for x in range(n):
                if m[x] > 0:
                    d[x] = float(d[x]) / m[x]
    else:
        m = float(max([x[1] for x in histogram]))
        if m > 0:
            for x in range(len(histogram)):
                histogram[x] = (histogram[x][0], histogram[x][1] / m)


def fill(iterator, bins):
    """fill a histogram from bins. 

    The values are given by an iterator so that the histogram
    can be built on the fly.

    Description:

    Count the number of times values from array a fall into
    numerical ranges defined by bins.  Range x is given by
    bins[x] <= range_x < bins[x+1] where x =0,N and N is the
    length of the bins array.  The last range is given by
    bins[N] <= range_N < infinity.  Values less than bins[0] are
    not included in the histogram.

    Arguments:
       iterator -- The iterator.
       bins -- 1D array.  Defines the ranges of values to use during
       histogramming.

    Returns:
    1D array.  Each value represents the occurences for a given
    bin (range) of values.
    """

    h = numpy.zeros(len(bins))

    # bins might not be uniform, so do a bisect
    for value in iterator:
        i = bisect.bisect_left(value)
        h[i] += 1

    return h


def fillHistograms(infile, columns, bins):
    """fill several histograms from several columns in a file.

    The histograms are built on the fly.

    Description:

    Count the number of times values from array a fall into
    numerical ranges defined by bins.  Range x is given by
    bins[x] <= range_x < bins[x+1] where x =0,N and N is the
    length of the bins array.  The last range is given by
    bins[N] <= range_N < infinity.  Values less than bins[0] are
    not included in the histogram.

    Arguments:
       file -- The input file.
       columns -- columns to use
       bins -- a list of 1D arrays.  Defines the ranges of values to use during
       histogramming.

    Returns:
    a list of 1D arrays.  Each value represents the occurences for a given
    bin (range) of values.

    WARNING: missing value in columns are ignored
    """

    assert(len(bins) == len(columns))

    hh = [numpy.zeros(len(bins[x]), numpy.float) for x in range(len(columns))]

    for line in infile:
        if line[0] == "#":
            continue
        data = line[:-1].split()
        # bins might not be uniform, so do a bisect
        for x, y in enumerate(columns):
            try:
                v = float(data[y])
            except IndexError:
                continue
            i = bisect.bisect(bins[x], v) - 1
            if i >= 0:
                hh[x][i] += 1

    return hh
