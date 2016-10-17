'''Intervals.py - Utility functions for working with intervals
===========================================================

This module contains utility functions for working intervals
or list of intervals.

An interval is a tuple of a start and end coordinate in python's
0-based, half-open notation such as::

    (12, 20)

An interval list is simply a list of such intervals.

The majority of the functions in this module take one or more lists of
intervals and return one or more new lists of intervals.

Reference
---------

'''


def getLength(intervals):
    """return sum of intervals lengths.

    >>> getLength([(10,20), (30,40)])
    20

    """
    return sum([x[1] - x[0] for x in intervals])


def combine(intervals):
    """combine overlapping and adjacent intervals.

    >>> combine([(10,20), (30,40)])
    [(10, 20), (30, 40)]
    >>> combine([(10,20), (20,40)])
    [(10, 40)]
    >>> combine([(10,20), (15,40)])
    [(10, 40)]
    """
    if not intervals:
        return []

    new_intervals = []

    intervals.sort()
    first_from, last_to = intervals[0]

    for this_from, this_to in intervals[1:]:
        if this_from > last_to:
            new_intervals.append((first_from, last_to))
            first_from, last_to = this_from, this_to
            continue

        if last_to < this_to:
            last_to = this_to

    new_intervals.append((first_from, last_to))

    return new_intervals


def prune(intervals, first=None, last=None):
    """truncates all intervals that are extending beyond first or last.

    Empty intervals are deleted.
    """
    new_intervals = []
    for start, end in intervals:
        if end <= first or start >= end:
            continue
        new_intervals.append((max(first, start), min(end, last)))
    return new_intervals


def complement(intervals, first=None, last=None):
    """complement a list of intervals with intervals not in list.

    >>> complement([(10,20), (15,40)])
    []
    >>> complement([(10,20), (30,40)])
    [(20, 30)]
    >>> complement([(10,20), (30,40)], first=5)
    [(5, 10), (20, 30)]

    Arguments
    ---------
    intervals : list
        List of intervals
    first : int
        First position. If given, the interval from `first` to
        the first position in `intervals` is added.
    last : int
        Last position. If given, the interval from the last position
        in `intervals` to `last` is added.

    Returns
    -------
    intervals : list
        A new list of intervals
    """

    if len(intervals) == 0:
        if first is not None and last is not None:
            return [(first, last)]
        else:
            return []

    new_intervals = []

    intervals.sort()
    last_from, last_to = intervals[0]

    if first is not None and first < last_from:
        new_intervals.append((first, last_from))

    for this_from, this_to in intervals:
        if this_from > last_to:
            new_intervals.append((last_to, this_from))

        last_from = this_from
        last_to = max(last_to, this_to)

    if last and last > last_to:
        new_intervals.append((last_to, last))

    return new_intervals


def addComplementIntervals(intervals, first=None, last=None):
    """complement a list of intervals with intervals not
    in list and return both.

    The resulting interval list is sorted.
    """
    return sorted(intervals + complement(intervals, first, last))


def joined_iterator(intervals1, intervals2):
    """iterate over the combination of two intervals.

    returns the truncated intervals delineating the
    ranges of overlap between intervals1 and intervals2.
    """

    intervals1.sort()
    intervals2.sort()

    x, y = 0, 0
    while x < len(intervals1) and y < len(intervals2):

        xfrom, xto = intervals1[x]
        yfrom, yto = intervals2[y]

        if xto <= yfrom:
            x += 1
        elif yto <= xfrom:
            y += 1
        else:
            mto = min(xto, yto)
            mfrom = max(xfrom, yfrom)

            yield (mfrom, mto)

            if xto < yto:
                x += 1
            elif yto < xto:
                y += 1
            else:
                x += 1
                y += 1


def intersect(intervals1, intervals2):
    """intersect two interval sets.

    Return a set of intervals that is spanned by intervals in
    both sets. Returns the union of the two intervals.
    """
    if not intervals1 or not intervals2:
        return []

    result = []
    for start, end in joined_iterator(intervals1, intervals2):

        overlap = end - start
        if overlap:
            result.append((start, end))

    return result


def truncate(intervals1, intervals2):
    """truncate intervals in intervals1 by intervals2

    Example: truncate( [(0,5)], [(0,3)] ) = [(3,5)]
    """
    if not intervals1:
        return []
    if not intervals2:
        return intervals1

    intervals2 = combine(intervals2)

    result = []
    intersection = intersect(intervals1, intervals2)

    if not intersection:
        return intervals1

    y = 0
    for start, end in intervals1:

        # retrieve all overlapping in intersection
        ovl = []
        while y < len(intersection):
            mask_start, mask_end = intersection[y]
            if mask_start >= end:
                break
            if mask_start >= start:
                ovl.append((mask_start, mask_end))
            y += 1

        # build truncated segments
        last = start
        for mask_start, mask_end in ovl:
            if mask_start - last > 0:
                result.append((last, mask_start))
            last = mask_end
        if end - last > 0:
            result.append((last, end))

    return result


def calculateOverlap(intervals1, intervals2):
    """calculate overlap between two list of intervals.

    The intervals within each set should not be overlapping.
    """

    if not intervals1 or not intervals2:
        return 0

    overlap = 0
    for start, end in joined_iterator(intervals1, intervals2):
        overlap += end - start

    return overlap


def fromArray(a):
    """get intervals from a binary array."""

    if len(a) == 0:
        return []

    was_a = a[0]
    if was_a:
        start = 0

    intervals = []
    for x in range(1, len(a)):
        is_a = a[x]
        if is_a and not was_a:
            start = x
            was_a = is_a
        elif not is_a and was_a:
            if start is not None:
                intervals.append((start, x))
            was_a = is_a
    if is_a:
        intervals.append((start, len(a)))

    return intervals


def combineAtDistance(intervals, min_distance):
    """combine a list intervals and merge those that are less than a
    certain distance apart.
    """

    if not intervals:
        return []

    # merge overlapping intervals
    intervals = combine(intervals)

    new_intervals = []
    intervals.sort()

    first_start, last_end = intervals[0]

    for this_start, this_end in intervals[1:]:

        if this_start - last_end >= min_distance:
            new_intervals.append((first_start, last_end))
            first_start = this_start

        last_end = this_end

    new_intervals.append((first_start, last_end))

    return new_intervals


def getIntersections(intervals):
    """return regions were two intervals are overlapping.

    """
    if not intervals:
        return []

    intervals.sort()
    max_to = intervals[0][1]

    all_sections = []
    sections = [intervals[0][0], intervals[0][1]]
    for this_from, this_to in intervals[1:]:

        # no overlap: write everything and reset
        if this_from > max_to:
            all_sections.append(sections)
            max_to = this_to
            sections = []

        max_to = max(max_to, this_to)

        sections.append(this_from)
        sections.append(this_to)

    all_sections.append(sections)

    new_intervals = []
    for sections in all_sections:
        sections.sort()
        last_x = sections[0]
        for x in sections[1:]:
            if last_x == x:
                continue
            new_intervals.append((last_x, x))
            last_x = x

    return new_intervals


def RemoveIntervalsContained(intervals):
    """remove intervals that are fully contained in another

    [(10, 100), (20, 50), (70, 120), (130, 200), (10, 50), (140, 210), (150, 200)]

    results:

    [(10, 100), (70, 120), (130, 200), (140, 210)]
    """
    if not intervals:
        return []

    new_intervals = []

    intervals.sort()
    last_from, last_to = intervals[0]

    for this_from, this_to in intervals[1:]:
        # this is larger:
        if this_from <= last_from and this_to >= last_to:
            last_from, last_to = this_from, this_to
            continue

        # last is larger
        if last_from <= this_from and last_to >= this_to:
            continue

        # no complete overlap
        new_intervals.append((last_from, last_to))

        last_from, last_to = this_from, this_to

    new_intervals.append((last_from, last_to))

    return new_intervals


def RemoveIntervalsSpanning(intervals):
    """remove intervals that are full covering
    another, i.e. always keep the smallest.

    [(10, 100), (20, 50), (70, 120), (40,80), (130, 200), (10, 50), (140, 210), (150, 200)]

    result:

    [(20, 50), (40, 80), (70, 120), (150, 200)]
    """

    if not intervals:
        return []

    intervals.sort()

    last_intervals = intervals

    while 1:

        new_intervals = []

        last_from, last_to = last_intervals[0]

        for this_from, this_to in last_intervals[1:]:
            # print last_from, last_to, this_from, this_to
            # this is larger:
            if this_from <= last_from and this_to >= last_to:
                continue

            # last is larger:
            if last_from <= this_from and last_to >= this_to:
                last_from, last_to = this_from, this_to
                continue

            # no complete overlap
            new_intervals.append((last_from, last_to))
            last_from, last_to = this_from, this_to

        new_intervals.append((last_from, last_to))

        if len(last_intervals) == len(new_intervals):
            break

        last_intervals = new_intervals

    return new_intervals


def ShortenIntervalsOverlap(intervals, to_remove):
    """shorten intervals, so that there is no
    overlap with another set of intervals.

    assumption: intervals are not overlapping

    """
    if not intervals:
        return []

    if not to_remove:
        return interalls

    new_intervals = []

    intervals.sort()
    to_remove.sort()

    current_to_remove = 0

    for this_from, this_to in intervals:

        for remove_from, remove_to in to_remove:
            # print this_from, this_to, remove_from, remove_to
            if remove_to < this_from:
                continue
            if remove_from > this_to:
                continue

            if remove_from <= this_from and remove_to >= this_to:
                this_from = remove_to
                break

            if this_from < remove_from:
                new_intervals.append((this_from, remove_from))
                # print "adding", this_from, remove_from

            this_from = max(this_from, remove_to)

            if this_to < this_from:
                break

        if this_to > this_from:
            # print "adding", this_from, this_to
            new_intervals.append((this_from, this_to))

    return new_intervals
