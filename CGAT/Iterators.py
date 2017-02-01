"""
Iterators.py - Iterator functions
=================================

A collection of general purpose iterators.

"""

import random


def sample(iterable, sample_size=None):
    """sample # copies from iterator without replacement.

    Stores a temporary copy of the items in iterable. The function has
    thus a possibly high memory footprint and long pre-processing time
    to yield the first element.

    If sample_size is not given, the iterator returns elements in
    random order (see random.shuffle())

    .. note::
        This snippet was downloaded from an unknown source.

    """

    saved = []
    for element in iterable:
        saved.append(element)

    random.shuffle(saved)
    if not sample_size:
        sample_size = len(saved)
    for x in range(sample_size):
        yield saved[x]


def group_by_distance(iterable, distance=1):
    """group integers into non-overlapping intervals that
    are at most *distance* apart.

    >>> list( group_by_distance( (1,1,2,4,5,7) ) )
    [(1, 3), (4, 6), (7, 8)]

    >>> list( group_by_distance( [] ) )
    []

    >>> list( group_by_distance( [3] ) )
    [(3, 4)]

    >>> list( group_by_distance( [3,2] ) )
    Traceback (most recent call last):
    ...
    ValueError: iterable is not sorted: 2 < 3

    .. note:: 
        This snippet was downloaded from an unknown source.

    """
    i = iter(iterable)
    end = None
    start = end = cur = next(i)

    for cur in i:
        if cur < end:
            raise ValueError("iterable is not sorted: %i < %i" % (cur, end))
        if cur - end > distance:
            yield (start, end + 1)
            start = cur
        end = cur
    yield (start, end + 1)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
