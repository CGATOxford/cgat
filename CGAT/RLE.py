"""RLE.py - a simple run length encoder
=======================================

:Tags: Python

Taken from: http://rosettacode.org/wiki/Run-length_encoding#Python
"""
from itertools import groupby
import array


def encode(input_array):
    """encode array or string.

    return tuples of (count, value).

    >>> encode(array.array( "i", (10,10,10,10,20,20,20,20) ) )
    [(4, 10), (4, 20)]

    >>> encode("aaaaahhhhhhmmmmmmmuiiiiiiiaaaaaa")
    [(5, 'a'), (6, 'h'), (7, 'm'), (1, 'u'), (7, 'i'), (6, 'a')]

    """
    return [(len(list(g)), k) for k, g in groupby(input_array)]


def decode(lst, typecode):
    """decode to array

    >>> decode( [(4, 10), (4, 20)], typecode="i" )
    array('i', [10, 10, 10, 10, 20, 20, 20, 20])

    >>> decode( [(5, 'a'), (6, 'h'), (7, 'm'), (1, 'u'), (7, 'i'), (6, 'a')], typecode="c" )
    array('c', 'aaaaahhhhhhmmmmmmmuiiiiiiiaaaaaa')

    """
    a = array.array(typecode)
    for n, c in lst:
        a.extend(array.array(typecode, (c,) * n))
    return a


def compress(input_string, bytes=1):
    """return compressed stream."""
    r = []
    for n, c in encode(input_array):
        while n > 255:
            n -= 255
            r.append(chr(255))
            r.append(ord(c))
        r.append(chr(n))
        r.append(ord(c))
    return "".join(r)


def uncompress(stream):

    n, r = None, []
    for c in stream:
        if n is None:
            n = c
        else:
            r.append(n, c)
            n = None
    return decode(r)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
