'''
AString.py - strings as arrays of characters
============================================

This module provides the :class:`AString` class to efficiently
represent long, chromosomal nucleotide sequences in memory.

Reference
---------

'''
from array import array


import sys

IS_PY3 = sys.version_info.major >= 3


class AString(array):
    """implementation of a string as an array.

    This class conserves memory as it uses only 1 byte per letter,
    while python strings use the machine word size for a letter.

    It adds a subset of the python string class such as upper() and
    lower() for convenience. Slicing and printing return strings.

    The :class:`AString` can be constructed by any iterable that is
    accepted by the constructor of :py:class:`array.array`.

    """

    def __new__(cls, *args):
        if IS_PY3:
            return array.__new__(cls, "b", *[x.encode("ascii") for x in args])
        else:
            return array.__new__(cls, "b", *args)

    def upper(self):
        """return upper case version."""
        if IS_PY3:
            return AString(str(self).upper().encode('ascii'))
        else:
            return AString(str(self).upper())

    def lower(self):
        """return lower case version."""
        if IS_PY3:
            return AString(str(self).lower().encode('ascii'))
        else:
            return AString(str(self).lower())

    def __getitem__(self, *args):
        """return slice as a string."""

        if IS_PY3:
            return array.__getitem__(self, *args).tostring().decode("ascii")
        else:
            return array.__getitem__(self, *args).tostring()

    def __setitem__(self, start, end, sub):
        """set slice start:end from a string sub."""
        return array.__setitem__(self,
                                 start, end,
                                 array("b", sub))

    def __str__(self):
        if IS_PY3:
            return self.tostring().decode("ascii")
        else:
            return self.tostring()
