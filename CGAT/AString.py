##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
AString.py - strings as arrays of characters
============================================

This module provides the :class:`AString` class to efficiently
represent long, chromosomal nucleotide sequences in memory.

Reference
---------

'''
from array import array


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
        return array.__new__(cls, "c", *args)

    def upper(self):
        """return upper case version."""
        return AString(self.tostring().upper())

    def lower(self):
        """return lower case version."""
        return AString(self.tostring().lower())

    def __getslice__(self, *args):
        """return slice as a string."""
        return array.__getslice__(self, *args).tostring()

    def __setslice__(self, start, end, sub):
        """set slice start:end from a string sub."""
        return array.__setslice__(self,
                                  start, end,
                                  array("c", sub))

    def __str__(self):
        return self.tostring()


