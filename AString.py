################################################################################
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
#################################################################################
'''
AString.py - a compact string of characters
===========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import array

##------------------------------------------------------------
class AString:
    """an array posing as a sequence.

    This class conserves memory as it uses only 1 byte per letter,
    while python strings use the machine word size for a letter.

    It exports a mixture of the methods in the python string and 
    python array classes.

    .. note::

       Using this class will incur a heavy penalty compared to
       using :class:array.array directly. Only use sparingly for
       heavy computations.
    """
    
    def __init__(self, *args):
        self.mArray = array.array("c", *args )

    def __getitem__(self, *args ):
        return self.mArray.__getitem__( *args )

    def upper( self ):
        """return upper case version."""
        return AString( self.mArray.tostring().upper() )
    
    def lower(self):
        """return lower case version."""
        return AString( self.mArray.tostring().lower() )

    def insert( self, *args ):
        self.mArray.insert( *args )
        
    def index( self, *args ):
        return self.mArray.index( *args )

    def count( self, *args ):
        return self.mArray.count( *args )

    def reverse( self, *args ):
        self.mArray.reverse( *args )

    def extend( self, *args ):
        self.mArray.extend( *args )

    def remove( self, *args ):
        self.mArray.remove( *args )

    def fromstring( self, *args ):
        self.mArray.fromstring( *args )

    def tostring( self, *args ):
        return self.mArray.tostring( *args )

    def __getslice__( self, *args):
        """return slice as a string."""
        return self.mArray.__getslice__(*args).tostring()

    def __setslice__( self, start, end, sub):
        """set slice start:end from a string sub."""
        return self.mArray.__setslice__(start, end,
                                        array.array("c", sub ))

    def __str__(self):
        return self.mArray.tostring()

    def __getattr__(self, name):
        return getattr( self.mArray, name)
