"""
Wiggle.py - Support for the `wiggle`_ format used by `ucsc`_.
=============================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

.. _wiggle: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _ucsc: ??
"""

# I had to re-implement the reader, so that random access would work.


import psyco_full
import numpy
from bx.misc.seekbzip2 import SeekableBzip2File
from bx import interval_index_file
import bx.wiggle



class WiggleIndexedAccess( interval_index_file.AbstractIndexedAccess ):
    """
    Indexed access to a wiggle file.
    """
    def get( self, src, start, end, dtype = numpy.float ):
        """return a list of numpy arrays.
        
        Each entry in the list is a tuple (pos, array),
        where pos refers to the residue number of the 
        first position in pos.

        dtype: numpy datatype

        """
        intersections = self.indexes.find( src, start, end ) 
        intersections.sort()
        result = []
        for istart, iend, ival in intersections:
            reader = self.get_at_offset( ival ) 
            xstart=max( istart,start)
            xend=min(end,iend)
            a = numpy.zeros( (min(end,iend)-xstart,), dtype )
            for chr, pos, value in reader:
                ## stay within current contiguous block
                if chr != src or pos >= xend: break
                if pos < start: continue
                a[pos-xstart] = value
            
            result.append( (xstart, a) )
        return result

    def read_at_current_offset( self, file, **kwargs ):
        """
        Read the MAF block at the current position in `file` and return an
        instance of `Alignment`.
        """
        return bx.wiggle.Reader( file, **kwargs )

class WiggleMultiIndexedAccess( interval_index_file.AbstractMultiIndexedAccess ):
    """
    Indexed access to multiple MAF files.
    """
    indexed_access_class = WiggleIndexedAccess

    def get( self, src, start, end, dtype = numpy.float ):
        blocks = []
        for index in self.indexes: 
            blocks.extend( index.get( src, start, end, dtype ) ) 
        return blocks 

