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
IndexedGenome.py - Wrappers for several interval indices
=========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import NCL as ncl
from bx.intervals.intersection import Intersecter, Interval

class IndexedGenome:
    '''Genome with indexed intervals.
    '''

    index_factory = ncl.NCL
    def __init__( self ): 
        self.mIndex = {}

    def add( self, contig, start, end, value ):
        
        if contig not in self.mIndex: 
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add( start,end,value )

    def __getitem__(self, args):
        '''return intervals overlapping with key.'''
        if args[0] not in self.mIndex: 
            raise KeyError("contig %s not in index" % args[0] )
        
        return self.mIndex[args[0]].find( args[1],args[2] )

    def contains( self, contig, start, end ):
        if contig not in self.mIndex: return False
        return len(list(self.mIndex[contig].find( start, end))) > 0

    def get(self, contig, start, end):
        '''return intervals overlapping with key.'''
        if contig not in self.mIndex: 
            raise KeyError("contig %s not in index" % contig )
        
        return self.mIndex[contig].find( start, end )

    def __len__(self):
        '''return number of contigs.'''
        return len(self.mIndex)

class Simple( IndexedGenome ):
    '''index intervals without storing a value.'''
    index_factory = ncl.NCLSimple
    
    def __init__(self, *args, **kwargs ):
        IndexedGenome.__init__(self, *args, **kwargs )

    def add( self, contig, start, end ):
        
        if contig not in self.mIndex: 
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add( start,end )


class Quicksect( IndexedGenome ):
    '''index intervals using quicksect. 
    
    Permits finding closest interval in case there is 
    no overlap.
    '''
    index_factory = Intersecter
    
    def __init__(self, *args, **kwargs ):
        IndexedGenome.__init__(self, *args, **kwargs )

    def add( self, contig, start, end, value ):
        
        if contig not in self.mIndex: 
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add_interval( Interval( start,end, value ) )
        
    def get(self, contig, start, end ):
        '''return intervals overlapping with key.'''
        if contig not in self.mIndex: 
            raise KeyError("contig %s not in index" % contig )
        
        return [(x.start, x.end, x.value) for x in self.mIndex[contig].find( start, end ) ]
    
    def before( self, contig, start, end, num_intervals = 1, max_dist = 2500 ):
        '''get closest interval before *start*.''' 
        if contig not in self.mIndex: 
            raise KeyError("contig %s not in index" % contig )
        return [(x.start, x.end, x.value) for x in self.mIndex[contig].before_interval( Interval( start,end), 
                                                                                        num_intervals = 1,
                                                                                        max_dist = max_dist ) ]


    def after( self, contig, start, end, num_intervals = 1, max_dist = 2500 ):
        '''get closest interval after *end*.''' 
        if contig not in self.mIndex: 
            raise KeyError("contig %s not in index" % contig )
        return [(x.start, x.end, x.value) for x in self.mIndex[contig].after_interval( Interval( start,end), 
                                                                                       num_intervals = 1,
                                                                                       max_dist = max_dist ) ]


