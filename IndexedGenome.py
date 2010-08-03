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
IndexedGenome.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import ncl

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
