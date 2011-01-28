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
Pipeline.py - Tools for ruffus pipelines
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''

import sys, re, string, glob, collections, copy

FILE_SEPARATOR="-"
TABLE_SEPARATOR="_"
AGGREGATE_PLACEHOLDER="agg"

def to_aggregate( x ):
    if x: return x
    else: return AGGREGATE_PLACEHOLDER

class Sample:
    '''a sample/track.
    '''

    attributes = ( "condition", "tissue", "replicate")

    def __init__(self, filename = None ):
        self.data = collections.OrderedDict()
        if filename: self.fromFile( filename )

    def clone( self ):
        return copy.deepcopy( self )

    def asFile( self ):
        return FILE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )
    
    def asTable( self):
        return TABLE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )

    def split( self, s, sep ):
        d = s.split(sep)
        assert len(d) == len(self.attributes)
        self.data = collections.OrderedDict( zip( self.attributes, d ) )
        
    def fromFile( self, fn ):
        self.split( fn, FILE_SEPARATOR )

    def fromTable( self, fn ):
        self.split( fn, TABLE_SEPARATOR )
    
    def isReplicate( self ):
        '''return true if track/sample is a replicate of an experiment.'''
        return self.replicate != None

    def asAggregate( self, *args ):
        '''return a new aggregate Sample.'''
        assert len(args) > 0, "requiring at least one attribute to aggregate with"
        n = copy.deepcopy( self )
        for x in self.attributes: 
            if x not in args:
                n.data[x] = None
        return n

    def __str__(self ):
        return self.asFile()

    def __repr__(self ):
        return self.asFile()
    
    def __eq__(self, other):
        return str(other) == str(self)
    
    def __hash__(self):
        '''return hash value.'''
        return hash(self.asFile())

    def __getattr__(self, key ):
        return self.data[key]

class Aggregate:
    
    def __init__(self, tracks, *args ):
        '''an aggregate of samples. 

        Aggregate tracks via attributes in *args.
        '''
        self.factory = tracks.factory
        self.aggregates = args
        assert len(args) > 0, "requiring at least one attribute to aggregate with"
        self.track2groups = collections.defaultdict( list )
        for track in tracks:
            key = []
            for arg in args: key.append( getattr( track, arg ) )
            key = track.asAggregate( *args )
            self.track2groups[key].append( track )

    def getTracks( self, key ):
        return self.track2groups[key]

    def __str__(self):
        x = []
        for key, v in self.track2groups.iteritems():
            x.append( "%s\t%s" % (str(key), "\t".join( map(str, v))) )
        return "\n".join(x)
    
    def __iter__(self):
        return self.track2groups.keys().__iter__()

    def keys(self):
        return self.track2groups.keys()

class Tracks:
    '''a collection of samples/tracks.'''
    factory = Sample
    
    def __init__(self, factory = Sample): 
        self.factory = factory
        self.tracks = []

    def loadFromDirectory( self, files, pattern ):
        '''load tracks from a list of files, applying pattern.'''
        tracks = []
        rx = re.compile(pattern)

        for f in files:
            tracks.append( self.factory( filename = rx.search( f ).groups()[0] ) )

        self.tracks = tracks
        return self

    def __iter__(self):
        return self.tracks.__iter__()

def getReplicateGroups( tracks ):
    '''return a dictionary mapping experiment to replicates.'''
    t = [ splitTrack( x ) for x in tracks ]
    t.sort()
    result = dict( [ (x, list(y)) for x,y in itertools.groupby( t, key = lambda x: (x[0],x[1]) ) ] )
    return result
