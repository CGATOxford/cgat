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
PipelineTracks.py - Track definition in pipelines
=================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline processes the data streams several experimental
data sources. Building a consistent pipeline is difficult
as the data-streams frequently cross. Each project brings with it different 
ways of combining experimental data files. For example, there might be

   * several biological or technical replicates, which need to be grouped by experiment.
   * several experiments in different tissues.
   * several experiments in different conditions, which need to be compared against a reference 
      condition, for example stimulated versus unstimulated.
   * control experiments which provide base lines, for example the ``input``
      tracks in ChIP-Seq experiments. These need to be matched with the correct experiment.

Furthermore, tracks need to be encoded for use in different places. For example,
they might be used as part of file names (no space), as SQL table names (no ``-+.``)
or within R (no ``_``).

This module provides some tools to map tracks to different representations and to
group them in flexible ways in order to provide convenient short-cuts in pipelines.

Code
----

'''

import sys, re, string, glob, collections, copy

# '-' as separator
FILE_SEPARATOR="-"

# no '_'
R_SEPARATOR="."

# no '-', use '_'
TABLE_SEPARATOR="_"

AGGREGATE_PLACEHOLDER="agg"

def to_aggregate( x ):
    if x: return x
    else: return AGGREGATE_PLACEHOLDER

def from_aggregate( x ):
    if x == AGGREGATE_PLACEHOLDER: return None
    else: return x

# Hacky, think about improvements:
# 1. use a named-tuple factory style approach?
# 2. read-only
class Sample(object):
    '''a sample/track.
    '''

    attributes = ( "experiment", )
    representation = "file"

    def __init__(self, filename = None ):
        collections.namedtuple.__init__(self)
        self.data = collections.OrderedDict( zip( self.attributes, [None] * len(self.attributes) ) )
        if filename: self.fromFile( filename )

    def clone( self ):
        return copy.deepcopy( self )

    def asFile( self ):
        return FILE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )
    
    def asTable( self):
        return TABLE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )

    def asR( self):
        return R_SEPARATOR.join( map( to_aggregate, self.data.values() ) )

    def split( self, s, sep ):
        d = map( from_aggregate, s.split(sep) )
        assert len(d) == len(self.attributes)
        self.data = collections.OrderedDict( zip( self.attributes, d ))
        
    def fromFile( self, fn ):
        self.split( fn, FILE_SEPARATOR )

    def fromTable( self, fn ):
        self.split( fn, TABLE_SEPARATOR )

    def fromR( self, fn ):
        self.split( fn, R_SEPARATOR )
    
    def isReplicate( self ):
        '''return true if track/sample is a replicate of an experiment.'''
        return self.replicate != None

    def asAggregate( self, *args ):
        '''return a new aggregate Sample.'''
        n = copy.deepcopy( self )
        for x in self.attributes: 
            if x not in args:
                n.data[x] = None
        return n

    def toLabels( self ):
        '''return attributes that this track is an aggregate of.'''
        return [ x for x in self.attributes if self.data[x] != None ]
        
    def __str__(self ):
        return self.__repr__()

    def __repr__(self):
        return self.asFile()

    def __eq__(self, other):
        return str(other) == str(self)
    
    def __hash__(self):
        '''return hash value.'''
        return hash(self.asFile())

    def __getattr__(self, key ):
        if key in self.attributes:
            return object.__getattribute__(self, "data")[key]
        else:
            return object.__getattribute__(self, key )

    def __setattr__(self, key, val ):
        if key in self.attributes:
            object.__getattribute__(self, "data" )[key] = val
        else:
            object.__setattr__( self, key, val)

    @classmethod
    def setDefault( cls, representation = None ):
        '''set default representation for tracks to *representation*.
        If *represenation* is None, the representation will be set to
        the library default (asFile()).
        '''
        if representation == None or representation == "asFile":
            cls.__repr__ = cls.asFile
        elif representation == "asTable":
            cls.__repr__ = cls.asTable
        elif representation == "asR":
            cls.__repr__ = cls.asR

class Sample3(Sample):
    '''a sample/track.
    '''
    attributes = ( "tissue", "condition", "replicate")

class Aggregate:
    
    def __init__(self, tracks, 
                 track = None,
                 labels = None,
                 ):
        '''an aggregate of samples. 

        Aggregate tracks.

        *label* aggregate via labels

        *track* aggregate via an aggregate track

        If neither *track* nor *labels* is given,
        all tracks will be aggregated.

        '''
        self.factory = tracks.factory

        if labels:
            self.aggregates = labels
        elif track:
            self.aggregates = track.toLabels()
        else:
            # aggregate all
            self.aggregates = []
            
        self.track2groups = {}
        for track in tracks:
            key = track.asAggregate( *self.aggregates )
            if key not in self.track2groups:
                self.track2groups[key] = []
            self.track2groups[key].append( track )

    def getTracks( self, pattern = None ):
        '''return all tracks within this aggregate.'''
        r = sum([x for x in self.track2groups.values()], [] )
        if pattern:
            return [ pattern % x for x in r ]
        else:
            return r

    def __str__(self):
        x = []
        for key, v in self.track2groups.iteritems():
            x.append( "%s\t%s" % (str(key), "\t".join( map(str, v))) )
        return "\n".join(x)

    def __getitem__(self, key):
        return self.track2groups[key]
    
    def __iter__(self):
        return self.track2groups.keys().__iter__()
    
    def __len__(self):
        return len(self.track2groups)

    def keys(self):
        return self.track2groups.keys()
    
    def iteritems( self ):
        return self.track2groups.iteritems()

class Tracks:
    '''a collection of samples/tracks.'''
    factory = Sample
    
    def __init__(self, factory = Sample): 
        self.factory = factory
        self.tracks = []

    def loadFromDirectory( self, files, pattern, exclude = None ):
        '''load tracks from a list of files, applying pattern.

        If set, exclude files matching patterns in *exclude*.
        '''
        tracks = []
        rx = re.compile(pattern)
        
        if exclude: to_exclude = [ re.compile(x) for x in exclude]

        for f in files:
            if exclude:
                skip = False
                for x in to_exclude:
                    if x.search( f ): 
                        skip = True
                        break
                if skip: continue

            tracks.append( self.factory( filename = rx.search( f ).groups()[0] ) )

        self.tracks = tracks
        return self

    def __iter__(self):
        return self.tracks.__iter__()

    def __len__(self):
        return len(self.tracks)

    def __iadd__(self, other ):
        assert self.factory == other.factory
        self.tracks.extend( other.tracks )
        return self

    def __add__(self, other ):
        assert self.factory == other.factory
        n = copy.deepcopy( self )
        n.tracks.extend( other.tracks )
        return n

    def __contains__(self, key ):
        # do parsing (i.e. filenames versus tablenames?)
        return key in self.tracks

    def getTracks( self, pattern = None ):
        '''return all tracks in set. '''
        if pattern:
            return [ pattern % x for x in self.tracks ]
        else:
            return self.tracks

def getSamplesInTrack( track,tracks ):
    '''return all tracks in *tracks* that constitute *track*.'''
    return Aggregate( tracks, track = track )[track]
