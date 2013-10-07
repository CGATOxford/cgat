'''
PipelineTracks.py - Definition of tracks in pipelines
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Motivation
----------

A pipeline typically processes the data streams from several experimental
data sources. These data streams are usually processed separately (processing,
quality control) and as aggregates. For example, consider the following 
experimental layout:

+------------------------------+--------------------------------------------------+
|*Filename*                    |*Content*                                         |
+------------------------------+--------------------------------------------------+
|liver-stimulated-R1           |liver, stimulated, replicate 1                    |
+------------------------------+--------------------------------------------------+
|liver-stimulated-R2           |liver, stimulated, replicate 2                    |
+------------------------------+--------------------------------------------------+
|liver-unstimulated-R1         |liver, unstimulated, replicate 1                  |
+------------------------------+--------------------------------------------------+
|liver-unstimulated-R2         |liver, unstimulated, replicate 2                  |
+------------------------------+--------------------------------------------------+
|heart-stimulated-R1           |heart, stimulated, replicate 1                    |
+------------------------------+--------------------------------------------------+
|heart-stimulated-R2           |heart, stimulated, replicate 2                    |
+------------------------------+--------------------------------------------------+
|heart-unstimulated-R1         |heart, unstimulated, replicate 1                  |
+------------------------------+--------------------------------------------------+
|heart-unstimulated-R2         |heart, unstimulated, replicate 2                  |
+------------------------------+--------------------------------------------------+

The experiment measured in two tissues with two conditions with two replicates each
giving eight data streams. During the analysis, the streams are merged in a variety 
of combinations:

    * unmerged for initial processing, QC, etc.
    * by replicates to assess reproducibility of measurements
    * by condition to assess the size of the response to the stimulus
    * by tissue to assess differences between tissue and address the biological question.

The crossing of data streams complicates the building of pipelines, especially
as no two experiments are the same. The :mod:`PipelineTracks` module assists
in controlling these data streams. This module provides some tools to map tracks to different 
representations and to group them in flexible ways in order to provide convenient short-cuts 
in pipelines.

There are three class within :mod:`PipelineTracks`: :class:`Sample`, :class:`Tracks` and :class:`Aggregate`.

A Track
+++++++

The basic atomic data structure is a :class:`Sample` or :term:`track`. A :term:`track` is a 
single measurement that can be combined with other tracks. A track identifier consists of a tuple of
attributes. Each track in and experimental design has the same number of labels in the same order. 
In the example above, there are three attributes: tissue, condition and replicate. Identifiers are thus 
``('liver', 'stimulated','R1')`` or ``('heart','unstimulated','R2')``.

The same track can be represented by different names depending on context, for example when 
it is used as a filename or a database table. As filename, the track ``('heart','unstimulated','R2')``
is rendered as ``heart-unstimulated-R2`` (avoiding spaces), while as a table, it reads
``heart-unstimulated-R2``, avoiding ``-+.``. The :class:`Sample` class provides convenience methods to 
convert names from  one context to another. 

Track containers
++++++++++++++++

A container of type :class:`Tracks` stores one or more objects of type :class:`Sample`.

Aggregates
++++++++++

Tracks can be combined into aggregates. Aggregation is indicated by the ``agg`` keyword.

For example, the ``liver-stimulated-agg`` aggregate combines the tracks ``liver-stimulated-R1``
and ``liver-stimulated-R2``. The aggregate ``agg-stimulated-agg`` combines all replicates and
all tissues (``liver-stimulated-R1``, ``liver-stimulated-R2``, ``heart-stimulated-R1``, 
``heart-stimulated-R2``)

Usage
-----

Defining tracks and aggregates
++++++++++++++++++++++++++++++

To use tracks, you need to first define a new :class:`Sample`. In the example above with the attributes 
tissue, condition and replicate, the :class:`Sample` could be::

   import PipelineTracks
   
   class MySample( PipelineTracks.Sample ):
        attributes = ( "tissue", "condition", "replicate" )

Once defined, you can add tracks to a :class:`tracks` container. For example::

   TRACKS = PipelineTracks.Tracks( MySample ).loadFromDirectory( glob.glob( "*.fastq.gz" ), 
                                                                 pattern = "(\S+).fastq.gz" )

will collect all files ending in ``.fastq.gz``. The track identifiers will be derived by removing the ``fastq.gz``
suffix. The variable ``TRACKS`` contains all the tracks derived from files ending in ``*.fastq.gz``::

   >>> print TRACKS
   [liver-stimulated-R2, heart-stimulated-R2, liver-stimulated-R1, liver-unstimulated-R1, heart-unstimulated-R2, heart-stimulated-R1, heart-unstimulated-R1, liver-unstimulated-R2]

To build aggregates, use :class:`PipelineTracks.Aggregate`. The following combines replicates for each experiment::

   EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )

Aggregates are simply containers of associated data sets. To get a list of experiments, type::

   >>> EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
   >>> print list(EXPERIMENT)
   [heart-stimulated-agg, heart-unstimulated-agg, liver-stimulated-agg, liver-unstimulated-agg]

or::

   >>> print EXPERIMENT.keys()
   [heart-stimulated-agg, heart-unstimulated-agg, liver-stimulated-agg, liver-unstimulated-agg]

To obtain all replicates in the experiment ``heart-stimulated``, use dictionary access::

   >>> print EXPERIMENTS['heart-stimulated-agg']
   [heart-stimulated-R2, heart-stimulated-R1]
   
The returned objects are tracks. To use a :term:`track` as a tablename or as a file, use data
access functions :meth:`Sample.asTable` or :meth:`Sample.asFile`, respectively::

   >>> print [x.asFile() for x in EXPERIMENTS['heart-stimulated-agg'] ]
   ['heart-stimulated-R2', 'heart-stimulated-R1']

   >>> print [str(x) for x in EXPERIMENTS['heart-stimulated-agg'] ]
   ['heart-stimulated-R2', 'heart-stimulated-R1']

   >>> print [x.asTable() for x in EXPERIMENTS['heart-stimulated-agg'] ]
   ['heart_stimulated_R2', 'heart_stimulated_R1']

Note how the ``-`` is converted to ``_`` as the former are illegal as SQL table names.

The default representation is file-based. By using the class method::

   MySample.setDefault( "asTable" )

the default representation can be changed for all tracks simultaneously.

You can have multiple aggregates. For example, some tasks might require all conditions or all 
tissues::

   CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
   TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

You can have several :class:`Tracks` within a directory. :class:`Tracks` are simply
containers and as such do not have any actions associated with them.

Using tracks in pipelines
+++++++++++++++++++++++++

Unfortunately, tracks and aggregates do not work yet directly as ruffus_ 
task lists. Instead, they need to be converted to files explicitely using 
list comprehensions.

If you wanted to process all tracks separately, use::

   @files( [ ("%s.fastq.gz" % x.asFile(),
               "%s.qc" % x.asFile()) for x in TRACKS ] )
   def performQC( infile, outfile ):
      ....

The above statement will create the following list of input/output files for the ``performQC`` task::

   [ ( "liver-stimulated-R1.fastq.gz", "liver-stimulated-R1.qc" )
     ( "liver-stimulated-R2.fastq.gz" , "liver-stimulated-R2.qc" ),
     ...
   ]

Using aggregates works similarly, though you will need to create the file
lists yourself using nested list comprehensions. The following creates
an analysis per experimemnt::

   @files( [( ([ "%s.fastq.gz" % y.asFile() for y in EXPERIMENTS[x]]), 
                     "%s.out" % x.asFile()) 
                     for x in EXPERIMENTS ] )
   def checkReproducibility( infiles, outfile ):
      ....

The above statement will create the following list of input/output files::

   [ ( ( "liver-stimulated-R1.fastq.gz", "liver-stimulated-R2.fastq.gz" ), "liver-stimulated-agg.out" ),
     ( ( "liver-unstimulated-R1.fastq.gz", "liver-unstimulated-R2.fastq.gz" ), "liver-unstimulated-agg.out" ),
     ( ( "heart-stimulated-R1.fastq.gz", "heart-stimulated-R2.fastq.gz" ), "heart-stimulated-agg.out" ),
     ( ( "heart-unstimulated-R1.fastq.gz", "heart-unstimulated-R2.fastq.gz" ), "heart-unstimulated-agg.out" ),
   ]

The above code makes sure that the file dependencies are observed. Thus, if ``heart-stimulated-R1.fastq.gz``
changes, only ``heart-stimulated-agg.out`` will be re-computed.

Tracks and aggregates can be used within a task. The following code will collect all replicates for 
the experiment ``liver-stimulated-agg``

    >>> track = TRACKS.factory( filename = "liver-stimulated-agg" )
    >>> replicates = PipelineTracks.getSamplesInTrack( track, TRACKS )
    >>> print replicates
    [liver-stimulated-R2, liver-stimulated-R1]

API
----

.. _ruffus: http://www.ruffus.org.uk/


'''

import sys
import re
import string
import glob
import collections
import copy

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
    '''a sample/track with one attribute called ``experiment``.
    '''

    attributes = ( "experiment", )
    representation = "file"

    def __init__(self, filename = None ):
        '''create a new Sample.

        If filename is given, the sample name will be derived from *filename*.
        '''
        
        collections.namedtuple.__init__(self)
        self.data = collections.OrderedDict( zip( self.attributes, [None] * len(self.attributes) ) )
        if filename: self.fromFile( filename )

    def clone( self ):
        '''return a copy of self.'''
        return copy.deepcopy( self )

    def asFile( self ):
        '''return sample as a filename'''
        return FILE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )
    
    def asTable( self):
        '''return sample as a tablename'''
        return TABLE_SEPARATOR.join( map( to_aggregate, self.data.values() ) )

    def asR( self):
        '''return sample as valid R label'''
        return R_SEPARATOR.join( map( to_aggregate, self.data.values() ) )

    def fromFile( self, fn ):
        '''build sample from filename *fn*'''
        self._split( fn, FILE_SEPARATOR )

    def fromTable( self, tn ):
        '''build sample from tablename *tn*'''
        self._split( tn, TABLE_SEPARATOR )

    def fromR( self, rn ):
        '''build sample from R name *rn*'''
        self._split( rn, R_SEPARATOR )

    def _split( self, s, sep ):
        if len(self.attributes) == 1:
            self.data = collections.OrderedDict( ((self.attributes,s),) )
        else:
            d = map( from_aggregate, s.split(sep) )
            if len(d) != len(self.attributes):
                raise ValueError("can not match %s (sep='%s') against attributes %s" % (s, sep, self.attributes))
            self.data = collections.OrderedDict( zip( self.attributes, d ))
    
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
    '''a sample/track with three attributes: tissue, condition and replicate.
    '''
    attributes = ( "tissue", "condition", "replicate" )

class Aggregate:
    
    def __init__(self, 
                 tracks, 
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
    '''a collection of tracks.'''
    factory = Sample
    
    def __init__(self, factory = Sample): 
        '''create a new container. 

        New tracks are derived using *factory*.
        '''

        self.factory = factory
        self.tracks = []

    def loadFromDirectory( self, files, pattern, exclude = None ):
        '''load tracks from a list of files, applying pattern.

        Pattern is a regular expression with at at least one 
        group, for example ``(.*).gz``.

        If set, exclude files matching regular expression in *exclude*.
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
        '''return true if *key* is in tracks.'''
        # do parsing (i.e. filenames versus tablenames?)
        return key in self.tracks

    def getTracks( self, pattern = None ):
        '''return all tracks in container.
        '''
        if pattern:
            return [ pattern % x for x in self.tracks ]
        else:
            return self.tracks

def getSamplesInTrack( track, tracks ):
    '''return all tracks in *tracks* that constitute *track*.'''
    return Aggregate( tracks, track = track )[track]
