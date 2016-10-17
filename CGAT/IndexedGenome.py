'''IndexedGenome.py - Random access to interval lists
==================================================

This module provides a consistent front-end to various interval containers.

Two implementations are available:

NCL
   Nested containment lists as described in
   http://bioinformatics.oxfordjournals.org/content/23/11/1386.short. The
   implemenation was taken from `pygr
   <http://code.google.com/p/pygr>`_.

quicksect
   Quicksect algorithm used in Galaxy, see `here
   <https://github.com/brentp/quicksect>`_.  This requires python.bx
   to be installed. The benefit of quicksect is that it allows also
   quick retrieval of intervals that are closest before or after an query.

The principal clas is :class:`IndexedGenome` which uses NCL and stores
a value associated with each interval. :class:`Quicksect` is equivalent
to :class:`IndexedGenome` but uses quicksect. The :class:`Simple` is a
light-weight version of :class:`IndexedGenome` that does not store a
value and thus preserves space.

The basic usage is::

   from IndexedGenome import IndexedGenome
   index = IndexedGenome()
   for contig, start, end, value in intervals:
      index.add(contig, start, end, value)

   print index.contains("chr1", 1000, 2000)
   print index.get("chr1", 10000, 20000)

The index is built in memory.

Reference
---------

'''
from CGAT import NCL as ncl
from bx.intervals.intersection import Intersecter, Interval


class IndexedGenome:

    '''Genome with indexed intervals.
    '''

    index_factory = ncl.NCL

    def __init__(self):
        self.mIndex = {}

    def add(self, contig, start, end, value):

        if contig not in self.mIndex:
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add(start, end, value)

    def __getitem__(self, args):
        '''return intervals overlapping with key.'''
        if args[0] not in self.mIndex:
            raise KeyError("contig %s not in index" % args[0])

        return self.mIndex[args[0]].find(args[1], args[2])

    def contains(self, contig, start, end):
        if contig not in self.mIndex:
            return False
        return len(list(self.mIndex[contig].find(start, end))) > 0

    def get(self, contig, start, end):
        '''return intervals overlapping with key.'''
        if contig not in self.mIndex:
            raise KeyError("contig %s not in index" % contig)

        return self.mIndex[contig].find(start, end)

    def __len__(self):
        '''return number of contigs.'''
        return len(self.mIndex)


class Simple(IndexedGenome):

    '''index intervals without storing a value.'''
    index_factory = ncl.NCLSimple

    def __init__(self, *args, **kwargs):
        IndexedGenome.__init__(self, *args, **kwargs)

    def add(self, contig, start, end):

        if contig not in self.mIndex:
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add(start, end)


class Quicksect(IndexedGenome):

    '''index intervals using quicksect.

    Permits finding closest interval in case there is
    no overlap.
    '''
    index_factory = Intersecter

    def __init__(self, *args, **kwargs):
        IndexedGenome.__init__(self, *args, **kwargs)

    def add(self, contig, start, end, value):

        if contig not in self.mIndex:
            self.mIndex[contig] = self.index_factory()
        self.mIndex[contig].add_interval(Interval(start, end, value))

    def get(self, contig, start, end):
        '''return intervals overlapping with key.'''
        if contig not in self.mIndex:
            raise KeyError("contig %s not in index" % contig)

        return [(x.start, x.end, x.value)
                for x in self.mIndex[contig].find(start, end)]

    def before(self, contig, start, end, num_intervals=1, max_dist=2500):
        '''get closest interval before *start*.'''
        if contig not in self.mIndex:
            raise KeyError("contig %s not in index" % contig)
        return [(x.start, x.end, x.value)
                for x in self.mIndex[contig].before_interval(
                    Interval(start, end),
                    num_intervals=1,
                    max_dist=max_dist)]

    def after(self, contig, start, end, num_intervals=1, max_dist=2500):
        '''get closest interval after *end*.'''
        if contig not in self.mIndex:
            raise KeyError("contig %s not in index" % contig)
        return [(x.start, x.end, x.value)
                for x in self.mIndex[contig].after_interval(
                    Interval(start, end),
                    num_intervals=1,
                    max_dist=max_dist)]
