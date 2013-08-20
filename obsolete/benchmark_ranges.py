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
"""

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

benchmarks various methods to do range queries.

The packages benchmarked are:

1. bx.python
2. rtree
3. ncl

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse, timeit, random, time
import CGAT.NCL as ncl

try:
    import rtree2
    HAVE_RTREE = True
except ImportError, msg:
    HAVE_RTREE = False

try:
    import bx.intervals.io
    import bx.intervals.intersection
    HAVE_BX = True
except ImportError, msg:
    HAVE_BX = False

try:
    from bme10.fast import *
    HAVE_BME = True
except ImportError, msg:
    HAVE_BME = False

try: 
    import quicksect
    HAVE_QUICKSECT = True
except ImportError, msg:
    HAVE_QUICKSECT = False

import CGAT.Experiment as E

def generate_uniform_segments( workspace_size, segment_size, num_segments ):
    """generate a uniformly covering set of segments without overlap

    If num_segments == 0, the segments are packed without gaps. Otherwise they are
    uniformly distributed.
    """

    if num_segments == 0:
        for x in xrange( 0, workspace_size, segment_size):
            yield (x, min(workspace_size,x+segment_size) )
    else:
        space = workspace_size - segment_size * num_segments
        increment = space // num_segments
        x = 0
        for y in xrange( num_segments):
            yield( x, min(workspace_size,x+segment_size) )
            x += segment_size + increment

def generate_overlapping_segments( workspace_size, segment_size, num_segments ):
    """generate a covering set of segments with overlap

    Generates a fully covering set of segments. A covering set of segments are built from segment_size.
    Later, segments of twice the size are added overlapping the previous set, and so on until
    the final segments covers the full workspace.
    """

    s = segment_size

    while s < workspace_size:
        for x in xrange( 0, workspace_size, s):
            yield (x, min(workspace_size,x+s) )
        s *= 2
        
    yield( 0, workspace_size )

def generate_random_segments( workspace_size, segment_size, num_segments ):
    """generate randomly positioned and sized segments."""
    
    for x in range(num_segments):
        start = random.randint( 0, workspace_size )
        end = random.randint( start, workspace_size )
        yield (start, end)

def sample_uniform_segments( workspace_size, segment_size ):
    """sample from a uniform set of covering segments without overlap
    """

    starts = list(range( 0, workspace_size, segment_size))
    random.shuffle( starts )
    
    for x in starts:
        yield (x, min(workspace_size,x+segment_size) )

def bx_create( segmenter ):
    
    intersector = bx.intervals.intersection.Intersecter()
    for start, end in segmenter:
        intersector.add_interval( bx.intervals.Interval(start,end) )
    return intersector

def bx_lookup( index, sampler ):
    for start, end in sampler:
        result = list(index.find( start, end ))
    return result

def ncl_create( segmenter ):
    index = ncl.IntervalDB()
    index.fromlist( [ (coords[0],coords[1],id) for id,coords in enumerate( segmenter ) ] )
    return index

def ncl_lookup( index, sampler ):
    for start, end in sampler:
        result = list(index.find_overlap( start, end ))
    return result

def nclpy_create( segmenter ):
    index = ncl.NCL()
    for id,coords in enumerate( segmenter ): index.add( coords[0],coords[1],id )
    return index

def nclpy_create_on_disk( segmenter ):
    index = ncl.NCL( filestem = "tmp", force=True )
    for id,coords in enumerate( segmenter ): index.add( coords[0],coords[1],id )
    return index

def nclpy_create_on_disk_and_flush( segmenter ):
    index = ncl.NCL( filestem = "tmp", force=True )
    for id,coords in enumerate( segmenter ): index.add( coords[0],coords[1],id )
    del index
    return ncl.NCL( filestem = "tmp" )

def nclpy_lookup( index, sampler ):
    for start, end in sampler:
        result = list(index.find( start, end ))
    return result

def bme_create( segmenter ):
    index = GFile( [ Interval( "1", start, end) for start, end in segmenter],
                   header = {CHROM:0,START:1,END:2} )
    return index

def bme_lookup( index, sampler ):
    for start, end in sampler:
        result = list(index.intersection( Interval( "1", start, end ) ))
    return result

def quicksect_create( segmenter ):

    IntervalNode = quicksect.IntervalNode
    Feature = quicksect.Feature
    x = segmenter.next()
    index = IntervalNode( Feature( x[0], x[1]) )
    for start, end in segmenter:
        index = index.insert( Feature( start, end) )

    return index

def quicksect_lookup( index, sampler ):

    for start, end in sampler:
        result = list(index.find( start, end ) )
    return result

class RTreeIndex:
    def __init__( self ):
        self.index = rtree.Rtree()

    def add( self, x, start, end ):
        self.index.add(x, (start, 0, end-1, 0))

    def find( self, start, end ):
        return self.index.intersection( (start, 0, end-1,0) )

def rtree_create( segmenter ):
    
    index = rtree.Rtree()
    x = 0
    for start, end in segmenter:
        x += 1
        index.add(x, (start, 0, end-1, 0) )

    return index

def rtree_lookup( index, sampler ):

    for start, end in sampler:
        result = list(index.intersection( (start, 0, end-1,0) ))
    return result

def rtree2_create( segmenter ):
    
    index = RTreeIndex()
    x = 0
    for start, end in segmenter:
        x += 1
        index.add(x, start, end )

    return index

def rtree2_lookup( index, sampler ):

    for start, end in sampler:
        result = list(index.find( start, end) )
    return result

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )    

    parser.add_option( "--segments", dest="segment_function", type="choice",
                       choices=("uniform", "overlapping", "random", ),
                      help="segment layout [default=%default]."  )

    parser.set_defaults(
        workspace_size = 10000000,
        segment_size = 1000,
        segment_sizes = (10, 100, 1000),
        num_segments = (10, 100, 1000, 10000),
        verify = False,
        segment_function = "uniform",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser )

    pairs = [ 
        ("ncl:lowlevel-api", ncl_create, ncl_lookup),
        ("ncl:memory", nclpy_create, nclpy_lookup), 
        ("ncl:disk:create", nclpy_create_on_disk, nclpy_lookup),
        # ("nclpy:disk:use", nclpy_create_on_disk_and_flush, nclpy_lookup), 
        ]

    if HAVE_BME:
        pairs.append( ("bme:memory", bme_create, bme_lookup), )

    if HAVE_RTREE:
        pairs.extend( [
                ("rtree", rtree_create, rtree_lookup),
                ("rtree wrapper", rtree2_create, rtree2_lookup), ] )

    if HAVE_BX:
        pairs.append( ("python.bx", bx_create, bx_lookup) )

    if HAVE_QUICKSECT:
        pairs.append( ("quicksect", quicksect_create, quicksect_lookup) )

    # pairs = ( ("quicksect", quicksect_create, quicksect_lookup), )

    # options.num_segments = [0]

    segment_function = { "uniform" : generate_uniform_segments,
                         "overlapping" : generate_overlapping_segments,
                         "random" : generate_random_segments,
                         }[options.segment_function]
    
    if options.verify:
        for segment_size in options.segment_sizes:
            segmenter = segment_function( options.workspace_size, segment_size, num_segments )
            sampler = sample_uniform_segments( options.workspace_size, segment_size )
        
            segments = list( segmenter)

            indices = []
            for name, create, lookup in pairs:
                indices.append( create( segments ) )

            for start,end in sampler:
                for index, v in zip( indices, pairs):
                    name, create, lookup = v
                    intervals = list(lookup( index, ((start,end),) ))
                    print name, intervals

    else:
        sys.stdout.write("name\tsize\tnsegs\ttcreate\ttlookup\n" )

        for num_segments in options.num_segments:

            for segment_size in options.segment_sizes:

                for name, create, lookup in pairs:
                    segmenter = segment_function( options.workspace_size, segment_size, num_segments )
                    sampler = sample_uniform_segments( options.workspace_size, segment_size )

                    t0 = time.time()
                    index = create(segmenter)
                    t1 = time.time()
                    lookup( index, sampler )
                    t2 = time.time()

                    sys.stdout.write( "%s\t%i\t%i\t%f\t%f\n" % (name, segment_size, num_segments, t1-t0, t2-t1 ) )
                    
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

