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
maq2psl.py - wonvert the output of maqview to psl format
========================================================

:Author: Andreas Heger
:Release: $Id: maq2psl.py 2781 2009-09-10 11:33:14Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Convert the output of maqview to psl format.

Note: the script assumes that the coordinates and the matches are sorted
correctly. Careful, as "seg1" < "seg10" < "seg9".

Usage
-----

Example::

   python maq2psl.py --help

Type::

   python maq2psl.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.Blat as Blat
import CGAT.Maq as Maq
import bisect
import itertools

class Segment: 
    def __init__(self, segment, pos, contig, lstart, rstart ):
        self.mSegment = segment
        self.start = pos
        self.contig = contig
        self.mLeftStart = lstart
        self.mRightStart = rstart
    def __str__(self):
        return "\t".join( map(str, ( self.mSegment,
                                     self.start,
                                     self.contig,
                                     self.mLeftStart,
                                     self.mRightStart )))

def iterator_segments( infile, segment_length ):
    for line in infile: 
        if line[0] == "#": continue
        if line.startswith( "segment\tpos" ): continue
        segment, pos, contig, lstart, rstart = line[:-1].split("\t")
        yield Segment( segment, int(pos), contig, int(lstart), int(rstart) )

def match_smaller( iter1, iter2, matchfun1, matchfun2 = None):
    """output pairs where iter1 <= iter2."""

    if matchfun2 == None: matchfun2 = matchfun1
    value2 = iter2.next()
    key2 = matchfun2(value2)
    l = None
    for value1 in iter1:
        key1 = matchfun1(value1)
        while key2 < key1: 
            yield( l, value2 )
            value2 = iter2.next()
            key2 = matchfun2(value2)
        l = value1
    
    while 1:
        yield( l, value2 )
        value2 = iter2.next()
        key2 = matchfun2(value2)

def matchby_comparison( iter1, iter2, matchfun1, matchfun2 = None ):
    """return pairs of values from two streams that are matched by match 
    function.

    The two input streams need to be sorted by increasing key value.
    """

    if matchfun2 == None: matchfun2 = matchfun1
    groups1 = itertools.groupby( iter1, matchfun1 )
    groups2 = itertools.groupby( iter2, matchfun2 )
    key2, value2 = groups2.next()
    last_key1, last_key2 = key1, key2
    for key1, value1 in groups1:
        assert last_key1 < key1, "input needs to be sorted by keys in ascending order."
        last_key1 = key1
        while key1 > key2:
            last_key2 =  key2
            key2, value2 = groups2.next()
            assert last_key2 < key2, "input needs to be sorted by keys in ascending order."
        if key1 < key2:
            continue
        yield( value1, value2 )


def matchby_sequence( iter1, iter2, matchfun1, matchfun2 = None ):
    """return pairs of values from two streams that are matched by match 
    function.

    The two input streams should be sorted in the same order, not necessarily
    increasing, and all componets should be present.
    """

    if matchfun2 == None: matchfun2 = matchfun1
    groups1 = itertools.groupby( iter1, matchfun1 )
    groups2 = itertools.groupby( iter2, matchfun2 )

    while 1:
        key1, value1 = groups1.next()
        key2, value2 = groups2.next()

        assert key1 == key2, "input needs to be sorted in the same order."
        yield( value1, value2 )

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: maq2psl.py 2781 2009-09-10 11:33:14Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-c", "--filename-coordinates", dest="filename_coordinates", type="string",
                      help="filename with coordinates."  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename pattern for additional data [%default].")

    parser.set_defaults(
        genome_file = "genome",
        filename_coordinates = None,
        segment_length = 32,
        )

    (options, args) = E.Start( parser )

    if options.genome_file:
        genome = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        genome = None

    ninput, noutput = 0, 0

    if options.filename_coordinates:
        segment_length = options.segment_length
        a = matchby_sequence(  iterator_segments( open( options.filename_coordinates, "r"), options.segment_length ),
                               Maq.iterator( options.stdin ),
                               lambda x: (x.mSegment),
                               lambda x: (x.contig) )
        
        for segments, maqs in a:
            pairs = match_smaller( segments, maqs,
                                   lambda x: x.start, 
                                   lambda x: x.start ) 

            for segment, maq in pairs:
                ninput += 1

                assert maq.start >= segment.start, "maq start < segment start: %i < %i" % (maq.start, segment.start)
                assert maq.start + maq.mLength <= segment.start + 2 * segment_length, "maq end > segment end: %i < %i" % (maq.start + maq.mLength, segment.start + 2 * segment_length)
        
                psl = Blat.Match()
                psl.fromMaq( maq )

                match_start = maq.start
                segment_start = segment.start
                contig, left_start, right_start = segment.contig, segment.mLeftStart, segment.mRightStart

                if options.loglevel >= 2:
                    options.stdlog.write("# mapping: name=%s, match_start=%i, segment=%s\n" % (maq.contig, match_start, str(segment)))

                # build positions of the two blocks
                left_size = segment_length - (match_start - segment_start)
                right_size = segment_length - left_size
                mapped1_start = left_start + match_start - segment_start
                mapped1_end   = left_start + segment_length
                mapped2_start = right_start
                mapped2_end   = right_start + right_size

                if options.loglevel >= 3:
                    options.stdlog.write("# mapped: match_start=%i, segment_start=%i, left_size=%i, right_size=%i, mapped1=(%i-%i), mapped2=(%i-%i)\n" %\
                                             (match_start, segment_start, left_size, right_size, mapped1_start, mapped1_end, mapped2_start, mapped2_end) )


                psl.mSbjctId = contig
                if genome: psl.mSbjctLength = genome.getLength( contig )
                psl.mSbjctFrom = mapped1_start
                psl.mSbjctTo = mapped2_end
                psl.mNBlocks = 2
                psl.mBlockSizes= [left_size, right_size]
                psl.mQueryBlockStarts = [0, left_size]
                psl.mSbjctBlockStarts = [mapped1_start, mapped2_start]
                psl.mSbjctNGapsCounts = 1
                psl.mSbjctNGapsBases = mapped2_start - mapped1_end

                options.stdout.write( str(psl) + "\n" )
                noutput += 1

    else:

        for maq in Maq.iterator( options.stdin ):
            ninput += 1

            psl = Blat.Match()
            psl.fromMaq( maq )

            options.stdout.write( str(psl) + "\n" )
            noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i\n" % (ninput, noutput) )

    E.Stop()
