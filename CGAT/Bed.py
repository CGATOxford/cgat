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
Bed.py - Tools for working with bed files
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import re
import numpy
import bisect
import itertools

import NCL as ncl

Headers = [
    "contig", "start", "end",
    "name", "score", "strand",
    "thinkStart", "thickEnd",
    "itemRGB", "blockCount",
    "blockSizes", "blockStarts" ]

class Bed(object):
    """an interval in bed format."""

    map_key2field = { 'name' : 0,
                      'score' : 1,
                      'strand' : 2,
                      'thickStart' : 3,
                      'thickEnd' : 4,
                      'itemRGB' : 5,
                      'blockCount': 6,
                      'blockSizes': 7,
                      'blockStarts': 8 }

    def __init__(self):
        '''empty constructor.'''
        self.contig = None
        self.start = 0
        self.end = 0
        self.fields = []
        self.track = None
        
    def __str__(self):
        return "\t".join( (self.contig, str(self.start), str(self.end) ) + tuple(map(str, self.fields)))

    def fromGTF( self, gff, is_gtf = False, name=None ):
        """fill from gtf formatted entry."""
        self.contig, self.start, self.end = gff.contig, gff.start, gff.end
        try:
            self.fields = [getattr( gff, name), 
                            [gff.score,0][gff.score == None], 
                            gff.strand ]
        except AttributeError:
            self.fields = [gff[name],
                            [gff.score,0][gff.score == None], 
                            gff.strand ]
            
    def __contains__(self, key ):
        return self.map_key2field[key] < len(self.fields)

    def __getitem__(self, key):
        return self.fields[self.map_key2field[key]]

    def __setitem__(self, key, value):
        self.fields[self.map_key2field[key]] = value

    def __getattr__(self, key ):
        try: 
            return self.fields[self.map_key2field[key]]
        except IndexError:
            return None

    @property
    def columns(self):
        '''return number of columns in bed-entry.'''
        return 3 + len(self.fields)

class Track(object):
    '''bed track information.'''
    def __init__(self, line ):
        r= re.compile('([^ =]+) *= *("[^"]*"|[^ ]*)')

        self._d = {}
        for k, v in r.findall(line[:-1]):
            if v[:1]=='"':
                self._d[k] = v[1:-1]
            else:
                self._d[k] = v
            
        self._line = line[:-1]
    def __str__(self):
        return self._line

    def __getitem__(self, key): return self._d[key]
    def __setitem__(self, key,val): self._d[key] = val
        
def iterator( infile ):
    """iterate over a bed formatted file.

    The iterator is :term:`track` aware.

    This iterator yields :class:`Bed` objects. 
    """

    track = None
    for line in infile:
        if line.startswith("track"):
            track = Track( line )
            continue
        # ignore comments
        if line.startswith("#"): continue
        # ignore empty lines (in order to parse pseudo bed files)
        if line.strip() == "": continue

        b = Bed()
        # split at tab (Bed standard, do not split at space as this will split the name field)
        data = line[:-1].split("\t")
        try:
            b.contig, b.start, b.end = data[0], int(data[1]), int(data[2])
        except IndexError:
            raise ValueError("parsing error in line '%s'" % line[:-1])
        b.fields = data[3:]
        b.track = track
        yield b

# for compatibility, remove
def bed_iterator( infile ):
    return iterator( infile )

def setName( iterator ):
    '''yield bed entries in which name is set if unset.
    '''
    for i, bed in enumerate(iterator):
        if "name" not in bed:
            bed.name = str(i)
        yield bed

def grouped_iterator( iterator ):
    '''yield bed results grouped by track.'''
    return itertools.groupby( iterator, lambda x: x.track )

def blocked_iterator( iterator ):
    '''yield blocked bed results.'''

    last_id = None
    blocks = []

    def _update( bed, blocks ):
        blocks.sort()
        bed.start, bed.end = blocks[0][0], blocks[-1][1]
        s = bed.start
        # hacky - needs be abstracted into Bed object
        bed.fields.extend( [""] * (9 - len(bed.fields)))
        bed.fields[3] = str(bed.start)
        bed.fields[4] = str(bed.end)
        bed.fields[5] = 0
        bed.fields[6] = len(blocks)
        bed.fields[7] = ",".join( [str(y-x) for x,y in blocks ])
        bed.fields[8] = ",".join( [str(x-s) for x,y in blocks])
        return bed

    last_bed = None
    for bed in iterator:
        if last_id != bed.name:
            if last_id: yield _update(last_bed, blocks)
            blocks = []
            last_id = bed.name
        last_bed = bed
        blocks.append( (bed.start, bed.end) )
    
    yield _update(bed, blocks)

def readAndIndex( infile, with_values = False, per_track = False ):
    """read and index a bed formatted file in ``infile``.

    If ``with_values`` is set, the original bed entry will be kept for
    further reference. Otherwise only the intervals will be indexed and
    any additional fields in the bed entry will be ignored.

    The default is to use all intervals. If per_track is set,
    separate indices will be build for each track.
    """

    if with_values:
        idx_factory = ncl.NCL
    else:
        idx_factory = ncl.NCLSimple
        
    def _build( iter ):
        idx = {}
        for e in iter:
            if e.contig not in idx: idx[e.contig] = idx_factory()
            try:
                if with_values:
                    idx[e.contig].add( e.start,e.end,e )
                else:
                    idx[e.contig].add( e.start,e.end )
            except ValueError:
                # ignore zero-length intervals
                pass
        return idx

    if per_track:
        indices = {}
        for track, beds in grouped_iterator( bed_iterator( infile ) ):
            if track == None:
                return _build(beds)
            else:
                indices[track["name"]] = _build( beds )
        return indices
    else:
        return _build( bed_iterator( infile ) )

def binIntervals( iterator, num_bins = 5, method = "equal-bases", bin_edges = None ):
    '''merge adjacent bins by score.

    Several merging methods are possible:

    equal-bases
       merge intervals such that each bin contains the equal number of bases

    equal-intervals
       merge intervals such that each bin contains the equal number intervals

    This options requires the fifth field (score) of the bed input 
    file to be present.

    bins should be non-overlapping.

    returns a list of intervals (:class:`Bed`) and the bin_edges
    '''

    data = []
    beds = list(iterator)

    for bed in beds: bed.fields[1] = float( bed.fields[1])

    if bin_edges == None:
        if method == "equal-bases":
            data = numpy.array( sorted( [(x.fields[1], x.end-x.start) for x in beds ] ) )
        elif method == "equal-intervals":
            data = numpy.array( sorted( [(x.fields[1], 1) for x in beds ] ) )
        elif method == "equal-range":
            vals = [x.fields[1] for x in beds ]
            mi, ma = min( vals ), max( vals )
            increment = float(ma - mi) / num_bins
            bin_edges = numpy.arange( mi, ma, increment)
        else:
            raise ValueError("unknown method %s to compute bins, supply bin_edges" % method )

    if bin_edges == None:
        sums = data[:,1].cumsum( axis=0 )
        total = float(sums[-1])
        increment = float( total / num_bins )
        bin_edges = [data[0][0]]
        occupancies = []
        threshold = increment
        occ = 0
        for v, s in zip( data, sums ):
            occ += v[1]
            if s > threshold:
                bin_edges.append( v[0] )
                threshold += increment
                occupancies.append( occ )
                occ = 0

        bin_edges.append( data[-1][0] + 1 )

    beds.sort( key=lambda x: (x.contig,x.start) )
    last_contig, start, end, last_name = None, None, None, None
    new_beds = []
    for bed in beds:
        name = bisect.bisect_right(bin_edges, bed.fields[1])-1
        contig = bed.contig
        if name != last_name or last_contig != contig:
            if last_name != None:
                b = Bed()
                b.contig, b.start, b.end, b.fields = last_contig, start, end, [last_name]
                new_beds.append(b)
            start = bed.start
            last_name = name
            last_contig = contig
            
        end = bed.end

    if last_name != None:
        b = Bed()
        b.contig, b.start, b.end, b.fields = contig, start, end, [last_name]
        new_beds.append(b)
        
    return new_beds, bin_edges
    
def merge( iterator ):
    '''merge overlapping intervals.

    returns a list of merged intervals.
    '''

    beds = list(iterator)
    if len(beds) == 0: return []

    beds.sort( key = lambda x: (x.contig, x.start) )

    def iterate_chunks( beds ):
        
        last = beds[0]
        to_join = [last]
        end = last.end
        contig = last.contig

        for this in beds[1:]:
            d = this.start - end
            if this.contig != contig or d >= 0:
                yield to_join, end
                contig = this.contig
                to_join = []
                end = this.end

            end = max(end, this.end)
            to_join.append( this )
        
        yield to_join, end
        raise StopIteration

    n = []
    for to_join, end in iterate_chunks(beds):

        y = Bed()
        y.contig = to_join[0].contig
        y.start = to_join[0].start
        y.end = end
        n.append( y )

    return n
