import re
import ncl
import numpy
import bisect
import itertools

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
        self.contig = None
        self.start = 0
        self.end = 0
        self.mFields = []
        self.mTrack = None
        
    def __str__(self):
        return "\t".join( (self.contig, str(self.start), str(self.end) ) + tuple(map(str, self.mFields)))

    def fromGFF( self, gff, is_gtf = False ):
        """fill from gff formatted entry."""
        self.contig, self.start, self.end = gff.contig, gff.start, gff.end
        if is_gtf: self.mFields = [gff.gene_id]

    def __contains__(self, key ):
        return self.map_key2field[key] < len(self.mFields)

    def __getitem__(self, key):
        return self.mFields[self.map_key2field[key]]
        
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
    """iterate for parsing a bed formatted file.

    This iterator yields :class:`Bed` objects.
    """

    track = None
    for line in infile:
        if line.startswith("track"):
            track = Track( line )
            continue
        if line.startswith("#"): continue

        b = Bed()
        data = line[:-1].split()
        try:
            b.contig, b.start, b.end = data[0], int(data[1]), int(data[2])
        except IndexError:
            raise ValueError("parsing error in line '%s'" % line[:-1])
        b.mFields = data[3:]
        b.mTrack = track
        yield b

# for compatibility, remove
def bed_iterator( infile ):
    return iterator( infile )

def grouped_iterator( iterator ):
    '''return bed results grouped by track.'''
    return itertools.groupby( iterator, lambda x: x.mTrack )

def readAndIndex( infile, with_values = False, per_track = False ):
    """read and index a bed formatted file in *infile*.

    If *with_values* is set, the original bed entry will be kept for
    further referenc. Otherwise only the intervals will be indexed and
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

    Several merging methods are possible.
    equal-bases, equal-intervals.
    
    This options requires the fifth field of the bed input file to be
    present.

    bins should be non-overlapping.'''

    data = []
    beds = list(iterator)

    for bed in beds: bed.mFields[1] = float( bed.mFields[1])

    if bin_edges == None:
        if method == "equal-bases":
            data = numpy.array( sorted( [(x.mFields[1], x.end-x.start) for x in beds ] ) )
        elif method == "equal-intervals":
            data = numpy.array( sorted( [(x.mFields[1], 1) for x in beds ] ) )
        elif method == "equal-range":
            vals = [x.mFields[1] for x in beds ]
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
        name = bisect.bisect_right(bin_edges, bed.mFields[1])-1
        contig = bed.contig
        if name != last_name or last_contig != contig:
            if last_name != None:
                b = Bed()
                b.contig, b.start, b.end, b.mFields = last_contig, start, end, [last_name]
                new_beds.append(b)
            start = bed.start
            last_name = name
            last_contig = contig
            
        end = bed.end

    if last_name != None:
        b = Bed()
        b.contig, b.start, b.end, b.mFields = contig, start, end, [last_name]
        new_beds.append(b)
        
    return new_beds, bin_edges
    
def merge( iterator ):
    '''merge bed segments.'''

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
                yield to_join
                contig = this.contig
                to_join = []
                end = this.end

            end = max(end, this.end)
            to_join.append( this )
        
        yield to_join
        raise StopIteration

    n = []
    for to_join in iterate_chunks(beds):

        y = Bed()
        y.contig = to_join[0].contig
        y.start = to_join[0].start
        y.end = to_join[-1].end
        n.append( y )

    return n
