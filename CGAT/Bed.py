'''Bed.py - Tools for working with bed files
=========================================

This module contains methods for working with :term:`bed`
formatted files.

.. note::
   Another way to access the information in :term:`bed` formatted
   files is through pysam_.

The principal class is :class:`Bed` to represent :term:`bed` formatted
entries.  The method :func:`iterate` iterates over a bed file and is
aware of UCSC track information that might be embedded in the
file. Additional functions can process intervals (:func:`merge`,
:func:`binIntervals`, :func:`setName`, etc).

The method :func:`readAndIndex` can build an in-memory index of a bed-file
for quick cross-referencing.

Reference
---------

'''
import re
import numpy
import bisect
import itertools

from CGAT import NCL as ncl
from CGAT import IOTools as IOTools

Headers = [
    "contig", "start", "end",
    "name", "score", "strand",
    "thinkStart", "thickEnd",
    "itemRGB", "blockCount",
    "blockSizes", "blockStarts"]


class Bed(object):
    """an interval in bed format.

    Coordinates are represented as 0-based, half-open intervals.

    Fields in the record can be accessed as attributes or through
    a dictionary type access::

       print b.contig()
       print b["contig"]

    Bed-formatted records can have a variable number of columuns
    with a minimum of 3. Accessing an optional attribute that is not present
    will raise an IndexError.

    Attributes
    ----------
    contig : string
       Chromosome/contig.
    start : int
       Start position of the interval.
    end : int
       End position of the interval.
    name : string
       Name of the interval (optional).
    score : float
       Score associated with interval (optional).
    strand : char
       Strand of the interval (optional).
    thickStart
    thickEnd
    itemRGB
    blockCount : int
       Number of blocks for bed intervals spanning multiple blocks (BED12).
    blockSizes : string
       Comma-separated list of sizes of the blocks (BED12).
    blockStarts : string
       Comma-separated list of start positions of the blocks (BED12).
    """

    map_key2field = {'name': 0,
                     'score': 1,
                     'strand': 2,
                     'thickStart': 3,
                     'thickEnd': 4,
                     'itemRGB': 5,
                     'blockCount': 6,
                     'blockSizes': 7,
                     'blockStarts': 8}

    default_value = "."

    def __init__(self):
        self.contig = None
        self.start = 0
        self.end = 0
        self.fields = []
        self.track = None

    def __str__(self):
        return "\t".join((self.contig, str(self.start),
                          str(self.end)) + tuple(map(str, self.fields)))

    def copy(self):
        '''Returns a new bed object that is a copy of this one'''

        new_entry = Bed()
        new_entry.__dict__ = self.__dict__.copy()
        return new_entry

    def fromGTF(self, gff, is_gtf=False, name=None):
        """fill fields from gtf formatted entry

        Arguments
        ---------
        gff : a gff entry.
           The object should contain the fields ``contig``,
           ``start`` and ``end`` in 0-based, half-open coordinates.
        name : bool
           If given, attempt to set the name atttribute of the interval
           by this attribute of the `gff` object such as ``gene_id`` or
           ``transcript_id``.
        """

        self.contig, self.start, self.end = gff.contig, gff.start, gff.end
        try:
            self.fields = [getattr(gff, name),
                           [gff.score, 0][gff.score is None],
                           gff.strand]
        except AttributeError:
            self.fields = [gff[name],
                           [gff.score, 0][gff.score is None],
                           gff.strand]

    def toIntervals(self):
        """return intervals for BED12 entries.

        If the entry is not BED12, the whole region will be returned.

        Returns
        -------
        intervals : list
           A list of tuples (start,end) with the block coordinates in
           the Bed entry.

        """

        if self.columns >= 12:

            blockStarts = list(map(int, self.blockStarts.split(",")))
            starts = [self.start + blockStart for blockStart in blockStarts]
            blockLengths = list(map(int, self.blockSizes.split(",")))

            assert (len(blockStarts), len(blockLengths)) == \
                (int(self.blockCount), int(self.blockCount)), \
                "Malformed Bed12 entry:\n%s"

            ends = [start + length for start, length
                    in zip(starts, blockLengths)]

            assert ends[-1] == self.end, \
                "Malformed Bed12 entry:\n%s" % str(self)

            return list(zip(starts, ends))
        else:
            return [(self.start, self.end)]

    def fromIntervals(self, intervals):
        """Fill co-ordinates from list of intervals.

        If multiple intervals are provided and entry is BED12 then the
        blocks are automatically set.

        Arguments
        ---------
        intervals : list
           List of tuples (start, end) with block coordinates.
        """

        intervals = sorted(intervals)
        self.start = intervals[0][0]
        self.end = intervals[-1][1]

        if self.columns >= 12:
            self["thickStart"] = self.start
            self["thickEnd"] = self.end

            blockStarts = [interval[0] - self.start for interval in intervals]
            blockSizes = [end - start for start, end in intervals]

            blockCount = len(intervals)

            self["blockStarts"] = ",".join(map(str, blockStarts))
            self["blockSizes"] = ",".join(map(str, blockSizes))
            self["blockCount"] = str(blockCount)

        else:
            if len(intervals) > 1:
                raise ValueError(
                    "Multiple intervals provided to non-bed12 entry")

    def __contains__(self, key):
        return self.map_key2field[key] < len(self.fields)

    def compare(self, other):
        a = (self.contig, self.start, self.end)
        b = (other.contig, other.start, other.end)
        return (a > b) - (a < b)

    def __lt__(self, other):
        return self.__richcmp__(other, 0)

    def __richcmp__(self, other, op):
        retval = self.compare(other)
        return ((op == 0 and retval < 0) or
                (op == 1 and retval <= 0) or
                (op == 2 and retval == 0) or
                (op == 3 and retval != 0) or
                (op == 4 and retval > 0) or
                (op == 5 and retval >= 0))

    def __getitem__(self, key):
        return self.fields[self.map_key2field[key]]

    def __setitem__(self, key, value):
        try:
            position = self.map_key2field[key]
        except IndexError:
            raise IndexError("Unknown key: %s" % key)

        try:
            self.fields[position] = value
        except IndexError:

            self.fields.extend([self.default_value] *
                               (position - len(self.fields) + 1))

            self.fields[position] = value

    def __getattr__(self, key):
        try:
            return self.fields[self.map_key2field[key]]
        except IndexError:
            return None

    @property
    def columns(self):
        '''return number of columns in bed-entry.'''
        return 3 + len(self.fields)


class Track(object):
    """Bed track information."""

    def __init__(self, line):
        r = re.compile('([^ =]+) *= *("[^"]*"|[^ ]*)')

        self._d = {}
        for k, v in r.findall(line[:-1]):
            if v[:1] == '"':
                self._d[k] = v[1:-1]
            else:
                self._d[k] = v

        self._line = line[:-1]

    def __str__(self):
        return self._line

    def __getitem__(self, key):
        return self._d[key]

    def __setitem__(self, key, val):
        self._d[key] = val


def iterator(infile):
    """iterate over a :term:`bed` formatted file.

    Comments and empty lines are ignored. The iterator is
    :term:`track` aware and will set the ``track`` attribute for the
    Bed objects it yields.

    Arguments
    ---------
    infile : File

    Yields
    ------
    bed
       :class:`Bed` object

    """

    track = None
    for line in infile:
        if line.startswith("track"):
            track = Track(line)
            continue
        # ignore comments
        if line.startswith("#"):
            continue
        # ignore empty lines (in order to parse pseudo bed files)
        if line.strip() == "":
            continue

        b = Bed()
        # split at tab (Bed standard, do not split at space as this will split
        # the name field)
        data = line[:-1].split("\t")
        try:
            b.contig, b.start, b.end = data[0], int(data[1]), int(data[2])
        except IndexError:
            raise ValueError("parsing error in line '%s'" % line[:-1])
        b.fields = data[3:]
        b.track = track
        yield b


def bed_iterator(infile):
    """Deprecated, use :func:`iterator`."""
    return iterator(infile)


def setName(iterator):
    """yield bed entries in which name is set to the record number if
    unset.

    Yields
    ------
    bed
       :class:`Bed` object
    """
    for i, bed in enumerate(iterator):
        if "name" not in bed:
            bed.name = str(i)
        yield bed


def grouped_iterator(iterator):
    """yield bed results grouped by track.

    Note that the iterator supplied needs to be sorted by the track
    attribute. This is usually the case in :term:`bed` formatted
    files.

    Yields
    ------
    bed
       :class:`Bed` object

    """
    return itertools.groupby(iterator, lambda x: x.track)


def blocked_iterator(iterator):
    '''yield blocked bed results.

    Intervals with the same name are merged into a single entry. This
    method can be used to convert BED6 formatted entries to
    BED12. Note that the input iterator needs to be sorted by bed
    name.

    Yields
    ------
    bed
       :class:`Bed` object

    '''

    last_id = None
    blocks = []

    def _update(bed, blocks):
        blocks.sort()
        bed.start, bed.end = blocks[0][0], blocks[-1][1]
        s = bed.start
        # hacky - needs be abstracted into Bed object
        bed.fields.extend([""] * (9 - len(bed.fields)))
        bed.fields[3] = str(bed.start)
        bed.fields[4] = str(bed.end)
        bed.fields[5] = 0
        bed.fields[6] = len(blocks)
        bed.fields[7] = ",".join([str(y - x) for x, y in blocks])
        bed.fields[8] = ",".join([str(x - s) for x, y in blocks])
        return bed

    last_bed = None
    for bed in iterator:
        if last_id != bed.name:
            if last_id:
                yield _update(last_bed, blocks)
            blocks = []
            last_id = bed.name
        last_bed = bed
        blocks.append((bed.start, bed.end))

    yield _update(bed, blocks)


def readAndIndex(infile, with_values=False, per_track=False):
    """read and index a bed formatted file in ``infile``.

    The index is not strand-aware.

    Arguments
    ---------
    infile : File
       File object to read from.
    with_values : bool
       If True, store the actual bed entry. Otherwise, just the
       intervals are recorded and any additional fields will be ignored.
    per_track : bool
       If True build indices per track.

    Returns
    -------
    index : dict
       A dictionary of nested containment lists (:term:`NCL`). Each
       key is a contig. If `per_track` is set, the dictionary has an
       additional first level for the track.
    """

    if with_values:
        idx_factory = ncl.NCL
    else:
        idx_factory = ncl.NCLSimple

    def _build(iter):
        idx = {}
        for e in iter:
            if e.contig not in idx:
                idx[e.contig] = idx_factory()
            try:
                if with_values == "name":
                    idx[e.contig].add(e.start, e.end, e.name)
                elif with_values is True:
                    idx[e.contig].add(e.start, e.end, e)
                else:
                    idx[e.contig].add(e.start, e.end)
            except ValueError:
                # ignore zero-length intervals
                pass
        return idx

    if per_track:
        indices = {}
        for track, beds in grouped_iterator(iterator(infile)):
            if track is None:
                return _build(beds)
            else:
                indices[track["name"]] = _build(beds)
        return indices
    else:
        return _build(iterator(infile))


def binIntervals(iterator, num_bins=5, method="equal-bases", bin_edges=None):
    """merge adjacent intervals by the score attribute.

    This method takes all the intervals in the collection builds a histogram
    of all the scores in the collection. The partition into the bins can use
    one of the following merging methods:

    equal-bases
       merge intervals such that each bin contains the equal number of bases

    equal-intervals
       merge intervals such that each bin contains the equal number intervals

    This method requires the fifth field (score) of the bed input file
    to be present.

    Arguments
    ---------
    iterator
       Iterator yielding bed intervals
    num_bins : int
       Number of bins to create in the histogram
    method : string
       Binning method
    bin_edges : list
       List of bin edges. These take precedence over `method`.

    Returns
    -------
    intervals : list
       list of intervals (:class:`Bed`)
    bin_edges : list
       list of bin edges

    """

    data = []
    beds = list(iterator)

    for bed in beds:
        bed.fields[1] = float(bed.fields[1])

    if bin_edges is None:
        if method == "equal-bases":
            data = numpy.array(
                sorted([(x.fields[1], x.end - x.start) for x in beds]))
        elif method == "equal-intervals":
            data = numpy.array(sorted([(x.fields[1], 1) for x in beds]))
        elif method == "equal-range":
            vals = [x.fields[1] for x in beds]
            mi, ma = min(vals), max(vals)
            increment = float(ma - mi) / num_bins
            bin_edges = numpy.arange(mi, ma, increment)
        else:
            raise ValueError(
                "unknown method %s to compute bins, supply bin_edges" % method)

    if bin_edges is None:
        sums = data[:, 1].cumsum(axis=0)
        total = float(sums[-1])
        increment = float(total / num_bins)
        bin_edges = [data[0][0]]
        occupancies = []
        threshold = increment
        occ = 0
        for v, s in zip(data, sums):
            occ += v[1]
            if s > threshold:
                bin_edges.append(v[0])
                threshold += increment
                occupancies.append(occ)
                occ = 0

        bin_edges.append(data[-1][0] + 1)

    beds.sort(key=lambda x: (x.contig, x.start))
    last_contig, start, end, last_name = None, None, None, None
    new_beds = []
    for bed in beds:
        name = bisect.bisect_right(bin_edges, bed.fields[1]) - 1
        contig = bed.contig
        if name != last_name or last_contig != contig:
            if last_name is not None:
                b = Bed()
                b.contig, b.start, b.end, b.fields = last_contig, start, end, [
                    last_name]
                new_beds.append(b)
            start = bed.start
            last_name = name
            last_contig = contig

        end = bed.end

    if last_name is not None:
        b = Bed()
        b.contig, b.start, b.end, b.fields = contig, start, end, [last_name]
        new_beds.append(b)

    return new_beds, bin_edges


def merge(iterator):
    '''merge overlapping intervals and returns a list of merged intervals.
    '''

    beds = list(iterator)
    if len(beds) == 0:
        return []

    beds.sort(key=lambda x: (x.contig, x.start))

    def iterate_chunks(beds):

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
            to_join.append(this)

        yield to_join, end
        raise StopIteration

    n = []
    for to_join, end in iterate_chunks(beds):

        y = Bed()
        y.contig = to_join[0].contig
        y.start = to_join[0].start
        y.end = end
        n.append(y)

    return n


def getNumColumns(filename):
    '''return number of fields in bed-file by looking at the first
    entry.

    Returns
    -------
    ncolumns : int
       The number of columns. If the file is empty, 0 is returned.
    '''
    with IOTools.openFile(filename) as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            return len(line[:-1].split("\t"))
    return 0
