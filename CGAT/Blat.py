'''Blat.py - tools for working with PSL formatted files and data
=============================================================

This module provides a class to parse :term:`PSL` formatted
files such as those output by the BLAT tool.

This module defines the :class:`Blat.Match` class representing a
single entry and a series of iterators to iterate of :term:`PSL`
formatted files (:func:`iterator`, :func:`iterator_target_overlap`,
...).

Reference
---------

'''
import copy
import string
import collections

try:
    import alignlib_lite
except ImportError:
    pass

from CGAT import Components as Components
from CGAT import Experiment as E


class Error(Exception):
    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self._message)

    def _get_message(self):
        return self._message

    def _set_message(self, message):
        self._message = message
    message = property(_get_message, _set_message)


class ParsingError(Error):

    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line=None):
        if line:
            self.message = message + " at line:" + line
        else:
            self.message = message

HEADER = """psLayout version 3

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q       Q       Q       Q       T       T       T       T       block   blockSizes      qStarts  tStarts
        match                   count   bases   count   bases           name    size    start   end     name    size    start   end     count
---------------------------------------------------------------------------------------------------------------------------------------------------------------"""


class Match:

    """a :term:`psl` formatted alignment.

    Block coordinates are on the forward strand for target and on
    the forward/reverse strand for the query depending on the strand.

    The fields mQueryFrom/To and mSbjctFrom/To are always on the forward
    strand.
    """

    def __init__(self):
        self.mNMatches = 0
        self.mNMismatches = 0
        self.mNRepMatches = 0
        self.mNns = 0
        self.mQueryNGapsCounts = 0
        self.mQueryNGapsBases = 0
        self.mSbjctNGapsCounts = 0
        self.mSbjctNGapsBases = 0
        self.strand = "+"
        self.mQueryId = ""
        self.mQueryLength = 0
        self.mQueryFrom = 0
        self.mQueryTo = 0
        self.mSbjctId = ""
        self.mSbjctLength = 0
        self.mSbjctFrom = 0
        self.mSbjctTo = 0
        self.mNBlocks = 0
        self.mBlockSizes = []
        self.mQueryBlockStarts = []
        self.mSbjctBlockStarts = []

    def convertCoordinates(self):
        """convert coordinates.

        This rescales the block positions so that they start at 0 and converts
        the query to forward and the sbjct to forward/reverse coordinates.

        About the psl psl format from the manual at
        http://genome.ucsc.edu/google/goldenPath/help/pslSpec.html

        ::
           In general the coordinates in psl files are "zero based
           half open." The first base in a sequence is numbered zero
           rather than one.  When representing a range the end
           coordinate is not included in the range. Thus the first 100
           bases of a sequence are represented as 0-100, and the
           second 100 bases are represented as 100-200.

           There is a another little unusual feature in the .psl
           format. It has to do with how coordinates are handled on
           the negative strand.  In the qStart/qEnd fields the
           coordinates are where it matches from the point of view of
           the forward strand (even when the match is on the reverse
           strand). However on the qStarts[] list, the coordinates are
           reversed.

        This class works in forward coordinates for the query and
        forward/reverse coordinates for the sbjct.

        For a negative strand match, the following is done:
           * invert mSbjctFrom and mSbjctTo with mSbjctLength
           * add block sizes to mQueryStarts and mSbjctStarts
           * invert mQueryStarts and mSbjctStarts
           * reverse blocksize, mQueryStarts and mSbjctStarts

        """

        block_sizes = self.mBlockSize
        query_starts = self.mQueryBlockStarts
        sbjct_starts = self.mSbjctBlockStarts

        if self.mSbjctStrand == "-":
            self.mSbjctFrom, self.mSbjctTo = self.mSbjctLength - \
                self.mSbjctTo, self.mSbjctLength - self.mSbjctFrom
            query_starts = [
                self.mQueryLength - x - y for x, y in zip(query_starts, block_sizes)]
            sbjct_starts = [
                self.mSbjctLength - x - y for x, y in zip(sbjct_starts, block_sizes)]
            query_starts.reverse()
            sbjct_starts.reverse()

        self.mQueryStarts = [x - self.mQueryFrom for x in query_starts]
        self.mSbjctStarts = [x - self.mSbjctFrom for x in sbjct_starts]

    def _switchQueryStrandFlag(self):
        '''switch the strand flag of the query
        '''

        if self.strand == "-":
            self.strand = "+"
        else:
            self.strand = "-"

    def switchTargetStrand(self):
        """switch the target strand.

        Use in cases in which a feature has been defined on the
        negative target strand with reverse coordinates. The result
        will be the same alignment using forward coordinates on the
        target.

        This method will also update the query strand and coordinates.
        """

        block_sizes = self.mBlockSizes
        query_starts = self.mQueryBlockStarts
        sbjct_starts = self.mSbjctBlockStarts

        # invert alignment
        self.mQueryBlockStarts = [
            self.mQueryLength - x - y
            for x, y in zip(query_starts, block_sizes)]
        self.mQueryBlockStarts.reverse()
        self.mSbjctBlockStarts = [
            self.mSbjctLength - x - y
            for x, y in zip(sbjct_starts, block_sizes)]
        self.mSbjctBlockStarts.reverse()
        self.mBlockSizes.reverse()

        # update target coordinates only
        self.mSbjctFrom, self.mSbjctTo = self.mSbjctLength - \
            self.mSbjctTo, self.mSbjctLength - self.mSbjctFrom
        self._switchQueryStrandFlag()

    def fromTable(self, data):

        try:
            (nmatches, nmismatches, nrepmatches, nns,
             query_ngaps_counts, query_ngaps_bases,
             sbjct_ngaps_counts, sbjct_ngaps_bases,
             strand,
             query_id, query_length, query_from, query_to,
             sbjct_id, sbjct_length, sbjct_from, sbjct_to,
             nblocks, block_sizes,
             query_block_starts, sbjct_block_starts) = data[:21]
        except ValueError:
            raise ParsingError("parsing error: %i fields" %
                               len(data), "\t".join(data))

        nmatches, nmismatches, nrepmatches, nns = list(map(
            int,
            (nmatches, nmismatches, nrepmatches, nns)))

        self.mNMatches = nmatches
        self.mNMismatches = nmismatches
        self.mNRepMatches = nrepmatches
        self.mNns = nns
        self.mQueryNGapsCounts = int(query_ngaps_counts)
        self.mQueryNGapsBases = int(query_ngaps_bases)
        self.mSbjctNGapsCounts = int(sbjct_ngaps_counts)
        self.mSbjctNGapsBases = int(sbjct_ngaps_bases)
        self.strand = strand
        self.mQueryId = query_id
        self.mQueryLength = int(query_length)
        self.mQueryFrom = int(query_from)
        self.mQueryTo = int(query_to)
        self.mSbjctId = sbjct_id
        self.mSbjctLength = int(sbjct_length)
        self.mSbjctFrom = int(sbjct_from)
        self.mSbjctTo = int(sbjct_to)
        self.mNBlocks = int(nblocks)
        self.mBlockSizes = list(map(int, block_sizes[:-1].split(",")))
        self.mQueryBlockStarts = list(
            map(int, query_block_starts[:-1].split(",")))
        self.mSbjctBlockStarts = list(
            map(int, sbjct_block_starts[:-1].split(",")))

        # this makes sure that the block positions are rescaled
        if self.mQueryLength != 0:
            self.mQueryCoverage = 100.0 * \
                (self.mNMismatches + self.mNMatches) / self.mQueryLength
        else:
            self.mQueryCoverage = 0

        if self.mSbjctLength != 0:
            self.mSbjctCoverage = 100.0 * \
                (self.mNMismatches + self.mNMatches) / self.mSbjctLength
        else:
            self.mSbjctCoverage = 0

        if nmatches + nmismatches > 0:
            self.mPid = 100.0 * float(nmatches) / (nmatches + nmismatches)
        else:
            self.mPid = 100.0

    def fromMaq(self, maq):
        """build BLAT entry from a MAQ match.

        see :class:`Maq.Match`.
        """
        start = maq.start - 1
        l = maq.mLength
        self.mQueryId = maq.mName
        self.mSbjctId = maq.contig
        self.mQueryFrom = 0
        self.mQueryTo = l
        self.mQueryLength = l
        self.strand = maq.strand
        self.mSbjctFrom = start
        self.mSbjctTo = start + l
        self.mNMismatches = maq.mNMismatches
        self.mNMatches = l - self.mNMismatches
        self.mBlockSizes = [l]
        self.mQueryBlockStarts = [0]
        self.mSbjctBlockStarts = [start]
        self.mNBlocks = 1

    def getBlocks(self):
        """return a list of aligned blocks."""
        return list(zip(self.mQueryBlockStarts,
                        self.mSbjctBlockStarts,
                        self.mBlockSizes))

    def __str__(self):
        return "\t".join(map(str, (
            self.mNMatches,
            self.mNMismatches,
            self.mNRepMatches,
            self.mNns,
            self.mQueryNGapsCounts,
            self.mQueryNGapsBases,
            self.mSbjctNGapsCounts,
            self.mSbjctNGapsBases,
            self.strand,
            self.mQueryId,
            self.mQueryLength,
            self.mQueryFrom,
            self.mQueryTo,
            self.mSbjctId,
            self.mSbjctLength,
            self.mSbjctFrom,
            self.mSbjctTo,
            self.mNBlocks,
            ",".join(map(str, self.mBlockSizes)) + ",",
            ",".join(map(str, self.mQueryBlockStarts)) + ",",
            ",".join(map(str, self.mSbjctBlockStarts)) + ",")))

    def getMapQuery2Target(self):
        """return a map between query to target.

        If the strand is "-", the coordinates for query are on
        the negative strand.
        """

        map_query2target = alignlib_lite.py_makeAlignmentBlocks()

        f = alignlib_lite.py_AlignmentFormatBlat(
            "%i\t%i\t%i\t%i\t%s\t%s\t%s\n" % (
                min(self.mQueryBlockStarts),
                max(self.mQueryBlockStarts),
                min(self.mSbjctBlockStarts),
                max(self.mSbjctBlockStarts),
                ",".join([str(x) for x in self.mQueryBlockStarts]) + ",",
                ",".join([str(x) for x in self.mSbjctBlockStarts]) + ",",
                ",".join([str(x) for x in self.mBlockSizes]) + ","))
        f.copy(map_query2target)

        return map_query2target

    def getMapTarget2Query(self):
        """return a map between target to query.

        If the strand is "-", the coordinates for query are on
        the negative strand.
        """

        map_target2query = alignlib_lite.py_makeAlignmentBlocks()

        f = alignlib_lite.py_AlignmentFormatBlat(
            "%i\t%i\t%i\t%i\t%s\t%s\t%s\n" % (
                min(self.mSbjctBlockStarts),
                max(self.mSbjctBlockStarts),
                min(self.mQueryBlockStarts),
                max(self.mQueryBlockStarts),
                ",".join([str(x) for x in self.mSbjctBlockStarts]) + ",",
                ",".join([str(x) for x in self.mQueryBlockStarts]) + ",",
                ",".join([str(x) for x in self.mBlockSizes]) + ","))
        f.copy(map_target2query)
        return map_target2query

    def copy(self):
        return copy.copy(self)

    def fromMap(self, map_query2target, use_strand=None):
        """return a map between query to target."""

        self.mNMatches = map_query2target.getNumAligned()
        f = str(alignlib_lite.py_AlignmentFormatBlat(map_query2target))

        self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo, \
            self.mQueryBlockStarts, self.mSbjctBlockStarts, self.mBlockSizes = f.split(
                "\t")

        if self.mBlockSizes:
            self.mBlockSizes = list(map(int, self.mBlockSizes[:-1].split(",")))
            self.mQueryBlockStarts = list(map(
                int, self.mQueryBlockStarts[:-1].split(",")))
            self.mSbjctBlockStarts = list(map(
                int, self.mSbjctBlockStarts[:-1].split(",")))
        else:
            self.mBlockSizes = []
            self.mQueryBlockStarts = []
            self.mSbjctBlockStarts = []

        self.mNBlocks = len(self.mBlockSizes)

        self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo = \
            list(map(int,
                     (self.mQueryFrom, self.mQueryTo,
                      self.mSbjctFrom, self.mSbjctTo)))

        # queryfrom and queryto are always forward strand coordinates
        if use_strand and self.strand == "-":
            self.mQueryFrom, self.mQueryTo = self.mQueryLength - \
                self.mQueryTo, self.mQueryLength - self.mQueryFrom

    def fromPair(self,
                 query_start, query_size, query_strand, query_seq,
                 target_start, target_size, target_strand, target_seq):
        '''fill from two aligned sequences.

        Note that sequences are case-sensitive.'''

        self.mQueryLength = query_size
        self.mSbjctLength = target_size

        map_query2target = alignlib_lite.py_makeAlignmentBlocks()

        assert len(query_seq) == len(target_seq)

        x, y = query_start, target_start
        nmatches, nmismatches = 0, 0
        for q, t in zip(query_seq, target_seq):
            tq, tt = q != "-", t != "-"
            if tq and tt:
                map_query2target.addPair(x, y)
                if q == t:
                    nmatches += 1
                else:
                    nmismatches += 1

            if tq:
                x += 1
            if tt:
                y += 1

        self.mNMatches, self.mNMismatches = nmatches, nmismatches
        self.strand = query_strand
        # the following call will set query_from, query_to for the forward strand
        # though block coordinates might be on the negative strand
        self.fromMap(map_query2target, use_strand=True)

        # if target is on negative strand, swop strands
        if target_strand == "-":
            # swap target strand - this will also swap the query strand
            self.switchTargetStrand()

    def iterator_introns(self):

        b = self.getBlocks()
        last_tend = b[0][1] + b[0][2]
        for qstart, tstart, size in b[1:]:
            yield (last_tend, tstart)
            last_tend = tstart + size
        raise StopIteration

    def iterator_exons(self):

        for qstart, tstart, size in self.getBlocks():
            yield tstart, tstart + size
        raise StopIteration

    def iterator_introns(self):

        b = self.getBlocks()
        last_tend = b[0][1] + b[0][2]
        for qstart, tstart, size in b[1:]:
            yield (last_tend, tstart)
            last_tend = tstart + size
        raise StopIteration

    def iterator_query_exons(self):

        for qstart, tstart, size in self.getBlocks():
            yield qstart, qstart + size
        raise StopIteration

    def iterator_sbjct_exons(self):

        for qstart, tstart, size in self.getBlocks():
            yield tstart, tstart + size
        raise StopIteration

    def iterator_query_introns(self):

        b = self.getBlocks()
        last_qend = b[0][0] + b[0][2]
        for qstart, tstart, size in b[1:]:
            yield (last_qend, qstart)
            last_qend = qstart + size
        raise StopIteration

    def getHeader(self):

        return HEADER

    def getHeaders(self):
        return ("match", "mismatch", "repeats", "Ns",
                "qGapCount", "qGapBases", "tGapCount", "tGapBases", "strand",
                "qName", "qSize", "qStart", "qEnd",
                "tName", "tSize", "tStart", "tEnd",
                "blockCount", "blockSizes",
                "qStarts", "tStarts")


class MatchPSLX(Match):

    def __init__(self):
        Match.__init__(self)
        self.mQuerySequence = []
        self.mSbjctSequence = []

    def __str__(self):
        return "\t".join((Match.__str__(self),
                          ",".join(self.mQuerySequence) + ",",
                          ",".join(self.mSbjctSequence) + ","))

    def fromTable(self, data):

        Match.fromTable(self, data)

        try:
            query_sequence, sbjct_sequence = data[21:23]

        except ValueError:
            raise ParsingError("parsing error", "\t".join(data))

        self.mQuerySequence = query_sequence[:-1].split(",")
        self.mSbjctSequence = sbjct_sequence[:-1].split(",")

    def getHeaders(self):
        return Match.getHeaders(self) + ("qSequence", "tSequence")

    def fromPSL(self, other, query_sequence, sbjct_sequence):
        """fill entry from a psl match.

        sequences are on forward strand starting at
        query_from and sbjct_from, respectively.
        """

        (self.mNMatches,
         self.mNMismatches,
         self.mNRepMatches,
         self.mNns,
         self.mQueryNGapsCounts,
         self.mQueryNGapsBases,
         self.mSbjctNGapsCounts,
         self.mSbjctNGapsBases,
         self.strand,
         self.mQueryId,
         self.mQueryLength,
         self.mQueryFrom,
         self.mQueryTo,
         self.mSbjctId,
         self.mSbjctLength,
         self.mSbjctFrom,
         self.mSbjctTo,
         self.mNBlocks,
         self.mBlockSizes,
         self.mQueryBlockStarts,
         self.mSbjctBlockStarts) = (other.mNMatches,
                                    other.mNMismatches,
                                    other.mNRepMatches,
                                    other.mNns,
                                    other.mQueryNGapsCounts,
                                    other.mQueryNGapsBases,
                                    other.mSbjctNGapsCounts,
                                    other.mSbjctNGapsBases,
                                    other.strand,
                                    other.mQueryId,
                                    other.mQueryLength,
                                    other.mQueryFrom,
                                    other.mQueryTo,
                                    other.mSbjctId,
                                    other.mSbjctLength,
                                    other.mSbjctFrom,
                                    other.mSbjctTo,
                                    other.mNBlocks,
                                    other.mBlockSizes,
                                    other.mQueryBlockStarts,
                                    other.mSbjctBlockStarts)

        if self.strand == "-":
            query_sequence = query_sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

        self.mQuerySequence = []
        self.mSbjctSequence = []
        offset_query = min(self.mQueryBlockStarts)
        offset_sbjct = min(self.mSbjctBlockStarts)

        for start_query, start_sbjct, size in other.getBlocks():
            self.mQuerySequence.append(
                query_sequence[start_query - offset_query:start_query - offset_query + size])
            self.mSbjctSequence.append(
                sbjct_sequence[start_sbjct - offset_sbjct:start_sbjct - offset_sbjct + size])


def _iterate(infile):
    """iterator over psl output.

    The header is optional.
    """

    format = "ps3"

    while 1:
        line = infile.readline()
        if not line:
            raise StopIteration

        if line[0] == "#":
            continue

        if line.startswith("match"):
            continue

        if line.startswith("psLayout version 3"):
            format = "ps3"
            for x in range(4):
                infile.readline()
            continue

        match = Match()
        match.fromTable(line[:-1].split("\t"))
        yield match


class BlatIterator:

    def __init__(self, f, *args, **kwargs):
        self.mIterator = _iterate(f)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return next(self.mIterator)
        except StopIteration:
            return None

    def next(self):
        return self.__next__()


def iterator2(infile):
    """iterate over the contents of a psl file.
    """

    while 1:
        line = infile.readline()
        if not line:
            raise StopIteration
        if line[0] == "#":
            continue
        if line.startswith("match"):
            continue
        if line.startswith("psLayout version 3"):
            for x in range(4):
                infile.readline()
            continue
        match = Match()
        # match.fromTable( line[:-1].split() )
        yield match


def iterator(infile):
    """iterate over the contents of a psl file.
    """
    while 1:
        line = infile.readline()
        if not line:
            raise StopIteration
        if line[0] == "#":
            continue
        if line.startswith("match"):
            continue
        if line.startswith("psLayout version 3"):
            for x in range(4):
                infile.readline()
            continue
        match = Match()
        match.fromTable(line[:-1].split())
        yield match


def iterator_pslx(infile):
    """iterate over the contents of a pslx file.
    """
    while 1:
        line = infile.readline()
        if not line:
            raise StopIteration
        if line[0] == "#":
            continue
        if line.startswith("match"):
            continue
        if line.startswith("psLayout version 3"):
            for x in range(4):
                infile.readline()
            continue
        match = MatchPSLX()
        match.fromTable(line[:-1].split())
        yield match


def iterator_target_overlap(infile, merge_distance):
    '''iterate over psl formatted infile and return
    blocks of target overlapping alignments.'''

    last_sbjct_id = None
    start = None
    end = None
    matches = []
    processed_contigs = set()

    while 1:

        match = next(infile)

        if match is None:
            break

        if match.mSbjctId != last_sbjct_id or match.mSbjctFrom >= (end + merge_distance):
            if last_sbjct_id:
                yield matches
            matches = []
            if last_sbjct_id != match.mSbjctId and match.mSbjctId in processed_contigs:
                raise ValueError("input not sorted by target (contig,start): already encountered %s\n%s" %
                                 (match.mSbjctId, str(match)))

            processed_contigs.add(last_sbjct_id)

            last_sbjct_id = match.mSbjctId
            start = match.mSbjctFrom
            end = match.mSbjctTo

        if match.mSbjctFrom < start:
            raise ValueError("input not sorted by target (contig,start): %i < %i\n%s" %
                             (match.mSbjctFrom, start, str(match)))

        end = max(match.mSbjctTo, end)
        matches.append(match)

    if last_sbjct_id:
        yield matches


def iterator_query_overlap(infile, merge_distance):
    '''iterate over psl formatted infile and return
    blocks of target overlapping alignments.'''

    last_query_id = None
    start, end = None, None
    matches = []
    processed_contigs = set()

    while 1:

        match = next(infile)

        if match is None:
            break

        if match.mQueryId != last_query_id or match.mQueryFrom >= (end + merge_distance):
            if last_query_id:
                yield matches
            matches = []

            if last_query_id != match.mQueryId and match.mQueryId in processed_contigs:
                raise ValueError("input not sorted by query (contig,start): already encountered %s\n%s" %
                                 (match.mQueryId, str(match)))

            processed_contigs.add(last_query_id)

            last_query_id = match.mQueryId
            start, end = match.mQueryFrom, match.mQueryTo

        if match.mQueryFrom < start:
            raise ValueError("input not sorted by query (contig,start): %i < %i\n%s" %
                             (match.mQueryFrom, start, str(match)))

        end = max(match.mQueryTo, end)
        matches.append(match)

    if last_query_id:
        yield matches


def iterator_test(infile, report_step=100000):
    '''only output parseable lines from infile.'''

    ninput, noutput, nerrors = 0, 0, 0

    while 1:
        try:
            x = next(infile)
        except ParsingError as msg:
            nerrors += 1
            ninput += 1
            E.warn(str(msg))
            continue
        except StopIteration:
            break

        if not x:
            break

        ninput += 1

        if ninput % report_step == 0:
            E.info("progress: ninput=%i, noutput=%i" % (ninput, noutput))

        yield x
        noutput += 1

    E.info("iterator_test: ninput=%i, noutput=%i, nerrors=%i" %
           (ninput, noutput, nerrors))


def iterator_per_query(iterator_psl):
    """iterate over the contents of a psl file per query
    """

    data = collections.defaultdict(list)

    for match in iterator_psl:
        data[match.mQueryId].append(match)

    for x in list(data.values()):
        yield x

    raise StopIteration

FIELDS = ("matches",
          "misMatches",
          "repMatches",
          "nCount",
          "qNumInsert",
          "qBaseInsert",
          "tNumInsert",
          "tBaseInsert",
          "strand",
          "qName",
          "qSize",
          "qStart",
          "qEnd",
          "tName",
          "tSize",
          "tStart",
          "tEnd",
          "blockCount",
          "blockSizes",
          "qStarts",
          "tStarts")


def addAlignments(matches, shift=0, by_query=False):
    """building a genome to query alignment for all matches

    The genome alignment is shifted by *shift*.
    """

    if by_query:
        for match in matches:
            if not hasattr(match, "mMapQuery2Target"):
                map_query2target = match.getMapQuery2Target()
                if shift:
                    map_query2target.moveAlignment(shift, 0)
                match.mMapQuery2Target = map_query2target
    else:
        for match in matches:
            if not hasattr(match, "mMapTarget2Query"):
                map_target2query = match.getMapTarget2Query()
                if shift:
                    map_target2query.moveAlignment(shift, 0)
                match.mMapTarget2Query = map_target2query


def getComponents(matches, max_distance=0, min_overlap=0, by_query=False):
    """return overlapping matches.

    max_distance
       allow reads to be joined if they are # residues apart.
       Adjacent reads are 1 residue apart, overlapping reads are 0 residues
       apart
    min_overlap
       require at least # residues to be overlapping

    """

    addAlignments(matches, by_query=by_query)

    components = Components.IComponents()

    for x in range(0, len(matches)):
        components.add(x, x)

    if min_overlap > 0 and max_distance > 0:
        raise ValueError(
            "both min_overlap (%i) and max_distance (%i) > 0" % (min_overlap, max_distance))

    if by_query:
        if min_overlap > 0:
            f = lambda x, y: alignlib_lite.py_getAlignmentOverlap(
                matches[x].mMapQuery2Target,
                matches[y].mMapQuery2Target,
                alignlib_lite.py_RR) >= min_overlap
        else:
            f = lambda x, y: alignlib_lite.py_getAlignmentShortestDistance(
                matches[x].mMapQuery2Target,
                matches[y].mMapQuery2Target,
                alignlib_lite.py_RR) <= max_distance
    else:
        if min_overlap > 0:
            f = lambda x, y: alignlib_lite.py_getAlignmentOverlap(
                matches[x].mMapTarget2Query,
                matches[y].mMapTarget2Query,
                alignlib_lite.py_RR) >= min_overlap
        else:
            f = lambda x, y: alignlib_lite.py_getAlignmentShortestDistance(
                matches[x].mMapTarget2Query,
                matches[y].mMapTarget2Query,
                alignlib_lite.py_RR) <= max_distance

    for x in range(len(matches)):
        for y in range(0, x):
            if f(x, y):
                components.add(x, y)

    return components.getComponents()
