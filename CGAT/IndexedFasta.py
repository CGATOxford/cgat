'''IndexedFasta.py - fast random access in fasta files
===================================================

This module provides fast random access to :term:`fasta` formatted
files that have been previously indexed. The indexing can be done
either through the samtools faidx tool (accessible through pysam_) or
using the in-house methods implemented in this module.

The main class is :class:`IndexedFasta`. This is a factory function
that provides transparent access to both samtools or CGAT indexed
fasta files.  The basic usage to retrieve the sequence spanning the
region chr12:10,000-10,100 is::

   from IndexedFasta import IndexedFasta
   fasta = IndexedFasta("hg19")
   fasta.getSequence("chr12", "+", 10000, 10100)

To index a file, use the :mod:`scripts/index_fasta` command line utility or the
:func:`createDatabase` function::

   > python index_fasta.py hg19 chr*.fa

This module has some useful utility functions:

:func:`splitFasta`
   split a :term:`fasta` formatted file into smaller pieces.

:func:`parseCoordinates`
   parse a coordinate string in various formats

but otherwise the module contains a multitude of additional functions that are
only of internal use.

Reference
----------

'''
import os
import sys
import array
import string
import re
import struct
import math
import tarfile
import random
import zlib
import gzip
import tempfile
import io
from CGAT import Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
from CGAT.AString import AString
import pysam
from future.moves import dbm
from six import StringIO

IS_PY3 = sys.version_info.major >= 3


class Uncompressor:

    def __init__(self, filename, unmangler):
        self.mFile = open(filename, "rb")
        self.mUnMangler = unmangler

    def read(self, block_size, indices, start, end):
        """read an uncompressed block from start:end.

        The compressed chunk starts at first_pos.
        NOTE: This is poorly implemented - do better.
        """

        # skip over uncompressed blocks
        d = int(math.floor(float(start) / block_size))
        r = start % block_size
        assert(d < len(indices))
        self.mFile.seek(indices[d])

        # read x bytes of compressed data, at least one full chunk.
        nchunks = int(math.ceil(float((r + end - start)) / block_size))

        fragments = []
        for x in range(d, d + nchunks):
            s = self.mFile.read(indices[x + 1] - indices[x])
            fragments.append(self.mUnMangler(s))
        u = "".join(fragments)

        assert len(u) >= end - start, \
            "fragment smaller than requested size: %i > %i-%i=%i" %\
            (len(u), end, start, end - start)

        return u[r:r + end - start]


def writeFragments(outfile_fasta,
                   outfile_index,
                   fragments,
                   mangler, size,
                   write_all=False):
    """write mangled fragments to *outfile_fasta* in chunks of *size*
    updating *outfile_index*.

    returns part of last fragment that has not been written and is
    less than *size* and the number of fragments output.

    If *write_all* is True, all of the fragments are written to
    the file and the last file position is added to *outfile_index*
    as well.
    """

    s = "".join(fragments)
    rest = len(s) % size
    if len(s) > size:
        for x in range(0, len(s) - rest, size):
            outfile_index.write("\t%i" % outfile_fasta.tell())
            outfile_fasta.write(mangler(s[x:x + size]))

    if rest:
        if write_all:
            outfile_index.write("\t%i" % outfile_fasta.tell())
            outfile_fasta.write(mangler(s[-rest:]))
            outfile_index.write("\t%i" % outfile_fasta.tell())
            return ""
        else:
            return s[-rest:]
    else:
        return ""


def gzip_mangler(s):

    xfile = StringIO()
    gzipfile = gzip.GzipFile(fileobj=xfile, mode="wb")
    gzipfile.write(s)
    gzipfile.close()

    m = xfile.getvalue()
    xfile.close()
    # write identifier

    return m


def gzip_demangler(s):
    if sys.version_info.major >= 3:
        gzipfile = io.TextIOWrapper(
            gzip.GzipFile(fileobj=io.BytesIO(s), mode="rb"),
            encoding="ascii")
    else:
        gzipfile = gzip.GzipFile(fileobj=io.BytesIO(s), mode="rb")
    m = gzipfile.readline()

    return m


class Translator:

    """translate a sequence."""

    def __init__(self):
        self.mRegEx = re.compile(" +")

    def __call__(self, sequence):
        return "".join(self.mMapScore2Char[self.mMapScore2Score[int(x)]]
                       for x in self.mRegEx.split(sequence.strip()))

    def translate(self, sequence):
        raise NotImplementedError("translate not implemented.")


class TranslatorPhred(Translator):

    """translate phred quality scores."""

    def __init__(self, *args, **kwargs):
        Translator.__init__(self, *args, **kwargs)
        self.mMapScore2Char = [chr(33 + x) for x in range(0, 93)]
        self.mMapScore2Score = list(range(0, 93))

    def translate(self, sequence):
        return array.array("I", (ord(x) - 33 for x in sequence))


class TranslatorSolexa(Translator):

    """translate solexa quality scores."""

    def __init__(self, *args, **kwargs):
        Translator.__init__(self, *args, **kwargs)
        self.mMapScore2Char = [chr(64 + x) for x in range(0, 128)]
        self.mMapScore2Score = [
            int(10.0 * math.log(1.0 + 10 ** (x / 10.0)) / math.log(10) + .499)
            for x in range(-64, 65)]

    def translate(self, sequence):
        raise NotImplementedError("translate not implemented.")
        # return array.array("i", (ord(x) - 64))


class TranslatorRange200(Translator):

    """translate pcap quality scores.

    For example for PCAP scores.

    These scores range from 0 to 100 and are the
    "a weighted sum of input base quality values
    (Huang and Madan 1999)

    The numerical values from 0 to 200 are stored
    as values form 33 to 233
    "
    """

    def __init__(self, *args, **kwargs):
        Translator.__init__(self, *args, **kwargs)
        self.mMapScore2Char = [chr(33 + x) for x in range(0, 200)]

    def __call__(self, sequence):
        try:
            return "".join(self.mMapScore2Char[int(x)]
                           for x in self.mRegEx.split(sequence.strip()))
        except ValueError as msg:
            raise ValueError(msg + " parsing error in fragment: %s" % sequence)

    def translate(self, sequence):
        return array.array("I", (ord(x) - 33 for x in sequence))


class TranslatorBytes(Translator):

    """output binary values as bytes permitting values from 0 to 255

    Note the resulting file will not be iterable as newline is not
    a record-separator any more.
    """

    def __init__(self, *args, **kwargs):
        Translator.__init__(self, *args, **kwargs)

    def __call__(self, sequence):
        try:
            return "".join(chr(int(x)) for x in
                           self.mRegEx.split(sequence.strip()))
        except ValueError as msg:
            print("parsing error in line: %s" % sequence)
            print("message=%s" % str(msg))
            return ""

    def translate(self, sequence):
        return array.array("I", (ord(x) for x in sequence))


class MultipleFastaIterator:

    def __init__(self,
                 filenames,
                 regex_identifier=None,
                 format="auto"):

        if isinstance(filenames, str):
            self.filenames = [filenames]
        else:
            self.filenames = filenames

        self.regexIdentifier = regex_identifier
        self.iterator = self._iterate()
        self.format = format

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return next(self.iterator)
        except StopIteration:
            return None

    def next(self):
        try:
            return next(self.iterator)
        except StopIteration:
            return None

    def _iterate(self):
        """iterate over muliple files."""

        def _iter(infile):

            identifier = None
            is_new = False

            for line in infile:
                if line.startswith("#"):
                    continue
                if line.startswith(">"):

                    if self.regexIdentifier:
                        try:
                            identifier = re.search(
                                self.regexIdentifier, line[1:-1]).groups()[0]
                        except AttributeError:
                            raise ValueError(
                                "could not parse identifier from line %s "
                                "- check the input" % line[1:-1])
                    else:
                        identifier = re.split("\s", line[1:-1])[0]
                    is_new = True
                else:
                    if not identifier:
                        raise ValueError(
                            "refusing to emit sequence without identifier "
                            "- check the input")
                    yield is_new, identifier, line.strip()
                    is_new = False

        for filename in self.filenames:
            if self.format == "tar.gz" or self.format == "tar" or \
               (self.format == "auto" and filename.endswith("tar.gz")):
                if filename == "-":
                    if IS_PY3:
                        tf = tarfile.open(fileobj=sys.stdin.buffer, mode="r|*")
                    else:
                        tf = tarfile.open(fileobj=sys.stdin, mode="r|*")
                else:
                    tf = tarfile.open(filename, mode="r")
                for f in tf:
                    b, ext = os.path.splitext(f.name)
                    if ext.lower() in (".fasta", ".fa"):
                        E.info("extracting %s" % f.name)
                        if sys.version_info.major >= 3:
                            infile = io.TextIOWrapper(
                                tf.extractfile(f),
                                encoding="ascii")
                        else:
                            infile = tf.extractfile(f)
                        for x in _iter(infile):
                            yield x
                    else:
                        E.info("skipping %s" % f.name)

                if tf != sys.stdin:
                    tf.close()
                continue
            elif self.format == "fasta.gz" or (self.format == "auto" and
                                               filename.endswith(".gz")):
                infile = IOTools.openFile(filename, "r")
            elif filename == "-":
                infile = sys.stdin
            else:
                infile = IOTools.openFile(filename, "r")

            for x in _iter(infile):
                yield x
            if filename != "-":
                infile.close()

        raise StopIteration


def createDatabase(db, iterator,
                   force=False,
                   synonyms=None,
                   compression=None,
                   random_access_points=None,
                   regex_identifier=None,
                   clean_sequence=False,
                   ignore_duplicates=False,
                   allow_duplicates=False,
                   translator=None):
    """index files in filenames to create database.

    Two new files are created - db.fasta and db_name.idx

    If compression is enabled, provide random access points
    every # bytes.

    Dictzip is treated as an uncompressed file.

    regex_identifier: pattern to extract identifier from description line.
    If None, the part until the first white-space character is used.

    translator: specify a translator
    """

    if db.endswith(".fasta"):
        db = db[:-len(".fasta")]

    if compression:
        if compression == "lzo":
            import lzo

            def lzo_mangler(s):
                return lzo.compress(s, 9)
            mangler = lzo_mangler
            db_name = db + ".lzo"
            write_chunks = True
        elif compression == "zlib":
            def zlib_mangler(s):
                return zlib.compress(s, 9)
            mangler = zlib_mangler
            db_name = db + ".zlib"
            write_chunks = True
        elif compression == "gzip":
            mangler = gzip_mangler
            db_name = db + ".gz"
            write_chunks = True
        elif compression == "dictzip":
            from . import dictzip

            def mangler(x):
                return x

            db_name = db + ".dz"
            write_chunks = False
        elif compression == "bzip2":
            import bz2

            def bzip_mangler(x):
                return bz2.compress(x, 9)

            mangler = bzip_mangler
            db_name = db + ".bz2"
            write_chunks = True
        elif compression == "debug":
            def mangler(x):
                return x
            db_name = db + ".debug"
            write_chunks = True
        elif compression == "rle":
            from . import RLE
            mangler = RLE.compress
            db_name = db + ".rle"
            write_chunks = True
        else:
            raise ValueError("unknown compression library: %s" % compression)

        index_name = db + ".cdx"

        if write_chunks and random_access_points is None \
           or random_access_points <= 0:
            raise ValueError("specify chunksize in --random-access-points")

    else:
        def mangler(x):
            return x
        db_name = db + ".fasta"
        write_chunks = False
        index_name = db + ".idx"

    if os.path.exists(db_name) and not force:
        raise ValueError("database %s already exists." % db_name)

    if os.path.exists(index_name) and not force:
        raise ValueError("database index %s already exists." % index_name)

    outfile_index = open(index_name, "w")
    if compression == "dictzip":
        if random_access_points is None or random_access_points <= 0:
            raise ValueError(
                "specify dictzip chunksize in --random-access-points")
        outfile_fasta = dictzip.open(
            db_name, "wb", buffersize=1000000, chunksize=random_access_points)
        compression = None
    else:
        outfile_fasta = open(db_name, "w")

    identifiers = {}
    lsequence = 0
    identifier_pos, sequence_pos = 0, 0

    if sys.version_info.major >= 3:
        translation = str.maketrans("xX", "nN")
    else:
        translation = string.maketrans("xX", "nN")

    fragments = []
    lfragment = 0

    last_identifier = None

    while 1:

        try:
            result = next(iterator)
        except StopIteration:
            break

        if not result:
            break

        is_new, identifier, fragment = result

        if is_new:
            # check for duplicate identifiers
            if identifier in identifiers:
                if ignore_duplicates:
                    raise ValueError("ignore duplicates not implemented")
                elif allow_duplicates:
                    # the current implementation will fail if the same
                    # identifiers
                    # are directly succeeding each other
                    # better: add return to iterator that indicates a new
                    # identifier
                    out_identifier = identifier + \
                        "_%i" % (identifiers[identifier])
                    identifiers[identifier] += 1
                    identifiers[out_identifier] = 1
                else:
                    raise ValueError("%s occurs more than once" %
                                     (identifier,))
            else:
                identifiers[identifier] = 1
                out_identifier = identifier

            if last_identifier:
                if write_chunks:
                    writeFragments(outfile_fasta, outfile_index,
                                   fragments, mangler,
                                   size=random_access_points,
                                   write_all=True)

                    fragments = []
                    lfragment = 0
                else:
                    outfile_fasta.write("\n")

                outfile_index.write("\t%i\n" % lsequence)

            identifier_pos = outfile_fasta.tell()
            outfile_fasta.write(mangler(">%s\n" % out_identifier))
            sequence_pos = outfile_fasta.tell()

            outfile_index.write("%s\t%i" % (out_identifier,
                                            identifier_pos))
            if write_chunks:
                outfile_index.write("\t%i" % random_access_points)
            else:
                outfile_index.write("\t%i" % sequence_pos)

            fragments = []
            lsequence = 0
            last_identifier = identifier

        if translator:
            s = translator(fragment)
        else:
            s = re.sub("\s", "", fragment.strip())
            if clean_sequence:
                s = s.translate(translation)

        lsequence += len(s)

        if write_chunks:
            fragments.append(s)
            lfragment += len(s)
            if lfragment > random_access_points:
                rest = writeFragments(outfile_fasta,
                                      outfile_index,
                                      fragments,
                                      mangler,
                                      size=random_access_points,
                                      write_all=False)
                fragments = [rest]
                lfragment = len(rest)
        else:
            outfile_fasta.write(mangler(s))

    if write_chunks:
        writeFragments(outfile_fasta, outfile_index, fragments, mangler,
                       size=random_access_points, write_all=True)
    else:
        outfile_fasta.write("\n")

    outfile_index.write("\t%i\n" % lsequence)

    # add synonyms for the table
    if synonyms:
        for key, vals in list(synonyms.items()):
            for val in vals:
                outfile_index.write("%s\t%s\n" % (key, val))

NAME_MAP = {
    'uncompressed': ('fasta', 'idx', False),
    'lzo': ('lzo',   'cdx', True),
    'dictzip': ('dz',    'idx', False),
    'zlib': ('zlib',  'cdx', True),
    'gzip': ('gz',  'cdx', True),
    'bzip2': ('bz2',   'cdx', True),
    'debug': ('debug', 'cdx', True),
}

PREFERENCES = (
    'uncompressed', 'lzo', 'dictzip', 'zlib', 'gzip', 'bzip2', 'debug')


class CGATIndexedFasta:

    """an indexed fasta file."""

    def __init__(self, dbname):

        if dbname.endswith(".fasta"):
            dbname = dbname[:-len(".fasta")]

        for x in PREFERENCES:
            d = "%s.%s" % (dbname, NAME_MAP[x][0])
            i = "%s.%s" % (dbname, NAME_MAP[x][1])
            if os.path.exists(d) and os.path.exists(i):
                self.mMethod = x
                self.mDbname = d
                self.mNameIndex = i
                self.mNoSeek = NAME_MAP[x][2]
                break
        else:
            raise KeyError("unknown database %s" % dbname)

        self.mIsLoaded = False
        self.mSynonyms = {}
        self.mConverter = None
        self.mIndex = {}
        self.mTranslator = None

    def __len__(self):
        """return the number of sequences in fasta file."""
        if not self.mIsLoaded:
            self._loadIndex()
        return len(self.mIndex)

    def __contains__(self, contig):
        if not self.mIsLoaded:
            self._loadIndex()
        return contig in self.mIndex or contig in self.mSynonyms

    def __getitem__(self, key):
        """return full length sequence."""
        return self.getSequence(key, "+", 0, 0, as_array=True)

    def _loadIndex(self, compress=False):
        """load complete index into memory.

        if compress is set to true, the index will not be loaded,
        but a compressed index will be created instead.
        """
        if self.mMethod == "uncompressed":
            self.mDatabaseFile = open(self.mDbname, "r")
        elif self.mMethod == "dictzip":
            from . import dictzip
            self.mDatabaseFile = dictzip.GzipFile(self.mDbname)
        elif self.mMethod == "lzo":
            import lzo
            self.mDatabaseFile = Uncompressor(self.mDbname, lzo.decompress)
        elif self.mMethod == "gzip":
            self.mDatabaseFile = Uncompressor(self.mDbname, gzip_demangler)
        elif self.mMethod == "zlib":
            self.mDatabaseFile = Uncompressor(self.mDbname, zlib.decompress)
        elif self.mMethod == "bzip2":
            import bz2
            self.mDatabaseFile = Uncompressor(self.mDbname, bz2.decompress)
        elif self.mMethod == "debug":
            self.mDatabaseFile = Uncompressor(
                self.mDbname + ".debug", lambda x: x)

        filename_index = self.mNameIndex + ".dbm"

        if compress:
            # if os.path.exists(filename_index):
            #     raise OSError("file %s already exists" % filename_index)
            self.mIndex = dbm.open(filename_index, "n")
        elif os.path.exists(filename_index):
            self.mIndex = dbm.open(filename_index, "r")
            if len(self.mIndex) == 0:
                raise ValueError(
                    "length of index is 0, possibly a python version problem")
            self.mIsLoaded = True
            return
        else:
            self.mIndex = {}

        for line in open(self.mNameIndex, "r"):

            data = line[:-1].split("\t")

            if len(data) == 2:
                # ignore synonyms of non-existent contigs
                identifier = data[1]
                if data[0] not in self.mIndex:
                    continue
                self.mSynonyms[identifier] = data[0]
            else:
                # index with random access points
                if len(data) > 4:
                    (identifier, pos_id, block_size, lsequence) = data[
                        0], int(data[1]), int(data[2]), int(data[-1])
                    points = list(map(int, data[3:-1]))
                    self.mIndex[identifier] = (
                        pos_id, block_size, lsequence, points)
                else:
                    (identifier, pos_id, pos_seq, lsequence) = data[
                        0], int(data[1]), int(data[2]), int(data[-1])
                    self.mIndex[identifier] = struct.pack(
                        "QQi", pos_id, pos_seq, lsequence)

        self._addSynonyms()
        self.mIsLoaded = True

    def _addSynonyms(self):
        '''add synonyms to indices.
        '''

        # Treat common cases of naming incompatibilites like chr1 = 1.
        # Truncate or add known prefixes
        def _add(src, target):
            for p in ("chr", "contig", "scaffold", "Chr"):
                if src.startswith(p):
                    k = src[len(p):]
                    if k not in self.mIndex:
                        self.mSynonyms[k] = target
                    # add lower/upper-case version
                    k = src[0].upper() + src[1:]
                    if k not in self.mIndex:
                        self.mSynonyms[k] = target
                    k = src[0].lower() + src[1:]
                    if k not in self.mIndex:
                        self.mSynonyms[k] = target
                    break
            else:
                for p in ("chr", "contig", "scaffold"):
                    k = "%s%s" % (p, src)
                    if k not in self.mIndex:
                        self.mSynonyms[k] = target

        k = list(self.mSynonyms.items())

        # fix the ambiguity between chrMT and chrM between UCSC and ENSEMBL
        if "chrM" in self.mIndex and "MT" not in self.mIndex:
            self.mSynonyms["MT"] = "chrM"
        elif "chrM" not in self.mSynonyms and "MT" in self.mSynonyms:
            self.mSynonyms["chrM"] = "MT"

        # this does the same thing for "chrMito" used in the yeast
        # genome
        if "chrM" in self.mIndex and "chrMito" not in self.mIndex:
            self.mSynonyms["chrMito"] = "chrM"

        for key, val in k:
            _add(key, val)

        # add pointers to self
        for key in list(self.mIndex.keys()):
            _add(key, key)

    def setTranslator(self, translator=None):
        """set the :class:`Translator` to use."""
        self.mTranslator = translator

    def getDatabaseName(self):
        """returns the name of the database."""
        return self.mDbname

    def getToken(self, contig):
        """check if token is in index."""
        if not self.mIsLoaded:
            self._loadIndex()

        # First Convert the contig to string to compare string with strings
        contig = str(contig)

        if contig in self.mSynonyms:
            contig = self.mSynonyms[contig]

        if contig not in self.mIndex:
            raise KeyError("%s not in index" % contig)

        return contig

    def getLength(self, contig):
        """return sequence length for sbjct_token."""
        if not self.mIsLoaded:
            self._loadIndex()
        data = self.mIndex[self.getToken(contig)]
        try:
            pos_id, pos_seq, lcontig = struct.unpack("QQi", data)
        except struct.error:
            pos_id, pos_seq, lcontig, points = data
        return lcontig

    def getLengths(self):
        """return all sequence lengths."""
        if not self.mIsLoaded:
            self._loadIndex()
        results = []
        for contig in list(self.mIndex.values()):
            data = self.mIndex[self.getToken(contig)]
            try:
                pos_id, pos_seq, lcontig = struct.unpack("QQi", data)
            except struct.error:
                pos_id, pos_seq, lcontig, points = data
            results.append(lcontig)
        return results

    def compressIndex(self):
        """compress index.
        Creates a database interface to an index.
        """
        self._loadIndex(compress=True)

    def getContigs(self):
        """return a list of contigs (no synonyms)."""
        if not self.mIsLoaded:
            self._loadIndex()
        return list(self.mIndex.keys())

    def getContigSizes(self, with_synonyms=True):
        """return hash with contig sizes including synonyms."""
        if not self.mIsLoaded:
            self._loadIndex()

        contig_sizes = {}
        for key, val in list(self.mIndex.items()):
            contig_sizes[key] = self.getLength(key)

        if with_synonyms:
            for key, val in list(self.mSynonyms.items()):
                contig_sizes[key] = self.getLength(val)

        return contig_sizes

    def setConverter(self, converter):
        """set converter from coordinate system to 0-based, both strand,
        open/closed coordinate system.

        """
        self.mConverter = converter

    def getSequence(self,
                    contig,
                    strand="+",
                    start=0,
                    end=0,
                    converter=None,
                    as_array=False):
        """get a genomic fragment.

        A genomic fragment is identified by the coordinates
        contig, strand, start, end.

        The converter function supplied translated these coordinates
        into 0-based coordinates. By default, start and end are assumed
        to be pythonic coordinates and are forward/reverse coordinates.

        If as_array is set to true, return the AString object. This might
        be beneficial for large sequence chunks. If as_array is set to False,
        return a python string.
        """

        contig = self.getToken(contig)

        data = self.mIndex[contig]
        # dummy is
        # -> pos_seq for seekable streams
        # -> block_size for unseekable streams
        try:
            pos_id, dummy, lsequence = struct.unpack("QQi", data)
        except (struct.error, TypeError):
            pos_id, dummy, lsequence, points = data

        pos_seq = dummy
        block_size = dummy

        if end == 0:
            end = lsequence

        if end > lsequence:
            raise ValueError(
                "3' coordinate on %s out of bounds: %i > %i" %
                (contig, end, lsequence))

        if start < 0:
            raise ValueError(
                "5' coordinate on %s out of bounds: %i < 0" % (contig, start))

        if converter:
            first_pos, last_pos = converter(start, end,
                                            str(strand) in ("+", "1"),
                                            lsequence)
        elif self.mConverter:
            first_pos, last_pos = self.mConverter(start, end,
                                                  str(strand) in ("+", "1"),
                                                  lsequence)
        else:
            first_pos, last_pos = start, end
            if str(strand) in ("-", "0", "-1"):
                first_pos, last_pos = lsequence - \
                    last_pos, lsequence - first_pos

        if first_pos == last_pos:
            return ""

        assert first_pos < last_pos, \
            "first position %i is larger than last position %i " % \
            (first_pos, last_pos)

        p = AString()

        if self.mNoSeek:
            # read directly from position
            p.fromstring(
                self.mDatabaseFile.read(block_size, data[3],
                                        first_pos, last_pos))
        else:
            first_pos += pos_seq
            last_pos += pos_seq

            self.mDatabaseFile.seek(first_pos)
            p.fromstring(self.mDatabaseFile.read(last_pos - first_pos))

        if str(strand) in ("-", "0", "-1"):
            p = AString(Genomics.complement(str(p)))

        if self.mTranslator:
            return self.mTranslator.translate(p)
        elif as_array:
            return p
        else:
            if IS_PY3:
                return p.tostring().decode("ascii")
            else:
                return p.tostring()

    def getRandomCoordinates(self, size):
        """returns coordinates for a random fragment of size #.

        The coordinates are forward/reverse.

        Default sampling mode:

        Each residue has the same probability of being
        in a fragment. Thus, the fragment can be smaller than
        size due to contig boundaries.
        """
        if not self.mIsLoaded:
            self._loadIndex()

        token = random.choice(list(self.mIndex.keys()))
        strand = random.choice(("+", "-"))
        data = self.mIndex[token]
        try:
            pos_id, pos_seq, lcontig = struct.unpack("QQi", data)
        except (struct.error, TypeError):
            pos_id, pos_seq, lcontig, points = data

        start, end = 0, 0
        if size >= lcontig:
            end = lcontig
        else:
            while end - start == 0:
                rpos = random.randint(0, lcontig)
                if random.choice(("True", "False")):
                    start = rpos
                    end = min(rpos + size, lcontig)
                else:
                    start = max(0, rpos - size)
                    end = rpos

        return token, strand, start, end


class PysamIndexedFasta(CGATIndexedFasta):

    '''interface a  pysam/samtools indexed fasta file with the
    CGATIndexedFasta API.'''

    def __init__(self, dbname):

        # open database file and truncate
        if os.path.exists(dbname) and dbname.endswith(".fa"):
            self.mDatabaseFile = pysam.Fastafile(dbname)
            dbname = dbname[:-len(".fa")]
        elif os.path.exists(dbname) and dbname.endswith(".fasta"):
            self.mDatabaseFile = pysam.Fastafile(dbname)
            dbname = dbname[:-len(".fasta")]
        elif os.path.exists(dbname + ".fa"):
            self.mDatabaseFile = pysam.Fastafile(dbname + ".fa")
        elif os.path.exists(dbname + ".fasta"):
            self.mDatabaseFile = pysam.Fastafile(dbname + ".fasta")

        if os.path.exists(dbname + ".fai"):
            self.mNameIndex = dbname + ".fai"
        elif os.path.exists(dbname + ".fa.fai"):
            self.mNameIndex = dbname + ".fa.fai"
        elif os.path.exists(dbname + ".fasta.fai"):
            self.mNameIndex = dbname + ".fasta.fai"
        else:
            raise ValueError("PysamIndexedFasta requires pre-built index")

        self.mMethod = "faidx"
        self.mDbname = dbname
        self.mNoSeek = False
        self.mIsLoaded = False
        self.mSynonyms = {}
        self.mConverter = None
        self.mIndex = {}
        self.mTranslator = None

    def _loadIndex(self, compress=False):
        '''load index into memory.'''

        with open(self.mNameIndex, "r") as infile:

            for line in infile:
                contig, lsequence, offset, line_blen, line_len = line[
                    :-1].split()
                self.mIndex[contig] = struct.pack(
                    "QQi", int(offset), int(offset), int(lsequence))

        self.mIsLoaded = True
        self._addSynonyms()

    def getSequence(self,
                    contig,
                    strand="+",
                    start=0,
                    end=0,
                    converter=None,
                    as_array=False):

        contig = self.getToken(contig)

        data = self.mIndex[contig]
        try:
            pos_id, pos_seq, lsequence = struct.unpack("QQi", data)
        except struct.error:
            pos_id, pos_seq, lsequence, points = data

        # convert to 0-based positive strand coordinates
        if converter:
            first_pos, last_pos = converter(start, end,
                                            str(strand) in ("+", "1"),
                                            lsequence)
        elif self.mConverter:
            first_pos, last_pos = self.mConverter(start, end,
                                                  str(strand) in ("+", "1"),
                                                  lsequence)
        else:
            first_pos, last_pos = start, end
            if str(strand) in ("-", "0", "-1"):
                first_pos, last_pos = lsequence - \
                    last_pos, lsequence - first_pos

        sequence = self.mDatabaseFile.fetch(contig, first_pos, last_pos)
        if str(strand) in ("-", "0", "-1"):
            try:
                # works in py2 only
                sequence = string.translate(
                    sequence[::-1],
                    string.maketrans("ACGTacgtNn", "TGCAtgcaNn"))
            except AttributeError:
                # works in py3 only
                sequence = str(sequence[::-1]).translate(
                    str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))

        return sequence


def IndexedFasta(dbname, *args, **kwargs):
    '''factory function for IndexedFasta objects.'''

    if (os.path.exists(dbname) or os.path.exists(dbname + ".fa")) \
            and (os.path.exists(dbname + ".fai") or
                 os.path.exists(dbname + ".fa.fai")):
        return PysamIndexedFasta(dbname, *args, **kwargs)
    else:
        return CGATIndexedFasta(dbname, *args, **kwargs)


def _one_forward_closed(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    x -= 1
    if not c:
        x, y = l - y, l - x
    return x, y


def _zero_forward_closed(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    y += 1
    if not c:
        x, y = l - y, l - x
    return x, y


def _one_both_closed(x, y, c=None, l=None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x - 1, y


def _zero_both_closed(x, y, c=None, l=None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x, y + 1


def _one_forward_open(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    x -= 1
    y -= 1
    if not c:
        x, y = l - y, l - x
    return x, y


def _zero_forward_open(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    if not c:
        x, y = l - y, l - x
    return x, y


def _one_both_open(x, y, c=None, l=None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x - 1, y - 1


def _zero_both_open(x, y, c=None, l=None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.

    Parameters are from, to, is_positive_strand, length of contig.
    """
    return x, y


def getConverter(format):
    """return a converter function for converting various
    coordinate schemes into 0-based, both strand, closed-open ranges.

    converter functions have the parameters
    x, y, s, l: with x and y the coordinates of
    a sequence fragment, s the strand (True is positive)
    and l being the length of the contig.

    Format is a "-" separated combination of the keywords
    "one", "zero", "forward", "both", "open", "closed":

    zero/one: zero or one-based coordinates
    forward/both: forward coordinates or forward/reverse coordinates
    open/closed: half-open intervals (pythonic) or closed intervals

    """

    data = set(format.split("-"))

    if "one" in data:
        if "forward" in data:
            if "closed" in data:
                return _one_forward_closed
            else:
                return _one_forward_open
        else:
            if "closed" in data:
                return _one_both_closed
            else:
                return _one_both_open
    else:
        if "forward" in data:
            if "closed" in data:
                return _zero_forward_closed
            else:
                return _zero_forward_open
        else:
            if "closed" in data:
                return _zero_both_closed
            else:
                return _zero_both_open


def benchmarkRandomFragment(fasta, size):
    """returns a random fragment of size."""

    contig, strand, start, end = fasta.getRandomCoordinates(size)
    s = fasta.getSequence(contig, strand, start, end)
    return s


def verify(fasta1, fasta2, num_iterations, fragment_size,
           stdout=sys.stdout, quiet=False):
    """verify two databases.

    Get segment from fasta1 and check for presence in fasta2.
    """
    if not quiet:
        stdout.write(
            "verifying %s and %s using %i random segments of length %i\n" %
            (fasta1.getDatabaseName(),
             fasta2.getDatabaseName(),
             num_iterations,
             fragment_size))
        stdout.flush()
    nerrors = 0
    for x in range(num_iterations):
        contig, strand, start, end = fasta1.getRandomCoordinates(fragment_size)
        s1 = fasta1.getSequence(contig, strand, start, end)
        s2 = fasta2.getSequence(contig, strand, start, end)
        if s1 != s2:
            if not quiet:
                stdout.write("discordant segment: %s:%s:%i:%i\n%s\n%s\n" %
                             (contig, strand, start, end, s1, s2))
            nerrors += 1
    return nerrors


def splitFasta(infile, chunk_size, dir="/tmp", pattern=None):
    """split a fasta file into a subset of files.

    If pattern is not given, random file names are chosen.
    """

    n = 0
    chunk = 0

    def _getFilename(chunk):
        if pattern:
            outname = pattern % chunk
            outfile = os.open(outname, "w")
        else:
            (outfile, outname) = tempfile.mkstemp(dir=dir)
        return (outfile, outname)

    outfile, outname = _getFilename(chunk)
    filenames = [outname]
    noutput = 0
    for line in infile:
        if line[0] == "#":
            continue
        if line[0] == ">":
            n += 1
            if n > chunk_size:
                os.close(outfile)
                n = 1
                outfile, outname = _getFilename(chunk)
                filenames.append(outname)
            noutput += 1

        os.write(outfile, line)

    os.close(outfile)

    if noutput == 0:
        os.remove(outname)
        return []

    return filenames


def parseCoordinates(s):
    '''parse a coordinate string.

    The coordinate string can be various formats, such as
    ``chr1:+:10:1000``, ``chr1:10..1000``.

    Returns
    -------
    contig : string
        The chromosome/contig.
    strand : char
        Strand. If not present, set to "+".
    start : int
        Start of interval
    end : int
        End of interval. If not present, set to start + 1.
    '''

    if ":" in s:
        d = s.strip().split(":")
        if len(d) == 4:
            contig, strand, start, end = d
        elif len(d) == 2:
            contig = d[0]
            strand = "+"
            if ".." in d[1]:
                start, end = d[1].split("..")
            else:
                start = d[1]
        else:
            raise ValueError("format not recognized")
    else:
        contig = s
        strand = "+"
        start = "0"
        # full sequence
        end = "0"

    start = int(re.sub(",", "", start))
    if end:
        end = int(re.sub(",", "", end))
    else:
        end = start + 1

    return contig, strand, start, end
