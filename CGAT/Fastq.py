'''Fastq.py - methods for dealing with fastq files
===============================================

This module provides an iterator of :term:`fastq` formatted files
(:func:`iterate`). Additional iterators allow guessing of the quality
score format (:func:`iterate_guess`) or converting them
(:func:`iterate_convert`) while iterating through a file.

:func:`guessFormat` inspects a fastq file to guess the quality score format
and :func:`getOffset` returns the numeric offset for quality score conversion
for a particular quality score format.

.. note::
   Another way to access the information in :term:`fastq` formatted
   files is through pysam_.

Reference
---------

'''

import string

from math import log

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

RANGES = {
    'sanger': (33, 75),
    'illumina-1.8': (33, 79),
    'solexa': (59, 106),
    'phred64': (64, 106),
}


class Record:
    """A record representing a :term:`fastq` formatted record.

    Attributes
    ----------
    identifier : string
       Sequence identifier
    seq : string
       Sequence
    quals : string
       String representation of quality scores.
    format : string
       Quality score format. Can be one of ``sanger``,
       ``illumina-1.8``, ``solexa`` or ``phred64``.

    """

    def __init__(self, identifier, seq, quals, format=None):
        self.identifier, self.seq, self.quals, format = (
            identifier, seq, quals, format)
        self.format = None

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals)

    def guessFormat(self):
        '''return quality score format -
        might return several if ambiguous.'''

        c = [ord(x) for x in self.quals]
        mi, ma = min(c), max(c)
        r = []
        for format, v in RANGES.items():
            m1, m2 = v
            if mi >= m1 and ma <= m2:
                r.append(format)
        return r

    def guessDataType(self):
        '''return the datatype. This is done by inspecting the
        sequence for basecalls/colorspace ints'''

        datatypes = set()

        # don't check the first character as it may be a basecall even with cs
        if sum([x in ["0", "1", "2", "3", "."]
                for x in self.seq[1:]]) == len(self.seq[1:]):
            datatypes.add("colorspace")
        elif sum([x.upper() in string.ascii_uppercase
                  for x in self.seq]) == len(self.seq):
            datatypes.add("basecalls")

        return datatypes

    def trim(self, trim3, trim5=0):
        """remove nucleotides/quality scores from the 3' and 5' ends."""
        self.seq = self.seq[trim5:-trim3]
        self.quals = self.quals[trim5:-trim3]

    def trim5(self, trim5=0):
        """remove nucleotides/quality scores from the 5' ends."""
        self.seq = self.seq[trim5:]
        self.quals = self.quals[trim5:]

    def toPhred(self):
        '''return qualities as a list of phred-scores.'''
        assert self.format is not None, "format needs to be set for conversion"
        if self.format == "sanger":
            return [ord(x) - 33 for x in self.quals]
        elif self.format == "illumina-1.8":
            return [ord(x) - 33 for x in self.quals]
        elif self.format == "solexa":
            # from -5 to 40 (i.e., can be negative)
            log10x = log(10.0) + .499
            return [int(10.0 * log(1.0 + 10 ** (ord(x) / 10.0), 10) / log10x)
                    for x in self.quals]
        elif self.format == "phred64":
            return [ord(x) - 64 for x in self.quals]

    def fromPhred(self, quals, format):
        '''set qualities from a list of phred-scores.'''
        self.format = format
        # -1 for color space fastq file
        assert len(quals) == len(self.seq) or len(quals) == len(self.seq) - 1
        if self.format == "sanger":
            self.quals = "".join([chr(33 + x) for x in quals])
        elif self.format == "illumina-1.8":
            return [chr(33 + x) for x in quals]
        elif self.format == "solexa":
            log10x = log(10.0, 10) / 10.0
            q = [int(10.0 * (log(10 ** (x * log10x) - 1.0, 10)))
                 for x in quals]
            self.quals = "".join([chr(64 + x) for x in q])
        elif self.format == "phred64":
            self.quals = "".join([chr(64 + x) for x in quals])
        elif self.format == "integer":
            self.quals = " ".join(map(str, quals))


def iterate(infile):
    '''iterate over contents of fastq file.'''

    while 1:
        line1 = infile.readline()
        if not line1:
            break
        if not line1.startswith('@'):
            raise ValueError("parsing error: expected '@' in line %s" % line1)
        line2 = infile.readline()
        line3 = infile.readline()
        if not line3.startswith('+'):
            raise ValueError("parsing error: expected '+' in line %s" % line3)
        line4 = infile.readline()
        # incomplete entry
        if not line4:
            raise ValueError("incomplete entry for %s" % line1)

        yield Record(line1[1:-1], line2[:-1], line4[:-1])


def iterate_guess(infile, max_tries=10000, guess=None):
    '''iterate over contents of fastq file.

    Guess quality format by looking at the first `max_tries` entries and
    then subsequently setting the quality score format for each entry.

    Arguments
    ---------
    infile : File
       File or file-like object to iterate over
    max_tries : int
       Number of records to examine for guessing the quality score
       format.
    guess : string
       Default format. This format will be chosen in the quality
       score format is ambiguous. The method checks if the `guess`
       is compatible with the records read so far.

    Yields
    ------
    fastq
        An object of type :class:`Record`.

    Raises
    ------
    ValueError
        If the ranges of the fastq records are not compatible,
        are incompatible with guess or are ambiguous.

    '''
    quals = set(RANGES.keys())
    cache = []
    myiter = iterate(infile)
    lengths = []
    for c, record in enumerate(myiter):
        quals.intersection_update(set(record.guessFormat()))
        if len(quals) == 0:
            raise ValueError("could not guess format - ranges incompatible.")
        if len(quals) == 1:
            break
        cache.append(record)
        lengths.append(len(record.seq))
        if c > max_tries:
            break

    if len(quals) == 1:
        ref_format = list(quals)[0]
    elif guess in quals:
        E.warn("multiple input formats possible: %s. Continuing with %s" %
               (", ".join(quals), guess))
        ref_format = guess
    elif quals.issubset(set(["solexa", "phred64"])):
        # guessFormat will call phred64 reads as phred64 AND solexa
        # if both still remain after max_tries, assume phred64
        ref_format = "phred64"
    else:
        raise ValueError(
            "could not guess format - could be one of %s." % str(quals))

    for r in cache:
        r.format = ref_format
        yield r

    for r in myiter:
        r.format = ref_format
        yield r


def iterate_convert(infile, format, max_tries=10000, guess=None):
    '''iterate over contents of fastq file.

    The quality score format is guessed and all subsequent records
    are converted to `format`.

    Arguments
    ---------
    infile : File
       File or file-like object to iterate over
    format : string
       Quality score format to convert all records into.
    max_tries : int
       Number of records to examine for guessing the quality score
       format.
    guess : string
       Default format. This format will be chosen in the quality
       score format is ambiguous. The method checks if the `guess`
       is compatible with the records read so far.

    Yields
    ------
    fastq
        An object of type :class:`Record`.

    Raises
    ------
    ValueError
        If the ranges of the fastq records are not compatible,
        are incompatible with guess or are ambiguous.


    '''

    quals = set(RANGES.keys())
    cache = []
    myiter = iterate(infile)
    for c, record in enumerate(myiter):
        quals.intersection_update(set(record.guessFormat()))

        if len(quals) == 0:
            raise ValueError("could not guess format - ranges incompatible.")
        if len(quals) == 1:
            cache.append(record)
            break
        cache.append(record)

        if c > max_tries:
            break

    if len(quals) == 1:
        ref_format = list(quals)[0]
    elif quals.issubset(set(["solexa", "phred64"])):
        # guessFormat will call phred64 reads as phred64 AND solexa
        # if both still remain after max_tries, assume phred64
        ref_format = "phred64"
    elif guess in quals:
        E.warn("multiple input formats possible: %s. Continuing with %s" %
               (", ".join(quals), guess))
        ref_format = guess
    else:
        raise ValueError(
            "could not guess format - could be one of %s. "
            "If you know the format use the --format option" % str(quals))

    for r in cache:
        r.format = ref_format
        r.fromPhred(r.toPhred(), format)
        yield r

    for r in myiter:
        r.format = ref_format
        r.fromPhred(r.toPhred(), format)
        yield r


def guessFormat(infile, max_lines=10000, raises=True):
    '''guess format of FASTQ File.

    Arguments
    ---------
    infile : File
       File or file-like object to iterate over
    max_lines : int
       Number of lines to examine for guessing the quality score
       format.
    raises : bool
       Raise ValueError if format is ambiguous

    Returns
    -------
    formats : list
       list of quality score formats compatible with the file

    Raises
    ------
    ValueError
        If the ranges of the fastq records are not compatible.

    '''

    quals = set(RANGES.keys())
    myiter = iterate(infile)
    for c, record in enumerate(myiter):
        quals.intersection_update(set(record.guessFormat()))
        if len(quals) == 0:
            raise ValueError("could not guess format - ranges incompatible.")
        if len(quals) == 1:
            break
        if c > max_lines:
            break

    if len(quals) == 1:
        return list(quals)[0]
    elif raises is False:
        if quals.issubset(set(["solexa", "phred64"])):
            # guessFormat will call phred64 reads as phred64 AND solexa
            # if both still remain after max_lines, return phred64
            return "phred64"
        else:
            return quals
    else:
        raise ValueError(
            "could not guess format - could be one of %s." % str(quals))


def guessDataType(infile, max_lines=10000, raises=True):
    '''guess datatype of FASTQ File from [colourspace, basecalls]

    Arguments
    ---------
    infile : File
    File or file-like object to iterate over
    max_lines : int
       Number of lines to examine for guessing the datatype
    raises : bool
       Raise ValueError if format is ambiguous

    Returns
    -------
    formats : list
       list of datatypes compatible with the file (should only ever be one!)

    Raises
    ------
    ValueError
        If the ranges of the fastq records are not compatible.

    '''

    datatype = set(("colorspace", "basecalls"))

    myiter = iterate(infile)
    for c, record in enumerate(myiter):
        datatype.intersection_update(record.guessDataType())

        if len(datatype) == 0:
            raise ValueError("could not guess datatype")
        if len(datatype) == 1:
            break
        if c > max_lines:
            break

    if len(datatype) == 1:
        return list(datatype)[0]
    else:
        raise ValueError(
            "could not guess datatype - could be one of %s." % str(datatype))


def getOffset(format, raises=True):
    '''returns the ASCII offset for a certain format.

    If `raises` is set a ValueError is raised if there is not a single
    offset. Otherwise, a minimum offset is returned.

    Returns
    -------
    offset : int
       The quality score offset

    '''
    if type(format) in (set, list, tuple):
        offsets = set([RANGES[f][0] for f in format])
        if len(offsets) == 1:
            return list(offsets)[0]
        elif raises is False:
            return min(offsets)
        else:
            raise ValueError(
                "inconsistent offsets for multiple formats: %s" % offsets)

    return RANGES[format][0]


def getReadLength(filename):
    '''return readlength from a fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, the result
    will be inaccurate.

    Returns
    -------
    read_length : int
    '''

    with IOTools.openFile(filename) as infile:
        record = next(iterate(infile))
        return len(record.seq)
