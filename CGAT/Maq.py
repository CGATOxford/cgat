'''
Maq.py - Working with Maq formatted files
=========================================

This module provides a parser for the output of the maq mapper.
As this mapper is out of use, the contents here are of historical
interest only.

Reference
---------

'''


class Error(Exception):

    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)

    def _get_message(self, message):
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


class Match:

    """a maq match.
    """

    def __init__(self):
        self.mName = None
        self.contig = None
        self.start = 0
        self.strand = "+"
        self.mPairInsert = 0
        self.mPairFlag = 0
        self.mQuality = 0
        self.mQualitySingle = 0
        self.mQualityAlternative = 0
        self.mNMismatches = 0
        self.mMismatchQuality = 0
        self.mNMismatch0 = 0
        self.mNMismatch1 = 0
        self.mLength = 0

    def __str__(self):
        return "\t".join(map(str, (self.mName,
                                   self.contig,
                                   self.start + 1,
                                   self.strand,
                                   self.mPairInsert,
                                   self.mPairFlag,
                                   self.mQuality,
                                   self.mQualitySingle,
                                   self.mQualityAlternative,
                                   self.mNMismatches,
                                   self.mMismatchQuality,
                                   self.mNMismatch0,
                                   self.mNMismatch1,
                                   self.mLength)))

    def fromTable(self, data):

        if len(data) == 14:
            (self.mName,
             self.contig,
             self.start,
             self.strand,
             self.mPairInsert,
             self.mPairFlag,
             self.mQuality,
             self.mQualitySingle,
             self.mQualityAlternative,
             self.mNMismatches,
             self.mMismatchQuality,
             self.mNMismatches0,
             self.mNMismatches1,
             self.mLength) = data

        (self.start,
         self.mPairInsert,
         self.mPairFlag,
         self.mQuality,
         self.mQualitySingle,
         self.mQualityAlternative,
         self.mNMismatches,
         self.mMismatchQuality,
         self.mNMismatches0,
         self.mNMismatches1,
         self.mLength) = list(map(int, (
             self.start,
             self.mPairInsert,
             self.mPairFlag,
             self.mQuality,
             self.mQualitySingle,
             self.mQualityAlternative,
             self.mNMismatches,
             self.mMismatchQuality,
             self.mNMismatches0,
             self.mNMismatches1,
             self.mLength)))
        # maq starts counting at 1
        self.start -= 1


def iterator(infile):
    """iterate over the contents of a maq file.
    """
    while 1:
        line = infile.readline()
        if not line:
            raise StopIteration
        if line[0] == "#":
            continue
        if line.startswith("read\tcontig"):
            continue
        match = Match()
        match.fromTable(line[:-1].split())
        yield match
