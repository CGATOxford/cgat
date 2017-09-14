"""SequencePairProperties.py - Computing metrics for aligned sequences
======================================================================

This module provides methods for extracting and reporting sequence
properties of aligned nucleotide sequences such as percent identity,
substitution rate, etc. Usage is the same as
:mod:`SequencePairProperties`.


Reference
---------

"""
import numpy
from CGAT import Mali as Mali
from CGAT import WrapperCodeML as WrapperCodeML


class SequencePairProperties:

    mPseudoCounts = 1

    def __init__(self):

        self.mLength = 0
        self.mNCodons = 0

    def addProperties(self, other):
        """add properties."""

        self.mLength += other.mLength
        self.mNCodons += other.mNCodons

    def updateProperties(self):
        pass

    def loadPair(self, seq1, seq2):
        """load sequence properties from a pair."""

    def __str__(self):
        return "\t".join(self.getFields())

    def __call__(self, seq1, seq2):
        self.loadPair(seq1.upper(), seq2.upper())

    def getFields(self):
        self.updateProperties()
        return []

    def getHeaders(self):
        return []

    def getHeader(self):
        return "\t".join(self.getHeaders)


class SequencePairPropertiesDistance(SequencePairProperties):

    """base class for distance estimators."""

    def __init__(self, *args, **kwargs):
        SequencePairProperties.__init__(self, *args, **kwargs)


class SequencePairPropertiesBaseML (SequencePairPropertiesDistance):

    """Counts for nucleic acid sequences.

    The first characters are ACGT.
    """

    mFormat = "%6.4f"

    def __init__(self, options, *args, **kwargs):
        SequencePairPropertiesDistance.__init__(self, *args, **kwargs)

        self.mBaseml = WrapperCodeML.BaseML()
        self.mBaseml.SetOptions(options)

        if options.loglevel >= 3:
            self.mDump = True
            self.mTest = True
        else:
            self.mDump = False
            self.mTest = False

    def loadPair(self, seq1, seq2):

        temp_mali = Mali.Mali()
        temp_mali.addSequence("seq1", 0, len(seq1), seq1)
        temp_mali.addSequence("seq2", 0, len(seq2), seq2)

        try:
            self.mResult = self.mBaseml.Run(temp_mali,
                                            tree="(seq1,seq2);",
                                            dump=self.mDump,
                                            test=self.mTest)
        except WrapperCodeML.UsageError:
            self.mResult = None

    def getHeaders(self):
        return ["distance", "lnl", "kappa", "converged"]

    def __str__(self):

        if not self.mResult:
            return "\t".join(["na"] * 4)

        if self.mResult.mKappa != "na":
            o_kappa = self.mFormat % self.mResult.mKappa
        else:
            o_kappa = self.mResult.mKappa

        if self.mResult.mCheckConvergence:
            o_converged = "0"
        else:
            o_converged = "1"

        return "\t".join((
            self.mFormat % self.mResult.mDistanceMatrix["seq1"]["seq2"],
            self.mFormat % self.mResult.mLogLikelihood,
            o_kappa,
            o_converged))


class SequencePairPropertiesCountsNa (SequencePairProperties):

    """Counts for nucleic acid sequences.

    The first characters are ACGT.
    """

    mAlphabet = "ACGT"
    mGapChar = "-"

    def __init__(self, *args, **kwargs):
        SequencePairProperties.__init__(self, *args, **kwargs)
        self.mNIdentical = 0
        self.mNAligned = 0
        self.mNDifferent = 0
        self.mNTransitions = 0
        self.mNTransversions = 0
        self.mNUnaligned1 = 0
        self.mNUnaligned2 = 0
        self.mMapChar2Pos = {}
        self.mPercentGC = 0
        self.mMatrix = []

    def buildSubstitutionMatrix(self, seq1, seq2, alphabet):
        """given a pair of sequences, calculate
        a substitution matrix for the given alphabet.
        """
        if len(seq1) != len(seq2):
            raise ValueError("two sequences of unequal length submitted")

        nchars = len(alphabet)
        matrix = numpy.zeros((nchars, nchars), numpy.int)
        map_char2pos = {}
        for x in alphabet:
            map_char2pos[x] = len(map_char2pos)

        for x in range(len(seq1)):
            try:
                matrix[map_char2pos[seq1[x]], map_char2pos[seq2[x]]] += 1
            except KeyError:
                continue

        return matrix

    def loadPair(self, seq1, seq2):
        """load sequence properties from a pair.
        """
        SequencePairProperties.loadPair(self, seq1, seq2)

        alphabet = self.mAlphabet + self.mGapChar

        map_char2pos = {}
        for x in alphabet:
            map_char2pos[x] = len(map_char2pos)

        # build coordinates for various substitution subsets
        transitions, transversions = [], []
        for x in ("AG", "GA", "CT", "TC"):
            transitions.append((map_char2pos[x[0]], map_char2pos[x[1]]))

        for x in ("AT", "TA", "GT", "TG", "GC", "CG", "AC", "CA"):
            transversions.append((map_char2pos[x[0]], map_char2pos[x[1]]))

        matrix = self.buildSubstitutionMatrix(seq1, seq2, alphabet)
        matrix_acgt = matrix[0:4, 0:4]

        self.mMatrix = matrix
        self.mMapChar2Pos = map_char2pos
        self.mNAligned = numpy.sum(numpy.ravel(matrix_acgt))
        self.mNIdentical = numpy.sum(numpy.trace(matrix_acgt))
        self.mNTransitions = numpy.sum([matrix[x] for x in transitions])
        self.mNTransversions = numpy.sum([matrix[x] for x in transversions])
        self.mNDifferent = self.mNAligned - self.mNIdentical
        self.mNUnaligned1 = numpy.sum(numpy.ravel(matrix[0:4, 4]))
        self.mNUnaligned2 = numpy.sum(numpy.ravel(matrix[4, 0:4]))

        cp = self.mMapChar2Pos['C']
        gp = self.mMapChar2Pos['G']

        # sum all rows and columns that have a least one G or C
        # and remove those that have two in order to not double count
        gc = numpy.sum(
            self.mMatrix[0:4, cp] +
            self.mMatrix[0:4, gp] +
            self.mMatrix[cp, 0:4] +
            self.mMatrix[gp, 0:4]) \
            - self.mMatrix[cp, cp] \
            - self.mMatrix[gp, gp] \
            - self.mMatrix[cp, gp] \
            - self.mMatrix[gp, cp]

        try:
            self.mPercentGC = "%5.2f" % (100.0 * float(gc) / self.mNAligned)
        except ZeroDivisionError:
            self.mPercentGC = "na"

    def __str__(self):

        return "\t".join(map(str, (self.mNIdentical, self.mNAligned,
                                   self.mNDifferent,
                                   self.mNTransitions, self.mNTransversions,
                                   self.mNUnaligned1, self.mNUnaligned2,
                                   self.mPercentGC)))

    def getHeaders(self):
        return ["identical", "aligned", "different", "transitions",
                "transversions", "unaligned1", "unaligned2", "pgc"]


class SequencePairPropertiesCountsCodons(SequencePairPropertiesCountsNa):

    """the first characters are ACGT."""

    def __init__(self):
        SequencePairPropertiesCountsNa.__init__(self)
        self.mNNonSynonymous = 0
        self.mNSynonymous = 0

    def __str__(self):

        return "\t".join((SequencePairPropertiesCountsNa.__str__(self),
                          str(self.mNNonSynonymous),
                          str(self.mNSynonymous)))

    def getHeaders(self):
        return ["nonsyn", "nsyn"]


class SequencePairPropertiesPID(SequencePairPropertiesDistance):
    """Percent identity.

    The percent identity is the ratio of the number of identical
    residues divided by the number of aligned residues.

    """
    mGapChars = ("-", "."),

    def __init__(self, *args, **kwargs):
        SequencePairPropertiesDistance.__init__(self, *args, **kwargs)
        self.mPID = 0

    def loadPair(self, seq1, seq2):

        nidentical, naligned = 0, 0
        for c1, c2 in zip(seq1.upper(), seq2.upper()):
            if c1 in self.mGapChars or c2 in self.mGapChars:
                continue
            naligned += 1
            if c1 == c2:
                nidentical += 1

        self.mPID = 100.0 * float(nidentical) / naligned

    def getHeaders(self):
        return ["distance", ]

    def __str__(self):
        return self.mFormat % self.mPID
