"""SequenceProperties.py - Computing metrics of nucleotide sequences
=================================================================

This module provides methods for extracting and reporting sequence
properties of nucleotide sequences such as the composition, length,
etc.

The classes provide the algorithms to provide the property. They will
store the latest result for output. Thus, processing is a two-step
procedure::

   from SequenceProperties import SequencePropertiesLength
   from SequenceProperties import SequencePropertiesNA

   counters = [SequencePropertiesLength(), SequencePropertiesNA()]

   # output column headers
   headers = sum(c.getHeaders() for c in counters]
   print "\t".join(headers)

   for sequence in sequences:
      # load sequence in each counter
      for c in counters:
          c.loadSequence(sequence)
      # output results
      print "\t".join(map(str, counters))

This design is useful to compute multiple properties while iterating
only once over an input file and output a single, multi-column table.

.. note::
    While useful and in working order, the design of the classes is
    cumbersome.

Reference
---------

"""

import re
import math
import hashlib
import base64
import itertools
import six

from CGAT import Genomics as Genomics

import Bio.Alphabet.IUPAC


class SequenceProperties(object):
    """Base class.

    This class is the base class for SequenceProperty objects. Derived
    classes need to overload most of its methods.
    """

    def __init__(self):

        self.mLength = 0
        self.mSeqType = "na"
        self.mSequence = ""

    def addProperties(self, other):
        """add properties."""
        self.mLength += other.mLength

    def updateProperties(self):
        pass

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""

        self.mSeqType = seqtype
        self.mSequence = re.sub("[ -.]", "", sequence).upper()
        self.mLength = len(self.mSequence)

    def __str__(self):
        return "\t".join(self.getFields())

    def getFields(self):
        self.updateProperties()
        return []

    def getHeaders(self):
        return []


class SequencePropertiesSequence(SequenceProperties):
    """Add properties: the actual sequence.

    sequence
       The sequence

    This class outputs the actual sequence supplied.
    """

    def __init__(self):
        SequenceProperties.__init__(self)

    def addProperties(self, other):
        """add properties."""
        SequenceProperties.addProperties(self, other)

    def updateProperties(self):
        SequenceProperties.updateProperties(self)

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)
        # self.mSequence = sequence

    def getFields(self):
        fields = SequenceProperties.getFields(self)
        return fields + [self.mSequence, ]

    def getHeaders(self):
        headers = SequenceProperties.getHeaders(self)
        return headers + ["sequence"]


class SequencePropertiesHid(SequenceProperties):
    """Add properties: a hash of sequence

    hid
        Hash identifier of a sequence

    The hash is computed using the md5 algorithm and the resulting
    byte sequence is then translated into printable characters.
    """

    def __init__(self):
        SequenceProperties.__init__(self)
        self.mHid = ""

    def loadSequence(self, sequence, seqtype="na"):
        """load hid sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)

        # do the encryption
        h = hashlib.md5(six.b(sequence)).digest()
        # map to printable letters: hid has length 22, so the padded '=' are
        # truncated. You have to add them, if you ever want to decode,
        # but who would do such a thing :=)
        r = base64.encodestring(h)[0:22].decode("ascii")

        # finally substitute some characters:
        # '/' for '_', so we have legal file names
        # '[' for '+' and ']' for '=' for internet-applications
        hid = r.replace('/', '_').replace('+', '[').replace('=', ']')

        self.mHid = hid

    def getFields(self):
        fields = SequenceProperties.getFields(self)
        return fields + [self.mHid, ]

    def getHeaders(self):
        headers = SequenceProperties.getHeaders(self)
        return headers + ["hid", ]


class SequencePropertiesLength(SequenceProperties):
    """Add properties: sequence length and number of codons

    length
        Sequence length
    ncodons
        Length in codons

    The number of codons is 0 for an amino-acid sequence.
    """

    def __init__(self):
        SequenceProperties.__init__(self)
        self.mLength = 0
        self.mNCodons = 0

    def addProperties(self, other):
        """add properties."""
        self.mLength += other.mLength
        self.mNCodons += other.mNCodons

    def updateProperties(self):
        SequenceProperties.updateProperties(self)

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)
        self.mLength = len(sequence)
        if self.mSeqType == "na":
            self.mNCodons = len(sequence) / 3
        else:
            self.mNCodons = 0

    def getFields(self):
        fields = SequenceProperties.getFields(self)
        return fields + ["%i" % self.mLength,
                         "%i" % self.mNCodons]

    def getHeaders(self):
        headers = SequenceProperties.getHeaders(self)
        return headers + ["length", "ncodons"]


class SequencePropertiesNA(SequenceProperties):
    """Add properties: nucleotide composition

    nUnk
        Number of unknown residues
    nA, nC, nG, nT, nGC, nAT
        Nucleotide counts
    pA, pC, pG, pT, pGC, pAT
        Nucleotide frequencies

    """

    def __init__(self, reference_usage=[]):
        SequenceProperties.__init__(self)
        self.mCountsGC = 0
        self.mCountsAT = 0
        self.mCountsOthers = 0
        # counts of nucleotides
        self.mCountsNA = {}
        self.mAlphabet = Bio.Alphabet.IUPAC.unambiguous_dna.letters + "N"
        for x in self.mAlphabet:
            self.mCountsNA[x] = 0

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)
        self.mCountsGC += other.mCountsGC
        self.mCountsAT += other.mCountsAT
        self.mCountsOthers += other.mCountsOthers
        for na, count in list(other.mCountsNA.items()):
            self.mCountsNA[na] += count

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)
        # counts of nucleotides
        self.mCountsNA = {}
        for x in self.mAlphabet:
            self.mCountsNA[x] = 0

        for na in sequence.upper():
            if na in ('G', 'C'):
                self.mCountsGC += 1
            elif na in ('A', 'T'):
                self.mCountsAT += 1
            try:
                self.mCountsNA[na] += 1
            except KeyError:
                self.mCountsOthers += 1

    def getFields(self):
        fields = SequenceProperties.getFields(self)
        # Counts
        fields.append("%i" % self.mCountsOthers)

        t = 0
        for x in self.mAlphabet:
            fields.append("%i" % self.mCountsNA[x])
            t += self.mCountsNA[x]
        fields.append("%i" % self.mCountsGC)
        fields.append("%i" % self.mCountsAT)
        # Percentages
        if t == 0:
            for x in self.mAlphabet:
                fields.append("na")
        else:
            for x in self.mAlphabet:
                fields.append("%f" % (float(self.mCountsNA[x]) / t))

        t2 = float(self.mCountsGC + self.mCountsAT)
        if t2 == 0:
            fields.append("na")
            fields.append("na")
        else:
            fields.append("%f" % (float(self.mCountsGC) / t2))
            fields.append("%f" % (float(self.mCountsAT) / t2))

        return fields

    def getHeaders(self):
        headers = SequenceProperties.getHeaders(self)
        # Counts
        headers.append("nUnk")
        for x in self.mAlphabet:
            headers.append("n%s" % x)
        headers.append("nGC")
        headers.append("nAT")
        # Percentages
        for x in self.mAlphabet:
            headers.append("p%s" % x)
        headers.append("pGC")
        headers.append("pAT")

        return headers


class SequencePropertiesDN(SequenceProperties):
    """Add Properties : dinucleotide counts

    nAA, nAC, ...
        Dinucleotide counts
    mCountsOthers
        Unknown dinucleotides
    """

    def __init__(self, reference_usage=[]):

        SequenceProperties.__init__(self)
        self.mCountsDinuc = {}
        self.mCountsOthers = 0
        self.mAlphabet = Bio.Alphabet.IUPAC.unambiguous_dna.letters
        for dinucleotide in itertools.product(self.mAlphabet, repeat=2):
            self.mCountsDinuc["".join(dinucleotide)] = 0

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)
        self.mCountsOthers += other.mCountsOthers
        for dn, count in list(other.mCountsDinuc.items()):
            self.mCountsDinuc[dn] += count

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)

        # IMS: generator rather than list
        # to save memory for might be a neater way of doing this.
        for dinuc in (sequence[x - 2:x]
                      for x in range(2, len(sequence) + 1, 1)):
            try:
                self.mCountsDinuc[dinuc] += 1
            except KeyError:
                self.mCountsOthers += 1

    def getFields(self):

        fields = SequenceProperties.getFields(self)
        for x in itertools.product(self.mAlphabet, repeat=2):
            fields.append("%i" % self.mCountsDinuc["".join(x)])
        fields.append("%s" % self.mCountsOthers)
        return fields

    def getHeaders(self):

        headers = SequenceProperties.getHeaders(self)
        for dinucleotide in itertools.product(self.mAlphabet, repeat=2):
            headers.append("".join(dinucleotide))
        headers.append("mCountsOthers")
        return headers


class SequencePropertiesCpg(SequencePropertiesNA, SequencePropertiesDN):
    """Add Properties : CpG density and observed / expected.

    CpG_count
        Number of CpG in sequence
    CpG_density
        CpG density, number of CpG divided by 2 * sequence length
    CpG_ObsExp
        Ratio of observed to expected number of CpG. The latter is
        calculated as the product of nC * nG. The ratio is normalized
        by the sequence length.  Set to 0 if no ``C`` or ``G`` in
        sequence.

    """

    def __init__(self, reference_usage=[]):

        SequencePropertiesNA.__init__(self)
        SequencePropertiesDN.__init__(self)
        self.mCpG_density = 0
        self.mCpG_ObsExp = 0

    def addProperties(self, other):
        SequencePropertiesDN.addProperties(self, other)
        SequencePropertiesNA.addProperties(self, other)
        self.mCpG_density += other.mCpG_density
        self.mCpG_ObsExp += other.mCpG_ObsExp

    def loadSequence(self, sequence, seqtype="na"):
        SequencePropertiesNA.loadSequence(self, sequence, seqtype)
        SequencePropertiesDN.loadSequence(self, sequence, seqtype)

        if self.mLength > 0:
            self.mCpG_density = (float(self.mCountsDinuc["CG"]) /
                                 (float(self.mLength) / 2.0))
            if (float(self.mCountsNA["C"] * self.mCountsNA["G"])) / self.mLength > 0:
                self.mCpG_ObsExp = ((float(self.mCountsDinuc["CG"])) /
                                    ((float(self.mCountsNA["C"]) *
                                      float(self.mCountsNA["G"])) /
                                     float(self.mLength)))
            else:
                self.mCpG_ObsExp = 0.0
        else:
            self.mCpG_density = 0.0
            self.mCpG_ObsExp = 0.0

    def getFields(self):
        fields = ["%s" % self.mCountsDinuc["CG"]]
        fields.append("%6.4f" % self.mCpG_density)
        fields.append("%6.4f" % self.mCpG_ObsExp)
        return fields

    def getHeaders(self):
        headers = ["CpG_count"]
        headers.append("CpG_density")
        headers.append("CpG_ObsExp")
        return headers


class SequencePropertiesGaps(SequenceProperties):
    """Add Properties : number of gaps in a sequence

    Gaps are identified by unknown characters (``[XN]``)

    ngaps
         Number of gap characters in sequnce
    nseq_regions
         Number of ungapped regions
    ngap_regions
         Number of gapped regions
    """

    gap_chars = 'xXnN'

    def __init__(self,  gap_chars='xXnN', *args, **kwargs):
        SequenceProperties.__init__(self, *args, **kwargs)
        self.gap_chars = gap_chars

        self.ngaps = 0
        self.nseq_regions = 0
        self.ngap_regions = 0

    def getFields(self):

        fields = SequenceProperties.getFields(self)
        return fields + ["%i" % self.ngaps,
                         "%i" % self.nseq_regions,
                         "%i" % self.ngap_regions]

    def getHeaders(self):

        headers = SequenceProperties.getHeaders(self)
        return headers + ["ngaps", "nseq_regions", "ngap_regions"]

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequenceProperties.loadSequence(self, sequence, seqtype)

        gap_chars = self.gap_chars
        ngaps = 0
        was_gap = not (sequence[0] in gap_chars)

        ngap_regions = 0
        nseq_regions = 0

        x = 0
        for x, c in enumerate(sequence):

            is_gap = c in gap_chars
            if is_gap:
                ngaps += 1
                if not was_gap:
                    ngap_regions += 1
            else:
                if was_gap:
                    nseq_regions += 1
            was_gap = is_gap

        self.ngaps = ngaps
        self.ngap_regions = ngap_regions
        self.nseq_regions = nseq_regions

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)

        self.ngaps += other.ngaps
        self.nseq_regions += other.nseq_regions
        self.ngap_regions += other.ngap_regions


class SequencePropertiesDegeneracy (SequencePropertiesLength):
    """Add properties : codon degeneracy

    nstops
        Number of stop codons
    nsites1d
        Number of non-degenerate sites
    nsites2d, nsites3d, nsites4d
        Number 2-fold, 3-fold, 4-fold degenerate sites
    ngc
        Number of positions containing either G or C
    ngc3
        Number of 3rd codon position containing G or C
    ngc3
        Number of non-degenerate 3rd codon position containing G or C
    n2gc3, n3gc3, n4gc3
        Number of 2-fold, 3-fold, 4-fold degenerate 3rd codon positions
        containing G or C
    pgc
        Percentage of positions containing either G or C
    pgc3
        Percentage of 3rd codon position containing G or C
    pgc3
        Percentage of non-degenerate 3rd codon position containing G or C
    p2gc3, p3gc3, p4gc3
        Percentage of 2-fold, 3-fold, 4-fold degenerate 3rd codon positions
        containing G or C

    The degeneracies for amino acids are::

       2: MW are non-degenerate.
       9: EDKNQHCYF are 2-fold degenerate.
       1: I is 3-fold degenerate
       5: VGATP are 4-fold degenerate.
       3: RLS are 2-fold and four-fold degenerate.
          Depending on the first two codons, the codons are counted
          as two or four-fold degenerate codons. This is encoded
          in the file Genomics.py.

    The number of degenerate sites is computed across all
    codon positions.

    """

    def __init__(self):

        SequencePropertiesLength.__init__(self)

        self.mNGC = 0
        self.mNSites1D, self.mNSites2D, self.mNSites3D, self.mNSites4D = (
            0, 0, 0, 0)
        self.mNGC3 = 0
        self.mN2DGC3 = 0
        self.mN3DGC3 = 0
        self.mN4DGC3 = 0

        self.mNStopCodons = 0

        # setup counting arrays
        # nucleotide counts for each position (is not a sum of the counts
        # per degenerate site, as the codon might be intelligible, e.g. GNN).
        self.mCounts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0}]

        # nucleotide counts for each position per degeneracy
        self.mCountsDegeneracy = []

        for x in (0, 1, 2):
            xx = []
            for y in range(5):
                yy = {}
                for z in Bio.Alphabet.IUPAC.extended_dna.letters:
                    yy[z] = 0
                xx.append(yy)
            self.mCountsDegeneracy.append(xx)

    def addProperties(self, other):

        SequencePropertiesLength.addProperties(self, other)
        self.mNStopCodons += other.mNStopCodons

        for x in (0, 1, 2):
            for y in range(5):
                for z in Bio.Alphabet.IUPAC.extended_dna.letters:
                    self.mCountsDegeneracy[x][y][
                        z] += other.mCountsDegeneracy[x][y][z]

        for x in (0, 1, 2):
            for y, count in list(other.mCounts[x].items()):
                self.mCounts[x][y] += count

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""

        SequencePropertiesLength.loadSequence(self, sequence, seqtype)
        if len(sequence) % 3:
            raise ValueError(
                '''sequence length is not a multiple of 3 (length=%i)''' %
                (len(sequence)))

        # uppercase all letters
        sequence = sequence.upper()

        self.mNStopCodons = 0

        # setup counting arrays
        # nucleotide counts for each position (is not a sum of the counts
        # per degenerate site, as the codon might be intelligible, e.g. GNN).
        self.mCounts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0}]

        # nucleotide counts for each position per degeneracy
        self.mCountsDegeneracy = []

        for x in (0, 1, 2):
            xx = []
            for y in range(5):
                yy = {}
                for z in Bio.Alphabet.IUPAC.extended_dna.letters:
                    yy[z] = 0
                xx.append(yy)
            self.mCountsDegeneracy.append(xx)

        # use generator rather than list to save memory
        for codon in (sequence[x:x + 3] for x in range(0, len(sequence), 3)):

            for x in (0, 1, 2):
                self.mCounts[x][codon[x]] += 1

            if Genomics.IsStopCodon(codon):
                self.mNStopCodons += 1
                continue

            try:
                aa, deg1, deg2, deg3 = Genomics.GetDegeneracy(codon)
                degrees = (deg1, deg2, deg3)
                for x in range(len(degrees)):
                    self.mCountsDegeneracy[x][degrees[x]][codon[x]] += 1

            except KeyError:
                pass

    def updateProperties(self):
        """update fields from counts."""

        SequencePropertiesLength.updateProperties(self)
        self.mNGC = 0
        self.mNSites1D, self.mNSites2D, self.mNSites3D, self.mNSites4D = (
            0,) * 4
        self.mNSitesD3 = 0

        for x in (0, 1, 2):
            self.mNGC += self.mCounts[x]['C'] + self.mCounts[x]['G']
            self.mNSites1D += sum(self.mCountsDegeneracy[x][1].values())
            self.mNSites2D += sum(self.mCountsDegeneracy[x][2].values())
            self.mNSites3D += sum(self.mCountsDegeneracy[x][3].values())
            self.mNSites4D += sum(self.mCountsDegeneracy[x][4].values())

        self.mNGC3 = self.mCounts[2]['C'] + self.mCounts[2]['G']
        self.mN2DGC3 = self.mCountsDegeneracy[2][2][
            'C'] + self.mCountsDegeneracy[2][2]['G']
        self.mN3DGC3 = self.mCountsDegeneracy[2][3][
            'C'] + self.mCountsDegeneracy[2][3]['G']
        self.mN4DGC3 = self.mCountsDegeneracy[2][4][
            'C'] + self.mCountsDegeneracy[2][4]['G']

        # count all degenerate sites at third codon position
        self.mNSitesD3 = sum(self.mCountsDegeneracy[2][2].values()) +\
            sum(self.mCountsDegeneracy[2][3].values()) +\
            sum(self.mCountsDegeneracy[2][4].values())

        # number of GCs at degenerate third codon positions
        self.mNDGC3 = self.mN2DGC3 + self.mN3DGC3 + self.mN4DGC3

        if self.mLength > 0:
            self.mPNGC = "%5.4f" % (float(self.mNGC) / float(self.mLength))
        else:
            self.mPNGC = "na"

        if self.mNCodons > 0:
            self.mPNGC3 = "%5.4f" % (float(self.mNGC3) / float(self.mNCodons))
        else:
            self.mPNGC3 = "na"

        if self.mNSites2D > 0:
            self.mPN2DGC3 = "%5.4f" % (
                float(self.mN2DGC3) / float(self.mNSites2D))
        else:
            self.mPN2DGC3 = "na"

        if self.mNSites3D > 0:
            self.mPN3DGC3 = "%5.4f" % (
                float(self.mN3DGC3) / float(self.mNSites3D))
        else:
            self.mPN3DGC3 = "na"

        if self.mNSites4D > 0:
            self.mPN4DGC3 = "%5.4f" % (
                float(self.mN4DGC3) / float(self.mNSites4D))
        else:
            self.mPN4DGC3 = "na"

        if self.mNSitesD3 > 0:
            self.mPNDGC3 = "%5.4f" % (
                float(self.mNDGC3) / float(self.mNSitesD3))
        else:
            self.mPNDGC3 = "na"

    def getFields(self):
        fields = SequencePropertiesLength.getFields(self)
        return fields + ["%i" % self.mNStopCodons,
                         "%i" % self.mNSites1D,
                         "%i" % self.mNSites2D,
                         "%i" % self.mNSites3D,
                         "%i" % self.mNSites4D,
                         "%i" % self.mNSitesD3,
                         "%i" % self.mNGC,
                         "%i" % self.mNGC3,
                         "%i" % self.mNDGC3,
                         "%i" % self.mN2DGC3,
                         "%i" % self.mN3DGC3,
                         "%i" % self.mN4DGC3,
                         "%s" % self.mPNGC,
                         "%s" % self.mPNGC3,
                         "%s" % self.mPNDGC3,
                         "%s" % self.mPN2DGC3,
                         "%s" % self.mPN3DGC3,
                         "%s" % self.mPN4DGC3]

    def getHeaders(self):
        headers = SequencePropertiesLength.getHeaders(self)
        return headers + ["nstops",
                          "nsites1d",
                          "nsites2d",
                          "nsites3d",
                          "nsites4d",
                          "nsitesd3",
                          "ngc",
                          "ngc3",
                          "ndgc3",
                          "n2dgc3",
                          "n3dgc3",
                          "n4dgc3",
                          "pgc",
                          "pgc3",
                          "pdgc3",
                          "p2dgc3",
                          "p3dgc3",
                          "p4dgc3",
                          ]


class SequencePropertiesAA(SequenceProperties):
    """Add Properties : amino acid composition of nucleotide sequence.

    The codons in the nucleotide sequence are translated into amino
    acids before counting. The nucleotide sequence must be a multiple
    of 3.

    nA, nC, nD, ...
        Amino acid counts.
    pA, pC, pD, ...
        Amino acid frequencies.

    """

    def __init__(self, reference_usage=[]):

        SequenceProperties.__init__(self)

        self.mReferenceUsage = reference_usage

        # counts of amino acids
        self.mCountsAA = {}
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            self.mCountsAA[x] = 0

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)

        for aa, count in list(other.mCountsAA.items()):
            self.mCountsAA[aa] += count

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""

        SequenceProperties.loadSequence(self, sequence, seqtype)

        if len(sequence) % 3:
            raise ValueError(
                '''sequence length is not a multiple of 3 (length=%i)''' %
                (len(sequence)))

        # counts of amino acids
        self.mCountsAA = {}

        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            self.mCountsAA[x] = 0

        for codon in (sequence[x:x + 3] for x in range(0, len(sequence), 3)):
            aa = Genomics.MapCodon2AA(codon)
            self.mCountsAA[aa] += 1

    def getFields(self):

        fields = SequenceProperties.getFields(self)
        t = 0
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            fields.append("%i" % self.mCountsAA[x])
            t += self.mCountsAA[x]
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            fields.append("%f" % (float(self.mCountsAA[x]) / t))
        return fields

    def getHeaders(self):
        '''Return list of data headers'''
        headers = SequenceProperties.getHeaders(self)
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            headers.append("n%s" % x)
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            headers.append("p%s" % x)
        return headers


class SequencePropertiesAminoAcids(SequenceProperties):
    """Add Properties : amino acid composition

    nA, nC, nD, ...
        Amino acid counts.
    pA, pC, pD, ...
        Amino acid frequencies.
    """

    def __init__(self, reference_usage=[]):

        SequenceProperties.__init__(self)

        # counts of amino acids
        self.mCountsAA = {}
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            self.mCountsAA[x] = 0
        self.mOtherCounts = 0

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)

        for aa, count in list(other.mCountsAA.items()):
            self.mCountsAA[aa] += count

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence. """

        SequenceProperties.loadSequence(self, sequence, seqtype)

        # set to zero
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            self.mCountsAA[x] = 0
        self.mOtherCounts = 0

        for aa in sequence:
            if aa == "-":
                continue
            try:
                self.mCountsAA[aa] += 1
            except KeyError:
                self.mOtherCounts += 1

    def getFields(self):

        fields = SequenceProperties.getFields(self)

        t = 0

        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            fields.append("%i" % self.mCountsAA[x])
            t += self.mCountsAA[x]

        if t > 0:
            for x in Bio.Alphabet.IUPAC.extended_protein.letters:
                fields.append("%f" % (float(self.mCountsAA[x]) / t))
        else:
            for x in Bio.Alphabet.IUPAC.extended_protein.letters:
                fields.append("0")

        return fields

    def getHeaders(self):

        fields = SequenceProperties.getHeaders(self)
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            fields.append("n%s" % x)
        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            fields.append("p%s" % x)

        return fields


class SequencePropertiesCodons(SequencePropertiesLength):
    """Add Properties : codon frequencies

    nAAA, nAAC, ...
         Codon counts
    pAAA, pAAC, ...
         Codon frequencies
    """

    def __init__(self):

        SequencePropertiesLength.__init__(self)

        self.mCodonCounts = {}
        for c in list(Genomics.GeneticCodeAA.keys()):
            self.mCodonCounts[c] = 0

    def addProperties(self, other):
        SequencePropertiesLength.addProperties(self, other)

        for codon, count in list(other.mCodonCounts.items()):
            if codon not in self.mCodonCounts:
                self.mCodonCounts[codon] = 0
            self.mCodonCounts[codon] += count

    def updateProperties(self):

        SequencePropertiesLength.updateProperties(self)

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""

        if len(sequence) % 3:
            raise ValueError(
                '''sequence length is not a multiple of 3
                (length=%i)''' % (len(sequence)))

        SequencePropertiesLength.loadSequence(self, sequence, seqtype)

        # uppercase all letters and count codons
        self.mCodonCounts = Genomics.CountCodons(sequence.upper())

    def getFields(self):

        fields = SequencePropertiesLength.getFields(self)

        keys = list(self.mCodonCounts.keys())
        keys.sort()
        t = 0
        for key in keys:
            t += self.mCodonCounts[key]
            fields.append("%i" % self.mCodonCounts[key])
        t = float(t)
        for key in keys:
            fields.append("%f" % (float(self.mCodonCounts[key]) / t))

        return fields

    def getHeaders(self):

        fields = SequencePropertiesLength.getHeaders(self)

        keys = list(self.mCodonCounts.keys())
        keys.sort()
        for key in keys:
            fields.append("n%s" % key)
        for key in keys:
            fields.append("p%s" % key)

        return fields


class SequencePropertiesCodonUsage(SequencePropertiesCodons):
    """Add properties : Codon usage

    The codon frequency is the ratio of the number of times a
    particular codon is used for a particular amino acid, didived the
    number of times that particular amino acid appears in the
    sequence. A ratio of 1.0 means that this particular codon is
    always used to encode its amino acid, while a frequency of 0.5
    means it is used 50% of the times.

    rAAA, rAAC, ...
        Codon frequencies.

    """

    def __init__(self):

        SequencePropertiesCodons.__init__(self)

        self.mCodonFrequencies = {}

        for codon, aa in list(Genomics.GeneticCodeAA.items()):
            self.mCodonFrequencies[codon] = 0
        for codon in Genomics.StopCodons:
            self.mCodonFrequencies[codon] = 0

    def addProperties(self, other):

        SequencePropertiesCodons.addProperties(self, other)

    def updateProperties(self):

        SequencePropertiesCodons.updateProperties(self)

        self.mCodonFrequencies = Genomics.CalculateCodonFrequenciesFromCounts(
            self.mCodonCounts)

    def getFields(self):

        fields = SequenceProperties.getFields(self)

        keys = list(self.mCodonFrequencies.keys())
        keys.sort()
        for key in keys:
            fields.append("%f" % self.mCodonFrequencies[key])

        return fields

    def getHeaders(self):

        fields = SequenceProperties.getHeaders(self)

        keys = list(self.mCodonFrequencies.keys())
        keys.sort()
        for key in keys:
            fields.append("r%s" % key)

        return fields


class SequencePropertiesCodonTranslator(SequencePropertiesCodonUsage):
    """Add properties : codon sequence is translated into frequencies.

    tsequence
        comma separated list of codon frequencies. The frequencies are
        in percentages.
    """

    def __init__(self):

        SequencePropertiesCodonUsage.__init__(self)

        self.mSequence = ""

    def addProperties(self, other):

        self.mSequence += other.mSequence

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""
        SequencePropertiesCodons.loadSequence(self, sequence, seqtype)

        self.mSequence = sequence

    def getHeaders(self):

        fields = SequenceProperties.getHeaders(self)
        fields.append("tsequence")
        return fields

    def getFields(self):

        fields = SequenceProperties.getFields(self)

        data = []
        for x in range(0, len(self.mSequence), 3):
            codon = self.mSequence[x:x + 3]
            if codon in self.mCodonFrequencies:
                v = self.mCodonFrequencies[codon]
            else:
                v = 0.0

            data.append("%i" % (v * 100))

        fields.append(",".join(data))

        return fields


class SequencePropertiesBias(SequencePropertiesCodons):
    """Add properties : bias measures of codon sequence.

    This class outputs metrics showing how biased the codon usage in a
    particular sequence is compared to a reference codon usage.  The
    reference codon usage is given as a dictionary of codon
    frequencies and multiple dictionaries can be given to compute the
    bias against multiple codon usages.

    entropy
        Entropy of the sequence.
    ml0, ml1, ...
        Message length of sequence compared to reference codon
        usages.
    relml0, relml1, ...
        Relative message length of sequence compared to reference codon
        usages. The relative message length is the message lenght divided
        by the number of codons.
    relentropy0, relentropy1, ...
        Relative entropy of sequence compared to reference codon usages.
        Also called conditional entropy or encoding cost.
    kl0, kl1, ...
        Kullback-Leibler Divergence (relative entropy) of sequence compared
        to reference codon usages.

    Arguments
    ---------

    reference_usage : list
        A list of codon frequency tables. The bias will be computed
        against each.
    pseudocounts : int
        Pseudo-counts to add
    """

    def __init__(self, reference_usage=[], pseudocounts=0):

        SequencePropertiesCodons.__init__(self)
        self.mReferenceUsage = reference_usage
        self.mEntropy = 0
        self.mPseudoCounts = pseudocounts

    def updateProperties(self):

        SequenceProperties.updateProperties(self)
        self.mProperties = []
        self.mEntropy = self.getEntropy()
        for usage in self.mReferenceUsage:
            self.mProperties.append({'ml': self.getMessageLength(usage),
                                     'entropy': self.getEntropy(usage),
                                     'kl': self.getKL(usage)})

    def getMessageLength(self, usage):
        """return message length of a sequence
        in terms of its reference usage."""
        ml = 0
        for codon, count in list(self.mCodonCounts.items()):
            ml -= float(count) * math.log(usage[codon])
        return ml

    def getEntropy(self, usage=None):
        """return entropy of a source in terms of a reference usage.
        Also called conditional entropy or encoding cost.

        Note that here I compute the sum over 20 entropies,
        one for each amino acid.

        If not given, calculate entropy.
        """

        e = 0
        freqs = Genomics.CalculateCodonFrequenciesFromCounts(
            self.mCodonCounts, self.mPseudoCounts)
        if usage is None:
            usage = freqs
        for codon, count in list(self.mCodonCounts.items()):
            e -= freqs[codon] * math.log(usage[codon])
        return e

    def getKL(self, usage):
        """return Kullback-Leibler Divergence (relative entropy) of sequences with
        respect to reference codon usage.
        """
        e = 0
        freqs = Genomics.CalculateCodonFrequenciesFromCounts(
            self.mCodonCounts, self.mPseudoCounts)
        for codon, count in list(self.mCodonCounts.items()):
            e += usage[codon] * math.log(usage[codon] / freqs[codon])
        return e

    def getFields(self):

        fields = SequenceProperties.getFields(self)
        fields.append("%f" % self.mEntropy)
        for x in self.mProperties:
            fields.append("%f" % x['ml'])
            fields.append("%f" % (x['ml'] / self.mNCodons))
            fields.append("%f" % x['entropy'])
            fields.append("%f" % x['kl'])
        return fields

    def getHeaders(self):

        headers = SequenceProperties.getHeaders(self)
        headers.append("entropy")
        for x in range(len(self.mReferenceUsage)):
            headers.append("ml%i" % x)
            headers.append("relml%i" % x)
            headers.append("relentropy%i" % x)
            headers.append("kl%i" % x)
        return headers


class SequencePropertiesCounts(SequenceProperties):
    """Add Properties : Residue counts against arbirtrary alphabet

    nUnk
        Number of unknown residues
    nA, nB, ...
        Character counts
    pA, pB, ...
        Character frequencies

    Arguments
    ---------
    alphabet : string
        List of characters in alphabet

    """

    def __init__(self, alphabet):

        SequenceProperties.__init__(self)
        self.mAlphabet = alphabet
        self.mCounts = {}
        self.mCountsOthers = 0

        for x in self.mAlphabet:
            self.mCounts[x] = 0

    def addProperties(self, other):
        SequenceProperties.addProperties(self, other)
        for na, count in list(other.mCounts.items()):
            self.mCounts[na] += count
        self.mCountsOthers += self.mCountsOthers

    def loadSequence(self, sequence, seqtype="na"):
        """load sequence properties from a sequence."""

        SequenceProperties.loadSequence(self, sequence, seqtype)

        # counts of amino acids
        self.mCounts = {}
        for x in self.mAlphabet:
            self.mCounts[x] = 0
        self.mCountsOthers = 0

        for na in sequence.upper():
            try:
                self.mCounts[na] += 1
            except KeyError:
                self.mCountsOthers += 1

    def getFields(self):
        fields = SequenceProperties.getFields(self)
        fields.append("%i" % self.mCountsOthers)
        t = 0

        for x in self.mAlphabet:
            fields.append("%i" % self.mCounts[x])
            t += self.mCounts[x]

        if t == 0:
            for x in self.mAlphabet:
                fields.append("na")
        else:
            for x in self.mAlphabet:
                fields.append("%f" % (float(self.mCounts[x]) / t))

        return fields

    def getHeaders(self):
        headers = SequenceProperties.getHeaders(self)
        headers.append("nUnk")
        headers.extend(["n%s" % x for x in self.mAlphabet])
        headers.extend(["p%s" % x for x in self.mAlphabet])
        return headers


class SequencePropertiesEntropy(SequencePropertiesCounts):
    """Add properties : Entropy of a sequence

    entropy
        Entropy of the sequence

    Arguments
    ---------
    alphabet : string
        List of characters in alphabet
    pseudocounts : int
        Pseudo-counts to add

    """

    def __init__(self, alphabet, pseudocounts=0):

        SequencePropertiesCounts.__init__(self, alphabet=alphabet)
        self.mPseudoCounts = pseudocounts
        self.mEntropy = None

    def addProperties(self, other):
        SequencePropertiesCounts.addProperties(self, other)

    def updateProperties(self):

        SequencePropertiesCounts.updateProperties(self)

        self.mProperties = []

        self.mEntropy = self.getEntropy()

    def getEntropy(self, usage=None):
        """return entropy of a source in terms of a reference usage.

        Also called conditional entropy or encoding cost.
        """

        e = None

        v = list(self.mCounts.values())
        total = sum(v) + len(v) * self.mPseudoCounts
        pc = self.mPseudoCounts
        if total != 0:
            e = 0
            frequencies = [float(x + pc) / total for x in v]
            for f in frequencies:
                if f > 0:
                    e -= f * math.log(f)
        return e

    def getFields(self):

        fields = SequenceProperties.getFields(self)
        if self.mEntropy is not None:
            fields.append("%f" % self.mEntropy)
        else:
            fields.append("na")

        return fields

    def getHeaders(self):

        fields = SequenceProperties.getHeaders(self)
        fields.append("entropy")
        return fields
