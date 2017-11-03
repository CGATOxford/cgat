"""
Genomics.py - Tools for working with genomic data
=================================================

:Tags: Python

Reference
---------

"""

import numpy
import os
import string
import re
import hashlib
import base64
import tempfile
import sys
from functools import reduce

try:
    import alignlib_lite
except ImportError:
    pass

from CGAT import AString as AString

global_last_filename_genome = None
global_forward_sequences = {}
global_last_sbjct_token = None

if sys.version_info.major >= 3:
    global_translator = str.maketrans("ACGTacgt", "TGCAtgca")
else:
    global_translator = string.maketrans("ACGTacgt", "TGCAtgca")


def complement(s):
    """reverse complement a sequence.

    >>> complement("ACATACATACTA")
    'TAGTATGTATGT'

    Returns
    -------
    string
    """
    return s[::-1].translate(global_translator)


def GetHID(sequence):
    """returns a hash value for a sequence.

    The hash value is computed using md5 and converted
    into printable characters.

    >>> GetHID("ACATACATACTA")
    'trcPGx9VNT36XMlG0XvUBQ'

    Returns
    -------
        A hash value
    """

    # do the encryption
    h = hashlib.md5(sequence).digest()

    # map to printable letters: hid has length 22, so the padded '='
    # are truncated.
    r = base64.encodestring(h)[0:22]

    # finally substitute some characters:
    # '/' for '_', so we have legal file names
    # '[' for '+' and ']' for '=' for internet-applications

    hid = string.replace(r, '/', '_')
    hid = string.replace(hid, '+', '[')
    hid = string.replace(hid, '=', ']')

    return hid


def String2Location(s):
    """convert a string to location information.

    >>> String2Location("chr1:12:15")
    ('chr1', '+', 12, 15)

    Returns
    -------
    contig : string
    strand : string
    start : int
    end : int
    """

    data = re.split("[:]+", s)
    if len(data) == 3:
        return data[0], "+", int(data[1]), int(data[2])
    elif len(data) == 4:
        return data[0], data[1], int(data[2]), int(data[3])
    else:
        raise ValueError("unknown format %s" % (s))


def readContigSizes(infile):
    """read sizes of contigs from file.

    Arguments
    ---------
    infile : string
        Filename of :term:`tsv` separated file.

    Returns
    -------
    dict

    """
    sizes = {}
    for line in infile:
        if line[0] == "#":
            continue
        sbjct_token, size = line[:-1].split("\t")[:2]
        sizes[sbjct_token] = int(size)

    return sizes


def forceForwardCoordinates(start, end, strand, length):
    """return forward coordinates.

    If strand is negative, the coordinates in a and b
    will be converted. If they are on the positive
    strand, they will be returned as is.

    Arguments
    ---------
    start : int
        Start coordinate
    end : int
        End coordinate
    strand : string
        Strand of interval. The values of "-", "0", "-1" indicate
        a negative strand.
    length : int
        Length of chromosome.
    """
    if strand in ("-", "0", "-1", 0):
        return length - end, length - start
    else:
        return start, end


def CountGeneFeatures(first_position,
                      alignment,
                      genomic_sequence=None,
                      border_stop_codon=0,
                      stop_codons=("TAG", "TAA", "TGA")):
    """calculate number of genomic features in a peptide to genome
    alignment.

    Note that codons can be split, for example::

        S 0 2 5 0 2 I 0 17541 3 0 2 S 1 2 5 0 2 I 0 27979 3 0 2 S 1 2

    Arguments
    ---------
    first_position : int
         Start of alignment on genome.
    alignment : string
         Alignment in :term:`CIGAR` format, for example from exonerate_.
    genomic_sequence : string
         Genomic sequence for alignment
    border_stop_codon : int
         Number of codons that are ignored at the edges of match
         regions.  border_stop_codon should be divisible by three.
    stop_codons : list
         List of stop codons

    Returns
    -------
    nintrons : int
        Number of introns
    nframeshifts : int
        Number of frameshifts in aligment.
    ngaps : int
        Number of gaps in aligment.
    nsplit : int
        Number of codons split by introns in alignment.
    nstopcodons : int
        Number of stop codons in alignment.
    disruptions : list
       List of disruptions in the prediction. Each
       disruption is a tuple of ( "stop|frameshift", position in protein,
       position in cds, position on genomic sequence).

    """

    current_pos_genome = first_position
    current_pos_cds = 0

    nstopcodons = 0
    nintrons = 0
    nframeshifts = 0
    ngaps = 0
    nsplits = 0

    disruptions = []

    # position of nucleotides for split codons
    partial_codon = []

    first_state = True

    for state, l_protein, l_genome in alignment:

        if state in ("G", "M", "P"):

            # reset split codon. Do a sanity check, but ignore first
            # split codon in the case of exon alignments
            if len(partial_codon) % 3 != 0 and not first_state:
                raise ValueError(
                    "split codon was not multiple of three: "
                    "partial_codon=%s, alignment=%s" %
                    (partial_codon, alignment))
                # print alignment
                # pass

            partial_codon = []
            first_state = False

        if state in ("M", "G"):
            # check for stop codons, ignore the first
            # border_stop_codon nucleotides note: this should be a
            # multiple of three
            if genomic_sequence:

                for x in range(border_stop_codon,
                               l_genome - border_stop_codon,
                               3):
                    y = current_pos_cds + x
                    z = current_pos_genome + x
                    codon = genomic_sequence[z:z + 3].upper()

                    if codon in stop_codons:
                        nstopcodons += 1
                        disruptions.append(("stop",
                                            y, y + 3,
                                            z, z + 3))

        if state == "I":
            nintrons += 1

        elif state == "F":
            nframeshifts += 1
            disruptions.append((
                "frameshift",
                current_pos_cds, current_pos_cds + l_genome,
                current_pos_genome, current_pos_genome + l_genome))

        elif state == "G":
            ngaps += 1

        elif state == "P":
            ngaps += 1

        elif state == "S":

            # I used to ignore alignments that start with a split
            # codon as they might have been incomplet. However, the
            # code below seems to work for split codens if not
            # first_state:
            nsplits += 1

            # check for stop-codons in split codons as well:
            for x in range(l_genome):
                partial_codon.append(
                    (current_pos_cds + x, current_pos_genome + x))

            if border_stop_codon < 3 and genomic_sequence and \
               len(partial_codon) % 3 == 0:

                for x in range(0, len(partial_codon), 3):
                    codon = "".join(genomic_sequence[c[1]]
                                    for c in partial_codon[x:x + 3]).upper()
                    if codon in stop_codons:
                        nstopcodons += 1
                        disruptions.append(
                            ("split-stop",
                             partial_codon[x][
                                 0], partial_codon[x + 2][0],
                             partial_codon[x][1], partial_codon[x + 2][1]))
        # advance in cds
        if state in ("M", "G", "F", "S"):
            current_pos_cds += l_genome

        current_pos_genome += l_genome

    return nintrons, nframeshifts, ngaps, nsplits, nstopcodons, disruptions


def Alignment2String(alignment):
    """convert a tuple alignment to an alignment string.
    """
    return " ".join([" ".join(list(map(str, x))) for x in alignment])


def String2Alignment(source):
    """convert an alignment string to a tuple alignment.
    """

    d = source.split(" ")
    ali = []
    if len(d) < 3:
        return ali

    for x in range(0, len(d), 3):
        ali.append((d[x], int(d[x + 1]), int(d[x + 2])))

    return ali


def GetAlignmentLength(alignment):
    """return Alignment length"""

    q, s = 0, 0
    for state, l_query, l_sbjct in alignment:
        q += l_query
        s += l_sbjct

    return q, s


def Alignment2ExonBoundaries(alignment,
                             query_from=0,
                             sbjct_from=0,
                             add_stop_codon=1):
    """extract exon coordinates from a peptide2genome alignment.

    Arguments
    ---------
    aligment : list
        List of tuples of the alignment in CIGAR format.
    query_from : int
        Start position of alignment on peptide sequence.
    sbjct_from : int
        Start position of alignment on nucleotide sequence.
    add_stop_codon : int
        Add final stop codon to exon boundaries.

    Returns
    -------
    exons : list
        A list of exons. Each exon is a tuple of (query_from,
        query_pos, frame, sbjct_from, sbjct_pos, ali)

    """

    exons = []

    # count in nucleotides for query
    query_from *= 3
    query_pos = query_from
    sbjct_pos = sbjct_from
    frame = 0
    ali = []

    p_as_intron = True
    if p_as_intron:
        ali_states = ("S", "M", "G", "F")
    else:
        ali_states = ("S", "M", "G", "F", "P")

    for state, l_query, l_sbjct in alignment:

        # count as nucleotides
        l_query *= 3

        if state in ali_states:
            ali.append((state, l_query, l_sbjct))

        if state == "M":
            query_pos += l_query

        elif state == "G":
            query_pos += l_query

        elif state == "P":
            # treat as intron
            if p_as_intron:
                frame = query_from % 3
                if frame != 0:
                    frame = 3 - frame
                exons.append(
                    (query_from, query_pos, frame, sbjct_from, sbjct_pos, ali))
                ali = []
                query_from = query_pos
                sbjct_from = sbjct_pos + l_sbjct

            query_pos += l_query

        elif state == "S":
            query_pos += l_sbjct

        elif state == "I":
            pass

        elif state == "5":
            # frame is frame for this exon
            frame = query_from % 3
            if frame != 0:
                frame = 3 - frame
            exons.append(
                (query_from, query_pos, frame, sbjct_from, sbjct_pos, ali))
            ali = []
            query_pos += l_query

        elif state == "3":
            query_pos += l_query
            query_from = query_pos
            sbjct_from = sbjct_pos + l_sbjct

        sbjct_pos += l_sbjct

    # add three for the stop codon:
    if add_stop_codon:
        query_pos += 3

    frame = query_from % 3
    if frame != 0:
        frame = 3 - frame

    exons.append((query_from, query_pos, frame, sbjct_from, sbjct_pos, ali))

    return exons


def RemoveFrameShiftsFromAlignment(row_ali, col_ali, gap_char="-"):
    """remove frame shifts in an alignment.

    Frameshifts are gaps are 1, 2, 4, or 5 residues long.

    >>> RemoveFrameShiftsFromAlignment("ABC-EFG", "AB-DEFG")
    ('ABEFG', 'ABEFG')

    Arguments
    ---------
    row_ali : string
        Alignment string of row.
    col_ali : string
        Alignment string of column.
    gap_char : string
        Gap character to identify aligments.

    Returns
    -------
    new_row_ali : string
        New alignment string for row
    new_col_ali : string
        New aligment string for column
    """

    match_string = "[^%s]%s+[^%s]" %\
        (gap_char, gap_char, gap_char)

    positions = []
    for x in re.finditer(match_string, row_ali):
        positions.append((x.start() + 1, x.end() - 1))

    for x in re.finditer(match_string, col_ali):
        positions.append((x.start() + 1, x.end() - 1))

    positions.sort()
    # cut and paste new alignment
    new_row_ali = []
    new_col_ali = []
    a = 0
    for first, last in positions:
        if (last - first) % 3:
            new_row_ali.append(row_ali[a:first])
            new_col_ali.append(col_ali[a:first])
            a = last

    new_row_ali.append(row_ali[a:])
    new_col_ali.append(col_ali[a:])

    return "".join(new_row_ali), "".join(new_col_ali)


def IsStopCodon(codon, stop_codons=("TAG", "TAA", "TGA")):
    """return True if codon is a known stop codon.
    """
    return codon in stop_codons


def MaskStopCodons(sequence, stop_codons=("TAG", "TAA", "TGA")):
    """mask stop codons in a sequence.

    Stop codons are masked with ``NNN``.

    Arguments
    ---------
    sequence : string
        Nucleotide sequence to mask.
    stop_codons : string
        List of known stop codons.

    Returns
    -------
    masked_sequence : string
    """
    codons = []

    for x in range(0, len(sequence), 3):
        codon = sequence[x:x + 3]
        if codon in stop_codons:
            codon = "NNN"
        codons.append(codon)
    return codons.join("")


def Alignment2DNA(alignment, query_from=0, sbjct_from=0):
    """convert a peptide2genome alignment to a nucleotide2nucleotide
    alignment.

    Instead of peptide coordinates, the alignment will be
    in codon coordinates.

    Arguments
    ---------
    aligment : list
        List of tuples of the alignment in CIGAR format.
    query_from : int
        Start position of alignment on peptide sequence.
    sbjct_from : int
        Start position of alignment on nucleotide sequence.

    Returns
    -------
    alignment : object
       The alignment as an alignlib.AlignmentVector object.
    """

    map_query2sbjct = alignlib_lite.py_makeAlignmentVector()

    # count in nucleotides for query
    query_pos = query_from * 3
    sbjct_pos = sbjct_from

    for state, l_query, l_sbjct in alignment:

        # count as nucleotides
        l_query *= 3

        if state in ("A", "B", "C"):

            if state in ("A"):
                l_query = 0
            elif state in ("B"):
                l_query = 1
            elif state in ("C"):
                l_query = 2

        elif state in ("a", "b", "c"):

            if state in ("a"):
                l_query = 0
            elif state in ("b"):
                l_query = 2
            elif state in ("c"):
                l_query = 1

        elif state == "S":
            l_query = l_sbjct

        if l_query > 0 and l_sbjct > 0:
            alignlib_lite.addDiagonal2Alignment(map_query2sbjct,
                                                query_pos, query_pos +
                                                l_query,
                                                sbjct_pos - query_pos)

        query_pos += l_query
        sbjct_pos += l_sbjct

    return map_query2sbjct


Wobble = {
    "A": "GCW",
    "C": "TGY",
    "D": "GAY",
    "E": "GAR",
    "F": "TTY",
    "G": "GGW",
    "H": "CAY",
    "I": "ATW",
    "K": "AAR",
    "L": "WTW",
    "M": "ATG",
    "N": "AAY",
    "P": "CCW",
    "Q": "CAR",
    "R": "WGW",
    "S": "WWW",
    "T": "ACW",
    'U': "TGA",        # selenocysteine
    "V": "GTW",
    "W": "TGG",
    "X": "WWW",        # X can be stop codon or masked sequence.
    "Y": "TAW",
    "*": "NNN",
    "?": "WWW",        # ? can be stop codon or masked sequence.
    "-": "---",
}

GeneticCode = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "X",  # stop
    "TAG": "X",  # stop
    "TGT": "C",
    "TGC": "C",
    "TGA": "X",  # stop
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

GeneticCodeSeleno = GeneticCode.copy()
GeneticCodeSeleno["TGA"] = "U"

GeneticCodeAA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

StopCodons = ("TAG", "TAA", "TGA")

GeneticCodeDegenerate = {
    "TC": "S",
    "CT": "L",
    "CC": "P",
    "CG": "R",
    "AC": "T",
    "GT": "V",
    "GC": "A",
    "GG": "G",
}


Degeneracy = {
    "TTT": ("F", 1, 1, 2),
    "TTC": ("F", 1, 1, 2),
    "TTA": ("L", 2, 1, 2),
    "TTG": ("L", 2, 1, 2),
    "TCT": ("S", 1, 1, 4),
    "TCC": ("S", 1, 1, 4),
    "TCA": ("S", 1, 1, 4),
    "TCG": ("S", 1, 1, 4),
    "TAT": ("Y", 1, 1, 2),
    "TAC": ("Y", 1, 1, 2),
    "TGT": ("C", 1, 1, 2),
    "TGC": ("C", 1, 1, 2),
    "TGG": ("W", 1, 1, 1),
    "CTT": ("L", 1, 1, 4),
    "CTC": ("L", 1, 1, 4),
    "CTA": ("L", 2, 1, 4),
    "CTG": ("L", 2, 1, 4),
    "CCT": ("P", 1, 1, 4),
    "CCC": ("P", 1, 1, 4),
    "CCA": ("P", 1, 1, 4),
    "CCG": ("P", 1, 1, 4),
    "CAT": ("H", 1, 1, 2),
    "CAC": ("H", 1, 1, 2),
    "CAA": ("Q", 1, 1, 2),
    "CAG": ("Q", 1, 1, 2),
    "CGT": ("R", 1, 1, 4),
    "CGC": ("R", 1, 1, 4),
    "CGA": ("R", 2, 1, 4),
    "CGG": ("R", 2, 1, 4),
    "ATT": ("I", 1, 1, 3),
    "ATC": ("I", 1, 1, 3),
    "ATA": ("I", 1, 1, 3),
    "ATG": ("M", 1, 1, 1),
    "ACT": ("T", 1, 1, 4),
    "ACC": ("T", 1, 1, 4),
    "ACA": ("T", 1, 1, 4),
    "ACG": ("T", 1, 1, 4),
    "AAT": ("N", 1, 1, 2),
    "AAC": ("N", 1, 1, 2),
    "AAA": ("K", 1, 1, 2),
    "AAG": ("K", 1, 1, 2),
    "AGT": ("S", 1, 1, 2),
    "AGC": ("S", 1, 1, 2),
    "AGA": ("R", 2, 1, 2),
    "AGG": ("R", 2, 1, 2),
    "GTT": ("V", 1, 1, 4),
    "GTC": ("V", 1, 1, 4),
    "GTA": ("V", 1, 1, 4),
    "GTG": ("V", 1, 1, 4),
    "GCT": ("A", 1, 1, 4),
    "GCC": ("A", 1, 1, 4),
    "GCA": ("A", 1, 1, 4),
    "GCG": ("A", 1, 1, 4),
    "GAT": ("D", 1, 1, 2),
    "GAC": ("D", 1, 1, 2),
    "GAA": ("E", 1, 1, 2),
    "GAG": ("E", 1, 1, 2),
    "GGT": ("G", 1, 1, 4),
    "GGC": ("G", 1, 1, 4),
    "GGA": ("G", 1, 1, 4),
    "GGG": ("G", 1, 1, 4),
}

GAP_CHAR = "-"

DegeneracyAA = {
    "F": 2,
    "L": 6,
    "S": 4,
    "Y": 2,
    "C": 2,
    "W": 1,
    "P": 4,
    "H": 2,
    "Q": 2,
    "R": 4,
    "I": 3,
    "M": 1,
    "T": 4,
    "N": 2,
    "K": 2,
    "S": 2,
    "R": 2,
    "V": 4,
    "A": 4,
    "D": 2,
    "E": 2,
    "G": 4,
    "X": 0,
    "-": 0,
    ".": 0,
}

AMBIGUOUS_CODES_NA = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'N': 'ACGT',
    'X': 'ACGT',
    'U': 'T',
    'R': 'AG',
    'Y': 'CT',
    'M': 'AC',
    'K': 'GT',
    'S': 'CG',
    'W': 'AT',
    'H': 'ACT',
    'B': 'CGT',
    'V': 'ACG',
    'D': 'AGT'}

REVERSE_AMBIGUOUS_NA = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AA': 'A', 'CC': 'C', 'GG': 'G', 'TT': 'T', 'UU': 'U',
    'AG': 'R', 'GA': 'R',
    'CT': 'Y', 'TC': 'Y',
    'AC': 'M', 'CA': 'M',
    'GT': 'K', 'TG': 'K',
    'CG': 'S', 'GC': 'S',
    'AT': 'W', 'TA': 'W',
    'ACT': 'H', 'ATC': 'H', 'CAT': 'H', 'CTA': 'H', 'TAC': 'H', 'TCA': 'H',
    'CGT': 'B', 'CTG': 'B', 'GCT': 'B', 'GTC': 'B', 'TCG': 'B', 'TGC': 'B',
    'ACG': 'V', 'AGC': 'V', 'CAG': 'V', 'CGA': 'V', 'GAC': 'V', 'GCA': 'V',
    'AGT': 'D', 'ATG': 'D', 'GAT': 'D', 'GTA': 'D', 'TAG': 'D', 'TGA': 'D',
}

ENCODE_GENOTYPE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AA': 'A', 'CC': 'C', 'GG': 'G', 'TT': 'T', 'UU': 'U',
    'AG': 'r', 'GA': 'R',
    'CT': 'y', 'TC': 'Y',
    'AC': 'm', 'CA': 'M',
    'GT': 'k', 'TG': 'K',
    'CG': 's', 'GC': 'S',
    'AT': 'w', 'TA': 'W',
}

DECODE_GENOTYPE = {
    'A': 'AA',
    'C': 'CC',
    'G': 'GG',
    'T': 'TT',
    'r': 'AG', 'R': 'AG',
    'y': 'CT', 'Y': 'CT',
    'm': 'AC', 'M': 'AC',
    'k': 'GT', 'K': 'GT',
    's': 'CG', 'S': 'CG',
    'w': 'AT', 'W': 'AT',
}


def encodeGenotype(code):
    '''encode genotypes like GG, GA into a one-letter code.
    The returned code is lower case if code[0] < code[1], otherwise
    it is uppercase.
    '''
    return ENCODE_GENOTYPE[code.upper()]


def decodeGenotype(code):
    '''decode single letter genotypes like m, M into two letters.
    This is the reverse operation to :meth:`encodeGenotype`.
    '''
    return DECODE_GENOTYPE[code]


def resolveAmbiguousNA(code):
    '''resolve ambiguous nucleic acid letters.
    '''
    return AMBIGUOUS_CODES_NA[code.upper()]


def resolveReverseAmbiguousNA(genotype):
    '''map a genotype to a single letter amino acid amiguous code,
    for example, CT -> Y.
    '''
    return REVERSE_AMBIGUOUS_NA[genotype.upper()]


def GetMapAA2Codons():
    """returns a map of amino acids to codons

    No stop codons.
    ."""
    map_aa2codons = {}

    for codon, aa in list(GeneticCodeAA.items()):
        if aa not in map_aa2codons:
            map_aa2codons[aa] = []
        map_aa2codons[aa].append(codon)
    return map_aa2codons


def GetDegeneracy(codon):
    return Degeneracy[codon.upper()]


def IsStopCodon(codon, stop_codons=("TAG", "TAA", "TGA")):
    return codon in stop_codons


def MapCodon2AA(codon, is_seleno=False, ignore_n=True):
    '''map a codon to an amino acid using the standard translation
    tables

    The mapping returns gaps as gaps and will return an amino acid
    for incomplete codons if there is unambiguous mapping.

    If ``is_seleno`` is set, the codon is translated for a selenoprotein.

    If ``ignore_n`` is set, codons with ``n`` are returned
    as ``?`` in order to distinguish them from stop codons.

    Amino acids are returned as upper-case letters.
    '''

    codon = codon.upper()

    codon = re.sub("[.-]", "", codon)
    if len(codon) == 0:
        return GAP_CHAR

    if is_seleno:
        code = GeneticCodeSeleno
    else:
        code = GeneticCode

    if codon in code:
        return code[codon]
    elif codon == "---":
        return "-"
    elif codon[:2] in GeneticCodeDegenerate:
        # check for four-fold degenerate codons,
        # if they can be mapped (ENSEMBL does it).
        return GeneticCodeDegenerate[codon[:2]]
    elif ignore_n and ("n" in codon or "N" in codon):
        return "?"
    else:
        return "X"


def Protein2Wobble(s):

    c = []
    for x in s:
        c.append(Wobble[x])
    return "".join(c)


def Alignment2PeptideAlignment(alignment,
                               query_from=0,
                               sbjct_from=0,
                               genomic_sequence=None):
    """convert a Peptide2DNA aligment to a Peptide2Peptide alignment.

    How to handle frameshifts?
    """

    map_query2sbjct = alignlib_lite.py_makeAlignmentVector()

    query_pos = query_from
    sbjct_pos = 0
    sbjct_genome_pos = sbjct_from
    sbjct_residues = []
    codon = ""

    for state, l_query, l_sbjct in alignment:

        query_increment = 0
        sbjct_increment = 0

        if state == "M":

            query_increment = l_query
            sbjct_increment = l_sbjct / 3
            if genomic_sequence:
                codon = genomic_sequence[
                    sbjct_genome_pos:sbjct_genome_pos + l_sbjct]

        elif state == "S":
            if l_query:
                sbjct_increment = 1
                query_increment = 1

            if genomic_sequence:
                codon += genomic_sequence[sbjct_genome_pos:
                                          sbjct_genome_pos + l_sbjct]

        elif state == "G":
            query_increment = l_query
            sbjct_increment = l_sbjct / 3
            if genomic_sequence:
                codon += genomic_sequence[sbjct_genome_pos:
                                          sbjct_genome_pos + l_sbjct]

        elif state == "P":
            # only increment query, sbjct does not advance.
            query_increment = l_query

        if query_increment and sbjct_increment:
            alignlib_lite.py_addDiagonal2Alignment(map_query2sbjct,
                                                   query_pos, query_pos +
                                                   query_increment,
                                                   sbjct_pos - query_pos)

        if sbjct_increment and genomic_sequence:
            for x in range(0, len(codon), 3):
                sbjct_residues.append(MapCodon2AA(codon[x:x + 3]))
            codon = ""

        query_pos += query_increment
        sbjct_pos += sbjct_increment

        sbjct_genome_pos += l_sbjct

    return map_query2sbjct, "".join(sbjct_residues)


def translate(sequence,
              is_seleno=False,
              prefer_lowercase=True,
              ignore_n=False,
              ):
    '''convert DNA sequence to a peptide sequence

    If ``is_seleno`` is set, "TGA" codons are treated as
    encoding for selenocysteine.

    If ``ignore_n`` is set, codons with ``n`` are returned
    as ``?`` in order to distinguish them from stop codons.

    '''
    residues = []

    for x in range(0, len(sequence), 3):
        if prefer_lowercase:
            is_lower = False
            for c in sequence[x:x + 3]:
                is_lower = is_lower or c in "acgtnx"
        else:
            is_lower = True
            for c in sequence[x:x + 3]:
                is_lower = is_lower and c in "acgtnx"

        if is_lower:
            residues.append(MapCodon2AA(sequence[x:x + 3].upper(),
                                        is_seleno=is_seleno,
                                        ignore_n=ignore_n).lower())
        else:
            residues.append(MapCodon2AA(sequence[x:x + 3].upper(),
                                        is_seleno=is_seleno,
                                        ignore_n=ignore_n).upper())

    return "".join(residues)


def TranslateDNA2Protein(*args, **kwargs):
    """convert a DNA sequence to a peptide sequence.
    keep case.

    deprecated - use :meth:`translate` instead.
    """
    return translate(*args, **kwargs)


def Alignment2CDNA(alignment,
                   query_from=0,
                   sbjct_from=0,
                   genome=None,
                   remove_frameshifts=0):
    """build cDNA sequence from genomic fragment and
    return alignment of query to it.
    """

    fragments = []
    sbjct_pos = 0
    map_query2sbjct = alignlib_lite.py_makeAlignmentVector()

    # count in nucleotides for query
    query_pos = query_from * 3
    sbjct_pos = sbjct_from
    # position in cDNA
    cdna_pos = 0
    for state, l_query, l_sbjct in alignment:

        # count as nucleotides
        l_query *= 3

        keep = False

        if state == "M":
            keep = True
        elif state == "S":
            l_query = l_sbjct
            keep = True
        elif state == "F" and not remove_frameshifts:
            keep = True
        elif state == "G":
            if l_sbjct > 0:
                keep = True
        elif state == "P":
            keep = False

        if keep:
            if genome:
                fragments.append(genome[sbjct_pos:sbjct_pos + l_sbjct])

            if l_query > 0 and l_sbjct > 0:
                alignlib_lite.py_addDiagonal2Alignment(map_query2sbjct,
                                                       query_pos,
                                                       query_pos + l_query,
                                                       cdna_pos - query_pos)
            cdna_pos += l_sbjct

        query_pos += l_query
        sbjct_pos += l_sbjct

    return map_query2sbjct, fragments.join("")


def GetExon(exons, first_aa):

    first_na = (first_aa - 1) * 3
    combined_na = max(exons[0].mPeptideFrom, first_na)

    # add not covered positions
    g = exons[0].mGenomeFrom + max(0, combined_na - exons[0].mPeptideFrom)

    # add introns
    e = 0
    while exons[e].mPeptideTo < combined_na:
        e += 1
        g += exons[e].mGenomeFrom - exons[e - 1].mGenomeTo

    return g, e, combined_na / 3 + 1


def Exons2Alignment(exons):
    """build an cigar alignment string from a list of exons.
    """

    exons.sort()

    phase = 0
    alignment = []
    last_to = 0
    nexon = 0

    for this_from, this_to in exons:

        nexon += 1

        # print alignment
        # print this_from, this_to

        lintron = this_from - last_to
        lexon = this_to - this_from

        if lintron < 0:
            raise ValueError("error: invalid intron length: %i" % lintron)

        if last_to:

            if (phase != 0):
                missing = 3 - phase

                alignment.append(("S", 0, phase))

            if lintron >= 4:
                # deal with real introns:
                # at least 4 nucleotides for splice signal
                alignment.append(("5", 0, 2))
                alignment.append(("I", 0, lintron - 4))
                alignment.append(("3", 0, 2))

            else:
                # if phase is 0, regard as frameshift
                if phase == 0:
                    alignment.append(("F", 0, lintron))
                else:
                    # deal with frameshift introns
                    # add dummy entries for splice signals
                    alignment.append(("5", 0, 0))
                    alignment.append(("I", 0, lintron))
                    alignment.append(("3", 0, 0))

            if phase != 0:
                # deal with very small exons

                missing = 3 - phase

                # deal with small exons that can not
                # accomodate a full codon
                if missing > lexon:
                    alignment.append(("S", 1, missing))
                else:
                    alignment.append(("S", 1, missing))
                    phase = 0
                    this_from += missing

                # raise ValueError ("expecting a split codon,
                # but small exon %i-%i can not accomodate phase
                # %i in exons %s " % (this_from, this_to, phase, str(exons)))

        l = this_to - this_from

        if l < 0:
            # note: sometimes the last exon just contains
            # a split codon and the stop codon.
            # Thus: do not raise an error if it is the last exon
            if nexon == len(exons):
                break
            else:
                pass
                # raise ValueError ("error: negative length
                # for aligned residues at exon %i-%i in exons
                # %s" % (this_from, this_to, str(exons)))

        phase = l % 3

        ncodons = (l - phase) / 3
        alignment.append(("M", ncodons, ncodons * 3))

        last_to = this_to

    return alignment


def AlignmentProtein2CDNA(src, exons1=None, exons2=None):
    """convert a peptide alignment to a nucleotide
    alignment.

    multiplies coordinates with 3.
    Insert introns.

    Note: alignment starts at 1
    """
    result = src.getNew()

    # enter correct exons (for example, if alignment is not fully covering or
    # gene structure is incomplete)
    gx, exon_id1, first_x = GetExon(exons1, src.getRowFrom())
    gy, exon_id2, first_y = GetExon(exons2, src.getColFrom())

    # print "x=", gx, exon_id1, first_x
    # print "y=", gy, exon_id2, first_y

    last_x = first_x - 1
    last_y = first_y - 1

    exon1_to = exons1[exon_id1].mGenomeTo
    exon2_to = exons2[exon_id2].mGenomeTo

    # only iterate in alignment
    for x in range(first_x,
                   min(src.getRowTo() + 1, exons1[-1].mPeptideTo / 3)):

        y = src.mapRowToCol(x)
        if not y:
            continue
        if y < first_y:
            continue

        gx += (x - last_x - 1) * 3
        gy += (y - last_y - 1) * 3

        # print "# x=", x, "last_x=", last_x, "y=", y, "last_y=", last_y

        for i in range(0, 3):

            # print i, "\t", gx, exons1[exon_id1].mGenomeTo, exon_id2, "\t",
            # gy, exon_id2, exons2[exon_id2].mGenomeTo

            if gx >= exon1_to:
                exon_id1 += 1
                if exon_id1 < len(exons1):
                    gx += exons1[exon_id1].mGenomeFrom - exon1_to
                    exon1_to = exons1[exon_id1].mGenomeTo

            if gy >= exon2_to:
                exon_id2 += 1
                if exon_id2 < len(exons2):
                    gy += exons2[exon_id2].mGenomeFrom - exon2_to
                    exon2_to = exons2[exon_id2].mGenomeTo

            gx += 1
            gy += 1
            result.addPairExplicit(gx, gy, 0)

        last_x = x
        last_y = y

    return result


def GetDegenerateSites(seq1, seq2,
                       degeneracy=4,
                       position=3):
    """returns two new sequenes containing only degenerate sites.

    Only unmutated positions are counted.
    """

    new_seq1 = []
    new_seq2 = []
    for x in range(0, len(seq1), 3):

        c1 = seq1[x:x + 3]
        c2 = seq2[x:x + 3]

        if c1 in GeneticCodeAA and c2 in GeneticCodeAA:

            if GeneticCodeAA[c1] == GeneticCodeAA[c2]:

                if Degeneracy[c1][position] == degeneracy \
                   and Degeneracy[c2][position] == degeneracy:
                    new_seq1.append(c1[position - 1])
                    new_seq2.append(c2[position - 1])

    return "".join(new_seq1), "".join(new_seq2)


class SequencePairInfo:

    """the first characters are ACGT."""

    def __init__(self):
        self.mNIdentical = 0
        self.mNAligned = 0
        self.mNDifferent = 0
        self.mNTransitions = 0
        self.mNTransversions = 0
        self.mNUnaligned1 = 0
        self.mNUnaligned2 = 0
        self.mAlphabet = 0
        self.mMapChar2Pos = {}
        self.mMatrix = []

    def getGCContent(self):
        """return GC content."""
        if not self.mNAligned:
            raise ValueError("no data for calculating GC content")

        cp = self.mMapChar2Pos['C']
        gp = self.mMapChar2Pos['G']
        # two row and two columns, subtract to avoid double counting
        gc = numpy.sum(self.mMatrix[0:4, cp] +
                       self.mMatrix[0:4, gp] +
                       self.mMatrix[cp, 0:4] +
                       self.mMatrix[gp, 0:4]) \
            - self.mMatrix[cp, cp] \
            - self.mMatrix[gp, gp] \
            - self.mMatrix[cp, gp] \
            - self.mMatrix[gp, cp]

        return float(gc) / 2.0 / self.mNAligned

    def __str__(self):
        return "\t".join(map(str, (self.mNIdentical,
                                   self.mNAligned, self.mNDifferent,
                                   self.mNTransitions, self.mNTransversions,
                                   self.mNUnaligned1, self.mNUnaligned2)))

    def getHeader(self):
        return "\t".join(("identical", "aligned", "different",
                          "transitions", "transversions",
                          "unaligned1", "unaligned2"))


class SequencePairInfoCodons(SequencePairInfo):

    """the first characters are ACGT."""

    def __init__(self):
        SequencePairInfo.__init__(self)
        self.mNNonSynonymous = 0
        self.mNSynonymous = 0

    def __str__(self):

        return "\t".join((SequencePairInfo.__str__(self),
                          str(self.mNNonSynonymous),
                          str(self.mNSynonymous)))

    def getHeader(self):
        return "\t".join((SequencePairInfo.getHeader(self), "nnonsyn", "nsyn"))


def AlignedPair2SubstitutionMatrix(seq1, seq2, alphabet):
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


def CalculatePairIndices(seq1, seq2, gap_char="-", with_codons=False):
    """returns number of idential and transitions/transversions substitutions
    in the alignment.

    If with-codons = True, synonymous and nonsynonymous changes will
    be recorded as well. The routine assumes no frame-shifts and will
    count more than one change as non-synonymous.

    """
    alphabet = "ACGT" + gap_char

    map_char2pos = {}
    for x in alphabet:
        map_char2pos[x] = len(map_char2pos)

    # build coordinates for various substitution subsets
    transitions, transversions = [], []
    for x in ("AG", "GA", "CT", "TC"):
        transitions.append((map_char2pos[x[0]], map_char2pos[x[1]]))

    for x in ("AT", "TA", "GT", "TG", "GC", "CG", "AC", "CA"):
        transversions.append((map_char2pos[x[0]], map_char2pos[x[1]]))

    matrix = AlignedPair2SubstitutionMatrix(seq1, seq2, alphabet)
    matrix_acgt = matrix[0:4, 0:4]

    if with_codons:
        result = SequencePairInfoCodons()
    else:
        result = SequencePairInfo()

    result.mMatrix = matrix
    result.mMapChar2Pos = map_char2pos
    result.mNAligned = numpy.sum(numpy.ravel(matrix_acgt))
    result.mNIdentical = numpy.sum(numpy.trace(matrix_acgt))
    result.mNTransitions = numpy.sum([matrix[x] for x in transitions])
    result.mNTransversions = numpy.sum([matrix[x] for x in transversions])
    result.mNDifferent = result.mNAligned - result.mNIdentical
    result.mNUnaligned1 = numpy.sum(numpy.ravel(matrix[0:4, 4]))
    result.mNUnaligned2 = numpy.sum(numpy.ravel(matrix[4, 0:4]))

    if with_codons:
        nsyn, nnon = 0, 0
        pairs = list(zip(seq1, seq2))
        for x in range(len(pairs)):
            a, b = pairs[x]
            if a != b:

                l = (x // 3) * 3
                c1 = MapCodon2AA(seq1[l:l + 3])
                c2 = MapCodon2AA(seq2[l:l + 3])

                if c1 == GAP_CHAR or c2 == GAP_CHAR:
                    continue

                # print x, a, b, l, c1, c2, seq1[l:l+3], seq2[l:l+3], c1 == c2
                if c1 == c2:
                    nsyn += 1
                else:
                    nnon += 1

        result.mNSynonymous = nsyn
        result.mNNonSynonymous = nnon

    return result


def makeSubstitutionMatrix(type="EMBOSS"):
    """make alignator with DNA substitution matrix.

    EMBOSS style matrix:
    identity = 5
    mismatch = -4
    gop = -16
    gep = -4

    ClustalW style matrix:
    match = 1 mismatch = 0
    gop = -10 gep = -0.1
    """

    matrix = []

    if type == "emboss":
        match = 5
        mismatch = -4
        gop = -16.0
        gep = -4.0
    elif type == "blastn":
        match = 1
        mismatch = -3
        gop = -5.0
        gep = -2.0
    elif type == "clustal":
        match = 1
        mismatch = 0
        gop = -10.0
        gep = -0.1
    elif type == "iub":
        match = 1.9
        mismatch = 0
        gop = -10.0
        gep = -0.1
    elif type == "trans":
        match = 1
        mismatch = -5
        gop = -10.0
        gep = -0.1
    elif type == "backtrans":
        # mismatches have to be penalized - there should be none.
        match = 1
        mismatch = -10
        gop = -10.0
        gep = -10.0
    else:
        raise ValueError("unknown MATRIX")

    # build empty matrix
    for m in range(21):
        matrix.append([mismatch] * 21)

    # A <-> A
    matrix[0][0] = match
    # C <-> C
    matrix[1][1] = match
    # G <-> G
    matrix[5][5] = match
    # T <-> T
    matrix[16][16] = match

    if type == "trans":
        # A <-> G
        matrix[0][5] = matrix[5][0] = match
        # C <-> T
        matrix[1][16] = matrix[16][1] = match
    elif type == "backtrans":
        # AGCT <-> W
        for x in (0, 1, 5, 16):
            matrix[x][18] = matrix[18][x] = match
        # AG <-> V
        for x in (0, 5):
            matrix[x][17] = matrix[17][x] = match
        # CT <-> Y
        for x in (1, 16):
            matrix[x][19] = matrix[19][x] = match
        # N <-> N
        for x in (0, 1, 5, 11, 16, 17, 19):
            matrix[x][11] = matrix[11][x] = match

    handle_tmpfile, filename_tmpfile = tempfile.mkstemp()
    for m in matrix:
        os.write(handle_tmpfile, "\t".join(list(map(str, m))) + "\n")
    os.close(handle_tmpfile)

    smatrix = alignlib_lite.py_readSubstitutionMatrixAA(filename_tmpfile)
    os.remove(filename_tmpfile)
    return smatrix, gop, gep


def CalculateRCSUValuesFromCounts(counts, pseudo_counts=0):
    """calculate RCSU values for codons.

    RCSU = relative frequency / uniform frequency
    """


def CalculateCodonFrequenciesFromCounts(counts, pseudo_counts=0):
    """calculate codon frequencies from codon counts per amino acid.
    pseudo_counts are added if desired.
    """
    map_aa2codons = GetMapAA2Codons()
    weights = {}

    for aa, codons in list(map_aa2codons.items()):

        total_counts = reduce(
            lambda x, y: x + y, (counts[x] + pseudo_counts for x in codons))
        if total_counts > 0:
            for x in codons:
                weights[x] = float(counts[x] + pseudo_counts) / total_counts
        else:
            for x in codons:
                weights[x] = 0.0

    for codon in StopCodons:
        weights[codon] = 0.0

    return weights


def CalculateCAIWeightsFromCounts(counts, pseudo_counts=0):
    """calculate CAI weights from codon counts.
    pseudo_counts are added if desired.
    """
    map_aa2codons = GetMapAA2Codons()
    weights = {}

    for aa, codons in list(map_aa2codons.items()):

        max_counts = max(counts[x] for x in codons) + pseudo_counts

        for x in codons:
            weights[x] = float(counts[x] + pseudo_counts) / float(max_counts)

    for codon in StopCodons:
        weights[codon] = 0.0

    return weights


def IsJunk(contig):
    """returns true, if contigs is likely to be junk.

    This is done by name matching. Junk contigs contain either
    one of the following:

    random, unknown, chrU, chU.
    """

    c = contig.lower()
    negative_list = "random", "unknown", "chru", "chu"

    for x in negative_list:
        if re.search(x, c):
            return True

    return False


def CountCodons(sequence):
    """count the codons in a sequence."""
    if len(sequence) % 3 != 0:
        raise "sequence is not multiple of 3: %i" % len(sequence)

    counts = {}
    for c in list(GeneticCodeAA.keys()):
        counts[c] = 0
    for codon in [sequence[x:x + 3] for x in range(0, len(sequence), 3)]:
        try:
            counts[codon] += 1
        except KeyError:
            # skip stop-codons
            pass

    return counts


def GetUniformCodonUsage():
    """get list of frequencies for codons expected for uniform codon usage."""

    frequencies = {}

    aas = {}

    for codon, aa in list(GeneticCodeAA.items()):
        if aa not in aas:
            aas[aa] = 0
        aas[aa] += 1

    for codon, aa in list(GeneticCodeAA.items()):
        frequencies[codon] = 1.0 / aas[aa]

    return frequencies


def GetBiasedCodonUsage(bias=1.0):
    """get list of frequencies for codons according to some bias.

    The first codon for each aa is the most biased, all others are less biased.

    The ratio determines the relative bias between the first and all other
    codons. 0.0 is no bias, 1.0 is complete bias.
    """

    frequencies = {}

    aas = {}

    for codon, aa in list(GeneticCodeAA.items()):
        if aa not in aas:
            aas[aa] = 0
        aas[aa] += 1

    for codon, aa in list(GeneticCodeAA.items()):
        frequencies[codon] = 1.0 / aas[aa]

    return frequencies


def IsNegativeStrand(strand):
    return str(strand) in ("-", "0", "-1")


def IsPositiveStrand(strand):
    return not IsNegativeStrand(strand)


def convertStrand(strand):
    """convert various strand notations into [+-.].
    """
    s = str(strand)
    if s in ("-", "0", "-1"):
        return "-"
    elif s in ("+", "1"):
        return "+"
    else:
        return "."

INTRON_TYPES = (("U2-GT/AG", "GT", "AG"),
                ("U2-nc-GC/AG", "GC", "AG"),
                ("U12-AT/AC", "AT", "AC"))


def GetIntronType(sequence, both_strands=False):
    """return intron type for an intronic sequence.

    If both_strands is True, both strands are checked.
    """

    if both_strands is False:
        for name, prime5, prime3 in INTRON_TYPES:
            if sequence[:len(prime5)].upper() == prime5 and \
                    sequence[-len(prime3):].upper() == prime3:
                return name, prime5, prime3
        else:
            return "unknown", sequence[:5], sequence[-5:]
    else:
        r = complement(sequence)

        for name, prime5, prime3 in INTRON_TYPES:
            if sequence[:len(prime5)].upper() == prime5 and \
                    sequence[-len(prime3):].upper() == prime3:
                return name, prime5, prime3
            elif r[:len(prime5)].upper() == prime5 and \
                    r[-len(prime3):].upper() == prime3:
                return name, prime3, prime5
        else:
            return "unknown", sequence[:5], sequence[-5:]


def printPrettyAlignment(seq1, *args):
    """print a pretty alignment."""

    for other in args:
        assert len(seq1) == len(other)

    seqrow, otherrows = [], []
    nother = len(args)
    for x in range(nother):
        otherrows.extend([[], []])

    for p in range(len(seq1)):
        if p > 0 and p % 3 == 0:
            seqrow.append(" ")
            for x in range(2 * nother):
                otherrows[x] += " "
            if len(seqrow) > 120:
                print("".join(seqrow))
                for x in range(2 * nother):
                    print("".join(otherrows[x]))

                print()
                seqrow, otherrows = [], []
                nother = len(args)
                for x in range(nother):
                    otherrows.extend([[], []])

        l1 = len(seq1[p])
        c1 = seq1[p]
        max_len = max([len(x[p]) for x in args])
        seqrow.append(c1 + " " * (max_len - 11))

        for x in range(nother):

            c2 = args[x][p]
            l2 = len(c2)
            if c1 == c2:
                code = "|"
            elif l1 < l2:
                code = "-"
            elif l2 > l2:
                code = "+"
            else:
                code = ":"
            otherrows[x * 2].append(code + " " * (max_len - len(code)))

            if l2 == 0:
                c2, l2 = "-", 1
            otherrows[x * 2 + 1].append(c2 + " " * (max_len - l2))

    print("".join(seqrow))
    for x in range(2 * nother):
        print("".join(otherrows[x]))


def ReadPeptideSequences(infile, filter=None, as_array=False,
                         regex_identifier=None):
    """read peptide sequence from fasta infile.
    """
    sequences = ParseFasta2Hash(
        infile, filter, regex_identifier=regex_identifier)

    if not as_array:
        for k in list(sequences.keys()):
            sequences[k] = sequences[k][:]
    return sequences


def ParseFasta2Hash(infile, filter=None, regex_identifier=None):
    """read fasta formatted sequences file and build a hash.

    Keys are all characters before the first whitespace in the
    description line.

    Previously, if the key contained a ":", everything before the ":"
    was removed.  This is not true any more.

    Use array for higher space efficiency.

    If regex_identifier is given, this is used to extract the identifier
    from the fasta description line.

    """
    parsed = {}
    key = None
    p = AString.AString()
    if regex_identifier:
        rx = regex_identifier
    else:
        rx = re.compile("^(\S+)")

    for line in infile:
        if line[0] == "#":
            continue
        if not line:
            continue

        if line[0] == ">":
            if key:
                if not filter or key in filter:
                    parsed[key] = p

            key = rx.search(line[1:-1]).groups()[0]
            p = AString.AString()
            continue

        p.extend(AString.AString(re.sub("\s", "", line[:-1])))

    if not filter or key in filter:
        if key:
            parsed[key] = p
    return parsed
