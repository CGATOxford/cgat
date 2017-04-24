"""MAST.py - Parser for MAST files
==================================

:Tags: Python

As of biopython 1.5.6, the MAST parser is broken.
"""

import re
import math
import collections

import numpy


class MAST:
    # number of matches

    def __init__(self):
        self.motifs = []
        self.matches = []


class Motif:
    pass


class Match:

    """a MAST entry.
    """
    header = "\t".join(("id", "description", "length", "pvalue", "evalue",
                        "nmatches", "diagram", "start", "end", "positions", "motifs", "directions"))

    def __init__(self):
        pass

    def __str__(self):
        return "\t".join(map(str,
                             (self.id,
                              self.description,
                              self.length,
                              self.pvalue,
                              self.evalue,
                              self.nmotifs,
                              self.diagram,
                              self.start,
                              self.end,
                              ",".join(map(str, self.positions)),
                              ",".join(map(str, self.motifs)),
                              ",".join(self.directions))))


def parse(infile):
    '''parse verbose MAST output.'''

    line = None
    for line in infile:
        if line.startswith("DATABASE AND MOTIFS"):
            break

    mast = MAST()
    # empty input
    if not line:
        return mast

    keep = False
    for line in infile:
        l = line.strip()
        if l.startswith("MOTIF WIDTH"):
            keep = True
            continue
        if not keep:
            continue
        if keep:
            if l.startswith("-"):
                continue
            if not l:
                break
            motif = Motif()
            motif.id, motif.length, motif.sequence = re.split("\s+", l)
            motif.length = int(motif.length)
            mast.motifs.append(motif)

    for line in infile:
        if line.startswith("SECTION III"):
            break

    assert line, "could not find tag 'SECTION III'"

    def __myiter(infile):
        keep = False
        do_yield = False
        for line in infile:
            if line.startswith("*"):
                continue
            if keep:
                l = line[:-1].strip()
                if do_yield:
                    # collect continuation of DIAGRAM line
                    if len(line) > 2:
                        block[-1] += l
                        continue
                    else:
                        yield block
                        keep = False
                        do_yield = False

                block.append(l)
                # start end of block (DIAGRAM might be multiline)
                if l.startswith("DIAGRAM"):
                    do_yield = True
            else:
                if not re.match("[\d\s]", line):
                    # start of new block
                    keep = True
                    block = [line[:-1].strip()]

    for block in __myiter(infile):
        assert len(block) == 4, "invalid block: %s" % "".join(block)
        match = Match()
        match.id = block[0]
        match.description = block[1]
        match.length, match.pvalue, match.evalue = re.match(
            "LENGTH =\s*(\d+)\s*COMBINED P-VALUE =\s*(\S+)\s*E-VALUE =\s*(\S+)", block[2]).groups()
        match.length, match.pvalue, match.evalue = \
            int(match.length), float(match.pvalue), float(match.evalue)
        match.diagram = re.match("DIAGRAM:\s+(\S+)", block[3]).groups()[0]
        parts = match.diagram.split("_")
        match.nmotifs = len(parts) // 2
        match.positions, match.directions, match.motifs = [], [], []
        match.start = ""
        match.end = ""
        if match.nmotifs > 0:
            pos = 0
            last_end = 0
            try:
                for part in parts:
                    is_motif = part[0] in "[<"
                    if is_motif:
                        code = part[1:-1]
                        match.positions.append(pos)
                        match.directions.append(code[0])
                        motif_id = int(code[1:])
                        match.motifs.append(motif_id)
                        pos += mast.motifs[motif_id - 1].length
                        last_end = pos
                    else:
                        pos += int(part)
            except ValueError:
                raise ValueError(
                    "ParsingError for motif diagram %s" % match.diagram)

            match.start = match.positions[0]
            assert last_end > 0
            match.end = last_end

        mast.matches.append(match)

    return mast


def frequencies2logodds(counts, background_frequencies=None):
    '''write a motif from *counts* to outfile.

    Counts should be a numpy matrix with *nalphabet* columns and *motif_width* rows.
    '''

    motif_width, nalphabet = counts.shape

    if background_frequencies is None:
        # use uniform background
        f = 1.0 / nalphabet
        background_frequencies = [f] * nalphabet

    assert len(background_frequencies) == nalphabet
    background_freqs = numpy.array(background_frequencies)
    logodds_matrix = numpy.ones((motif_width, nalphabet), dtype=numpy.float)

    for x in range(motif_width):
        # divide by backgound_freqs
        # divide by column total
        # take log
        logodds_matrix[x] = numpy.log(
            counts[x] / counts[x].sum() / background_freqs)

    return logodds_matrix


def writeMast(outfile, logodds_matrix, alphabet):
    '''output logodds matrix in MAST format.'''

    motif_width, nalphabet, = logodds_matrix.shape
    assert len(alphabet) == nalphabet

    outfile.write("ALPHABET= %s\n" % alphabet)
    outfile.write("log-odds matrix: alength= %i w= %i\n" %
                  (nalphabet, motif_width))
    for row in logodds_matrix:
        outfile.write(" ".join(
            ["%5.3f" % x for x in row]) + "\n")

    outfile.write("\n")


def writeTomTom(outfile, counts_matrix, header=False):
    '''output counts matrix in tomtom format.

    output counts with columns as motif positions
    and rows as alphabet.
    '''

    if header:
        outfile.write(
            '--------------------------------------------------------------------------------\n')
        outfile.write('        Motif 1 position-specific probability matrix\n')
        outfile.write(
            '--------------------------------------------------------------------------------\n')
        outfile.write(
            "letter-probability matrix: alength= %i w= %i nsites= 18 E= 1.0\n" % counts_matrix.shape)

    for row in numpy.transpose(counts_matrix):
        outfile.write(" ".join(
            ["%6i" % x for x in row]) + "\n")

    outfile.write("\n")


def sequences2motif(outfile, sequences, background_frequencies=None, format="MAST"):
    '''write a motif defined by a collection of sequences to outfile.'''

    # collect all letters
    counts = collections.defaultdict(int)
    nsequences = len(sequences)
    seqs = [re.sub("\s", "", x.upper()) for x in sequences]
    for s in seqs:
        for c in s:
            counts[c] += 1

    lengths = [len(x) for x in seqs]
    motif_width = min(lengths)
    assert max(lengths) == motif_width, "seqs have to have the same length"
    alphabet = "".join(sorted(counts.keys()))
    lalphabet = len(alphabet)
    char2index = dict([(y, x) for x, y in enumerate(alphabet)])

    if background_frequencies is None:
        # use uniform background
        f = 1.0 / lalphabet
        background_frequencies = [f] * lalphabet

    if format == "MAST":
        outfile.write("ALPHABET= %s\n" % alphabet)
        outfile.write("log-odds matrix: alength= %i w= %i\n" %
                      (lalphabet, motif_width))
        for r in range(0, motif_width):
            # use a pseudocount of 1
            row_counts = [1] * lalphabet
            for s in seqs:
                row_counts[char2index[s[r]]] += 1
            outfile.write(" ".join(
                ["%5.3f" % (math.log(float(row_counts[x]) / nsequences / background_frequencies[x]))
                 for x in range(lalphabet)]) + "\n")
    elif format == "probability":
        for r in range(0, motif_width):
            # use a pseudocount of 1
            row_counts = [1] * lalphabet
            for s in seqs:
                row_counts[char2index[s[r]]] += 1
            outfile.write(" ".join(
                ["%5.3f" % (float(row_counts[x]) / nsequences)
                 for x in range(lalphabet)]) + "\n")

    elif format == "TOMTOM":
        counts = numpy.zeros((motif_width, lalphabet))
        for x in range(0, motif_width):
            for s in seqs:
                counts[x, char2index[s[x]]] += 1

        for y in range(lalphabet):
            for x in range(motif_width):
                outfile.write(" %6i" % counts[x, y])
            outfile.write("\n")
