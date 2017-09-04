'''
Motifs.py - 
======================================================

:Tags: Python

Code
----

'''
import collections
from CGAT import Genomics as Genomics
from CGAT import FastaIterator as FastaIterator


def countMotifs(infile, motifs):
    '''find regular expression *motifs* in
    sequences within fasta formatted *infile*.
    '''

    it = FastaIterator.FastaIterator(infile)
    positions = []
    while 1:
        try:
            seq = next(it)
        except StopIteration:
            break
        if not seq:
            break

        rseq = Genomics.complement(seq.sequence)
        lsequence = len(seq.sequence)
        pos = []
        for motif, pattern in motifs:

            for x in pattern.finditer(seq.sequence):
                pos.append((motif, "+", x.start(), x.end()))
            for x in pattern.finditer(rseq):
                pos.append(
                    (motif, "-", lsequence - x.end(), lsequence - x.start()))

        positions.append((seq.title, pos))

    return positions


def getCounts(matches):
    '''count numbers of motifs.'''
    counts = collections.defaultdict(int)
    for seq, pos in matches:
        for motif, strand, start, end in pos:
            counts[motif] += 1
    return counts


def getOccurances(matches):
    '''count numbers of motifs, but only once per sequence'''
    counts = collections.defaultdict(int)
    for seq, pos in matches:
        motifs = set([x[0] for x in pos])
        for m in motifs:
            counts[m] += 1
    return counts

iupacdict = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'M': 'AC',
    'R': 'AG',
    'W': 'AT',
    'S': 'CG',
    'Y': 'CT',
    'K': 'GT',
    'V': 'ACG',
    'H': 'ACT',
    'D': 'AGT',
    'B': 'CGT',
    'X': 'ACGT',
    'N': 'ACGT'}

regexdict = dict(((x[1], x[0]) for x in iupacdict.items()))


def iupac2regex(pattern):
    '''convert iupac to regex pattern'''


def regex2iupac(pattern):
    '''convert regex to iupac pattern'''

    def _split(p):
        g = p.split("[")
        for c in g[0]:
            yield c
        for x in g[1:]:
            a, b = x.split("]")
            yield "".join(sorted(a))
            for c in b:
                yield c

    a = []
    for x in _split(pattern):
        a.append(regexdict[x])
    return "".join(a)
