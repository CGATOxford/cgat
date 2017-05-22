"""
Bioprospector.py - Parser for Bioprospector output
==================================================

:Tags: Python
"""

import sys
import re
import collections

import weblogolib
import corebio.seq

Match = collections.namedtuple("Match",
                               '''id, site, strand, start, end,
                                  width1, width2, sequence''')

Motif = collections.namedtuple("Motif",
                               '''id, pattern, width1, width2,
                                min_gap, max_gap,
                                score, sites, matches''')


def parse(infile):
    '''parse Bioprospector output.'''

    motifs = []

    def iterate_blocks(infile):
        lines = []
        keep = False
        for line in infile:
            if line.startswith("Motif #"):
                if lines:
                    yield lines
                keep = True
                lines = []

            if line.startswith("*****"):
                continue

            if keep:
                lines.append(line)
        if lines:
            yield(lines)

    for lines in iterate_blocks(infile):
        motif_id, pattern = re.match(
            "Motif \#(\d+): \((.+)\)", lines[0]).groups()
        widths, gaps, score, sites = re.match(
            "Width \((.*)\); Gap \[(.*)\]; MotifScore (.*); Sites (.*)", lines[1]).groups()
        width1, width2 = list(map(int, widths.split(",")))
        min_gap, max_gap = list(map(int, gaps.split(",")))

        sites = int(sites)
        start = end = 3
        while not lines[end].startswith(">"):
            end += 1

        block = lines[start:end]

        def convert(coords, width):
            if coords[0] == "f":
                strand, start = "+", int(coords[1]) - 1
            else:
                strand, start = "-", int(coords[1]) - width
            return strand, start

        sequences = []
        for x in range(end, len(lines), 2):
            if lines[x].strip() == "":
                break
            assert lines[x].startswith(">"), lines[x]
            id, lseq, site_id, coords = re.match(
                ">(\S+).*len\s(\d+).*site \#(\d+)\s*(.*)", lines[x][:-1]).groups()
            lseq, site_id = list(map(int, (lseq, site_id)))
            coords = re.split("\s", coords)
            if len(coords) == 2:
                strand1, start1 = convert(coords[:2], width1)
                strand = strand1
                start, end = start, start + width1
            elif len(coords) == 4:
                strand1, start1 = convert(coords[:2], width1)
                strand2, start2 = convert(coords[2:], width2)
                start, end = min(start1, start2), max(
                    start1 + width1, start2 + width2)
                strand = strand1 + strand2

            seq = lines[x + 1][:-1]
            sequences.append(Match._make((id, site_id,
                                          strand,
                                          start, end,
                                          width1, width2,
                                          seq)))

        motif = Motif._make((motif_id, pattern,
                             width1, width2,
                             min_gap, max_gap,
                             score,
                             sites, sequences))
        motifs.append(motif)

    return motifs


def build_logo(sequences, outfilename):

    seqs = corebio.seq.SeqList(alphabet=corebio.seq.dna_alphabet)
    for sequence in sequences:
        seqs.append(corebio.seq.dna(re.sub("\s", "-", sequence)))

    data = weblogolib.LogoData.from_seqs(seqs)
    options = weblogolib.LogoOptions()
    options.color_scheme = weblogolib.classic

    options.title = 'A Logo Title'
    format = weblogolib.LogoFormat(data, options)
    fout = open(outfilename, 'w')
    weblogolib.png_formatter(data, format, fout)

if __name__ == "__main__":

    results = parse(sys.stdin)

    outname = "test%i.eps"

    from . import MAST
    for x, motifs in enumerate(results):
        build_logo([y.sequence for y in motifs.matches], outname % x)
        print(motifs.score, motifs.sites)

        outfile = open("test%i.motif" % x, "w")
        MAST.sequences2motif(outfile,
                             [y.sequence for y in motifs.matches])
        outfile.close()

        outfile = open("test%i.tomtom" % x, "w")
        MAST.sequences2motif(outfile,
                             [y.sequence for y in motifs.matches],
                             format="TOMTOM")
        outfile.close()
