"""
Glam2.py - Parser for MAST files.
=================================

:Tags: Python

As of biopython 1.5.6, the MAST parser is broken.

"""

import re


class Match:
    # number of matches

    def __init__(self):
        self.score = 0
        self.columns = 0
        self.sequences = 0


class Glam2:
    # number of matches

    def __init__(self):
        self.motifs = []
        self.version = ""
        self.sequences = 0


def parse(infile):
    '''parse Glam2 output.'''

    g = Glam2()
    for line in infile:
        if line.startswith("Version"):
            g.version = re.match("Version (\d+)", line).groups()[0]
        elif line.startswith("Sequences"):
            g.sequences = re.match("Sequences: (\d+)", line).groups()[0]
        elif line.startswith("Score"):
            m = Match()
            d = re.match(
                "Score:\s*(\S+)\s*Columns:\s*(\d+)\s*Sequences:\s*(\d+)", line).groups()
            m.score, m.columns, m.sequences = float(d[0]), int(d[1]), int(d[2])
            g.motifs.append(m)

    return g
