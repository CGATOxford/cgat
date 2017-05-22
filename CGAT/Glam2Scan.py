"""
Glam2Scan.py - Parser for MAST files
====================================

:Tags: Python

As of biopython 1.5.6, the MAST parser is broken.

"""

import re


class Glam2Scan:
    # number of matches

    def __init__(self):
        self.motifs = []
        self.matches = []


class Motif:
    pass


class Match:

    """a Glam2Scan entry.
    """
    header = "\t".join(("id", "start", "end", "strand", "score"))

    def __init__(self):
        pass

    def __str__(self):
        return "\t".join(map(str,
                             (self.id,
                              self.start,
                              self.end,
                              self.strand,
                              self.sequence,
                              self.score)))


def parse(infile):
    '''parse Glam2Scan output.'''

    x = -1
    for line in infile:
        if line.startswith("GLAM2SCAN"):
            x = 4
        x -= 1
        if x == 0:
            break

    if not line:
        raise ValueError("could not find tag 'GLAM2SCAN'")

    gl = Glam2Scan()
    keep = False
    for line in infile:
        if line.strip() == "":
            continue
        if line.startswith(" "):
            continue
        m = Match()
        d = re.split("\s+", line[:-1])
        if len(d) != 6:
            break
        m.id, m.start, m.sequence, m.end, m.strand, m.score =\
            d[0], int(d[1]) - 1, d[2], int(d[3]), d[4], float(d[5])
        gl.matches.append(m)

    return gl
