'''
WrapperZinba.py - utility functions for zinba output
====================================================

:Tags: Python

Purpose
-------

Usage
-----

Documentation
-------------

Code
----

'''

import collections

ZinbaPeak = collections.namedtuple(
    "ZinbaPeak",
    "contig unrefined_start unrefined_end strand posterior "
    "summit height refined_start refined_end median fdr")


def iteratePeaks(infile):
    '''iterate of zinba peaks in infile.'''

    for line in infile:

        if line.startswith("#"):
            continue
        if line.startswith("PEAKID\tChrom"):
            continue
        # skip empty lines
        if line.startswith("\n"):
            continue

        data = line[:-1].split("\t")

        if len(data) != 12:
            raise ValueError("could not parse line %s" % line)

        # I assume these are 1-based coordinates
        data[2] = max(int(data[2]) - 1, 0)
        # end
        data[3] = int(data[3])
        # posterior
        data[5] = float(data[5])
        # summit
        data[6] = max(int(data[6]) - 1, 0)
        # height
        data[7] = int(data[7])
        # refined_start
        data[8] = max(int(data[8]) - 1, 0)
        # end
        data[9] = int(data[9])
        # median
        data[10] = int(data[10])
        # qvalue
        data[11] = float(data[11])

        yield ZinbaPeak._make(data[1:])
