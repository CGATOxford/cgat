'''
MACS.py - Parser for MACS output
================================

:Tags: Python

The :mod:`Pipeline` module contains various utility functions
for parsing MACS output.

API
----

'''

import collections
import math
import csv

import CGAT.CSV as CSV


MacsPeak = collections.namedtuple(
    "MacsPeak", "contig start end length summit tags pvalue fold fdr")


def iterateMacsPeaks(infile):
    '''iterate over peaks.xls file and return parsed data.
    The fdr is converted from percent to values between 0 and 1.
    '''

    for line in infile:
        if line.startswith("#"):
            continue
        if line.startswith("chr\tstart"):
            continue
        # skip empty lines
        if line.startswith("\n"):
            continue

        data = line[:-1].split("\t")

        if len(data) == 9:
            # convert % fdr
            data[8] = float(data[8]) / 100.0
        elif len(data) == 8:
            # if no fdr given, set to 0
            # data.append( 0.0 )
            # Steve - I don't understand this so I'm commenting it out and
            # raising an error
            raise ValueError("FDR value not set line %s" % line)
        else:
            raise ValueError("could not parse line %s" % line)

        # these are 1-based coordinates
        # macs can have negative start coordinates
        # start
        data[1] = max(int(data[1]) - 1, 0)
        # end
        data[2] = int(data[2])
        # length
        data[3] = int(data[3])
        # summit
        data[4] = int(data[4])
        # ntags
        data[5] = int(data[5])
        # -10log10(pvalue)
        data[6] = float(data[6])
        # fold
        data[7] = float(data[7])

        yield MacsPeak._make(data)

Macs2Peak = collections.namedtuple(
    "Macs2Peak",
    "contig start end length pileup "
    "pvalue fold qvalue "
    "name")


def iterateMacs2Peaks(infile):
    '''iterate over peaks.xls file and return parsed data.

    pvalues and fdr are converted to values between 0 and 1
    from their -log10 values.
    '''

    for row in csv.DictReader(CSV.CommentStripper(infile),
                              dialect='excel-tab'):
        # these are 1-based coordinates
        # macs can have negative start coordinates
        # start
        try:
            yield Macs2Peak._make(
                (row['chr'],
                 max(int(row['start']) - 1, 0),
                 int(row['end']),
                 int(row['length']),
                 float(row['pileup']),
                 math.pow(10, -float(row['-log10(pvalue)'])),
                 float(row['fold_enrichment']),
                 math.pow(10, -float(row['-log10(qvalue)'])),
                 row['name']))
        except KeyError as msg:
            raise KeyError("%s: %s" % (msg, row))
