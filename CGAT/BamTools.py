'''
BamTools - utilities for working with BAM files
===============================================

'''

import pysam


def isPaired(bamfile, alignments=1000):
    '''check if a *bamfile* contains paired end data

    The method reads at most the first *alignments* and returns
    True if any of the alignments are paired.
    '''

    samfile = pysam.Samfile(bamfile)
    n = 0
    for read in samfile:
        if read.is_paired:
            break
        n += 1
        if n == alignments:
            break

    samfile.close()

    return n != alignments
