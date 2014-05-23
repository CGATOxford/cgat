'''
BamTools - utilities for working with BAM files
===============================================

'''

import numpy
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


def estimateInsertSizeDistribution(bamfile,
                                   alignments=1000):
    '''estimate insert size from first alignments in bam file.

    returns mean and stddev of insert sizes.
    '''

    assert isPaired(bamfile), \
        'can only estimate insert size from' \
        'paired bam files'

    samfile = pysam.Samfile(bamfile)
    # only get positive to avoid double counting
    inserts = numpy.array(
        [read.tlen for read in samfile.head(alignments)
         if read.is_proper_pair and read.tlen > 0])
    return numpy.mean(inserts), numpy.std(inserts)


def estimateTagSize(bamfile,
                    alignments=10,
                    multiple="error"):
    '''estimate tag size from first alignments in file.'''
    samfile = pysam.Samfile(bamfile)
    sizes = [read.rlen for read in samfile.head(alignments)]
    mi, ma = min(sizes), max(sizes)

    if mi == 0 and ma == 0:
        sizes = [read.inferred_length for read in samfile.head(alignments)]
        # remove 0 sizes (unaligned reads?)
        sizes = [x for x in sizes if x > 0]
        mi, ma = min(sizes), max(sizes)

    if mi != ma:
        if multiple == "error":
            raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
        elif multiple == "mean":
            mi = int(sum(sizes) / len(sizes))

    return mi


def getNumberOfAlignments(bamfile):
    '''return number of alignments in bamfile.
    '''
    samfile = pysam.Samfile(bamfile)
    return samfile.mapped

