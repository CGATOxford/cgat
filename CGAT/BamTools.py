"""
BamTools - Utilities for working with BAM files
===============================================

This module brings together convenience function for working
with :term:`bam` formatted files.

"""

import numpy
import pysam


def isPaired(bamfile, alignments=1000):
    '''check if a `bamfile` contains paired end data

    The method reads at most the first *alignments* and returns
    True if any of the alignments are paired.
    '''

    samfile = pysam.AlignmentFile(bamfile)
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
                                   alignments=10000,
                                   n=10,
                                   method="picard",
                                   similarity_threshold=1.0,
                                   max_chunks=1000):
    '''estimate insert size from a subset of alignments in a bam file.

    Several methods are implemented.

    picard
        The method works analogous to picard by restricting the estimates
        to a core distribution. The core distribution is defined as all
        values that lie within n-times the median absolute deviation of
        the full data set.
    convergence
        The method works similar to ``picard``, but continues reading
        `alignments` until the mean and standard deviation stabilize.
        The values returned are the median mean and median standard
        deviation encountered.

    The method `convergence` is suited to RNA-seq data, as insert sizes
    fluctuate siginificantly depending on the current region
    being looked at.

    Only mapped and proper pairs are considered in the computation.

    Returns
    -------
    mean : float
       Mean of insert sizes.
    stddev : float
       Standard deviation of insert sizes.
    npairs : int
       Number of read pairs used for the estimation
    method : string
       Estimation method
    similarity_threshold : float
       Similarity threshold to apply.
    max_chunks : int
       Maximum number of chunks of size `alignments` to be used
       in the convergence method.

    '''

    assert isPaired(bamfile), \
        'can only estimate insert size from' \
        'paired bam files'

    samfile = pysam.AlignmentFile(bamfile)

    def get_core_distribution(inserts, n):
        # compute median absolute deviation
        raw_median = numpy.median(inserts)
        raw_median_dev = numpy.median(numpy.absolute(inserts - raw_median))

        # set thresholds
        threshold_min = max(0, raw_median - n * raw_median_dev)
        threshold_max = raw_median + n * raw_median_dev

        # define core distribution
        return inserts[numpy.logical_and(inserts >= threshold_min,
                                         inserts <= threshold_max)]

    if method == "picard":

        # only get first read in pair to avoid double counting
        inserts = numpy.array(
            [read.template_length for read in samfile.head(n=alignments)
             if read.is_proper_pair
             and not read.is_unmapped
             and not read.mate_is_unmapped
             and not read.is_read1
             and not read.is_duplicate
             and read.template_length > 0])
        core = get_core_distribution(inserts, n)

        return numpy.mean(core), numpy.std(core), len(inserts)

    elif method == "convergence":

        means, stds, counts = [], [], []
        last_mean = 0
        iteration = 0
        while iteration < max_chunks:

            inserts = numpy.array(
                [read.template_length for read in samfile.head(
                    n=alignments,
                    multiple_iterators=False)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)
            means.append(numpy.mean(core))
            stds.append(numpy.std(core))
            counts.append(len(inserts))
            mean_core = get_core_distribution(numpy.array(means), 2)
            mm = numpy.mean(mean_core)
            if abs(mm - last_mean) < similarity_threshold:
                break
            last_mean = mm

        return numpy.median(means), numpy.median(stds), sum(counts)
    else:
        raise ValueError("unknown method '%s'" % method)


def estimateTagSize(bamfile,
                    alignments=10,
                    multiple="error"):
    '''estimate tag/read size from first alignments in file.

    Arguments
    ---------
    bamfile : string
       Filename of :term:`bam` formatted file
    alignments : int
       Number of alignments to inspect
    multiple : string
       How to deal if there are multiple tag sizes present.
       ``error`` will raise a warning, ``mean`` will return the
       mean of the read lengths found. ``uniq`` will return a
       unique list of read sizes found. ``all`` will return all
       read sizes encountered.

    Returns
    -------
    size : int
       The read size (actual, mean or list of read sizes)

    Raises
    ------
    ValueError
       If there are multiple tag sizes present and `multiple` is set to
       `error`.

    '''
    samfile = pysam.AlignmentFile(bamfile)
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
        elif multiple == "uniq":
            mi = list(sorted(set(sizes)))
        elif multiple == "all":
            return sizes

    return mi


def getNumberOfAlignments(bamfile):
    '''return number of alignments in bamfile.
    '''
    samfile = pysam.AlignmentFile(bamfile)
    return samfile.mapped


def getNumReads(bamfile):
    '''count number of reads in bam file.

    This methods works through pysam.idxstats.

    Arguments
    ---------
    bamfile : string
        Filename of :term:`bam` formatted file. The file needs
        to be indexed.
    Returns
    -------
    nreads : int
        Number of reads
    '''

    lines = pysam.idxstats(bamfile).splitlines()

    try:
        nreads = sum(
            map(int, [x.split("\t")[2]
                      for x in lines if not x.startswith("#")]))

    except IndexError as msg:
        raise IndexError(
            "can't get number of reads from bamfile, msg=%s, data=%s" %
            (msg, lines))
    return nreads


def isStripped(bamfile):
    '''
    Check if the sequence is stripped in a bam file.
    '''
    bam = pysam.AlignmentFile(bamfile)
    read = bam.next()
    if read.seq is None:
        return True
    else:
        return False
