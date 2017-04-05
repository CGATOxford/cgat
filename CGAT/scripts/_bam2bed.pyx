from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport *

from libc.stdlib cimport abs

import collections, array, struct, itertools
import CGAT.Experiment as E

def merge_pairs(AlignmentFile input_samfile,
                outfile,
                min_insert_size = 0,
                max_insert_size = 400,
                bed_format = None ):
    '''merge paired ended data.

    For speed reasons, the aligned region is only approximated using
    the read length. Indels within reads are ignored. Thus bed coordinates
    might be a few residues off.

    The strand is always set to '+'.

    Pairs with a maximum insert size larger than *max_insert_size* are removed.

    If `bed_format` is a number, only the first x columns will be output.
    '''

    cdef int ninput = 0
    cdef int nremoved_unmapped = 0
    cdef int nremoved_insert = 0
    cdef int nremoved_contig = 0
    cdef int nremoved_unpaired = 0
    cdef int nremoved_coordinate = 0
    cdef int nremoved_take_only_second = 0
    cdef int noutput = 0
    cdef int flag
    cdef int isize

    cdef AlignedSegment read
    cdef int c_max_insert_size = max_insert_size
    cdef int c_min_insert_size = min_insert_size
    cdef int start, end, xstart, xend
    cdef int take_columns = 6

    # point to array of contig lengths
    cdef uint32_t *contig_sizes = input_samfile.header.target_len
    
    if bed_format != None:
        if bed_format < 3 or bed_format > 6: 
            raise ValueError("a bed file must have at least 3 and at most 6 columns")
        take_columns = bed_format

    for read in input_samfile:
        ninput += 1

        flag = read._delegate.core.flag 
        # remove unmapped reads
        if flag & 4:
            nremoved_unmapped += 1
            continue

        if read.pos < read.mpos:
            # lower coordinate than mate, ignore
            nremoved_coordinate += 1
            continue
        elif read.pos == read.mpos and flag & 64:
            # disambiguate, ignore first in pair
            nremoved_take_only_second += 1
            continue
        else:
            # taking the downstream pair allows to incl
            xstart = read.next_reference_start
            xend = read.reference_end
            if xstart < xend:
                start, end = xstart, xend
            else:
                start, end = xend, xstart

        # remove unpaired
        if not flag & 2:
            nremoved_unpaired += 1
            continue
            
        if read.tid != read.mrnm:
            nremoved_contig += 1
            continue

        # isize can be negative - depending on the pair orientation
        isize = abs(read.isize)
        if (c_max_insert_size and isize > c_max_insert_size) or \
           (c_min_insert_size and isize < c_min_insert_size):
            nremoved_insert += 1
            continue

        # truncate at contig end - overhanging reads might cause problems with chrM
        if end > contig_sizes[read.mrnm]:
            end = contig_sizes[read.mrnm]

        # count output pair as two so that it squares with ninput
        noutput += 2

        if take_columns == 3:
            outfile.write("%s\t%i\t%i\n" %
                          (input_samfile.getrname(read.tid),
                           start, end))
        elif take_columns == 4:
            outfile.write("%s\t%i\t%i\t%s\n" % 
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname))
        elif take_columns == 5:
            outfile.write("%s\t%i\t%i\t%s\t%i\n" % 
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname,
                           read.mapq,
                       ))
        else:
            # As we output the downstream read, reverse orientation
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            outfile.write("%s\t%i\t%i\t%s\t%i\t%s\n" % 
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname,
                           read.mapq,
                           strand
                          ))

    c = E.Counter()
    c.input = ninput
    c.removed_insert = nremoved_insert
    c.removed_contig = nremoved_contig
    c.removed_unmapped = nremoved_unmapped
    c.removed_unpaired = nremoved_unpaired
    c.removed_take_only_second = nremoved_take_only_second
    c.removed_coordinate = nremoved_coordinate
    c.output = noutput

    return c
