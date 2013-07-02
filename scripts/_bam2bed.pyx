from pysam.csamtools cimport *

import collections, array, struct, itertools
import CGAT.Experiment as E

def merge_pairs( Samfile input_samfile,
                 outfile,
                 min_insert_size = 0,
                 max_insert_size = 400 ):
    '''merge paired ended data.

    For speed reasons, the aligned region is only approximated using
    the read length. Indels within reads are ignored. Thus bed coordinates
    might be a few residues off.

    The strand is always set to '+'.

    Pairs with a maximum insert size larger than *max_insert_size* are removed.
    '''

    cdef int ninput = 0
    cdef int nremoved_unmapped = 0
    cdef int nremoved_insert = 0
    cdef int nremoved_contig = 0
    cdef int nremoved_unpaired = 0
    cdef int noutput = 0
    cdef int flag

    cdef AlignedRead read
    cdef int c_max_insert_size = max_insert_size
    cdef int c_min_insert_size = min_insert_size
    cdef int start, end

    for read in input_samfile:

        ninput += 1

        flag = read._delegate.core.flag 
        # remove unmapped reads
        if flag & 4:
            nremoved_unmapped += 1
            continue

        # remove unpaired
        if not flag & 2:
            nremoved_unpaired += 1
            continue
            
        # this is second pair of read - skip
        if flag & 128:
            continue

        if read.tid != read.mrnm:
            nremoved_contig += 1
            continue

        if (c_max_insert_size and read.isize > c_max_insert_size) or (c_min_insert_size and read.isize < c_min_insert_size) :
            nremoved_insert += 1
            continue

        if read.pos < read.mpos:
            start = read.pos
            end = read.mpos + read.rlen
        else:
            start = read.mpos
            end = read.pos + read.rlen

        # count output pair as two so that it squares with ninput
        noutput += 2

        outfile.write( "%s\t%i\t%i\t%s\t%i\t%s\n" % 
                       (input_samfile.getrname( read.tid ),
                        start, end,
                        read.qname,
                        read.mapq,
                        "+"
                        ) )

    c = E.Counter()
    c.input = ninput
    c.removed_insert = nremoved_insert
    c.removed_contig = nremoved_contig
    c.removed_unmapped = nremoved_unmapped
    c.removed_unpaired = nremoved_unpaired
    c.output = noutput

    return c
