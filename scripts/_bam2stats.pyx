#cimport csamtools

from pysam.csamtools cimport *

import collections, array, struct
import CGAT.Experiment as E

FLAGS = {
    1: 'paired',
    2: 'proper_pair',
    4: 'unmapped',
    8: 'mate_unmapped',
    16: 'reverse',
    32: 'mate_reverse',
    64: 'read1',
    128: 'read2',
    256: 'secondary',
    512: 'qc_fail',
    1024: 'duplicate'  }

def count( Samfile samfile,
           remove_rna,
           rna):
    '''

    '''
    cdef bint _remove_rna = remove_rna

    cdef AlignedRead read

    # counters
    cdef int ninput = 0
    cdef int nduplicates = 0
    cdef int nrna = 0
    cdef int nfiltered = 0

    cdef int max_hi = 0
    # count nh, nm tags
    nh_filtered, nm_filtered = collections.defaultdict( int ), collections.defaultdict( int )
    nh_all, nm_all = collections.defaultdict( int ), collections.defaultdict( int )
    mapq_filtered, mapq_all = collections.defaultdict( int ), collections.defaultdict( int )

    cdef int * flags_counts = <int*>calloc( len(FLAGS), sizeof(int) )

    # helper variables
    cdef int last_tid = -1
    cdef int last_pos = 0

    duplicates = collections.defaultdict( int )
    counts = collections.defaultdict( int )
    cdef int count = 0

    cdef int tid = -1
    cdef int flag

    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int32_t hi
    cdef int x
    cdef int lflags = len(FLAGS)
    cdef int f

    for read in samfile:

        flag = read._delegate.core.flag

        # is paired and read2
        #if flag & 1 and flag & 128:
            # ignore second reads to avoid double counting pairs
            # this needs to be done properly by counting in
            # pairs
        #    continue

        ninput += 1

        f = 1
        for x from 0 <= x < lflags:
            if flag & f: flags_counts[x] += 1
            f = f << 1

        # get maximum NI field
        v = bam_aux_get(read._delegate, 'HI')
        if v != NULL:
            hi = <int32_t>bam_aux2i(v)
            if hi > max_hi: max_hi = hi

        v = bam_aux_get(read._delegate, 'NH')
        if v != NULL:
            nh = <int32_t>bam_aux2i(v)
            nh_all[nh] += 1
        else:
            nh = -1

        v = bam_aux_get(read._delegate, 'NM')
        if v != NULL:
            nm = <int32_t>bam_aux2i(v)
            nm_all[nm] += 1
        else:
            nm = -1

        mapq_all[read.mapq] += 1

        # skip unmapped reads
        if read._delegate.core.flag & 4: continue

        if read.tid != last_tid:
            contig = samfile.getrname( read.rname )

        # note: does not take into account gaps within reads
        # or partial overlap.
        if rna and rna.contains( contig, read.pos, read.pos + read.qlen ):
            nrna += 1
            if _remove_rna: continue
        
        nfiltered += 1

        if nh >= 0: nh_filtered[nh] += 1
        if nm >= 0: nm_filtered[nm] += 1
        mapq_filtered[read.mapq] += 1

        # duplicate analysis - simply count per start position
        # ignoring sequence and strand
        if read.tid == last_tid and read.pos == last_pos:
            count += 1
            nduplicates += 1
            continue

        if count > 1:
            counts[count] += 1

        count = 1
        last_tid, last_pos = read.tid, read.pos

    c = E.Counter()
    
    c.input = ninput
    c.filtered = nfiltered
    c.rna = nrna
    c.duplicates = nduplicates

    # convert flags to labels
    t = {}
    f = 1
    for x from 0 <= x < lflags:
        t[FLAGS[f]] = flags_counts[x]
        f = f << 1

    return c, t, nh_filtered, nh_all, nm_filtered, nm_all, mapq_filtered, mapq_all, max_hi
