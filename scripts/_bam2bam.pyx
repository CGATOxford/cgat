from pysam.csamtools cimport *

import collections, array, struct, itertools
import CGAT.Experiment as E

cdef class SetNH:

    cdef iter
    cdef stack

    def __cinit__(self, iter ):
        self.iter = itertools.groupby( iter, lambda x: x.qname )
        self.stack = []

    def __iter__(self):
        return self

    def __next__(self): 
        """python version of next().
        """
        
        while 1:
            if self.stack:
                return self.stack.pop(0)
            else:
                key, x = self.iter.next()
                self.stack = list(x)
                nh = len(self.stack)
                for read in self.stack:
                    if not read.is_unmapped:
                        t = dict( read.tags )
                        t['NH'] = nh
                        read.tags = list(t.iteritems())

def filter_bam( Samfile input_samfile,
                Samfile output_samfile,
                Samfile reference_samfile,
                remove_nonunique = False,
                remove_contigs = None,
                remove_unmapped = False,
                remove_mismatches = False,
                colour_mismatches = False ):

    '''
    To conserve memory, the tid and NM flag from *transcripts_samfile*
    are packed into memory. As a consequence, this method requires 
    max(NM) < 256 (2^8) and max(tid) < 16777216 (2^24)

    If *remove_nonunique* is set, only uniquely matching reads will be output.

    If *remove_contigs* is given, contigs that are in remove_contigs will
    be removed. Note that this will only remove the alignment, but not
    all alignments of a particuluar read and the NH flag will *NOT* be
    updated.

    If *remove_unmapped* is given, unmapped reads will be removed.

    If *remove_mismatches* is set, only reads with number of mismatches
    better than in reference_samfile will be kept.

    If *colour_mismatches* is set, the ``CM`` tag will be used
    to count differences. By default, the ``NM`` tag is used.
    The tag that is used needs to present in both *input_samfile*
    and *reference_samfile*.
    '''

    cdef int ninput = 0
    cdef int nunmapped = 0
    cdef int nmismatches = 0
    cdef int noutput = 0
    cdef int nnonunique = 0
    cdef int nremoved_contigs = 0

    cdef bint c_remove_nonunique = remove_nonunique
    cdef bint c_remove_mismatches = remove_mismatches
    cdef bint c_remove_unmapped = remove_unmapped

    cdef int * remove_contig_tids 
    cdef int nremove_contig_tids = 0
    cdef AlignedRead read
    cdef AlignedRead match

    E.info( "building index" )
    # build index
    # this method will start indexing from the current file position
    # if you decide
    cdef int ret = 1#
    cdef int x
    cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t) )
    cdef uint64_t pos
    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int tid
    cdef int transript_tid
    cdef long val
    cdef bint skip = 0

    cdef int32_t read_mismatches = 0

    # decide which tag to use
    cdef char * nm_tag = 'NM' 
    cdef char * cm_tag = 'CM'
    cdef char * tag

    if colour_mismatches:
        tag = cm_tag
    else:
        tag = nm_tag

    if c_remove_mismatches:
        if not reference_samfile:
            raise ValueError("require another bam file for mismatch filtering" )

        # L = 1 byte (unsigned char)
        def _gen(): return array.array('B') 
        index = collections.defaultdict( _gen )

        while ret > 0:
            ret = samread( reference_samfile.samfile, b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4: continue
                qname = bam1_qname( b )
                tid = b.core.tid
                v = bam_aux_get(b, tag)
                nm = <int32_t>bam_aux2i(v)
                index[qname].append( nm )

    # setup list of contigs to remove:
    if remove_contigs:
        nremove_contig_tids = len(remove_contigs)
        remove_contig_tids = <int*>malloc( sizeof(int) * nremove_contig_tids )
        for x, rname in enumerate( remove_contigs):
            remove_contig_tids[x] = input_samfile.gettid( rname )

    bam_destroy1( b )
    
    E.info( "built index for %i reads" % len(index))
    
    E.info( "starting filtering" )

    for read in input_samfile:

        ninput += 1
        # if ninput > 10000: break

        # remove unmapped reads
        if read._delegate.core.flag & 4:
            if c_remove_unmapped:
                nunmapped += 1
            else:
                output_samfile.write( read )
            continue

        # remove non-unique alignments
        if c_remove_nonunique:
            v = bam_aux_get(read._delegate, 'NH')
            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh > 1: 
                    nnonunique += 1
                    continue

        # remove reads matched to certain contigs
        if nremove_contig_tids:
            skip = 0
            for x from 0 <= x < nremove_contig_tids:
                if remove_contig_tids[x] == read.tid:
                    skip = 1
                    break
            if skip: 
                nremoved_contigs += 1
                continue

        # remove reads in other bam file if their are better matching
        if c_remove_mismatches:
            if read.qname in index:
                # can compress index before, depends on
                # how many reads you expect test to filter
                nm = min(index[read.qname])
                read_mismatches = read.opt(tag)
                if nm > read_mismatches: 
                    nmismatches += 1
                    continue
                
        noutput += 1
        output_samfile.write( read )

    c = E.Counter()
    c.input = ninput
    c.removed_nonunique = nnonunique
    c.removed_contigs = nremoved_contigs
    c.output = noutput
    c.removed_unmapped = nunmapped
    c.removed_mismatches = nmismatches

    if nremove_contig_tids:
        free( remove_contig_tids )

    E.info( "filtering finished" )

    return c
