#cimport csamtools

from pysam.csamtools cimport *

import collections, array, struct, sys
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

cdef struct CountsType:
    int is_unmapped
    int mate_is_unmapped
    int is_paired
    int is_read1
    int is_read2
    int is_proper_pair
    int is_qcfail
    int is_duplicate

def count( Samfile samfile,
           remove_rna,
           rna,
           filename_fastq = None,
           outfile_details = None):
    '''

    '''
    cdef bint _remove_rna = remove_rna

    cdef AlignedRead read

    # counters
    cdef int ninput = 0
    cdef int nduplicates = 0
    # number of reads overlapping RNA (if rna != None)
    cdef int nrna = 0
    # number of reads not overlap RNA (if rna != None)
    cdef int n_norna = 0
    # number of reads present after filtering
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

    # detailed counting
    cdef FastqProxy fq
    cdef int64_t index, fastq_nreads
    cdef CountsType * fastq_counts
    cdef CountsType * fastq_count
    cdef char * read_name
    cdef int count_fastq = filename_fastq != None

    if count_fastq:
        E.info( "reading fastq file" )
        # Using a python dictionary here is due
        # for a large amount of memory usage.
        # Alternatives to dictionary
        # 1. POSIX hash tables (hsearch,...) or trees: very slow
        # 2. custom hash implementation: worth the effort?
        # 3. Sorted list and binary search: too slow for many lookups
        reads = {}
        fastqfile = Fastqfile( filename_fastq )
        fastq_nreads = 0
        for fq in fastqfile:
            reads[fq.name] = fastq_nreads
            fastq_nreads += 1
        
        E.info( "read %i reads" % fastq_nreads )

        E.info( "allocating %i bytes" % (fastq_nreads * sizeof(CountsType) ) )

        fastq_counts = <CountsType *>calloc( fastq_nreads, sizeof( CountsType ) )
        if fastq_counts == NULL:
            raise ValueError( "could not allocated memory: %i bytes" % (fastq_nreads * sizeof(CountsType) ))
 
    for read in samfile:

        if count_fastq:
            try:
                fastq_count = &fastq_counts[ reads[read.qname] ]
            except KeyError:
                continue
        
            if read.is_unmapped: fastq_count.is_unmapped += 1
            if read.mate_is_unmapped: fastq_count.mate_is_unmapped += 1
            if read.is_paired: fastq_count.is_paired += 1
            if read.is_read1: fastq_count.is_read1 += 1
            if read.is_read2: fastq_count.is_read2 += 1
            if read.is_proper_pair: fastq_count.is_proper_pair += 1
            if read.is_qcfail: fastq_count.is_qcfail += 1
            if read.is_duplicate: fastq_count.is_duplicate += 1
 
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
        if rna:
            if rna.contains( contig, read.pos, read.pos + read.alen ):
                nrna += 1
                if _remove_rna: continue
            else:
                n_norna += 1

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

    E.info( "finished computing counts" )

    counter = E.Counter()
    
    counter.input = ninput
    counter.filtered = nfiltered
    counter.rna = nrna
    counter.no_rna = n_norna
    counter.duplicates = nduplicates

    # convert flags to labels
    t = {}
    f = 1
    for x from 0 <= x < lflags:
        t[FLAGS[f]] = flags_counts[x]
        f = f << 1

    # count based on fastq data
    cdef int total_pairs = 0
    cdef int total_pair_is_unmapped = 0
    cdef int total_pair_is_proper_uniq = 0
    cdef int total_pair_is_proper_mmap = 0
    cdef int total_pair_is_proper_duplicate = 0
    cdef int total_pair_not_proper_uniq = 0
    cdef int total_pair_is_incomplete = 0
    cdef int total_pair_is_other = 0

    if count_fastq:
        E.info( "aggregating counts" )
        for qname, index in reads.items():
            fastq_count = &fastq_counts[index]

            # paired read counting
            if fastq_count.is_paired : 
                total_pairs += 1

                if fastq_count.is_unmapped == 2:
                    # an unmapped read pair
                    total_pair_is_unmapped += 1
                elif fastq_count.is_proper_pair == 2: 
                    # a unique proper pair
                    total_pair_is_proper_uniq +=1
                    # a duplicate unique proper pair
                    if fastq_count.is_duplicate == 2:
                        total_pair_is_proper_duplicate +=1
                elif fastq_count.is_proper_pair > 2:
                    # proper pairs that map to mulitple locations
                    total_pair_is_proper_mmap +=1
                elif fastq_count.is_read1 == 1 \
                        and fastq_count.is_read2 == 1 \
                        and fastq_count.is_proper_pair == 0:
                    # pair which map each read once, but not uniquely
                    total_pair_not_proper_uniq += 1
                elif (fastq_count.is_read1 == 1 and fastq_count.is_read2 == 0) \
                        or (fastq_count.is_read1 == 0 and fastq_count.is_read2 == 1):
                    # an incomplete pair - one read of a pair matches uniquel, but not the other
                    total_pair_is_incomplete += 1
                else:
                    total_pair_is_other += 1

        counter.total_pairs = total_pairs
        counter.total_pair_is_unmapped = total_pair_is_unmapped
        counter.total_pair_is_proper_uniq = total_pair_is_proper_uniq
        counter.total_pair_is_incomplete = total_pair_is_incomplete
        counter.total_pair_is_proper_duplicate = total_pair_is_proper_duplicate
        counter.total_pair_is_proper_mmap = total_pair_is_proper_mmap
        counter.total_pair_not_proper_uniq = total_pair_not_proper_uniq
        counter.total_pair_is_other = total_pair_is_other

        if outfile_details:
            if outfile_details != sys.stdout:
                # later: get access FILE * object
                outfile_details.write( "read\tis_unmapped\tmate_is_unmapped\tis_paired\tis_read1\tis_read2\tis_proper_pair\tis_qcfail\tis_duplicate\n" )
                for qname, index in reads.items():
                    fastq_count = &fastq_counts[index]

                    outfile_details.write( "%s\t%s\n" % (qname, "\t".join( \
                                map(str,
                                    (fastq_count.is_unmapped,
                                     fastq_count.mate_is_unmapped,
                                     fastq_count.is_paired,
                                     fastq_count.is_read1,
                                     fastq_count.is_read2,
                                     fastq_count.is_proper_pair,
                                     fastq_count.is_qcfail,
                                     fastq_count.is_duplicate
                                    ) )) ))

            else:
                # output to stdout much quicker
                printf("read\tis_unmapped\tmate_is_unmapped\tis_paired\tis_read1\tis_read2\tis_proper_pair\tis_qcfail\tis_duplicate\n" )
                for qname, index in reads.items():
                    fastq_count = &fastq_counts[index]
                    read_name = qname
                    printf( "%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
                            read_name,
                            fastq_count.is_unmapped,
                            fastq_count.mate_is_unmapped,
                            fastq_count.is_paired,
                            fastq_count.is_read1,
                            fastq_count.is_read2,
                            fastq_count.is_proper_pair,
                            fastq_count.is_qcfail,
                            fastq_count.is_duplicate )


    return counter, t, nh_filtered, nh_all, nm_filtered, nm_all, mapq_filtered, mapq_all, max_hi
