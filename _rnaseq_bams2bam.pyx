#cimport csamtools

from csamtools cimport *

import collections, array, struct
import Experiment as E

def filter( Samfile genome_samfile,
            Samfile transcripts_samfile,
            Samfile output_samfile,
            Samfile output_mismapped,
            transcripts,
            unique = False ):
    '''
    This method uses the flag NM.

    To conserve memory, the tid and NM flag from *transcripts_samfile*
    are packed into memory. As a consequence, this method requires 
    max(NM) < 256 (2^8) and max(tid) < 16777216 (2^24)
    '''
    cdef int ninput = 0
    cdef int nunmapped_genome = 0
    cdef int nunmapped_transcript = 0
    cdef int nmismapped = 0
    cdef int noutput = 0
    cdef int nnotfound = 0
    cdef int nlocation_ok = 0
    cdef int nnonunique = 0
    cdef bint c_unique = unique

    cdef AlignedRead read
    cdef AlignedRead match

    E.info( "building index" )
    # build index
    # this method will start indexing from the current file position
    # if you decide
    cdef int ret = 1
    cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t) )
    cdef uint64_t pos
    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int tid
    cdef int transript_tid
    cdef long val

    # L = 4 bytes
    def _gen(): return array.array('L') 
    index = collections.defaultdict( _gen )

    while ret > 0:
        ret = samread( transcripts_samfile.samfile, b)
        if ret > 0:
            # ignore unmapped reads
            if b.core.flag & 4: continue
            qname = bam1_qname( b )
            tid = b.core.tid
            v = bam_aux_get(b, 'NM')
            nm = <int32_t>bam_aux2i(v)
            index[qname].append( (tid << 8) | nm )

    bam_destroy1( b )
    
    E.info( "built index for %i reads" % len(index))
    
    E.info( "starting filtering" )

    for read in genome_samfile:
        ninput += 1
        # if ninput > 10000: break

        if c_unique:
            v = bam_aux_get(read._delegate, 'NH')
            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh > 1: 
                    nnonunique += 1
                    continue

        # is unmapped?
        if read._delegate.core.flag & 4:
            nunmapped_genome += 1
            noutput += 1
            output_samfile.write( read )
            continue

        # get transcripts that read matches to
        try:
            matches = index[read.qname]
        except KeyError:
            nnotfound += 1
            noutput += 1
            output_samfile.write( read )
            continue
        
        g_contig = genome_samfile.getrname( read.tid )

        # set mapped = True, if read is mapped to transcripts
        # set location_ok = True, if read matches in expected location
        # according to transcripts
        location_ok = False
        mapped = False
        read_mismatches = read.opt("NM")

        for val in matches:
            transcript_tid = val >> 8
            nm = val & 255

            # ignore reads that are mapped to transcripts with
            # more mismatches than the genomic location
            if nm > read_mismatches: continue

            mapped = True

            # find transcript
            transcript_id = transcripts_samfile._getrname( transcript_tid )
            gtfs = transcripts[transcript_id]
            t_contig, t_start, t_end = gtfs[0].contig, gtfs[0].start, gtfs[-1].end
            
            # does read map to genomic transcript location?
            if g_contig == t_contig and t_start <= read.pos <= t_end:
                location_ok = True
                break

        if location_ok:
            nlocation_ok += 1
            noutput += 1
            output_samfile.write( read )

        elif mapped:
            nmismapped += 1
            if output_mismapped:
                output_mismapped.write( read )

        else:
            nunmapped_transcript += 1
            noutput += 1
            output_samfile.write( read )

    c = E.Counter()
    c.input = ninput
    c.unmapped_genome = nunmapped_genome
    c.non_unique = nnonunique
    c.unmapped_transcript = nunmapped_transcript
    c.mismapped = nmismapped
    c.output = noutput
    c.notfound = nnotfound
    c.location_ok = nlocation_ok 

    E.info( "filtering finished" )

    return c
