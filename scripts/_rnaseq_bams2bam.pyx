from pysam.csamtools cimport *

import collections, array, struct, sys
import CGAT.Experiment as E



def filter( Samfile genome_samfile,
            Samfile output_samfile,
            Samfile output_mismapped,
            Samfile transcripts_samfile,
            Samfile junctions_samfile,
            transcripts,
            regions = None,
            unique = False,
            remove_contigs = None,
            colour_mismatches = False,
            ignore_mismatches = False,
            ignore_junctions = True,
            ignore_transcripts = False ):
    '''
    To conserve memory, the tid and NM flag from *transcripts_samfile*
    are packed into memory. As a consequence, this method requires 
    max(NM) < 256 (2^8) and max(tid) < 16777216 (2^24)

    If *unique* is set, only uniquely matching reads will be output.

    If *remove_contigs* is given, contigs that are in remove_contigs will
    be removed. Note that this will only remove the alignment, but not
    all alignments of a particuluar read and the NH flag will *NOT* be
    updated.

    If *colour_mismatches* is set, the ``CM`` tag will be used
    to count differences. By default, the ``NM`` tag is used.
    The tag that is used needs to present in both *transcripts_samfile*
    and *genome_samfile*.

    If *ignore_mismatches* is set, the number of mismatches is ignored.

    If *regions* is given, alignments overlapping regions will be removed.

    '''
    cdef int ninput = 0
    cdef int nunmapped_genome = 0
    cdef int nunmapped_transcript = 0
    cdef int nmismapped = 0
    cdef int noutput = 0
    cdef int nnotfound = 0
    cdef int nlocation_ok = 0
    cdef int nnonunique = 0
    cdef int ntested = 0
    cdef int nremoved_contigs = 0
    cdef int nremoved_junctions = 0
    cdef int nskipped_junctions = 0
    cdef int nadded_junctions = 0
    cdef int nremoved_regions = 0

    cdef bint c_unique = unique
    cdef bint c_test_mismatches = not ignore_mismatches
    cdef bint c_test_junctions = not ignore_junctions
    cdef bint c_test_transcripts = not ignore_transcripts
    cdef bint c_test_regions = regions
    cdef int * remove_contig_tids 
    cdef int nremove_contig_tids = 0
    cdef AlignedRead read
    cdef AlignedRead match
    
    # build index
    # this method will start indexing from the current file position
    # if you decide
    cdef int ret = 1
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
    cdef int * map_tid2tid

    cdef int32_t read_mismatches = 0

    # decide which tag to use
    cdef char * nm_tag = 'NM' 
    cdef char * cm_tag = 'CM'
    cdef char * tag

    if colour_mismatches:
        tag = cm_tag
    else:
        tag = nm_tag

    # set with junctions that are ignored
    skip_junctions = set()

    # build list of junctions
    if c_test_junctions:
        E.info( "building junction read index" )

        def _gen2(): return array.array('B') 
        junctions_index = collections.defaultdict( _gen2 )
        ret = 1
        while ret > 0:
            ret = samread( junctions_samfile.samfile, b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4: continue
                qname = bam1_qname( b )
                v = bam_aux_get(b, tag)
                nm = <int32_t>bam_aux2i(v)
                junctions_index[qname].append( nm )

        E.info( "built index for %i junction reads" % len(junctions_index))

        map_tid2tid = <int*>calloc(len(junctions_samfile.references), sizeof(int) )
        
        for x, contig_name in enumerate(junctions_samfile.references):
            map_tid2tid[x] = genome_samfile.gettid( contig_name )

    if c_test_transcripts:
        E.info( "building transcriptome read index" )
        # L = 4 bytes
        def _gen(): return array.array('L') 
        transcriptome_index = collections.defaultdict( _gen )
        ret = 1
        while ret > 0:
            ret = samread( transcripts_samfile.samfile, b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4: continue
                qname = bam1_qname( b )
                tid = b.core.tid
                v = bam_aux_get(b, tag)
                nm = <int32_t>bam_aux2i(v)
                transcriptome_index[qname].append( (tid << 8) | nm )

        E.info( "built index for %i transcriptome reads" % len(transcriptome_index))

    bam_destroy1( b )

    # setup list of contigs to remove:
    if remove_contigs:
        nremove_contig_tids = len(remove_contigs)
        remove_contig_tids = <int*>malloc( sizeof(int) * nremove_contig_tids )
        for x, rname in enumerate( remove_contigs):
            remove_contig_tids[x] = genome_samfile.gettid( rname )
    
    E.info( "starting filtering" )

    for read in genome_samfile:

        ninput += 1
        # if ninput > 10000: break

        # is unmapped?
        if read._delegate.core.flag & 4:
            nunmapped_genome += 1
            noutput += 1
            output_samfile.write( read )
            continue

        # optionally remove non-unique reads
        if c_unique:
            v = bam_aux_get(read._delegate, 'NH')
            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh > 1: 
                    nnonunique += 1
                    continue

        # optionally remove reads matched to certain contigs
        if nremove_contig_tids:
            skip = 0
            for x from 0 <= x < nremove_contig_tids:
                if remove_contig_tids[x] == read.tid:
                    nremoved_contigs += 1
                    skip = 1
                    break
            if skip: continue

        g_contig = genome_samfile.getrname( read.tid )
                
        # optionally remove reads mapped to certain regions
        if c_test_regions:
            intervals = regions.get( g_contig, read.pos, read.aend )
            skip = 0
            for start, end, value in intervals:
                if read.overlap( start, end ):
                    nremoved_regions += 1
                    skip = 1
                    break

            if skip: continue

        if c_test_junctions:
            if read.qname in junctions_index:
                # can compress index before, depends on
                # how many reads you expect to test
                nm = min(junctions_index[read.qname])
                read_mismatches = read.opt(tag)
                if nm > read_mismatches: 
                    nremoved_junctions += 1
                    continue
                else:
                    skip_junctions.add(read.qname)
                                        
        # set mapped = True, if read is mapped to transcripts
        # set location_ok = True, if read matches in expected location
        # according to transcripts
        location_ok = False
        mapped = False

        if c_test_transcripts:
            # get transcripts that read matches to
            try:
                matches = transcriptome_index[read.qname]
            except KeyError:
                nnotfound += 1
                matches = None
                
            if matches:

                if c_test_mismatches:
                    read_mismatches = read.opt(tag)

                for val in matches:
                    transcript_tid = val >> 8
                    # ignore reads that are mapped to transcripts with
                    # more mismatches than the genomic location
                    if c_test_mismatches:
                        nm = val & 255
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
            ntested += 1
            nlocation_ok += 1
            noutput += 1
            output_samfile.write( read )

        elif mapped:
            ntested += 1
            nmismapped += 1
            if output_mismapped:
                output_mismapped.write( read )

        else:
            nunmapped_transcript += 1
            noutput += 1
            output_samfile.write( read )

    # add junctions
    if c_test_junctions:
        E.info("adding junctions")
        junctions_samfile.reset()

        # rebuild contig map for junctions:
        if remove_contigs:
            remove_contig_tids = <int*>malloc( sizeof(int) * nremove_contig_tids )
            for x, rname in enumerate( remove_contigs):
                remove_contig_tids[x] = junctions_samfile.gettid( rname )
        
        for read in junctions_samfile:

            # optionally remove reads matched to certain contigs
            if nremove_contig_tids:
                skip = 0
                for x from 0 <= x < nremove_contig_tids:
                    if remove_contig_tids[x] == read.tid:
                        nremoved_contigs += 1
                        skip = 1
                    break
                if skip: 
                    nskipped_junctions += 1
                    continue

            if read.qname in skip_junctions:
                nskipped_junctions += 1
            else:
                nadded_junctions += 1
                noutput += 1
                # map tid from junction database to genome database
                read.tid = map_tid2tid[read.tid]
                output_samfile.write(read)
        
        free( map_tid2tid )

    c = E.Counter()
    c.input = ninput
    c.removed_nonunique = nnonunique
    c.removed_mismapped = nmismapped
    c.removed_contigs = nremoved_contigs
    c.removed_junctions = nremoved_junctions
    c.removed_regions = nremoved_regions
    c.skipped_junction_reads = len(skip_junctions)
    c.skipped_junctions = nskipped_junctions
    c.added_junctions = nadded_junctions
    c.output = noutput
    c.unmapped_genome = nunmapped_genome
    c.unmapped_transcript = nunmapped_transcript
    c.notfound = nnotfound
    c.location_ok = nlocation_ok 
    c.location_tested = ntested

    if nremove_contig_tids:
        free( remove_contig_tids )

    E.info( "filtering finished" )

    return c
