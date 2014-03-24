# cython: profile=True
#cimport csamtools
from pysam.csamtools cimport *

import collections, array, struct
import CGAT.Experiment as E

import os
import sys
import string
import re
import optparse
import math
import time
import tempfile
import subprocess
import types
import bisect
import array
import collections
import itertools
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Stats as Stats
import CGAT.SequenceProperties as SequenceProperties
import CGAT.Genomics as Genomics
import CGAT.Intervals as Intervals
import CGAT.IndexedGenome as IndexedGenome

import numpy
import pysam

class Counter:
    """
    This class does not remove small exons/introns,
    so beware ENSEMBL that stores frameshifts as mini-exons.
    """

    header = [ "contig", "strand" ]

    mMinIntronSize = 10

    def __init__(self, fasta = None, section = None, 
                 options = None, prefix = None):

        self.fasta = fasta
        self.section = section
        self.options = options
        
        if prefix:
            self.header = tuple(["%s%s" % (prefix, x) for x in self.header ])

        # if true, entry is skipped
        self.skip = False

        # counter 
        self.counter = E.Counter()

    def update(self, gffs):
        self.mGFFs = gffs
        self.skip = False
        self.count()

    def __str__(self):
        return "\t".join( (self.contig, self.strand) )

    def getHeader(self):
        if self.section:
            return "\t".join( ["%s_%s" % (self.section, x) for x in self.header] )
        else:
            return "\t".join( self.header )

    def count(self):
        self.contig = self.getContig()
        self.strand = self.getStrand()

    def getContig(self):
        contig = self.mGFFs[0].contig
        if self.fasta: contig = self.fasta.getToken( contig)
        return contig

    def getStrand(self):
        return self.mGFFs[0].strand

    def getGeneId( self ):
        return self.mGFFs[0].gene_id
    
    def getSequence(self, segments):
        """get sequence from a set of segments."""

        contig = self.getContig()
        strand = self.getStrand()
        
        s = []
        for start,end in segments:
            s.append(self.fasta.getSequence( contig, strand, start, end ))
            
        if Genomics.IsNegativeStrand( strand ):
            return( string.translate( "".join(s), string.maketrans("ACGTacgt", "TGCAtgca") )[::-1] )
        else:
            return "".join(s)

    def getExons( self ):
        """merge small introns into single exons. The following features are
        aggregated as exons: exon, CDS, UTR, UTR3, UTR5

        Exons are always sorted by coordinate, irrespective of strand.
        """
        ranges = GTF.asRanges( self.mGFFs, feature = ("exon", "CDS", "UTR", "UTR5", "UTR3" ) )
        assert len(ranges) > 0, "no exons in gene"
        return Intervals.combineAtDistance( ranges,
                                            self.mMinIntronSize )

    def getCDS( self ):
        """merge small introns into single exons. The following features are aggregated
        as exons: exon, CDS, UTR, UTR3, UTR5
        """
        ranges = GTF.asRanges( self.mGFFs, feature = ( "CDS", ) )
        return Intervals.combineAtDistance( ranges,
                                            self.mMinIntronSize )

    def getIntrons( self ):
        exons = self.getExons()
        assert len(exons) > 0, "no exons in gene %s" % self.getGeneId()
        introns = []
        last = exons[0][1]
        for e in exons[1:]:
            introns.append( (last, e[0] ) )
            last = e[1]
        return introns

    def getSegments( self ):
        if self.section == "exons":
            return self.getExons()
        elif self.section == "introns":
            return self.getIntrons()
        else:
            return self.getExons()

    def getJunctions(self, segments):
        '''return a list of junctions from segments.'''
        junctions = []
        last_end = segments[0][1]
        for start,end in segments[1:]:
            junctions.append((last_end,start))
            last_end = end
        return junctions

    def getUTRs(self):
        '''return a tuple with 5' and 3' UTRs.'''
        exons = self.getExons()
        cds = self.getCDS()
        utr3, utr5 = [], []
        if len(cds) == 0: return utr5, utr3
        strand = self.getStrand()
        cds_start, cds_end = cds[0][0], cds[-1][1]
        midpoint = ( cds_end - cds_start ) / 2 + cds_start
        for start, end in Intervals.truncate( exons, cds ):
            if end - start > 3:
                if start < midpoint:
                    utr5.append( (start,end) )
                else:
                    utr3.append( (start,end) )
        if strand == "-":
            utr5, utr3 = utr3, utr5
        return utr5, utr3
    
cimport numpy
DTYPE_INT = numpy.int
ctypedef numpy.int_t DTYPE_INT_t
DTYPE_FLOAT = numpy.float
ctypedef numpy.float_t DTYPE_FLOAT_t


cdef __add( numpy.ndarray[DTYPE_INT_t, ndim=1]counts, positions, int offset, int length ):
    '''add positions to counts vector'''
    cdef int p, pos
    for p in positions:
        pos = p - offset
        for pos from 0 <= pos < length:
            counts[pos] += 1

##-----------------------------------------------------------------------------------
class CounterReadCoverage(Counter):
    '''compute read coverage for all exons in a transcript. 

    Requires bam files to compute that coverage. Multiple bam
    files can be supplied, these will be summed up.

    Counts are separated into sense, antisense and any sense.
    '''
    
    header = ("length",) +\
        tuple( [ "%s_%s" % (x,y) for x,y in itertools.product( ("sense", "antisense", "anysense"),
                                                               ( ("pcovered", "nreads", ) + Stats.Summary().getHeaders() )) ] )
               
    # discard segments with size > max_length in order
    # to avoid out-of-memory
    max_length = 100000

    def __init__(self, bamfiles, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles

    def count(self):
        
        segments = self.getSegments()

        cdef AlignedRead read

        # remove segments with excessive length
        segments = [ x for x in segments if (x[1] - x[0]) < self.max_length ]

        if len(segments) == 0: 
            self.length = 0
            return

        cdef int length = sum( [x[1] - x[0] for x in segments ] )
        cdef numpy.ndarray[DTYPE_INT_t, ndim=1] counts_sense = numpy.zeros( length, dtype = numpy.int )
        cdef numpy.ndarray[DTYPE_INT_t, ndim=1] counts_antisense = numpy.zeros( length, dtype = numpy.int )
        cdef int p, pos, offset

        reads_sense, reads_antisense = set(), set()

        contig = self.getContig()
        if self.getStrand() == "+":
            is_reverse = False
        else:
            is_reverse = True

        l = 0
        for start, end in segments:

            offset = start - l
            for samfile in self.mBamFiles:
                for read in samfile.fetch( contig, start, end ):
                    # only count positions actually overlapping
                    positions = read.positions
                    if not positions: continue
                    if is_reverse == read.is_reverse:
                        for p in positions:
                            pos = p - offset
                            if 0 <= pos < length:
                                counts_sense[ pos ] += 1
                        reads_sense.add( read.qname )
                    else:
                        for p in positions:
                            pos = p - offset
                            if 0 <= pos < length:
                                counts_antisense[ pos ] += 1
                        reads_antisense.add( read.qname )

            l += end - start

        self.length = length
        self.nreads_sense = len(reads_sense)
        self.nreads_antisense = len(reads_antisense)
        self.nreads_anysense = self.nreads_sense + self.nreads_antisense

        counts_anysense = counts_sense + counts_antisense

        counts_sense = counts_sense[ counts_sense > 0]
        counts_antisense = counts_antisense[ counts_antisense > 0]
        counts_anysense = counts_anysense[ counts_anysense > 0]

        self.counts_sense = counts_sense
        self.counts_antisense = counts_antisense
        self.counts_anysense = counts_anysense
        
    def __str__(self):

        r = [ "%i" % self.length ]
        if self.length == 0:
            # set to 1 to permit division below
            self.length = 1

        for direction, counts, nreads in zip ( ("sense", "antisense", "anysense"),
                                               (self.counts_sense, self.counts_antisense, self.counts_anysense),
                                               (self.nreads_sense, self.nreads_antisense, self.nreads_anysense) ):
            r.append( "%5.2f" % (100.0 * len(counts) / self.length) )
            r.append( "%i" % (nreads) )
            r.append( str( Stats.Summary( counts, mode = "int" ) ) )

        return "\t".join( r )

##-----------------------------------------------------------------------------------
class CounterReadCounts(Counter):
    '''compute number of reads overlapping with exons.

    Requires bam files to compute that coverage. Multiple bam
    files can be supplied, these will be summed up.

    Both unique and non-unique counts are collected. Uniqueness
    is simply checked through alignment start position.

    Counts are separated into sense, antisense and any sense.

    If ``weight_multi_mapping`` is set, counts are weigthed by the NH flag.
    '''
    
    header = ( [ "%s_%s" % (x,y) for x,y in itertools.product( ( "sense", "antisense", "anysense"),
                                                                ( "unique_counts", "all_counts") ) ] )
    
    def __init__(self, bamfiles,
                 *args,
                 minimum_read_quality = 0,
                 weight_multi_mapping = False,
                 **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: 
            raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles
        self.weight_multi_mapping = weight_multi_mapping
        self.minimum_read_quality = minimum_read_quality

        if self.minimum_read_quality != 0:
            raise NotImplementedError('quality filtering not implemented')

    def count(self):
        '''count reads.'''
        # For performance reasons, two separate counting implementations.
        if self.weight_multi_mapping:
            self.countFloat()
        else:
            self.countInteger()

    def countFloat(self):
        '''count by weighting multi-mapping reads.'''
        
        segments = self.getSegments()
        contig = self.getContig()

        cdef float nsense_unique_counts = 0
        cdef float nsense_all_counts = 0
        cdef float nantisense_unique_counts = 0
        cdef float nantisense_all_counts = 0
        cdef float nanysense_unique_counts = 0
        cdef float nanysense_all_counts = 0
        cdef float weight = 0
        cdef int nh = 0
        cdef long last_any_pos = -1 
        cdef long last_sense_pos = -1
        cdef long last_anti_pos = -1

        if self.getStrand() == "+":
            is_reverse = False
        else:
            is_reverse = True
	    
        # count only once per read name
        counted = set()
	 
        for start, end in segments:
            for samfile in self.mBamFiles:
                last_any_pos = -1
                last_sense_pos = -1
                last_anti_pos = -1
                for read in samfile.fetch( contig, start, end ):
                    if not read.overlap( start, end ): continue
                    if read.qname in counted: continue
                    counted.add(read.qname)

                    try:
                        nh = read.opt( 'NH' )
                    except KeyError:
                        nh = 1
                    
                    weight = 1.0/nh
                        
                    nanysense_all_counts += weight
                    if last_any_pos != read.pos:
                        last_any_pos = read.pos
                        nanysense_unique_counts += weight
                        
                    if is_reverse == read.is_reverse:
                        nsense_all_counts += weight
                        if last_sense_pos != read.pos:
                            last_sense_pos = read.pos
                            nsense_unique_counts += weight
                    else:
                        nantisense_all_counts += weight
                        if last_anti_pos != read.pos:
                            last_anti_pos = read.pos
                            nantisense_unique_counts += weight

        self.result = (nsense_unique_counts,
                       nsense_all_counts,
                       nantisense_unique_counts,
                       nantisense_all_counts,
                       nanysense_unique_counts,
                       nanysense_all_counts )

    def countInteger(self):
        '''count all reads equally.'''

        segments = self.getSegments()
        contig = self.getContig()

        cdef int nsense_unique_counts = 0
        cdef int nsense_all_counts = 0
        cdef int nantisense_unique_counts = 0
        cdef int nantisense_all_counts = 0
        cdef int nanysense_unique_counts = 0
        cdef int nanysense_all_counts = 0
        cdef long last_any_pos = -1 
        cdef long last_sense_pos = -1
        cdef long last_anti_pos = -1

        if self.getStrand() == "+":
            is_reverse = False
        else:
            is_reverse = True

        # count only once per read name
        counted = set()

        for start, end in segments:
            for samfile in self.mBamFiles:
                last_any_pos = -1
                last_sense_pos = -1
                last_anti_pos = -1
                for read in samfile.fetch( contig, start, end ):
                    if not read.overlap( start, end ): continue
                    if read.qname in counted: continue
                    counted.add(read.qname)

                    nanysense_all_counts += 1
                    if last_any_pos != read.pos:
                        last_any_pos = read.pos
                        nanysense_unique_counts += 1
                        
                    if is_reverse == read.is_reverse:
                        nsense_all_counts += 1
                        if last_sense_pos != read.pos:
                            last_sense_pos = read.pos
                            nsense_unique_counts += 1
                    else:
                        nantisense_all_counts += 1
                        if last_anti_pos != read.pos:
                            last_anti_pos = read.pos
                            nantisense_unique_counts += 1

        self.result = (nsense_unique_counts,
                       nsense_all_counts,
                       nantisense_unique_counts,
                       nantisense_all_counts,
                       nanysense_unique_counts,
                       nanysense_all_counts )

    def __str__(self):
        return "\t".join( map(str, (self.result)))

##-----------------------------------------------------------------------------------
class CounterReadPairCountsFull(Counter):
    '''compute number of read pairs overlapping with exoIsoform
    Requires bam files to compute that coverage. Multiple bam
    files can be supplied, these will be summed up.

    A read pair is evaluated in the context of the intron-exon structure
    of a transcript.  The method computes the following four attributes 
    for each pair:

    1. Pair status

    proper
       The read is paired and proper and both parts align to the transcript.

    unmapped
       One read aligns to the transcript modele, the other is unmapped.

    outer
       One read aligns to the transcript models, the other is out of range.

    other
       Other configurations, such as three reads aligning to the same
       transcript, improper pairs, ...

    2. Direction status
    
    The direction status presents the direction of the reads
    with respect to the transcript. The four directions are
    FF, FR, RF, FF for a pair with F meaning forward and R
    meaning reverse.

    3. Overlap status

    Secondly, the method compiles the following attributes taking into
    account the transcript structure:

    exonic
       All aligned bases are within exons.
       
    extension
       All aligned bases are within introns, except some
       bases extend beyond the terminal exons.

    intronic
       All aligned bases are within introns.

    untronic
       Aligned bases are both in introns and exons.

    4. Splice status:

    The splice status denotes if the read is spliced and if so,
    if the splice juncions is consistent with the transcript 
    structure

    spliced
        The pair/read contains at least one splice junction and 
        the splice junction is part of the transcript structure.

        .. note::
            A read is also considered spliced if both reads in a pair 
            end and start on an exon boundary.
    
    unspliced
        The pair/read contains no splice juctions

    misspliced
        The pair/read contains splice junctions, but at least one
        splice junction not part of the transcript model.


    These three attributes are combined into counters.

    Both unique and non-unique counts are collected. Uniqueness
    is simply checked through alignment start position.

    If ``weight_multi_mapping`` is set, counts are weigthed by the NH
    flag.
    '''
    
    headers_status = ('proper', 'unmapped', 'outer', 'other')
    headers_direction = ('FF', 'FR', 'RF', 'RR')
    headers_exons = ('exonic', 'extension', 'intronic', 'untronic')
    headers_splicing = ('unspliced', 'spliced', 'misspliced')

    # minimum intron size - splicing only checked if gap within
    # read is larger than this.
    min_intron_size = 20
    
    def __init__(self, bamfiles, 
                 *args,
                 weight_multi_mapping = False,
                 minimum_read_quality = 0,
                 **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: 
            raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles
        self.weight_multi_mapping = weight_multi_mapping
        self.minimum_read_quality = minimum_read_quality
        self.header = [ '_'.join(x) 
                        for x in itertools.product( 
                                self.headers_status,
                                self.headers_direction,
                                self.headers_exons,
                                self.headers_splicing)] +\
            ['quality_pairs', 'quality_reads']

    def count(self):
        '''count all reads equally.'''

        # Obtain transcript information
        exons = self.getExons()
        introns = self.getIntrons()
        contig = self.getContig()
        junctions = self.getJunctions(exons)
        cdef bint is_reverse = False
        if self.getStrand() == "-":
            is_reverse = True

        cdef int min_intron_size = self.min_intron_size
        cdef float minimum_read_quality = self.minimum_read_quality
        cdef bint weight_multi_mapping = self.weight_multi_mapping

        # status variables
        cdef int npair_status = len(self.headers_status)
        cdef int ndirection_status = len(self.headers_direction)
        cdef int nexons_status = len(self.headers_exons)
        cdef int nspliced_status = len(self.headers_splicing)
        cdef int spliced_status = 0
        cdef int pair_status = 0
        cdef int direction_status = 0
        cdef int exons_status = 0
        cdef int quality_read_status = 0
        cdef int quality_pair_status = 0

        cdef int ncounters = \
            npair_status * ndirection_status *\
            nexons_status * nspliced_status
        
        # helper counting variables
        cdef bint is_rev1 = 0
        cdef bint is_rev2 = 0
        cdef int nbad_splice = 0
        cdef int ngood_splice = 0
        cdef int nbases_outside = 0
        cdef int nbases_block = 0
        cdef int nbases_exons = 0
        cdef int nbases_introns = 0
        cdef int nbases_total = 0
        cdef int nreads_in_pair_after = 0
        cdef int nreads_in_pair = 0
        cdef int current_exon = 0
        cdef int nexons = len(exons)
        cdef int nh = 0
        cdef float weight = 1.0
        cdef long current_exon_end = 0
        cdef long current_exon_start = 0
        cdef long exons_start = exons[0][0]
        cdef long exons_end = exons[-1][1]
        cdef long exon_start = 0
        cdef long exon_end = 0
        cdef long block_start = 0
        cdef long block_end = 0
        cdef long block_first_start = 0
        cdef long block_last_end = 0
        cdef long read_last_end = 0
        # TODO: define block struc 
        # maximum number of blocks is 100 
        # we expect few (<10)
        cdef long * block_starts = <long*>malloc(100 *sizeof(long))
        cdef long * block_ends = <long*>malloc(100 * sizeof(long))
        cdef int nblocks = 0
        cdef int current_block = 0
        cdef long max_start = 0
        cdef long min_end = 0

        if block_starts == NULL or block_ends == NULL:
            raise ValueError('could not allocated memory for blocks')

        cdef AlignedRead read1

        # define counters, add 1 for quality filtered reads
        counters = numpy.zeros(ncounters, dtype=numpy.float)
        counters.shape = (npair_status,
                          ndirection_status,
                          nexons_status,
                          nspliced_status)

        # retrieve all reads
        reads = []

        for samfile in self.mBamFiles:
            # make sure you get more than a proxy
            reads.extend(list(samfile.fetch(contig, 
                                            exons_start, 
                                            exons_end)))
            # sort by read name and position
            reads.sort(key=lambda x: (x.qname, x.pos))

            # group by read name
            for key, reads_in_pair in itertools.groupby(
                    reads,
                    key=lambda x: x.qname):
                
                reads_in_pair = list(reads_in_pair)
                pair_status = 0
                nreads_in_pair = len(reads_in_pair)
                if minimum_read_quality > 0:
                    reads_in_pair = [x for x in reads_in_pair \
                                     if x.mapq >= minimum_read_quality]
                    nreads_in_pair_after = len(reads_in_pair)
                    if nreads_in_pair != nreads_in_pair_after:
                        quality_read_status += nreads_in_pair - nreads_in_pair_after
                        quality_pair_status += 1
                    nreads_in_pair = nreads_in_pair_after
                    if nreads_in_pair == 0:
                        continue

                read1 = reads_in_pair[0]
                
                # Do I need to test for read1/read2?

                # test for read arrangement
                if nreads_in_pair == 2 and \
                   read1.is_proper_pair:
                    # proper pair, both parts in transcript model region
                    pair_status = 0
                    is_rev1 = reads_in_pair[0].is_reverse
                    is_rev2 = reads_in_pair[1].is_reverse

                elif nreads_in_pair == 1:
                    # only one read in transcript model
                    if read1.mate_is_unmapped:
                        # other read is unmapped
                        pair_status = 1
                        if read1.is_read1:
                            is_rev1 = read1.is_reverse
                            is_rev2 = read1.mate_is_reverse
                        else:
                            is_rev1 = read1.mate_is_reverse
                            is_rev2 = read1.is_reverse
                    elif reads_in_pair[0].is_proper_pair:
                        # outer: proper pair, but
                        # other read maps outside transcript model
                        pair_status = 2
                        is_rev1 = read1.is_reverse
                        is_rev2 = read1.mate_is_reverse
                    else:
                        # not a proper pair
                        pair_status = 3 
                    if read1.is_read1:
                        is_rev1 = read1.is_reverse
                        is_rev2 = read1.mate_is_reverse
                    else:
                        is_rev1 = read1.mate_is_reverse
                        is_rev2 = read1.is_reverse

                else:
                    # other configuration (> 2 reads)
                    # not proper pairs
                    pair_status = 3
                    is_rev1 = reads_in_pair[0].is_reverse
                    is_rev2 = reads_in_pair[0].mate_is_reverse

                # Iterate over blocks within reads and 
                # compute overlap with exons, introns, etc.
                # 
                # Blocks are already sorted by position
                # because reads are sorted by position
                # and blocks are always returned from
                # left-most coordinate
                ngood_splice, nbad_splice = 0, 0
                nbases_total = 0
                nbases_exons = 0
                nbases_introns = 0
                nbases_outside = 0
                nblocks = 0
                block_last_end = -1
                block_first_start = -1
                
                for read in reads_in_pair:
                    # end coordinate in genomic coordinates
                    blocks = read.blocks
                    read_last_end = -1

                    for block_start, block_end in blocks:
                        
                        # print '# a: block_start=', block_start, \
                        #         'block_end=', block_end, \
                        #         'block_first_start=', block_first_start, \
                        #         'block_last_end=', block_last_end, \
                        #         'nblocks=', nblocks
                        # check introns within read, not
                        # overall with blocks
                        if read_last_end > 0 and \
                           block_start - read_last_end >= min_intron_size:
                            if (read_last_end, block_start) in junctions:
                                ngood_splice += 1
                            else:
                                nbad_splice += 1
                        read_last_end = block_end

                        if block_last_end < block_start:
                            # new block, not overlapping with previous
                            if block_first_start >= 0:
                                nbases_total += block_last_end - block_first_start
                                block_starts[nblocks] = block_first_start
                                block_ends[nblocks] = block_last_end
                                nblocks += 1
                            block_first_start = block_start

                        block_last_end = block_end
                
                # print '# b: block_start=', block_start, \
                #             'block_end=', block_end, \
                #             'block_first_start=', block_first_start, \
                #             'block_last_end=', block_last_end, \
                #             'nblocks=', nblocks

                # close of loop
                nbases_total += block_last_end - block_first_start
                block_starts[nblocks] = block_first_start
                block_ends[nblocks] = block_last_end
                nblocks += 1

                #---------------------------------------------------
                # compute overlap with exons
                current_exon = 0
                current_block = 0
                block_start = block_starts[current_block]
                block_end = block_ends[current_block]

                while current_block < nblocks and block_start < exons_start:
                    nbases_outside += min(block_end, exons_start) - block_start
                    block_start = block_starts[current_block]
                    block_end = block_ends[current_block]
                    current_block += 1
                    
                while current_exon < nexons and current_block < nblocks:

                    exon_start, exon_end = exons[current_exon]
                    block_start = block_starts[current_block]
                    block_end = block_ends[current_block]
                    
                    if exon_end <= block_start:
                        current_exon += 1
                    elif block_end <= exon_start:
                        current_block += 1
                    else:
                        max_start = max(exon_start, block_start)
                        min_end = min(exon_end, block_end)
                        
                        nbases_exons += min_end - max_start

                        if exon_end < block_end:
                            current_exon += 1
                        elif block_end < exon_end:
                            current_block += 1
                        else:
                            current_exon += 1
                            current_block += 1

                # compute bases outside transcript
                while current_block < nblocks:
                    block_start = block_starts[current_block]
                    block_end = block_ends[current_block]
                    nbases_outside += block_end - max(block_start, exons_end)
                    current_block += 1

                #---------------------------------------------------
                # compute intron overlap
                nbases_introns = nbases_total - nbases_exons - nbases_outside

                #####################################################
                # sort out the splicing attribute
                if nbad_splice > 0:
                    # bad splice sites present
                    spliced_status = 2
                elif ngood_splice > 0:
                    # only good splice sites
                    spliced_status = 1
                else:
                    # no spliced reads
                    spliced_status = 0

                # sort out the direction attribute
                if is_rev1:
                    if is_rev2:   # RR
                        direction_status = 3
                    else:         # RF
                        direction_status = 2
                else:
                    if is_rev2:   # FR
                        direction_status = 1
                    else:         # FF
                        direction_status = 0

                # swap direction of reads according to transcript strand
                # RR = FF, FR = RF, etc.
                if is_reverse:
                    direction_status = 3 - direction_status
                    
                # print '#', exons
                # print '# bases total=', nbases_total, \
                #                       'exons=', nbases_exons,\
                #                       'introns=', nbases_introns, \
                #                       'outside=', nbases_outside, \
                #                       'nblocks=', nblocks 
                # print '# splice good=', ngood_splice, 'bad=', nbad_splice
                # print '# status=', pair_status, direction_status, exons_status, spliced_status
                # sort out the exon attribute
                if nbases_introns == 0:
                    if nbases_exons == nbases_total:
                        # only exonic
                        exons_status = 0
                    else:
                        # exonic + extension
                        exons_status = 1
                elif nbases_introns == nbases_total:
                    # only intronic
                    exons_status = 2
                else:
                    # other
                    exons_status = 3

                if weight_multi_mapping:
                    try:
                        nh = read.opt('NH')
                    except KeyError:
                        nh = 1
                    weight = 1.0/nh

                counters[pair_status][direction_status][exons_status][spliced_status] += weight

        free(block_starts)
        free(block_ends)

        self.counters = counters
        self.reads_below_quality = quality_read_status
        self.pairs_below_quality = quality_pair_status

    def __str__(self):
        return "\t".join(map(str, (self.counters.flat) )) +\
            "\t" +\
            "\t".join(map(str, (self.pairs_below_quality,
                                self.reads_below_quality)))

##-----------------------------------------------------------------------------------
class CounterReadPairCounts(CounterReadPairCountsFull):
    '''compute number of read pairs overlapping with exons.

    Requires bam files to compute that coverage. Multiple bam
    files can be supplied, these will be summed up.

    A read pair is evaluated in the context of the intron-exon structure
    of a transcript. 
    
    This method takes the output of the full counting results
    and compiles it into a subset of informative categories.
    '''
    
    def __init__(self, *args, library_type="unstranded", **kwargs ):
        CounterReadPairCountsFull.__init__(self, *args, **kwargs )
        self.header = ['counted_all',
                       'counted_spliced',
                       'counted_unspliced',
                       'sense_intronic',
                       'sense_inconsistent',
                       'sense_other',
                       'antisense',
                       'notproper',
                       'quality_pairs',
                       'quality_reads',
                       'total']
        self.library_type = library_type

    def count(self):
        '''count reads.'''
        CounterReadPairCountsFull.count(self)

        # some short-cuts
        # pair axis
        proper = 0
        # sense axis
        sense = 0
        antisense = None
        # exonic axis
        exonic = 0
        intronic = 2
        # splice axis
        spliced = 0
        unspliced = 1

        self.total_pairs = numpy.sum(self.counters.flat)

        # count everything that is not proper
        self.improper_pairs = \
            sum(numpy.sum(self.counters, axis=(1,2,3))[1:].flat)

        # only work with everything that is proper
        work = self.counters[proper]

        # set sense/antisense depending on library method
        if self.library_type == 'fr-unstranded':
            # sum over all read configurations
            # there is no antisense
            work = [numpy.sum(work, axis=0)] 
            sense, antisense = 0, None
        elif self.library_type == 'fr-secondstrand':
            # FR, RF
            sense, antisense = 1, 2
        elif self.library_type == 'fr-firststrand':
            # RF, FR
            sense, antisense = 2, 1
        else:
            raise NotImplementedError('unknown library type %s' %
                                      self.library_type)

        self.sense = numpy.sum(work[sense].flat)

        # exonic, proper, spliced              
        self.sense_proper_pairs_spliced = \
            work[sense][exonic][spliced]
        self.sense_proper_pairs_unspliced = \
            work[sense][exonic][unspliced]

        if antisense is None:
            self.antisense = 0
        else:
            # sum over direction axis
            self.antisense = \
                numpy.sum(work, axis=(1,2))[antisense]

        # intronic pairs that are sense
        self.sense_intronic_pairs = \
            numpy.sum(work[sense][intronic])

        self.sense_inconsistent_pairs = \
            numpy.sum(work[sense][intronic+1:].flat)
        
        # read pairs in sense direction, but that are not counted
        self.sense_other_pairs = self.sense - (
            self.sense_proper_pairs_spliced +
            self.sense_proper_pairs_unspliced +
            self.sense_intronic_pairs +
            self.sense_inconsistent_pairs)

        assert self.total_pairs == (
            self.sense_proper_pairs_spliced +
            self.sense_proper_pairs_unspliced +
            self.sense_intronic_pairs +
            self.sense_inconsistent_pairs +
            self.sense_other_pairs +
            self.antisense + 
            self.improper_pairs)

    def __str__(self):
        return "\t".join(map(str,(
            self.sense_proper_pairs_spliced +\
            self.sense_proper_pairs_unspliced,
            self.sense_proper_pairs_spliced,
            self.sense_proper_pairs_unspliced,
            self.sense_intronic_pairs,
            self.sense_inconsistent_pairs,
            self.sense_other_pairs,
            self.antisense,
            self.improper_pairs,
            self.pairs_below_quality,
            self.reads_below_quality,
            self.total_pairs)))
                              

