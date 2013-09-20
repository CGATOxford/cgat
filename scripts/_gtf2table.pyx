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

import bx
import bx.bbi.bigwig_file
import bx.intervals.io
import bx.intervals.intersection
import alignlib
import numpy
import CGAT.IndexedGenome as IndexedGenome
import pysam


class Counter:
    """
    This class does not remove small exons/introns,
    so beware ENSEMBL that stores frameshifts as mini-exons.
    """

    header = [ "contig", "strand" ]

    mMinIntronSize = 10

    def __init__(self, fasta = None, section = None, options = None, prefix = None):

        self.fasta = fasta
        self.section = section
        self.options = options
        
        if prefix:
            self.header = tuple(["%s%s" % (prefix, x) for x in self.header ])

        # if true, entry is skipped
        self.skip = False

        # counter 
        self.counter = E.Counter()

    def __call__(self, gffs):
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
        """merge small introns into single exons. The following features are aggregated
        as exons: exon, CDS, UTR, UTR3, UTR5
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
        assert len(exons) > 0, "no exons in gene"
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

    def getUTRs( self ):
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
            counts[ pos ] += 1

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
               
    # discard segments with size > mMaxLength in order
    # to avoid out-of-memory
    mMaxLength = 100000

    def __init__(self, bamfiles, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles

    def count(self):
        
        segments = self.getSegments()

        cdef AlignedRead read
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
            if end - start > self.mMaxLength: return []
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
    
    def __init__(self, bamfiles, *args, weight_multi_mapping = False, **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles
        self.weight_multi_mapping = weight_multi_mapping

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

        for start, end in segments:
            for samfile in self.mBamFiles:
                last_any_pos = -1
                last_sense_pos = -1
                last_anti_pos = -1
                for read in samfile.fetch( contig, start, end ):
                    if not read.overlap( start, end ): continue
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

        for start, end in segments:
            for samfile in self.mBamFiles:
                last_any_pos = -1
                last_sense_pos = -1
                last_anti_pos = -1
                for read in samfile.fetch( contig, start, end ):
                    if not read.overlap( start, end ): continue

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
