#cimport csamtools

from csamtools cimport *
import pysam

import collections, array, struct, sys, itertools
import Experiment as E
import Intervals, GTF
import numpy

CountResult = collections.namedtuple( "Counts", "upstream upstream_utr cds downstream_utr downstream" )

class RangeCounter:
    
    def __init__(self ):
        pass

    def setup( self, ranges ):
        self.counts = numpy.zeros( Intervals.getLength( ranges ) )
        
    def getCounts( self, contig, ranges, length = 0 ):
        '''count from a set of ranges.

        The ranges are scaled towards the length of the intervals that
        are being counted.
        '''

        self.setup( ranges )
        self.count( contig, ranges )
        
        if length > 0:
            lnormed = float(length)
            lunnormed = float(len(self.counts))
            if lunnormed == 0: return numpy.zeros(0)

            if lunnormed > lnormed:
                # compress by taking only counts at certain points within intervals
                take = numpy.unique( numpy.array( numpy.floor( 
                            numpy.arange( 0, lunnormed - 1, lunnormed / lnormed ) ), dtype = int ) )
            elif lunnormed == lnormed:
                # same size, no scaling
                take = numpy.arange( 0, lnormed, dtype = int )
            else:
                # expand by taking more counts
                take = numpy.array( numpy.arange(0, lunnormed, lunnormed / lnormed ), dtype = int )[:length]
                
            assert len(take) == length, "size of arrays unequal: %i != %i from %i" % (len(take), length, lunnormed)

            return self.counts[take]
        else:
            return self.counts

class RangeCounterBAM(RangeCounter):

    def __init__(self, Samfile samfile, int shift, *args, **kwargs ):
        RangeCounter.__init__(self, *args, **kwargs )
        self.samfile = samfile
        self.shift = shift

    def count(self, contig, ranges ):

        if len(ranges) == 0: return

        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end
        cdef int interval_width
        cdef int current_offset
        cdef Samfile samfile = self.samfile
        cdef AlignedRead read
        cdef int length

        current_offset = 0
        counts = self.counts
        length = len(counts)
        # for peak counting follow the MACS protocol:
        # see the function def __tags_call_peak in PeakDetect.py
        #
        # In words
        # Only take the start of reads (taking into account the strand)
        # for counting, extend reads by shift
        # on + strand shift tags upstream
        # i.e. look at the downstream window
        # note: filtering?
        # note: does this work with paired-end data?
        cdef int offset = self.shift // 2

        cdef int last_end = ranges[0][0]

        current_offset = 0

        for start, end in ranges:

            if offset > 0:
                # on the + strand, shift tags upstream
                xstart, xend = max(0, start - offset), max(0, end - offset)

                for read in samfile.fetch( contig, xstart, xend ):
                    if read.is_reverse: continue
                    rstart = max( 0, read.pos - xstart - offset)
                    rend = min( length, read.pos - xstart + offset) 
                    for i from rstart <= i < rend: counts[i] += 1

                # on the - strand, shift tags downstream
                xstart, xend = max(0, start + offset), max(0, end + offset)

                for read in samfile.fetch( contig, xstart, xend ):
                    if not read.is_reverse: continue
                    rstart = max( 0, read.pos + read.rlen - xstart - offset)
                    rend = min( length, read.pos + read.rlen - xstart + offset) 
                    for i from rstart <= i < rend: counts[i] += 1
            else:
                for read in samfile.fetch( contig, start, end ):
                    rstart = max( start, read.pos ) - start + current_offset
                    rend = min( end, read.aend) - start + current_offset
                    for i from rstart <= i < rend: counts[i] += 1

            current_offset += end - start

class RangeCounterBed(RangeCounter):

    def __init__(self, bedfile, *args, **kwargs ):
        RangeCounter.__init__(self, *args, **kwargs )
        self.bedfile = bedfile
        
    def count(self, contig, ranges ):
        
        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end, tstart, tend

        if len(ranges) == 0: return

        bedfile = self.bedfile
        counts = self.counts

        cdef int length = len(counts)
        cdef int current_offset = ranges[0][0]
        cdef int last_end = ranges[0][0]
        cdef int interval_width = 0

        for start, end in ranges:
            current_offset += start - last_end
            interval_width = end - start
            
            try:
                for bed in bedfile.fetch( contig, max(0,start), end, parser = pysam.asBed() ):
                    # truncate to range of interest
                    tstart, tend = max( bed.start, start), min( bed.end, end )
                    rstart = max( 0, tstart - current_offset )
                    rend = min( length, tend - current_offset )
                    for i from rstart <= i < rend: counts[i] += 1
            except ValueError:
                # contig not present
                pass

            last_end = end

class IntervalsCounter:
    
    format = "%i"

    def __init__(self, normalization = None, *args, **kwargs ):

        # normalization method to use before aggregation
        self.normalization = normalization        
        # lables of categories to count (introns, exons, ...)
        self.fields = []
        # number of items counted in each category (1012 introns, 114 exons, ...)
        self.counts = []
        # aggregate counts
        self.aggregate_counts = []

    def add( self, field, length ):
        self.counts.append( 0 )
        self.fields.append( field )
        if self.normalization in ("max", "sum", "total-max", "total-sum" ):
            self.dtype = numpy.float
        else:
            self.dtype = numpy.int

        self.aggregate_counts.append( numpy.zeros( length, dtype = self.dtype ) )

    def setNormalization( self, normalization ):
        self.normalization = normalization

        if self.normalization in ("max", "sum", "total-max", "total-sum" ):
            self.dtype = numpy.float
            self.format = "%6.4f"
        else:
            self.dtype = numpy.int
            self.format = "%i"

        for x,c in enumerate(self.aggregate_counts):
            self.aggregate_counts[x] = numpy.zeros( len(c), dtype = self.dtype ) 

    def aggregate( self, *counts ):

        if self.normalization == "total-max":
            norm = float(max( [max(c) for c in counts if len(c) > 0] ))
        elif self.normalization == "total-sum":
            norm = float(sum( [sum(c) for c in counts ] ))
        else:
            norm = 0

        for x, xx in enumerate( zip( self.aggregate_counts, counts)):
            agg, c = xx
            self.counts[x] += 1

            if len(c) == 0: continue
            
            if norm > 0:
                cc = c / norm
            elif self.normalization == "max":
                m = max(c)
                if m == 0: continue
                cc = c / m
            elif self.normalization == "sum":
                m = sum(c)
                if m == 0: continue
                cc = c / m
            else:
                cc = c
                
            agg += cc

    def buildMatrix( self, normalize = False ):
        '''build single matrix with all counts.
        
        cols = intervals counted (exons, introns) )
        rows = number of bins in intervals (every 10b / 10kb = 1000 bins)
        '''
        max_counts = max( [len(x) for x in self.aggregate_counts] )

        matrix = numpy.concatenate( list(itertools.chain.from_iterable( 
                [ (x, [0] * (max_counts - len(x))) for x in self.aggregate_counts ] )))
        
        matrix.shape = (len(self.aggregate_counts), max_counts )

        # normalize
        for x in range( len(self.counts)):
            matrix[:,x] /= self.counts[x]

        return matrix

    def writeMatrix( self, outfile ):
        '''write aggregate counts to *outfile*.

        Output format by default is a tab-separated table.
        '''
        outfile.write("bin\t%s\n" % "\t".join(self.fields) )        

        matrix = self.buildMatrix()

        for row, cols in enumerate(matrix.transpose()):
            outfile.write( "%i\t%s\n" % (row, "\t".join( [ self.format % x for x in cols ] ) ))

    def __str__(self):
        return "%s=%s" % (self.name, ",".join( [str(sum(x)) for x in self.aggregate_counts]) )

class GeneCounter( IntervalsCounter ):

    name = "geneprofile"
    
    def __init__(self, counter, 
                 int resolution_upstream,
                 int resolution_upstream_utr,
                 int resolution_cds,
                 int resolution_downstream_utr,
                 int resolution_downstream,
                 int extension_upstream = 0, 
                 int extension_downstream = 0,
                 *args,
                 **kwargs ):

        IntervalsCounter.__init__(self, *args, **kwargs )

        self.counter = counter
        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream 
        self.resolution_cds = resolution_cds
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        self.resolution_upstream_utr = resolution_upstream_utr
        self.resolution_downstream_utr = resolution_downstream_utr

        for field, length in zip( 
            ("upstream", "upstream_utr", "cds", "downstream_utr", "downstream"),
            (resolution_upstream,
             resolution_upstream_utr,
             resolution_cds,
             resolution_downstream_utr,
             resolution_downstream ) ):
            self.add( field, length )
        
    def count( self, gtf ):
        '''build ranges to be analyzed from a gene model.

        Returns a tuple with ranges for cds, upstream_utr, downstream_utr, upstream, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        exon_start, exon_end = exons[0][0], exons[-1][1]
        self.cds = GTF.asRanges( gtf, "CDS" )

        if len(self.cds) == 0: return 0

        cds_start, cds_end = self.cds[0][0], self.cds[-1][1]
        utrs = Intervals.truncate( exons, self.cds )
        self.upstream_utr = [ x for x in utrs if x[1] <= cds_start ]
        self.downstream_utr = [ x for x in utrs if x[0] >= cds_end ]
        self.upstream = [ ( exon_start - self.extension_upstream, exon_start ), ] 
        self.downstream = [ ( exon_end, exon_end + self.extension_downstream ), ]
        
        E.debug("counting cds" )
        self.counts_cds = self.counter.getCounts( contig, self.cds, self.resolution_cds  )
        E.debug("counting upstream_utr" )
        self.counts_upstream_utr = self.counter.getCounts( contig, self.upstream_utr, self.resolution_upstream_utr )
        E.debug("counting downstream_utr" )
        self.counts_downstream_utr = self.counter.getCounts( contig, self.downstream_utr, self.resolution_downstream_utr )
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts( contig, self.upstream, self.resolution_upstream ) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts( contig, self.downstream, self.resolution_downstream )
        E.debug("counting finished" )

        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_cds = self.counts_cds[::-1]
            self.counts_upstream_utr, self.counts_downstream_utr = self.counts_downstream_utr[::-1], self.counts_upstream_utr[::-1]
            self.counts_upstream, self.counts_downstream = self.counts_downstream[::-1], self.counts_upstream[::-1]

        self.aggregate( self.counts_upstream,
                        self.counts_upstream_utr,
                        self.counts_cds,
                        self.counts_downstream_utr,
                        self.counts_downstream )

        return 1

class TSSCounter( IntervalsCounter ):
    '''count profile at transcription start/end site.

    The transcription start site is the beginning of the first exon.
    The transcription termination site is the end of the last exon.

    '''

    name = "tssprofile"
    
    def __init__(self, counter, extension_out = 0, extension_in = 0, *args, **kwargs):
        IntervalsCounter.__init__(self, *args, **kwargs )

        self.counter = counter
        self.extension_out = extension_out
        self.extension_in = extension_in 

        for field, length in zip( ("tss", "tts"),
                                  ( extension_out + extension_in,
                                    extension_out + extension_in ) ):
            self.add( field, length )
        
    def count( self, gtf ):

        contig, strand = gtf[0].contig, gtf[0].strand
        exons = GTF.asRanges( gtf, "exon" )
        self.tss, self.tts = exons[0][0], exons[-1][1]

        # no max(0, ...) here as these ranges need to have always the same length
        self.tss_ranges = [ (self.tss - self.extension_out, 
                             self.tss + self.extension_in), ]
        self.tts_ranges = [ (self.tts - self.extension_in, 
                             self.tts + self.extension_out), ]
        E.debug( "tss=%s, tts=%s" % (self.tss_ranges, self.tts_ranges) )

        self.counts_tss = self.counter.getCounts( contig, self.tss_ranges )
        self.counts_tts = self.counter.getCounts( contig, self.tts_ranges )

         ## revert for negative strand
        if strand == "-":
            self.counts_tss, self.counts_tts = self.counts_tts[::-1], self.counts_tss[::-1]
            
        self.aggregate( self.counts_tss, self.counts_tts )

        return 1

def count( counters,
           gtf_iterator):
    '''
    '''

    c = E.Counter()
    counts = [0] * len(counters)

    iterations = 0

    for iteration, gtf in enumerate(gtf_iterator):
        c.input += 1
        for x, counter in enumerate(counters):
            counter.count( gtf )
            counts[x] += 1
            
        if iteration % 100 == 0:
            E.debug( "iteration %i: counts=%s" % (iteration, ",".join( map( str, counters) ) ))

    E.info( "counts: %s: %s" % (str(c), ",".join( map(str,counts)))) 
    return



        
    
        
