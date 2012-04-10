#cimport csamtools

from csamtools cimport *
import collections, array, struct
import Experiment as E
import numpy

PeakShapeResult = collections.namedtuple( "PeakShapeResult",
                                          "interval_width npeaks "
                                          "peak_center peak_width peak_height peak_relative_pos "
                                          "median closest_half_height furthest_halfheight "
                                          "bins counts" )
                                          

def count( samfiles,
           contig, 
           int start, 
           int end,
           bins,
           shift,
           float peak_ratio = 0.90,
           normalization = "none"):
    '''
    '''

    # collect pileup profile in region bounded by start and end.
    cdef int interval_width = end - start
    cdef int peak_nreads = 0
    cdef int npeaks = 0
    cdef int peak_center = 0
    cdef int i, x
    cdef int xstart, xend, rstart, rend
    cdef int nreads = 0
    cdef int * ccounts = <int*>calloc( interval_width, sizeof(int))

    # for peak counting follow the MACS protocol:
    # see the function def __tags_call_peak in PeakDetect.py
    #
    # In words
    # Only take the start of reads (taking into account the strand)
    # add d/2=shift to each side of peak and start accumulate counts.
    # for counting, extend reads by shift
    # on + strand shift tags upstream
    # i.e. look at the downstream window
    # note: filtering?
    # note: does this work with paired-end data?
    #cdef int offset = shift // 2

    for offset, samfile in zip(shift,samfiles):
    
        offset = offset // 2
        xstart, xend = max(0, start - offset), max(0, end - offset)

        for read in samfile.fetch( contig, xstart, xend ):
            if read.is_reverse: continue
            nreads += 1
            rstart = max( 0, read.pos - xstart - offset)
            rend = min( interval_width, read.pos - xstart + offset) 
            for i from rstart <= i < rend: ccounts[i] += 1

        # on the - strand, shift tags downstream
        xstart, xend = max(0, start + offset), max(0, end + offset)

        for read in samfile.fetch( contig, xstart, xend ):
            if not read.is_reverse: continue
            nreads += 1
            rstart = max( 0, read.pos + read.rlen - xstart - offset)
            rend = min( interval_width, read.pos + read.rlen - xstart + offset) 
            for i from rstart <= i < rend: ccounts[i] += 1

    counts = numpy.zeros( interval_width, dtype = numpy.int )
    for i from 0 <= i < interval_width: counts[i] = ccounts[i]
    peak_nreads = max(counts)
    peaks = numpy.array( range(0,interval_width) )[ counts >= peak_nreads ]
    peak_center = peaks[len(peaks) // 2] 

    # define peak height
    cdef int peak_height = numpy.ceil( peak_ratio * peak_nreads )
    
    peaks = numpy.array( range(0,interval_width) )[ counts >= peak_height ]
    npeaks = len(peaks)
    peak_start = peaks[0]
    peak_end = peaks[-1]
    cdef int peak_width = peak_end - peak_start

    # find size of peak
    # cdef int peak_start, peak_end
    # peak_start, peak_end = peak_center, peak_center

    # for peak_center >= i >= 0:
    #     if counts[i] < peak_height:
    #         peak_start = i
    #         break
    
    # for peak_center <= i < interval_width:
    #     if counts[i] < peak_height:
    #         peak_end = i
    #         break
    
    # closest and furthest distance of peak to half-height

    cdef int half_height = peak_height // 2
    cdef int left_first, left_last, right_first, right_last 
    left_first, left_last, right_first, right_last = peak_center, 0, interval_width, peak_center

    for i from 0 <= i < peak_center:
        if counts[i] >= half_height: 
            left_first = i
            break
        
    for i from peak_center > i >= 0:
        if counts[i] <= half_height: 
            left_last = i
            break
        
    for i from peak_center < i < interval_width:
        if counts[i] <= half_height: 
            right_first = i
            break

    for i from interval_width > i >= peak_center:
        if counts[i] >= half_height: 
            right_last = i
            break

    cdef int furthest_dist = max( peak_center - left_first, right_last - peak_center )
    cdef int closest_dist = min( peak_center - left_last, right_first - peak_center )

    # output bins
    nbins = len(bins) - 1
    hist = numpy.zeros( nbins, dtype = numpy.int )
    xstart = peak_center + bins[0] 

    for i from 0 <= i < nbins:
        xend = peak_center + bins[i+1]
        if xstart >= 0 and xend < interval_width:
            hist[i] = sum( counts[xstart:xend] ) 
        xstart = xend
        
    # debugging
    #for x,v in enumerate( hist ):
        # print x, "*" * v
    # for x,v in enumerate( counts ):
    #     print x, "*" * v
        
    result = PeakShapeResult._make( (interval_width, npeaks,
                                     start + peak_center, 
                                     peak_width, peak_nreads, 
                                     abs((interval_width // 2) - peak_center), 
                                     numpy.median(counts),
                                     closest_dist, furthest_dist,
                                     bins,
                                     hist ) )

    free( ccounts )
    
    return result

    
