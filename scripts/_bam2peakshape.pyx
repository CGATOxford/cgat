#cimport csamtools

from pysam.csamtools cimport *

import collections, array, struct
import CGAT.Experiment as E

import numpy

PeakShapeResult = collections.namedtuple( "PeakShapeResult",
                                          "interval_width npeaks "
                                          "peak_center peak_width peak_height peak_relative_pos "
                                          "nreads "
                                          "median closest_half_height furthest_halfheight "
                                          "bins counts" )

PeakShapeCounts = collections.namedtuple( "PeakShapeCounts",
                                          "nreads median counts" )

cdef coverageInInterval( Samfile samfile,
                         contig,
                         int start, 
                         int end,
                         int shift = 0):
    '''return coverage in window on *contig* bounded by *start* and *end*.

    Reads are optionally shifted by shift/2. Reads are shifted upstream
    on the forward strand and downstream on the reverse strand.

    returns a tuple:
       nreads = number of reads/tags counted
       counts = numpy array with reads per base
    '''
    
    # for peak counting follow the MACS protocol:
    # see the function def __tags_call_peak in PeakDetect.py
    #
    # In words
    # Only take the start of reads (taking into account the strand)
    # add d/2=shift to each side of peak and start accumulate counts.
    # for counting, extend reads by 2 * shift
    # on + strand shift tags upstream
    # i.e. look at the downstream window
    # note: filtering?
    # note: does this work with paired-end data?
    cdef int offset = shift // 2
    cdef int nreads = 0

    # collect pileup profile in region bounded by maxwindow_start and maxwindow_end.
    cdef int interval_width = end - start
    cdef int * ccounts = <int*>calloc( interval_width, sizeof(int))

    if interval_width <= 0: 
        return 0, numpy.zeros(0)

    if shift:
        xstart, xend = max(0, start - offset), max(0, end - offset)

        for read in samfile.fetch( contig, xstart, xend ):
            if read.is_reverse: continue
            nreads += 1
            # extend reads from upstream pos
            rstart = read.pos
            rend = rstart + 2 * offset 
            # truncate to interval
            rstart = max( 0, rstart - xstart )
            rend = min( interval_width, rend - xstart )
            for i from rstart <= i < rend: ccounts[i] += 1

        # on the - strand, shift tags downstream
        xstart, xend = start + offset, end + offset

        for read in samfile.fetch( contig, xstart, xend ):
            if not read.is_reverse: continue
            nreads += 1
            # shift and extend reads from downstream pos
            rend = read.aend
            rstart = rend - 2 * offset
            # truncate to interval
            rstart = max( 0, rstart - xstart )
            rend = min( interval_width, rend - xstart )

            for i from rstart <= i < rend: ccounts[i] += 1
            
    else:
        for read in samfile.fetch( contig, start, end ):
            nreads += 1
            # truncate to interval
            rstart = max( 0, read.pos - start )
            rend = min( interval_width, read.aend - start )
            for i from rstart <= i < rend: ccounts[i] += 1

    # transfer to numpy count object
    counts = numpy.zeros( interval_width, dtype = numpy.int )
    for i from 0 <= i < interval_width: counts[i] = ccounts[i]

    free( ccounts )

    return nreads, counts

def count( Samfile samfile,
           contig, 
           int start, 
           int end,
           bins,
           int window_size = 0,
           float peak_ratio = 0.90,
           int shift = 0,
           only_interval = False,
           centring_method = "reads" ):
    '''
    '''

    cdef int peak_nreads = 0
    cdef int npeaks = 0
    cdef int peak_center = 0
    cdef int i, x
    cdef int xstart, xend, rstart, rend

    assert start >= 0, "start < 0"

    # maximum extend of window = interval +- window_size
    cdef int maxwindow_start = max(0, start - window_size )
    cdef int maxwindow_end = end + window_size

    # bases added at right/left of interval
    cdef int offset_right = maxwindow_end - end
    cdef int offset_left = start - maxwindow_start

    cdef int interval_width = end - start

    # get counts in window
    nreads, counts_in_window = coverageInInterval( samfile, 
                                                   contig, 
                                                   maxwindow_start, 
                                                   maxwindow_end, 
                                                   shift )

    # counts only in interval - used to define peak center
    counts_in_interval = counts_in_window[offset_left:-offset_right]

    if len(counts_in_interval) == 0:
        raise ValueError( "empty interval: %i - %i for %s:%i-%i" % (offset_left, -offset_right, contig, start, end) )

    #################################################
    # compute peak shape parameters
    peak_nreads = max(counts_in_interval)
    peaks = numpy.array( range(0,interval_width) )[ counts_in_interval >= peak_nreads ]
    if centring_method == "reads":
        peak_center = peaks[len(peaks) // 2] 
    elif centring_method == "middle":
        peak_center = interval_width // 2
    else:
        raise ValueError( "unknown centring method '%s'" % centring_method)
 
    # define peak height
    cdef int peak_height = numpy.ceil( peak_ratio * peak_nreads )
    
    peaks = numpy.array( range(0,interval_width) )[ counts_in_interval >= peak_height ]
    npeaks = len(peaks)
    peak_start = peaks[0]
    peak_end = peaks[-1]
    cdef int peak_width = peak_end - peak_start

    # find size of peak
    # cdef int peak_start, peak_end
    # peak_start, peak_end = peak_center, peak_center

    # for peak_center >= i >= 0:
    #     if counts_in_interval[i] < peak_height:
    #         peak_start = i
    #         break
    
    # for peak_center <= i < interval_width:
    #     if counts_in_interval[i] < peak_height:
    #         peak_end = i
    #         break
    
    # closest and furthest distance of peak to half-height
    cdef int half_height = peak_height // 2
    cdef int left_first, left_last, right_first, right_last 
    left_first, left_last, right_first, right_last = peak_center, 0, interval_width, peak_center

    for i from 0 <= i < peak_center:
        if counts_in_interval[i] >= half_height: 
            left_first = i
            break
        
    for i from peak_center > i >= 0:
        if counts_in_interval[i] <= half_height: 
            left_last = i
            break
        
    for i from peak_center < i < interval_width:
        if counts_in_interval[i] <= half_height: 
            right_first = i
            break

    for i from interval_width > i >= peak_center:
        if counts_in_interval[i] >= half_height: 
            right_last = i
            break

    cdef int furthest_dist = max( peak_center - left_first, right_last - peak_center )
    cdef int closest_dist = min( peak_center - left_last, right_first - peak_center )

    #################################################
    # compute histogram
    nbins = len(bins) - 1
    hist = numpy.zeros( nbins, dtype = numpy.int )

    # decide in which region to count - interval or window
    if only_interval:
        counts = counts_in_interval
        # offset = peak
        offset = peak_center
    else:
        counts = counts_in_window
        # offset = peak to its corresponding location in
        # counts_in_window
        offset = peak_center + offset_left

    cdef int lcounts = len(counts)
        
    xstart = offset + bins[0] 
    for i from 0 <= i < nbins:
        xend = offset + bins[i+1]
        # only take complete bins
        if xstart >= 0 and xend < lcounts:
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
                                     nreads,
                                     numpy.median(counts),
                                     closest_dist, furthest_dist,
                                     bins,
                                     hist ) )

    
    return result

def countAroundPos( Samfile samfile,
                    contig, 
                    int pos,
                    shift,
                    bins ):
    '''count and bin in bins.
    '''
    cdef int xstart, xend, rstart, rend, i, offset

    nreads, counts = coverageInInterval( samfile, 
                                         contig, 
                                         max( 0, pos + bins[0]),
                                         pos + bins[-1],
                                         shift )
    
    nbins = len(bins) - 1
    hist = numpy.zeros( nbins, dtype = numpy.int )
    cdef int lcounts = len(counts)

    offset = -bins[0]
    xstart = offset + bins[0] 
    for i from 0 <= i < nbins:
        xend = offset + bins[i+1]
        # only take complete bins
        if xstart >= 0 and xend < lcounts:
            hist[i] = sum( counts[xstart:xend] ) 
        xstart = xend

    result = PeakShapeCounts._make( (nreads, 
                                     numpy.median(counts),
                                     hist ) )
    
    return result


    
