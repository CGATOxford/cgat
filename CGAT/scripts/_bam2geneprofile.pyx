#cimport csamtools

from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport *
import pysam

import collections, array, struct, sys, itertools
import CGAT.Experiment as E
import CGAT.Intervals as Intervals
import CGAT.GTF as GTF
import CGAT.Stats as Stats
import numpy

CountResult = collections.namedtuple(
    "Counts", "upstream upstream_utr cds downstream_utr downstream")

class RangeCounter:
    
    def __init__(self, countfiles, 
                 controlfiles=None, 
                 control_factor=None, 
                 *args, **kwargs ):

        self.countfiles = countfiles
        self.controlfiles = controlfiles
        self.control_factor = control_factor 
        if self.control_factor is None:
            if self.controlfiles is not None:
                # count number of tags in each file for normalization purposes
                self.norm_fg, self.norm_bg = [], []
                for f in self.countfiles:
                    E.debug("computing read total for %s" % f)
                    self.norm_fg.append(self.getTotal(f))
                for f in self.controlfiles:
                    E.debug("computing read total for %s" % f)
                    self.norm_bg.append(self.getTotal(f))
                fg = sum(self.norm_fg)
                bg = sum(self.norm_bg)
                if bg != 0:
                    self.control_factor = float(fg) / float(bg)
                else:
                    raise ValueError('no background data for normalization')
                E.info("normalization: fg=%i, bg=%f, factor=%f" %
                       (fg, bg, self.control_factor))

    def setup(self, ranges):
        self.counts = numpy.zeros(Intervals.getLength(ranges))
        if self.controlfiles:
            self.counts_bg = numpy.zeros(Intervals.getLength(ranges))

    def getTotal(self):
        raise NotImplementedError(
            'implementations of RangeCounter need to implement getTotal')

    def getCounts(self, contig, ranges, length=0):
        '''count from a set of ranges.

        The ranges are scaled towards the length of the intervals that
        are being counted.
        '''

        self.setup(ranges)
        self.count(self.counts, self.countfiles, contig, ranges)
        if self.controlfiles:
            self.count(self.counts_bg, self.controlfiles, contig, ranges)

        # subsample for length
        if length > 0:
            lnormed = float(length)
            lunnormed = float(len(self.counts))
            if lunnormed == 0: return numpy.zeros(0)
            # E.debug( "normed length = %i, unnormed length =%i " % (lnormed, lunnormed) )
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

            self.counts = self.counts[take]
            if self.controlfiles:
                self.counts_bg = self.counts_bg[take]

        if self.controlfiles:
            # scale background 
            self.counts += 1
            self.counts_bg += 1
            self.counts = numpy.array( self.counts, dtype=float) * \
                self.control_factor / self.counts_bg

        return self.counts

class RangeCounterBAM(RangeCounter):
    '''count densities using bam files.
    '''

    def getTotal(self, bamfile):
        '''return total number of tags in bedfile.'''
        return bamfile.mapped

    def __init__(self, *args, **kwargs ):
        '''
        :param samfiles: list of :term:`bam` formatted files
        :param shifts: list of read shifts
        :param extends: list of read extensions
        :param merge_pairs: merge read pairs for counting
        :param min_insert_size: remove paired reads with insert size below this threshold
        :param max_insert_size: remove paired reads with insert size above this threshold
        '''

        RangeCounter.__init__(self, *args, **kwargs )

    def count(self, counts, files, contig, ranges):

        if len(ranges) == 0: return

        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end
        cdef int interval_width
        cdef int current_offset
        cdef AlignedSegment read
        cdef int shift_extend
        cdef int work_offset
        cdef int pos
        cdef int length

        cdef AlignmentFile samfile

        for samfile in files:

            current_offset = 0

            for start, end in ranges:
                length = end - start

                # skip bamfiles where contig is not present
                if samfile.gettid(contig) < 0:
                    continue

                for read in samfile.fetch(contig, start, end):
                    rstart = max(start, read.pos) - start + current_offset
                    # skip unmapped reads that are assigend a position.
                    if read.aend is None:
                        continue
                    rend = min( end, read.aend) - start + current_offset
                    for i from rstart <= i < rend: counts[i] += 1

                current_offset += length

    def getTotal(self, samfile):
        '''return total number of mapped tags in samfile.'''
        return samfile.mapped


class RangeCounterBAMShift(RangeCounterBAM):
    '''count densities using bam files. 

    Before counting, reads are shifted and extended by a fixed amount.
    '''
    def __init__(self, countfiles, shifts, extends, 
                 controlfiles = None,
                 *args, **kwargs ):
        '''
        :param samfiles: list of :term:`bam` formatted files
        :param shifts: list of read shifts
        :param extends: list of read extensions
        '''

        RangeCounterBAM.__init__(self, countfiles,
                                 controlfiles = controlfiles, *args, **kwargs)

        self.shifts, self.extends = [], []

        for x in xrange(max( (len(shifts), len(extends), len(countfiles)) )):
            if x >= len(shifts):
                shift = 0
            else:
                shift = shifts[0]
            if x >= len(extends):
                extend = 0
            else:
                extend = extends[0]

            if shift > 0 and extend == 0:
                E.warn( "no extension given for shift > 0 - extension will be 2 * shift" )
                extend = shift * 2

            self.shifts.append( shift )
            self.extends.append( extend )

    def count(self, counts, files, contig, ranges ):

        if len(ranges) == 0: return

        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end
        cdef int interval_width
        cdef int current_offset
        cdef AlignedSegment read
        cdef int shift_extend
        cdef int work_offset
        cdef int pos
        cdef int length

        # shifting:
        # forward strand reads:
        #   - shift read upstream by shift
        #   - extend read upstream by extend 
        # reverse strand: 
        #   - shift read downstream by shift
        #   - extend read downstream by extend

        # There is a problem of discontinuity through shifting at
        # exon boundaries
        # 1. Density will accumulate at the ends as reads are shifted
        # up and downstream. However, shifted density will end up
        # in the introns and thus not counted. 
        # 2. The densities along exon boundaries will always be 
        # discontinuous.

        cdef int extend
        cdef int shift
        cdef AlignmentFile samfile

        for samfile, shift, extend in zip(files, self.shifts, self.extends):

            current_offset = 0
            shift_extend = shift + extend

            for start, end in ranges:
                length = end - start
                # collect reads including the regions left/right of interval
                xstart, xend = max(0, start - shift_extend), max(0, end + shift_extend)

                # skip bamfiles where contig is not present
                if samfile.gettid(contig) < 0:
                    continue

                for read in samfile.fetch(contig, xstart, xend):
                    if read.is_reverse: 
                        rstart = read.aend - start - shift_extend
                    else:
                        rstart = read.pos - start + shift

                    rend = min( length, rstart + extend ) + current_offset
                    rstart = max( 0, rstart ) + current_offset
                    for i from rstart <= i < rend: counts[i] += 1

                current_offset += length

class RangeCounterBAMMerge(RangeCounterBAM):
    '''count densities using bam files.

    Before counting, paired end reads are merged and counts
    are computed over the full merged range.

    Reads are not shifted.
    '''
    def __init__(self, countfiles,
                 merge_pairs, 
                 min_insert_size, max_insert_size, 
                 controlfiles = None,
                 *args, **kwargs ):
        '''
        :param samfiles: list of :term:`bam` formatted files
        :param min_insert_size: remove paired reads with insert size below this threshold
        :param max_insert_size: remove paired reads with insert size above this threshold
        '''

        RangeCounterBAM.__init__(self, countfiles,
                                 controlfiles = controlfiles,
                                 *args, **kwargs )
        self.merge_pairs = merge_pairs
        self.min_insert_size = min_insert_size
        self.max_insert_size = max_insert_size

    def count(self, counts, files, contig, ranges ):

        if len(ranges) == 0:
            return

        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end
        cdef int interval_width
        cdef int current_offset
        cdef AlignedSegment read
        cdef int work_offset
        cdef int pos
        cdef int length

        cdef AlignmentFile samfile
        cdef int min_insert_size = self.min_insert_size
        cdef int max_insert_size = self.max_insert_size

        for samfile in files:

            current_offset = 0

            for start, end in ranges:
                length = end - start
                
                xstart, xend = start, end

                # skip bamfiles where contig is not present
                if samfile.gettid(contig) < 0:
                    continue

                for read in samfile.fetch(contig, xstart, xend):
                    flag = read._delegate.core.flag 
                    # remove unmapped reads
                    if flag & 4:
                        continue
                    # remove unpaired
                    if not flag & 2:
                        continue
                    # this is second pair of read - skip to avoid double counting
                    if flag & 128:
                        continue
                    # remove reads on different contigs
                    if read.tid != read.mrnm:
                        continue
                    # remove if insert size too large
                    if (read.isize > max_insert_size) or \
                       (read.isize < min_insert_size):
                        continue
                    if read.pos < read.mpos:
                        rstart = max(start, read.pos)
                        rend = min(end, read.mpos + read.rlen)
                    else:
                        rstart = max(start, read.mpos)
                        rend = min(end, read.pos + read.rlen)

                    rstart += -start + current_offset
                    rend += -start + current_offset
                    for i from rstart <= i < rend:
                        counts[i] += 1

                current_offset += length

class RangeCounterBAMBaseAccuracy(RangeCounterBAM):
    '''count densities using bam files with base accuracy.
    '''
    def __init__(self, 
                 *args, **kwargs):
        '''
        '''

        RangeCounterBAM.__init__(self, *args, **kwargs )

    def count(self, counts, files, contig, ranges):

        if len(ranges) == 0:
            return

        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end
        cdef int interval_width
        cdef int current_offset
        cdef AlignedSegment read
        cdef int work_offset
        cdef int pos
        cdef int length
        
        cdef AlignmentFile samfile

        for samfile in files:

            current_offset = 0

            for start, end in ranges:
                length = end - start
                # skip bamfiles where contig is not present
                if samfile.gettid(contig) < 0:
                    continue

                for read in samfile.fetch(contig, start, end):
                    # note: not optimized for speed yet - uses a python list
                    positions = read.positions
                    for i in positions:
                        if i < start: continue
                        if i >= end: continue
                        try:
                            counts[i - start + current_offset] += 1
                        except IndexError:
                            print len(counts)
                            print "i=", i, "start=",start, "end=", end, "current_offset=", current_offset
                            print "positions=", positions
                            raise

                current_offset += length

class RangeCounterBed(RangeCounter):

    def __init__(self, *args, **kwargs):
        RangeCounter.__init__(self, *args, **kwargs)
                
    def getTotal(self, bedfile):
        '''return total number of tags in bedfile.'''
        cdef int total = 0
        for x in bedfile.fetch():
            total += 1
        return total

    def count(self, counts, files, contig, ranges):
        
        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int xstart, xend, rstart, rend, start, end

        if len(ranges) == 0: return

        cdef int length
        cdef int current_offset
        cdef int interval_width

        for bedfile in files:
            current_offset = 0

            for start, end in ranges:
                length = end - start
                try:
                    for bed in bedfile.fetch(contig, max(0, start), end,
                                             parser = pysam.asBed()):
                        # truncate to range of interest
                        rstart = max(0, bed.start - start) + current_offset
                        rend = min( length, bed.end - start) + current_offset
                        for i from rstart <= i < rend: counts[i] += 1
                except (ValueError, KeyError):
                    # contig not present
                    pass

                current_offset += length

        

class RangeCounterBigWig(RangeCounter):

    def __init__(self, *args, **kwargs ):
        RangeCounter.__init__(self, *args, **kwargs )
        
    def getTotal(self, filename):
        raise NotImplementedError(
            "getTotal not implemented for RangeCounterBigwig")

    def count(self, counts, files, contig, ranges ):
        
        # collect pileup profile in region bounded by start and end.
        cdef int i
        cdef int rstart, rend, start, end, tstart, tend

        if len(ranges) == 0:
            return

        cdef int length
        cdef int current_offset

        for wigfile in files:

            current_offset = 0
            for start, end in ranges:
                length = end - start

                values = wigfile.get(contig, start, end)
                if values == None:
                    continue

                for tstart, tend, value in values:
                    rstart = tstart - start + current_offset
                    rend = tend - start + current_offset
                    for i from rstart <= i < rend: counts[i] = value
                    
                current_offset += length

class IntervalsCounter:
    '''aggregate reads/intervals within a certain arrangement of
    intervals (sections). Sections are defined by a transcript-model 
    and might be configurations such as upstream/exon/downstream, 
    exon/intron, etc.

    Derived class should provide a ``self.count()`` method that accepts
    a list of intervals within a transcript and gene model,
    collects counts within some configuration based on the
    transcript structure and then calls ``self.aggregate()`` to
    aggregate these counts.
    
    The base class provides normalization and book keeping methods
    for the aggregation process.
    '''

    format = "%i"

    def __init__(self, counter, 
                 normalization = None, 
                 outfile_profiles = None,
                 *args, **kwargs ):
        '''counter is an object to provide counts in a list of regions
        scaled to a common size.
        '''

        self.counter = counter

        # normalization method to use before aggregation
        self.normalization = normalization        
        # lables of categories to count (introns, exons, ...)
        self.fields = []
        # number of items counted in each category (1012 introns, 114 exons, ...)
        self.counts = []
        # lengths of items counted in each category
        self.lengths = []
        # number of bins in this category
        self.nbins = []
        # aggregate counts
        self.aggregate_counts = []
        # number of aggregates skipped due to shape mismatch
        # (happens if region extends beyound start/end of contig
        self.nskipped = 0

        # for keeping totals without aggregation
        self.outfile_profiles = outfile_profiles

    def add(self, field, length):
        '''add a new region to counter.

        The region is named *field* and has length *length*.

        '''
        self.counts.append(0)
        self.lengths.append([])
        self.fields.append(field)
        self.nbins.append(length)
        if self.normalization in ("max", "sum", "total-max", "total-sum"):
            self.dtype = numpy.float
        else:
            self.dtype = numpy.int
        self.aggregate_counts.append(numpy.zeros(length, dtype=self.dtype))

    def setOutputProfiles(self, outfile_profiles):
        if outfile_profiles:
            self.outfile_profiles = outfile_profiles
            sum_counts = sum([len(x) for x in self.aggregate_counts])
            self.outfile_profiles.write(
                "\t%s\n" % ("\t".join(map(str, range(sum_counts)))))

    def getNumBins( self ):
        '''returns the number of bins in this counter.'''
        return sum(self.nbins )

    def closeOutputProfiles( self ):
        if self.outfile_profiles:
            self.outfile_profiles.close()

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

    def addName( self, name ):
        '''record name currently being processed.'''
        self.transcript_names.append( name )

    def addLengths( self, *lengths ):
        '''add interval lengths to this counter.'''

        assert len(lengths) == len(self.lengths)
        for idx, interv in enumerate(lengths):
            self.lengths[idx].append(sum([x[1]-x[0] for x in interv]))

    def aggregate(self, *counts):

        # save for outputting
        self.last_counts = counts

        # normalization across all regions
        if self.normalization == "total-max":
            norm = float(max([max(c) for c in counts if len(c) > 0]))
        elif self.normalization == "total-sum":
            norm = float(sum([sum(c) for c in counts] ))
        else:
            norm = 0

        for x, xx in enumerate( zip( self.aggregate_counts, counts)):
            agg, c = xx
            self.counts[x] += 1

            if len(c) == 0:
                continue
            
            if norm > 0:
                # apply uniform normalization
                cc = c / norm
            elif self.normalization == "max":
                # normalize by max count per region
                m = max(c)
                if m == 0: continue
                cc = c / m
            elif self.normalization == "sum":
                # normalize by sum counts per region
                m = sum(c)
                if m == 0: continue
                cc = c / m
            else:
                # no normalization
                cc = c

            try:
                agg += cc.astype(agg.dtype)
            except ValueError:
                self.nskipped += 1

    def getProfile(self,
                   normalize=None,
                   background_region_bins=10):
        '''return meta gene profile.

        apply various normalization methods.

        area
           Normalize by total area under curve of
           the complete profile.

        counts
           Normalize each segment individually by
           the number of regions included.
        
        background
           Apply background normalization to complete
           profile.

        '''
        if normalize is None or normalize == "none":
            profile = numpy.concatenate(self.aggregate_counts)
        
        elif normalize in ("area", "background"):
            profile = numpy.concatenate(self.aggregate_counts)
            profile = numpy.array(profile, dtype=numpy.float)
            
            if normalize == "area":
                background = profile.sum()
                
            elif normalize == "background":
                # take counts at either end of profile
                background_counts = numpy.concatenate( 
                    [profile[:background_region_bins],
                     profile[-background_region_bins:]])
                # compute mean()
                background = numpy.mean(background_counts)

            profile /= background

        elif normalize == "counts":
            profile = numpy.concatenate(
                [numpy.array(x, dtype=numpy.float) / y
                 for x, y in zip(self.aggregate_counts,
                                 self.counts)])

        return profile

    def writeLengthStats(self, outfile):
        '''output length stats to outfile.'''
        
        outfile.write("region\t%s\n" % "\t".join(Stats.Summary.fields))
        for field, l in zip(self.fields, self.lengths):
            outfile.write("%s\t%s\n" % (field, str(Stats.Summary(l))))

    def update(self, gtf):
        
        counted = self.count(gtf)

        if counted and self.outfile_profiles:
            name = gtf[0].transcript_id
            self.outfile_profiles.write("%s\t%s\n" % (name,
                                        "\t".join( [ "\t".join( map(str, x) ) for x in self.last_counts ] ) ) )
        
        return counted

    def __str__(self):
        return "%s=%s" % (self.name, ",".join( [str(sum(x)) for x in self.aggregate_counts]) )

class UnsegmentedCounter(IntervalsCounter):
    '''clone a segmented counter as a single
    unsegmented counter.
    '''
    
    name = None
    
    def __init__(self, counter,
                 *args,
                 **kwargs ):
        IntervalsCounter.__init__(self, counter, *args, **kwargs)
        self.name = counter.name
        self.add( 'all', counter.getNumBins() )


class SeparateExonCounter(IntervalsCounter):
    '''count reads over 1st, middle and last exons as well as upstream/
    downstream of genes/transcripts.  Only protein-coding genes (with a
    CDS) are counted
    '''
    name = "separateexonprofile"

    def __init__(self, counter, 
                 int resolution_upstream,
                 int resolution_first_exon,
                 int resolution_last_exon,
                 int resolution_exons,
                 int resolution_downstream,
                 int extension_upstream = 0,
                 int extension_downstream = 0,
                 int scale_flanks = 0,
                 *args,
                 **kwargs):

        IntervalsCounter.__init__(self, counter, *args, **kwargs)
	      
        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream
        self.resolution_first_exon = resolution_first_exon
        self.resolution_last_exon = resolution_last_exon
        self.resolution_exons = resolution_exons
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        self.scale_flanks = scale_flanks

        for field, length in zip(
                ("upstream", "first_exon", "exons", "last_exon", "downstream"),
                (resolution_upstream,
                 resolution_first_exon,
                 resolution_exons,
                 resolution_last_exon,
                 resolution_downstream)):
            self.add(field, length)


    def count(self, gtf):
        '''build ranges to be analyzed from a gene model.
        Returns a tuple with ranges for first, middle and last exons,
        upstream and downstream.  Does not count genes with 1 exon.
        '''
        
        contig = gtf[0].contig
        exons = GTF.asRanges(gtf, "exon")
        if len(exons) == 0:
            E.warn("no exons in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        if len(exons) == 2:
            first = [exons[0]]
            last = [exons[-1]]
            exons = [(0, 0)]

        elif len(exons) == 1:
            E.warn("only one exon in %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        else:
            first = [exons[0]]
            last = [exons[-1]]
            exons.remove(exons[0])
            exons.remove(exons[-1])

        first_start, first_end = first[0][0], first[0][1]
        last_start, last_end = last[0][0], last[0][1]
        exons_start, exons_end = exons[0][0], exons[-1][1]
        if self.scale_flanks > 0:
            self.extension_downstream = (last_end - first_start)*self.scale_flanks
            self.extension_upstream = (last_end - first_start)*self.scale_flanks
            E.dubg("scale flanks")
            
        if gtf[0].strand == "-":
            downstream = [(max(0, first_start - self.extension_downstream),
                           first_start),]
            upstream = [(last_end,
                         last_end + self.extension_upstream), ]
        else:
            upstream = [(max(0, first_start - self.extension_upstream),
	                 first_start), ]
            downstream = [(last_end,
	                   last_end + self.extension_downstream), ]

        E.debug("counting first exons")
        self.counts_first_exon = self.counter.getCounts(contig, first,
                                                        self.resolution_first_exon)

        E.debug("counting exons")
        self.counts_exons = self.counter.getCounts(contig,exons,
                                                   self.resolution_exons)

        E.debug("counting last exons")
        self.counts_last_exon = self.counter.getCounts(contig, last,
                                                       self.resolution_last_exon)

        E.debug("counting upstream")
        self.counts_upstream = self.counter.getCounts(contig, upstream,
                                                      self.resolution_upstream)

        E.debug("counting downstream")
        self.counts_downstream = self.counter.getCounts(contig, downstream,
                                                        self.resolution_downstream)
        
        E.debug("counting finished" )
        
        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_first_exon = self.counts_first_exon[::-1]
            self.counts_last_exon = self.counts_last_exon[::-1]
            self.counts_exons = self.counts_exons[::-1]
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]

        self.addLengths( upstream, first, exons, last, downstream )

        self.aggregate(self.counts_upstream,
                       self.counts_first_exon,
                       self.counts_exons,
                       self.counts_last_exon,
                       self.counts_downstream )

        return 1


class SeparateExonWithIntronCounter(IntervalsCounter):
    '''count reads over 1st, middle and last exons as well as introns, upstream/
    downstream of genes/transcripts.  Only protein-coding genes (with a
    CDS) are counted, introns are not split.
    '''
    name = "separateexonprofilewithintrons"

    def __init__(self, counter, 
                 int resolution_upstream,
                 int resolution_first_exon,
                 int resolution_last_exon,
                 int resolution_exons,
                 int resolution_introns,
                 int resolution_downstream,
                 int extension_upstream = 0,
                 int extension_downstream = 0,
                 int scale_flanks = 0,
                 *args,
                 **kwargs):

        IntervalsCounter.__init__(self, counter, *args, **kwargs)
	      
        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream
        self.resolution_first_exon = resolution_first_exon
        self.resolution_last_exon = resolution_last_exon
        self.resolution_exons = resolution_exons
        self.resolution_introns = resolution_introns
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        self.scale_flanks = scale_flanks

        for field, length in zip(
                ("upstream", "first_exon", "exons",
                 "introns", "last_exon", "downstream"),
                (resolution_upstream,
                 resolution_first_exon,
                 resolution_exons,
                 resolution_introns,
                 resolution_last_exon,
                 resolution_downstream)):
            self.add(field, length)

    def count(self, gtf):
        '''build ranges to be analyzed from a gene model.
        Returns a tuple with ranges for first, middle and last exons, introns,
        upstream and downstream.  Does not count genes with 1 exon.
        '''
        
        contig = gtf[0].contig
        exons = GTF.asRanges(gtf, "exon")
        introns = Intervals.complement(exons)

        if len(exons) == 0:
            E.warn("no exons in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        # Exclude intron-less genes
        if len(introns) == 0:
            E.warn("no introns in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0
	  
        if len(exons) == 2:
            first = [exons[0]]
            last = [exons[-1]]
            exons = [(0, 0)]

        elif len(exons) == 1:
            E.warn("only one exon in %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        else:
            first = [exons[0]]
            last = [exons[-1]]
            exons.remove(exons[0])
            exons.remove(exons[-1])

        first_start, first_end = first[0][0], first[0][1]
        last_start, last_end = last[0][0], last[0][1]
        exons_start, exons_end = exons[0][0], exons[-1][1]
        if self.scale_flanks > 0:
            self.extension_downstream = (last_end - first_start)*self.scale_flanks
            self.extension_upstream = (last_end - first_start)*self.scale_flanks
            E.dubg("scale flanks")
            
        if gtf[0].strand == "-":
            downstream = [(max(0, first_start - self.extension_downstream),
                           first_start),]
            upstream = [(last_end,
                         last_end + self.extension_upstream), ]
        else:
            upstream = [(max(0, first_start - self.extension_upstream),
	                 first_start), ]
            downstream = [(last_end,
	                   last_end + self.extension_downstream), ]

        E.debug("counting first exons")
        self.counts_first_exon = self.counter.getCounts(contig, first,
                                                        self.resolution_first_exon)

        E.debug("counting exons")
        self.counts_exons = self.counter.getCounts(contig, exons,
                                                   self.resolution_exons)

        E.debug("counting introns")
        self.counts_introns = self.counter.getCounts(contig, introns,
                                                     self.resolution_introns)

        E.debug("counting last exons")
        self.counts_last_exon = self.counter.getCounts(contig, last,
                                                       self.resolution_last_exon)

        E.debug("counting upstream")
        self.counts_upstream = self.counter.getCounts(contig, upstream,
                                                      self.resolution_upstream)

        E.debug("counting downstream")
        self.counts_downstream = self.counter.getCounts(contig, downstream,
                                                        self.resolution_downstream)
        
        E.debug("counting finished" )
        
        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_first_exon = self.counts_first_exon[::-1]
            self.counts_last_exon = self.counts_last_exon[::-1]
            self.counts_exons = self.counts_exons[::-1]
            self.counts_introns = self.counts_introns[::-1]
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]

        self.addLengths(upstream, first, exons, introns, last, downstream)

        self.aggregate(self.counts_upstream,
                       self.counts_first_exon,
                       self.counts_exons,
                       self.counts_introns,
                       self.counts_last_exon,
                       self.counts_downstream )

        return 1


class GeneCounter( IntervalsCounter ):
    '''count reads in exons and upstream/downstream of genes/transcripts.
    
    Only protein coding transcripts are counted (i.e. those with a CDS)
    '''
    
    name = "geneprofile"
    
    def __init__(self, counter, 
                 int resolution_upstream,
                 int resolution_exons,
                 int resolution_downstream,
                 int extension_upstream = 0, 
                 int extension_downstream = 0,
                 int scale_flanks = 0,
                 *args,
                 **kwargs ):

        IntervalsCounter.__init__(self, counter, *args, **kwargs )

        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream 
        self.resolution_exons = resolution_exons
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        self.scale_flanks = scale_flanks

        for field, length in zip( 
            ("upstream", "exons", "downstream"),
            (resolution_upstream,
             resolution_exons,
             resolution_downstream ) ):
            self.add( field, length )
        
    def count( self, gtf ):
        '''build ranges to be analyzed from a gene model.

        Returns a tuple with ranges for exons, upstream, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        if len(exons) == 0:
            E.warn( "no exons in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        exon_start, exon_end = exons[0][0], exons[-1][1]
        if self.scale_flanks > 0:
            self.extension_downstream = (exon_end - exon_start)*self.scale_flanks
            self.extension_upstream = (exon_end - exon_start)*self.scale_flanks
            E.debug("scale flanks")

        if gtf[0].strand == "-":
            downstream = [( max(0, exon_start - self.extension_downstream),
                            exon_start ),] 
            upstream = [( exon_end,
                          exon_end + self.extension_upstream),]
        else:
            upstream = [(max(0, exon_start - self.extension_upstream),
                         exon_start ),] 
            downstream = [(exon_end,
                           exon_end + self.extension_downstream ),]

        E.debug("counting exons" )
        self.counts_exons = self.counter.getCounts(contig,
                                                   exons,
                                                   self.resolution_exons)
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts(contig, 
                                                      upstream,
                                                      self.resolution_upstream) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts(contig,
                                                        downstream,
                                                        self.resolution_downstream)

        E.debug("counting finished" )

        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_exons = self.counts_exons[::-1]
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]
            
        self.addLengths( upstream, exons, downstream )

        self.aggregate( self.counts_upstream,
                        self.counts_exons,
                        self.counts_downstream )

        return 1


class GeneCounterWithIntrons( IntervalsCounter ):
    '''count reads in exons, in introns and the upstream/downstream of genes/transcripts.
    
    Assumes all the behavior of GeneCounter class, because the code is modified based on 
    GeneCounter Class. Class created by Tim on 16 May. 2013.
    
    Both  protein coding transcripts, and non coding transcripts are counted.
    (i.e. both those with a CDS, and those without a CDS are counted)
    '''
    
    name = "geneprofilewithintrons"
    
    def __init__(self, 
                 counter,
                 int resolution_upstream,
                 int resolution_exons,
                 int resolution_introns,  #Tim 16 May. 2013. for introns
                 int resolution_downstream,
                 int extension_upstream = 0, 
                 int extension_downstream = 0,
                 int scale_flanks = 0,
                 *args,
                 **kwargs ):

        IntervalsCounter.__init__(self, counter, *args, **kwargs )

        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream 
        self.resolution_exons = resolution_exons
        self.resolution_introns = resolution_introns   #Tim 16 May. 2013. for introns
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        self.scale_flanks = scale_flanks

        for field, length in zip( 
            ("upstream", "exons", "introns", "downstream"),  #Tim 16 May. 2013. for introns
            (resolution_upstream,
             resolution_exons,
             resolution_introns,  #Tim 16 May. 2013. for introns
             resolution_downstream ) ):
            self.add( field, length )
        
    def count( self, gtf ):
        '''build ranges to be analyzed from a gene model.

        Returns a tuple with ranges for upstream, exons, introns, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        introns = Intervals.complement(exons) #Tim 16 May. 2013. for introns
        
        if len(exons) == 0:
            E.warn( "no exons in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0
            
        # Tim 16 May. 2013. for introns. 
        # Genes have no introns are omitted for counting in all other regions. Because 
        # include them will distort the relative height of exons/introns in the profile.
        # In most gene profile plots, the relative height of different regions 
        # is very important info, thus omit the gene will ensure the correctness for the 
        # relative height across different regions.
        # Particularly true when user use normalization method: total-max.
        if len(introns) == 0:   
            E.debug( "Omitted this gene entirely because there are no introns in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0

        exon_start, exon_end = exons[0][0], exons[-1][1]
        if self.scale_flanks > 0:
            self.extension_downstream = (exon_end - exon_start)*self.scale_flanks
            self.extension_upstream = (exon_end - exon_start)*self.scale_flanks
            E.debug("scale flanks")

        if gtf[0].strand == "-":
            downstream = [ ( max(0, exon_start - self.extension_downstream), exon_start ), ] 
            upstream = [ ( exon_end, exon_end + self.extension_upstream ), ]
        else:
            upstream = [ ( max(0, exon_start - self.extension_upstream), exon_start ), ] 
            downstream = [ ( exon_end, exon_end + self.extension_downstream ), ]

        E.debug("counting exons" )
        self.counts_exons = self.counter.getCounts( contig, exons, self.resolution_exons  )
        E.debug("counting introns" ) #Tim 16 May. 2013. for introns. 
        self.counts_introns = self.counter.getCounts( contig, introns, self.resolution_introns  )        
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts( contig, upstream, self.resolution_upstream ) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts( contig, downstream, self.resolution_downstream )

        E.debug("counting finished" )

        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_exons = self.counts_exons[::-1]
            self.counts_introns = self.counts_introns[::-1]   #Tim 16 May. 2013. for introns. 
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]
            
        self.addLengths( upstream, exons, introns, downstream )   #Tim 16 May. 2013. for introns. 

        self.aggregate( self.counts_upstream,
                        self.counts_exons,
                        self.counts_introns,  #Tim 16 May. 2013. for introns. 
                        self.counts_downstream )

        return 1



class GeneCounterAbsoluteDistanceFromThreePrimeEnd( IntervalsCounter ):
    '''count reads in exons of genes/transcripts from the three prime
    polyA tail, without shrinking or lengthening the gene as it would
    in other genecounter mode in this script.  

    In other words, only count the reads that fall on the exons of the
    (virtual maximal) mRNA transcript.  Note that the distance is
    relative to TTS (three prime polyA tail) on the mNRA transcript,
    insteads of on the genomic assembly is used for the counting and
    plotting.
    
    For mRNA with multiple exons, the exons are first stiched together
    into one piece of mRNA (the virtual maximal transcript for each
    gene). Subsequently, the mRNA is being used for counting. This is
    the (only) proper way to avoid counting in introns, which screw up
    the actual genebody coverage profile. And this only works with
    --base-accuracy option, because otherwise spliced reads will be
    considered as covering the entire intron, or spliced out introns
    will be counted as covered by reads.  (this option imply the
    --base-accuracy option).
    
    Also count reads in introns (in exactly the same manner as if they
    are exons described above ) of genes/transcripts from the three
    prime polyA tail.
    
    Note:

    * Both protein coding transcripts and non coding
    transcripts are counted.  i.e. both those with a CDS, and those
    without a CDS are counted.
    
    * Only genes/transcripts with total length longer than the
    --extension-exons-absolute-distance-topolya is counted (i.e. those
    genes/transcripts shorter than this is omitted because they will
    only contribute to the last len(gene) base pair of the graph, thus
    end-up ploting a graph appears to have three prime bias even when
    there is no three prime bias at all.
    
    There is one alternative way of dealing with this issue: is to
    count the number of times a nucleotide is being visited, and
    normalize against the counts curve with the visits counting
    curve. However, under current architecture of the counter, it
    seems hard to support this method.  Assumeing you get enough
    number of genes longer than 4000b (in hg19) to plot the three
    prime bias to satisfactory quality, it seems not really worth the
    effort to implement the alternative way,

    '''
    
    name = "geneprofileabsolutedistancefromthreeprimeend"   #Tim 31th Aug 2013
    
    def __init__(self, 
                 counter,
                 int resolution_upstream,
                 int resolution_downstream,
                 int resolution_exons_absolute_distance_topolya,
                 int resolution_introns_absolute_distance_topolya,
                 # int resolution_exons_absolute_distance_tostartsite,   # Tim 31th Aug 2013: 
                                                                         # a possible feature for future,  if five prime bias is of your interest. 
                                                                         #(you need to create another class). It is not very difficult to derive from this class, but is not implemented yet
                                                                         # This future feature is slightly different the TSS profile already implemented, because in this future feature introns are skipped, 
                 # int resolution_introns_absolute_distance_tostartsite, # feature not implemented yet
                 
                 int extension_upstream = 0, 
                 int extension_downstream = 0,
                 int extension_exons_absolute_distance_topolya = 0,
                 int extension_introns_absolute_distance_topolya = 0,
                 # int extension_exons_absolute_distance_tostartsite = 0,   # feature not implemented yet
                 # int extension_introns_absolute_distance_tostartsite = 0, # feature not implemented yet
                 
                 int scale_flanks = 0,
                 *args,
                 **kwargs ):

        IntervalsCounter.__init__(self, counter, *args, **kwargs )

        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream 
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream
        
        self.resolution_exons_absolute_distance_topolya = resolution_exons_absolute_distance_topolya     #Tim 31th Aug 2013
        self.resolution_introns_absolute_distance_topolya = resolution_introns_absolute_distance_topolya #Tim 31th Aug 2013
        self.extension_exons_absolute_distance_topolya = extension_exons_absolute_distance_topolya       #Tim 31th Aug 2013
        self.extension_introns_absolute_distance_topolya = extension_introns_absolute_distance_topolya   #Tim 31th Aug 2013
        
        self.scale_flanks = scale_flanks

        for field, length in zip(
                ("upstream%dbp_zoomedTo%dbp"%(extension_upstream, resolution_upstream),
                 "exonsLast%dbp_zoomedTo%dbp"%(extension_exons_absolute_distance_topolya, resolution_exons_absolute_distance_topolya) ,
                 "intronsLast%dbp_zoomedTo%dbp"%(extension_introns_absolute_distance_topolya, resolution_introns_absolute_distance_topolya) ,
                 "downstream%dbp_zoomedTo%dbp"%(extension_downstream, resolution_downstream),
             ),
                (resolution_upstream,
                 resolution_exons_absolute_distance_topolya,
                 resolution_introns_absolute_distance_topolya,
                 resolution_downstream) 
        ):
            self.add(field, length)

            
    #Tim 31th Aug 2013: Important function in this class to stich together the exons into one complete of mRNA
    def __chopAllTranscriptToAFixedLengthForTheThreePrimeBiasCounting(self, gtf , exons, extension_threePrimeBiasNotZoomed ):
        '''
        works for both exons and introns, the variable name exons is only symolic here,
        it means either exons, or introns depending on the regions you passed in 
        during function call
        
        in order to use the gene2profile code to plot 3prime bias plot, 
        if a gene is longer than my plotting size (e.g. 5kb), because I 
        do not want to let the code to shrink the gene, I need to create 
        a "chopped gene" that only has the first several bp (e.g. 5kb) 
        closest to the 3prime, so the gene always appears to be of the length 
        "extension_threePrimeBiasNotZoomed", (e.g. 5kb)
        
        In this version, shorter genes are omitted at count() becase len([])==0.
        It is oimitted in the same mannar that was used to
        omit the genes without any exons, by returnning [], which can be interpreted
        by the calling function as omit, or take any other further actions.
       
        '''
        
        if len(exons)==0: return []  # this is to cater for the case where there's 
                                     # one exon, and no intron, but I am plotted\
                                     # unzoomed introns as well
        
        leftoverBasepairAllowance = extension_threePrimeBiasNotZoomed
        exons_chopped = []
        if gtf[0].strand == "-":
            for exon in exons:
                exonSize = exon[1]-exon[0]
                if leftoverBasepairAllowance > exonSize:
                    leftoverBasepairAllowance -= exonSize
                    exons_chopped.append(exon)            # append is the correct way, coz we process the first exon first. 
                                                            # so subsequent exons shall be appended so that the basepair ordering 
                                                            # is in reverse ordering alone the gene body, likewise, the exon 
                                                            # ordering is also the reverse ordering alone the gene body. In the 
                                                            # end, just before counts are aggregated, there will be code:
                                                            # if "-": counts_exons_NotZoomedAndChopped_ForThreePrimeBias[::-1] to flip 
                                                            # everything in base pair resolution.
                    E.debug(str(exons_chopped))
                else:
                    E.debug("Exon intervals preparation function: gene %s:%s, leftoverBasepairAllowance= %s " % (gtf[0].transcript_id, gtf[0].strand, leftoverBasepairAllowance)  )
                    E.debug("Exon intervals preparation function: The potential half exon: "+str(exon) )
                    exon_chopped = (exon[0], exon[0]+leftoverBasepairAllowance)
                    exons_chopped.append( exon_chopped )    # append is the correct way, same as above
                    exon_choppedSize = exon_chopped[1]-exon_chopped[0]
                    leftoverBasepairAllowance -= exon_choppedSize
                    E.debug( "Exon intervals preparation function: gene %s:%s, leftoverBasepairAllowance= %s " % (gtf[0].transcript_id, gtf[0].strand, leftoverBasepairAllowance)  )
                    E.debug( str(exons_chopped) )
                    break
                    # gene is longer than my window, thus chopped, 
                    # once in this section of code, exons_chopped is the final chopped gene, thus break from loop
                    # and later, exons_chopped is being return at the end of the function.

        elif gtf[0].strand == "+":
            for exon in exons[::-1]:
                exonSize = exon[1]-exon[0]
                if leftoverBasepairAllowance > exonSize:
                    leftoverBasepairAllowance -= exonSize
                    exons_chopped.insert(0, exon)#prepend is the correct way, coz we process the last exon first. so subsequent exons shall be prepended so that the basepair ordering is the same as exon ordering.
                    E.debug( str(exons_chopped) )
                else:
                    E.debug( "Exon intervals preparation function: gene %s:%s, leftoverBasepairAllowance= %s " % (gtf[0].transcript_id, gtf[0].strand, leftoverBasepairAllowance)  )
                    E.debug( "Exon intervals preparation function: The potential half exon: "+str(exon) )
                    exon_chopped = (exon[1]-leftoverBasepairAllowance, exon[1])
                    exons_chopped.insert(0,exon_chopped)#prepend is the correct way, same as above
                    exon_choppedSize = exon_chopped[1]-exon_chopped[0]
                    leftoverBasepairAllowance -= exon_choppedSize
                    E.debug( "Exon intervals preparation function: gene %s:%s, leftoverBasepairAllowance= %s " % (gtf[0].transcript_id, gtf[0].strand, leftoverBasepairAllowance)  )
                    E.debug( str(exons_chopped) )
                    break
                    # gene is longer than my window, thus chopped, 
                    # once in this section of code, exons_chopped is the final chopped gene, thus break from loop
                    # and later, exons_chopped is being return at the end of the function.
        else:
            E.error( "In exon intervals preparation function: strand is neither + nor - in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id) )
            # in this situation, gene direction column in the GTF file have some problems, 
            # exons_chopped is still [], and is being return 
            # at the end of the function. 
            # so that his componet is skipped by the counter
            
        #gene is nicely chopped, so that no alloance is left.
        if leftoverBasepairAllowance==0: 
            return exons_chopped
            
        # gene is shorter than my window, thus it is not chopped to any shorter at all.
        # you can return three different signals to the caller function, and they have different meanings:
        # if you return the gene, it will be expanded and counted, thus might not be your desired behaviour.
        # if you return None, it will be throwing an exception in count()
        # if you return [], thus len([])==0, the whole gene will be quitely omitted. 
        elif leftoverBasepairAllowance>0: 
            E.debug( "In exon intervals preparation function: gene %s:%s is shorter than the defined threeprimebiaswindow: %s ." % (gtf[0].gene_id, gtf[0].transcript_id, extension_threePrimeBiasNotZoomed) )            
            return []
        
        #leftoverBasepairAllowance is neither ==0 nor >0, there must be a bug in the code.
        else: 
            E.error( "In exon intervals preparation function: leftoverBasepairAllowance after chopping is neither ==0 nor >0 in gene %s:%s. It must be a bug" % (gtf[0].gene_id, gtf[0].transcript_id) )
            E.error( "In exon intervals preparation function: details:: before: s%, after: s%, strand: s%."%(str(exons), str(exons_chopped), gtf[0].strand)    )  
            
        #the __chopgenestoshowthreePrimeBiasNotZoomed subfunction ends here # this return below is never be used according to the logic above. i.e. unless it enters the else: warn due to a bug, this return statement is never executed.
        return
        
            
    def count( self, gtf ):
        '''build ranges to be analyzed from a gene model.

        Returns a tuple with ranges for upstream, exons, introns, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        # skip genes without exons
        if len(exons) == 0: return 0

        introns = Intervals.complement(exons) #Tim 31th Aug 2013
        
        if len(exons) == 0:
            E.warn( "In counter function: seems like a seriouse problem, there are no exons in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            return 0
            
        if len(introns) == 0:
            E.debug( "In counter function: no introns in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            # Tim 31th Aug 2013 :
            # do not warn here because there could be many single exon genes.
            # return 0 here because if no intron, we do not want to count this whole gene.
            # otherwise, the relative height among intron, exon, upstream, downstream will be wrong
            return 0
            
        exon_start, exon_end = exons[0][0], exons[-1][1]
        if self.scale_flanks > 0:
            self.extension_downstream = (exon_end - exon_start)*self.scale_flanks
            self.extension_upstream = (exon_end - exon_start)*self.scale_flanks
            E.debug("In counter function: scale flanks")

        if gtf[0].strand == "-":
            downstream = [ ( max(0, exon_start - self.extension_downstream), exon_start ), ] 
            upstream = [ ( exon_end, exon_end + self.extension_upstream ), ]
        else:
            upstream = [ ( max(0, exon_start - self.extension_upstream), exon_start ), ] 
            downstream = [ ( exon_end, exon_end + self.extension_downstream ), ]
        
        # within the __chopAllTranscriptToAFixedLengthForTheThreePrimeBiasCounting function, the flipping for genes on "-" strand is already taken care of. Tim 31th Aug 2013
        exons_NotZoomedAndChopped_ForThreePrimeBias   = self.__chopAllTranscriptToAFixedLengthForTheThreePrimeBiasCounting( gtf, exons,   self.extension_exons_absolute_distance_topolya   )
        introns_NotZoomedAndChopped_ForThreePrimeBias = self.__chopAllTranscriptToAFixedLengthForTheThreePrimeBiasCounting( gtf, introns, self.extension_introns_absolute_distance_topolya )
        if len(exons_NotZoomedAndChopped_ForThreePrimeBias) == 0:
            E.debug( "In counter function: after exons concat together, it is not long enough in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            # Tim 31th Aug 2013 :
            # do not warn here because there are many genes's total exon length might be short.
            # return 0 here because if no exons, we do not want to count this whole gene.
            # otherwise, the relative height of intron and exon will be wrong
            return 0        
        if len(introns_NotZoomedAndChopped_ForThreePrimeBias) == 0:
            E.debug( "In counter function: after introns concat together, it is not long enough in gene %s:%s" % (gtf[0].gene_id, gtf[0].transcript_id))
            # Tim 31th Aug 2013  :
            # do not warn here because there are many genes's total exon length might be short.
            # return 0 here because if no intron, we do not want to count this whole gene.
            # otherwise, the relative height of intron and exon will be wrong
            return 0
                        
            
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts( contig, upstream, self.resolution_upstream ) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts( contig, downstream, self.resolution_downstream )
        E.debug("counting exons NotZoomedAndChopped ForThreePrimeBias" )
        self.counts_exons_NotZoomedAndChopped_ForThreePrimeBias = self.counter.getCounts( contig, exons_NotZoomedAndChopped_ForThreePrimeBias, self.resolution_exons_absolute_distance_topolya )
        E.debug("counting introns NotZoomedAndChopped ForThreePrimeBias" )
        self.counts_introns_NotZoomedAndChopped_ForThreePrimeBias = self.counter.getCounts( contig, introns_NotZoomedAndChopped_ForThreePrimeBias, self.resolution_introns_absolute_distance_topolya )
        E.debug("counting finished" )

        
        if gtf[0].strand == "+":
            # nothings needs to be flipped for positive strand gene 
            pass 
        elif gtf[0].strand == "-":
            # flip for negative strand genes
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]
            #Tim 31th Aug 2013 
            self.counts_exons_NotZoomedAndChopped_ForThreePrimeBias   = self.counts_exons_NotZoomedAndChopped_ForThreePrimeBias[::-1]
            self.counts_introns_NotZoomedAndChopped_ForThreePrimeBias = self.counts_introns_NotZoomedAndChopped_ForThreePrimeBias[::-1]
          
        self.addLengths( upstream, 
                         exons_NotZoomedAndChopped_ForThreePrimeBias, 
                         introns_NotZoomedAndChopped_ForThreePrimeBias, 
                         downstream )                   #Tim 31th Aug 2013 
                
        self.aggregate( self.counts_upstream,
                        self.counts_exons_NotZoomedAndChopped_ForThreePrimeBias,
                        self.counts_introns_NotZoomedAndChopped_ForThreePrimeBias,                        
                        self.counts_downstream )        #Tim 31th Aug 2013 

        return 1

class UTRCounter( IntervalsCounter ):
    '''counts reads in 3'UTR and 5'UTR in addition
    to the extension outside.'''

    name = "utrprofile"
    
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

        IntervalsCounter.__init__(self, counter, *args, **kwargs )

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

        Returns a tuple with ranges for cds, upstream_utr,
        downstream_utr, upstream, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        # skip genes without exons
        if len(exons) == 0: return 0
        exon_start, exon_end = exons[0][0], exons[-1][1]
        self.cds = GTF.asRanges( gtf, "CDS" )
        # skip genes without CDS
        if len(self.cds) == 0: return 0

        cds_start, cds_end = self.cds[0][0], self.cds[-1][1]
        utrs = Intervals.truncate( exons, self.cds )

        if gtf[0].strand == "-":
            self.upstream_utr = [ x for x in utrs if x[0] >= cds_end ]
            self.downstream_utr = [ x for x in utrs if x[1] <= cds_start ]
            self.downstream = [ ( max(0, exon_start - self.extension_downstream), exon_start ), ] 
            self.upstream = [ ( exon_end, exon_end + self.extension_upstream ), ]
        else:
            self.upstream_utr = [ x for x in utrs if x[1] <= cds_start ]
            self.downstream_utr = [ x for x in utrs if x[0] >= cds_end ]
            self.upstream = [ ( max(0, exon_start - self.extension_upstream), exon_start ), ] 
            self.downstream = [ ( exon_end, exon_end + self.extension_downstream ), ]
        
        E.debug("counting cds" )
        self.counts_cds = self.counter.getCounts( contig, self.cds, self.resolution_cds  )
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts( contig, self.upstream, self.resolution_upstream ) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts( contig, self.downstream, self.resolution_downstream )

        if self.upstream_utr:
            E.debug("counting upstream_utr" )
            self.counts_upstream_utr = self.counter.getCounts( contig, self.upstream_utr, self.resolution_upstream_utr )
        else:
            # add pseudocounts for those genes without a UTR.
            self.counts_upstream_utr = numpy.zeros( self.resolution_upstream_utr )
            self.counts_upstream_utr += (self.counts_upstream[-1] + self.counts_cds[0]) // 2

        if self.downstream_utr:
            E.debug("counting downstream_utr" )
            self.counts_downstream_utr = self.counter.getCounts( contig, self.downstream_utr, self.resolution_downstream_utr )
        else:
            # add a pseudocount for those genes without a UTR.
            self.counts_downstream_utr = numpy.zeros( self.resolution_downstream_utr )
            self.counts_downstream_utr += (self.counts_downstream[0] + self.counts_cds[-1]) // 2

        E.debug("counting finished" )

        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_cds = self.counts_cds[::-1]
            self.counts_upstream_utr = self.counts_upstream_utr[::-1]
            self.counts_downstream_utr = self.counts_downstream_utr[::-1]
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]

        self.aggregate( self.counts_upstream,
                        self.counts_upstream_utr,
                        self.counts_cds,
                        self.counts_downstream_utr,
                        self.counts_downstream )

        return 1

class RegionCounter( GeneCounter ):
    '''count in complete region given by gene/transcript.'''

    name = "intervalprofile"
    
    def __init__(self, 
                 *args,
                 **kwargs ):

        GeneCounter.__init__(self, *args, **kwargs )

        # substitute field name 'exons' with 'interval'
        assert self.fields[1] == "exons"
        self.fields[1] = "interval"

    def count( self, gtf ):
        '''count intervals.'''

        g = GTF.Entry()
        g.copy( gtf[0] )
        g.start = gtf[0].start
        g.end = gtf[-1].end
        
        return GeneCounter.count( self, [ g ] )

class MidpointCounter( GeneCounter ):
    '''count in complete region given by gene/transcript.'''

    name = "midpointprofile"
    
    def __init__(self, counter, 
                 int resolution_upstream,
                 int resolution_downstream,
                 int extension_upstream = 0, 
                 int extension_downstream = 0,
                 int scale_flanks = 0,
                 *args,
                 **kwargs ):

        IntervalsCounter.__init__(self, counter, *args, **kwargs )

        self.extension_upstream = extension_upstream
        self.extension_downstream = extension_downstream 
        self.resolution_upstream = resolution_upstream
        self.resolution_downstream = resolution_downstream

        for field, length in zip( 
            ("upstream", "downstream"),
            (resolution_upstream,
             resolution_downstream ) ):
            self.add( field, length )
        
    def count( self, gtf ):
        '''build ranges to be analyzed from a gene model.

        Returns a tuple with ranges for exons, upstream, downstream.
        '''

        contig = gtf[0].contig 
        exons = GTF.asRanges( gtf, "exon" )
        # skip genes without exons
        if len(exons) == 0: return 0
        exon_start, exon_end = exons[0][0], exons[-1][1]
        midpoint = (exon_end - exon_start) // 2 + exon_start

        if gtf[0].strand == "-":
            downstream = [ ( max(0, midpoint - self.extension_downstream), midpoint ), ] 
            upstream = [ ( midpoint, midpoint + self.extension_upstream ), ]
        else:
            upstream = [ ( max(0, midpoint - self.extension_upstream), midpoint ), ] 
            downstream = [ ( midpoint, midpoint + self.extension_downstream ), ]
        
        E.debug("counting upstream" )
        self.counts_upstream = self.counter.getCounts( contig, upstream, self.resolution_upstream ) 
        E.debug("counting downstream" )
        self.counts_downstream = self.counter.getCounts( contig, downstream, self.resolution_downstream )

        E.debug("counting finished" )

        ## revert for negative strand
        if gtf[0].strand == "-":
            self.counts_upstream = self.counts_upstream[::-1]
            self.counts_downstream = self.counts_downstream[::-1]
            
        self.addLengths( upstream, downstream )

        self.aggregate( self.counts_upstream,
                        self.counts_downstream )

        return 1
        
class TSSCounter( IntervalsCounter ):
    '''count profile at transcription start/end site.

    The transcription start site is the beginning of the first exon.
    The transcription termination site is the end of the last exon.

    '''

    name = "tssprofile"
    
    def __init__(self, counter, extension_out = 0, extension_in = 0, *args, **kwargs):
        IntervalsCounter.__init__(self, counter, *args, **kwargs )

        self.extension_out = extension_out
        self.extension_in = extension_in 

        for field, length in zip( ("tss", "tts"),
                                  ( extension_out + extension_in,
                                    extension_out + extension_in ) ):
            self.add( field, length )
        
    def count( self, gtf ):

        contig, strand = gtf[0].contig, gtf[0].strand
        exons = GTF.asRanges( gtf, "exon" )
        # skip genes without exons
        if len(exons) == 0: return 0

        self.tss, self.tts = exons[0][0], exons[-1][1]

        # no max(0, ...) here as these ranges need to have always the same length
        if strand == "-":
            self.tts_ranges = [ (max(0, self.tss - self.extension_out), 
                                 self.tss + self.extension_in), ]
            self.tss_ranges = [ (max(0, self.tts - self.extension_in), 
                                 self.tts + self.extension_out), ]
        else:
            self.tss_ranges = [ (max(0, self.tss - self.extension_out), 
                                 self.tss + self.extension_in), ]
            self.tts_ranges = [ (max(0, self.tts - self.extension_in), 
                                 self.tts + self.extension_out), ]
        E.debug( "tss=%s, tts=%s" % (self.tss_ranges, self.tts_ranges) )

        self.counts_tss = self.counter.getCounts( contig, self.tss_ranges )
        self.counts_tts = self.counter.getCounts( contig, self.tts_ranges )

         ## revert for negative strand
        if strand == "-":
            self.counts_tss = self.counts_tss[::-1]
            self.counts_tts = self.counts_tts[::-1]
            
        self.aggregate( self.counts_tss, self.counts_tts )

        return 1

def countFromGTF(counters,
                 gtf_iterator):
    '''compute counts using counters for
    transcripts in gtf_iterator.
    '''

    c = E.Counter()
    counts = [0] * len(counters)

    E.info("starting counting" )
    names = []
    for iteration, gtf in enumerate(gtf_iterator):
        name = gtf[0].transcript_id
        E.debug( "processing %s" % (name))
        names.append(names)
        gtf.sort( key = lambda x: x.start )
        c.input += 1
        for x, counter in enumerate(counters):
            counter.update( gtf )
            counts[x] += 1

        if iteration % 100 == 0:
            E.debug( "iteration %i: counts=%s" % (iteration, ",".join( map( str, counters) ) ))

    E.info( "counts: %s: %s" % (str(c), ",".join( map(str,counts)))) 
    return names


def countFromCounts( counters,
                     all_counts ):
    '''collect counts from dataframe all_counts
    for counters.
    '''

    c = E.Counter()
    counts = [0] * len(counters)

    E.info("starting counting" )

    for iteration, x in enumerate(all_counts.iterrows()):
        c.input += 1
        name, read_counts = x
        for x, counter in enumerate(counters):
            counter.aggregate( numpy.array(read_counts) )
            counts[x] += 1

        if iteration % 100 == 0:
            E.debug( "iteration %i: counts=%s" % (iteration, ",".join( map( str, counters) ) ))

    E.info( "counts per counter: %s: %s" % (str(c), ",".join( map(str,counts)))) 
    return

        
    
        
