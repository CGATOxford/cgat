################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: snp2table.py 2861 2010-02-23 17:36:32Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
annotator_distance.py - statistical significance of distance between genomic segments
=====================================================================================

Purpose
-------

The script :file:`annotator_distance.py` computes the statististical significance 
of the association between segments on a genome.

The genome is arranged into ``workspaces``, contiguous blocks in the genome that are
labeled at their termini. An example of workspaces are intergenic spaces labeled ``3'`` 
or ``5'`` according to the adjacent gene.

Within the workspace are ``segments`` of interest, for example transcription factor binding
sites. The script counts the distance of a segment to the nearest ``workspace`` terminus. 
The counts are aggregated over all workspaces.

Next, the script randomly rearranges ``segments`` within a ``workspace`` in order to test the 
statistical significance of the observed counts.

The script implements various sampling and counting methods.

Usage
-----

Type::

   python annotator_distances.py --help

for command line usage.

Code
----
"""

import os
import sys
import optparse
import collections
import itertools
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.Intervals as Intervals
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import bx.intervals.io
import bx.intervals.intersection
import numpy
import random
import math
import array
import scipy
import scipy.stats

import matplotlib.pyplot as plt
from matplotlib import _pylab_helpers

## global functions, defined once for optimization purposes
normalize_transform = lambda x, y: numpy.array( x, float) / ( sum(x) + y)
cumulative_transform = lambda x, y: numpy.cumsum( numpy.array( x, float) / (sum(x) + y) )

def readWorkspace( infile, 
                   workspace_builder = "raw", 
                   label="none",
                   map_id2annotation = {} ):
    """read workspace from infile.

    A workspace is a collection of intervals with two labels associated
    to each interval, one for the 5' and one for the 3' end.

    Available workspace builders are:

    gff
       take a gff file. 

    gtf-intergenic
       build workspace from intergenic segments in a gtf file. 

    gtf-intronic
       build workspace from intronic segments in a gtf file

    gtf-genic
       the workspace is built from genes (first to last exon).

    Available labels are:

    none
       no labels are given to the ends of workspaces

    direction
       labels are given based on the 5'/3' end of the
       bounding exon

    annotation
       labels are given based on a gene2annotation map.

    returns a list of segments for each contig in a dictionary
    """

    if label == "none":
        label_f = lambda x,y: (("X",),("X",))
        info_f = lambda x: None
    elif label == "direction":
        label_f = lambda x,y: ( ( ("5","3")[x],), ( ("3","5")[y],) )
        info_f = lambda x: x.strand == "+"
    elif label == "annotation":
        label_f = lambda x,y: (map_id2annotation[x], map_id2annotation[y])
        info_f = lambda x: x.gene_id
        
    if workspace_builder == "gff":
        workspace = GTF.readAsIntervals( GFF.iterator( infile ) )

    elif workspace_builder == "gtf-intergenic":

        workspace = collections.defaultdict( list )
        # get all genes
        for e in GTF.merged_gene_iterator( GTF.iterator( infile ) ):
            workspace[e.contig].append( (e.start, e.end, info_f( e ) ) )

        # convert to intergenic regions.
        # overlapping genes are merged and the labels
        # of the right-most entry is retained
        for contig in workspace.keys():
            segs = workspace[contig]
            segs.sort()
            last = segs[0]
            new_segs = []
            for this in segs[1:]:
                if last[1] >= this[0]:
                    if this[1] > last[1]:
                        last = (last[0], this[1], this[2] )
                    continue
                assert last[1] < this[0], "this=%s, last=%s" % (this, last)

                new_segs.append( (last[1], this[0], 
                                  label_f( last[2], this[2] ) ) )
                last = this
            workspace[contig] = new_segs

    elif workspace_builder == "gtf-intronic":

        workspace = collections.defaultdict( list )
        
        # the current procedure will count nested genes
        # twice
        for ee in GTF.flat_gene_iterator( GTF.iterator( infile) ):
            
            exons = Intervals.combine( [ (e.start, e.end) for e in ee ] )
            introns = Intervals.complement( exons ) 

            r = ee[0]
            for start, end in introns:
                workspace[r.contig].append( (start, 
                                              end, 
                                              label_f( info_f( r ), info_f(r) ) 
                                              ) )
    elif workspace_builder == "gtf-genic":

        workspace = collections.defaultdict( list )
        
        # the current procedure will count nested genes
        # twice
        for ee in GTF.flat_gene_iterator( GTF.iterator( infile) ):

            exons = Intervals.combine( [ (e.start, e.end) for e in ee ] )
            start, end = exons[0][0], exons[-1][1]
            r = ee[0]
            workspace[r.contig].append( (start, 
                                          end, 
                                          label_f( info_f( r ), info_f(r) ) 
                                          ) )

    else:
        raise ValueError("unknown workspace_builder %s" % workspace_builder)

    return workspace

def readSegments( infile, indexed_workspace, 
                  truncate = False, 
                  format="gtf", 
                  keep_ambiguous = False,
                  remove_overhangs = False ):
    """read segments from infile.

    segments not overlapping with indexed_workspace are removed.

    If :attr: truncate is given, segments extending beyond the workspace
    are truncated.
    
    returns a list of segments for each contig in a dictionary
    """
    counter = E.Counter()

    segments = collections.defaultdict( list )

    def addSegment( contig, start, end, counter ):
        if contig in indexed_workspace:
            r = indexed_workspace[contig].find( start, end )
            if not r:
                counter.nskipped += 1
                return
            if len(r) > 1:
                counter.nambiguous += 1
                if not keep_ambiguous: 
                    return
            if truncate:
                for x in r:
                    wstart, wend = x.start, x.end
                    rstart, rend = max( start, wstart), min( end, wend)
                    if start < wstart or end > wend: counter.ntruncated += 1
                    segments[contig].append( (rstart, rend) )
                    counter.added += 1
            elif remove_overhangs:
                for x in r:
                    wstart, wend = x.start, x.end
                    rstart, rend = max( start, wstart), min( end, wend)
                    if start < wstart or end > wend: 
                        counter.overhangs +=1
                        break
                else:
                    segments[contig].append( (start, end) )
            else:
                segments[contig].append( (start, end) )
                counter.added += 1

            counter.nkept += 1
    
    if format == "gtf":
        gtf_iterator = GTF.flat_gene_iterator( GTF.iterator( infile ) )

        for gene in gtf_iterator:
            # get start and end ignoring introns
            # contig, start, end = gene[0].contig, min( [x.start for x in gene] ), max( [x.end for x in gene] )

            contig, coords = gene[0].contig, [(x.start, x.end) for x in gene ]
            counter.ninput += 1
            for start,end in coords:
                addSegment( contig, start, end, counter )

    elif format == "bed":
        bed_iterator = Bed.iterator( infile )
        for bed in bed_iterator:
            counter.ninput += 1
            addSegment( bed.contig, bed.start, bed.end, counter )

    E.info( "read segments: %s" % str(counter) )

    return segments

class Sampler(object):
    """base clase for objcects that create a sample of
    randomly arranged segments in a workspace.

    """

    def __init__(self, observed, work_start, work_end ):
        self.mObserved = observed

        self.mWorkStart, self.mWorkEnd = work_start, work_end

        self.mLengths = [x[1] - x[0] for x in observed ]
        self.mTotalLength = sum( self.mLengths )
        self.mFreeLength = work_end - work_start - self.mTotalLength

        assert self.mFreeLength >= 0, "negative length: workspace=(%i,%i) %i-%i<0, segments=%s, lengths=%s" % \
            (work_start, work_end, work_end-work_start, self.mTotalLength, self.mObserved, self.mLengths)

    def sample( self ):
        raise NotImplementedError( "define sample() in base classes" )
        
class SamplerPermutation(Sampler):
    """permute order of fragments and distribute randomly.

    The permutation works like this:
    1. Randomly permutate the order of segments
    2. Split the free space (:attr:mFreeSpace) within the workspace into n+1 randomly sized gaps
    3. Insert the gaps between permutated segments

    """
    def sample( self ):
        """return simulated fragments."""
        simulated = []
        # 1. permutate order of segments
        random.shuffle( self.mLengths )
        # 2. determine size of space between samples
        points = []
        for x in range(len(self.mLengths) + 1 ):
            points.append( random.randint( 0, self.mFreeLength ) )
        points.sort()
        # 3. move segments to appropriate place
        start = self.mWorkStart
        simulated = []
        last = 0
        for x in range(len(self.mLengths)):
            start += points[x] - last
            simulated.append( (start, start + self.mLengths[x]) )
            start += self.mLengths[x]
            last = points[x]

        assert start + (points[-1] - last) <= self.mWorkEnd, "start=%i, points[-1]=%i, work_end=%i" % \
            (start, points[-1] -last, self.mWorkEnd)

        return simulated

class SamplerBlocks(Sampler):
    """move blocks of fragments to take into account clustering."""
    def sample( self ):
        """return simulated fragments."""
        simulated = []
        raise NotImplementedError

class SamplerGaps(Sampler):
    """rearrange gaps within a block randomly. 

    This sampler will preserve same of the clustering structure of segments."""
    def __init__(self, *args, **kwargs):
        Sampler.__init__(self,*args, **kwargs)
        self.mGapLengths = [x[1] - x[0] for x in Intervals.complement( self.mObserved, self.mWorkStart, self.mWorkEnd )]

    def sample( self ):
        """return simulated fragments."""
        simulated = []
        gaps = self.mGapLengths
        random.shuffle( gaps )
        
        start = self.mWorkStart

        for x in range(len(self.mLengths)):
            start += gaps[x]
            simulated.append( (start, start + self.mLengths[x]) )
            start += self.mLengths[x]

        return simulated

class CountingResults(object):
    """a container for observed and simulated counts.
    """
    def __init__(self,labels):
        self.mLabels = labels

        self.mTransform = None
        self.mEnvelopes = {}
        self.mMedians = {}

        self.mObservedCounts = None
        self.mSimulatedCounts = None
        self.mStats = None

    def updateFDR( self, obs_pvalues, sim_pvalues ):
        """compute fdr stats with given counts.

        If obs_pvalues and sim_pvalues are given, computes the FDR (q-value) for the observed p-value.
        The q-value is the expected proportion of false positive observations at
        the observed p-value.

        qvalue = A / B
        A: average proportion of simulated data with P-Values < pvalue (expected false positive RATE)
        B: number of observed data with P-Values < pvalue (NUMBER of true positives)

        As there are several counters and labels, all observed and simulated pvalues
        are taken into account.

        The method needs to be called after :meth:update.
        """

        assert self.mStats != None, "updateFDR called before calling update."

        for label in self.mLabels:
            pvalue = self.mStats[label].pvalue
            a = scipy.stats.percentileofscore( sim_palues, pvalue ) / 100.0
            b = scipy.stats.percentileofscore( obs_pvalues, pvalue ) / 100.0 * len(obs_pvalues ) 
            if b >= 0:
                qvalue = min( 1.0, a / b )
            else:
                qvalue = 0

            self.mStats[label] = self.mStats[label]._replace( qvalue=qvalue)

    def update( self ):
        """update stats from given counts.
        """
        
        assert self.mObservedCounts != None, "update called without observed counts."
        assert self.mSimulatedCounts != None, "update called without simulated counts."

        self.mStats = {}
        
        cls = collections.namedtuple( "st", "observed expected ci95lower ci95upper pvalue qvalue" )

        for label in self.mLabels:
            obs = cumulative_transform(self.mObservedCounts[label], self.mObservedCounts.mOutOfBounds[label])
            pobs = findMedian( obs )
            medians = self.getMedians( label )
            medians.sort()
            pvalue = float(scipy.stats.percentileofscore( medians, pobs )) / 100.0
            self.mStats[label] = cls( 
                pobs,
                scipy.mean( medians ),
                scipy.stats.scoreatpercentile( medians, 5 ),
                scipy.stats.scoreatpercentile( medians, 95 ),
                pvalue, None )

    def getLabels( self ):
        return self.mLabels

    def getMedians( self, label ):
        """compute medians of all samples."""
        if label not in self.mMedians: 

            num_samples = len(self.mSimulatedCounts)

            medians = []

            for x in range(num_samples):
                data = self.mSimulatedCounts[x][label]
                threshold = self.mSimulatedCounts[x].mTotals[label] / 2
                t = 0
                for d in range(len(data)):
                    if t > threshold: break
                    t += data[d]
                medians.append( d )

            self.mMedians[label] = medians

        return self.mMedians[label]

    def getEnvelope( self, label, transform ):
        """compute envelope for label using transform.

        The envelope is the min, max and mean of the observed counts
        add a certain position.

        This function does a lazy evaluation. Pre-computed results are stored
        and returned if the same transform is applied.
        """
        if label in self.mEnvelopes and transform == self.mTransform:
            E.debug( "returning cached envelope for transform %s" % str(transform) )
            return self.mEnvelopes[label]

        E.debug( "computing new envelope for transform %s" % str(transform) )

        num_samples = len(self.mSimulatedCounts)

        mmin = numpy.array( transform( self.mSimulatedCounts[0][label], self.mSimulatedCounts[0].mOutOfBounds[label]), numpy.float )
        msum = numpy.array( transform( self.mSimulatedCounts[0][label], self.mSimulatedCounts[0].mOutOfBounds[label]), numpy.float )
        mmax = numpy.array( transform( self.mSimulatedCounts[0][label], self.mSimulatedCounts[0].mOutOfBounds[label]), numpy.float )

        for x in range(1,num_samples):
            v = transform( self.mSimulatedCounts[x][label], self.mSimulatedCounts[x].mOutOfBounds[label])
            mmin = numpy.minimum( mmin, v )
            mmax = numpy.maximum( mmax, v )
            msum = msum + v
                
        msum /= num_samples

        self.mTransform = transform

        self.mEnvelopes[label] = ( mmin, mmax, msum )

        return self.mEnvelopes[label]

class Counter( object ):
    """return object that counts segments in a workspace.

    A counter will implement an addCounts method that expects a sorted
    list of intervals within a region bounded by start,end.
    """

    ## python list is fastest for single value access, but requires a lot of
    ## memory. python array is a good compromise - slightly slower than python list
    ## but uses much less space. For range access, use numpy arrays.
    mBuildCounts = lambda self, num_bins, dtype: array.array( "I", [0] * num_bins )

    def __init__(self, labels, num_bins, resolution = 1, dtype = numpy.int8 ):

        self.mCounts = {}
        self.mTotals = {}
        # keep separate out-of-bounds counts in order to not interfere with
        # dtype
        self.mOutOfBounds = {}

        b = self.mBuildCounts
        for l in labels:
            self.mCounts[l] = self.mBuildCounts( num_bins, dtype )
            self.mTotals[l] = 0
            self.mOutOfBounds[l] = 0

        self.mNumBins = num_bins
        self.mResolution = resolution
        
    def __getitem__( self, key ):
        return self.mCounts[key]

    def getLabels( self ):
        return self.mCounts.keys()

    def resolve( self, value ):
        if self.mResolution > 1:
            return int(math.floor( float(value) / self.mResolution) )
        else:
            return value

class CounterTranscription( Counter ):
    """count transcription per base."""

    mName = "Transcription"

    ## numpy is fastest for counting with blocks of data
    mBuildCounts = lambda self, num_bins, dtype: numpy.zeros( num_bins, dtype )

    def addCounts( self, rr, start, end, left_labels, right_labels ):

        counts = self.mCounts
        totals = self.mTotals
        ofb = self.mOutOfBounds
        nbins = self.mNumBins
        resolve = self.resolve

        for istart, iend in rr:
            l = iend - istart
            dl = istart - start
            dr = end - iend
            l = self.resolve( l )
            if dl < dr:
                pos = self.resolve( dl )
                labels = left_labels
            elif dl > dr:
                pos = self.resolve( dr )
                labels = right_labels
            else:
                continue

            if pos >= nbins:
                for label in labels:
                    ofb[label] += l
                    totals[label] += l
            else:
                for label in labels:
                    counts[label][pos:pos+l] += 1
                    totals[label] += l

class CounterClosestDistance( Counter ):
    """count closest distance."""

    mName = "Closest distance"
    def addCounts( self, rr, start, end, left_labels, right_labels ):

        counts = self.mCounts
        totals = self.mTotals
        ofb = self.mOutOfBounds
        nbins = self.mNumBins
        resolve = self.resolve

        def __add( pos, labels ):
            if pos >= nbins:
                for label in labels:
                    self.mOutOfBounds[label] += 1
                    totals[label] += 1
            else:
                for label in labels:
                    counts[label][pos] += 1
                    totals[label] += 1

        pos = self.resolve( rr[0][0] - start )
        __add( pos, left_labels )
        pos = self.resolve( end - rr[-1][1] )
        __add( pos, right_labels )

class CounterAllDistances( Counter ):
    """count all distances."""

    mName = "All distances"
    def addCounts( self, rr, start, end, left_labels, right_labels ):

        counts = self.mCounts
        totals = self.mTotals
        ofb = self.mOutOfBounds
        nbins = self.mNumBins
        resolve = self.resolve

        for istart, iend in rr:
            dl = istart - start
            dr = end - iend
            if dl < dr:
                pos = resolve( dl )
                labels = left_labels
            elif dl > dr:
                pos = resolve( dr )
                labels = right_labels
            else: 
                continue

            if pos >= nbins:
                for label in labels:
                    ofb[label] += 1
                    totals[label] += 1
            else:
                for label in labels:
                    counts[label][pos] += 1
                    totals[label] += 1

def indexIntervals( intervals, with_values = False ):
    """index intervals using bx.
    """

    indexed = {}

    for contig, values in intervals.iteritems():
        intersector = bx.intervals.intersection.Intersecter()
        if with_values:
            for start, end, value in values:
                intersector.add_interval( bx.intervals.Interval(start,end,value=value) )
        else:
            for start, end in values:
                intersector.add_interval( bx.intervals.Interval(start,end) )
        indexed[contig] = intersector
    return indexed

def plotCounts( counter, options, transform = lambda x: x ):
    """create plots from counter."""

    num_bins = options.num_bins
    resolution = options.resolution
    bins = numpy.array( xrange(num_bins) ) * resolution
    for label in counter.getLabels():
        fig = plt.figure()

        if options.plot_samples:
            for x in range(options.num_samples):
                counts = transform( counter.mSimulatedCounts[x][label], counter.mSimulatedCounts[x].mOutOfBounds[label] )
                plt.plot( bins, counts / t, label = "sample_%i" % x )

        if options.plot_envelope:
            # counts per sample are in row
            mmin, mmax, mmean = counter.getEnvelope( label, transform )

            plt.plot( bins, mmin, label = "min" )
            plt.plot( bins, mmax, label = "max" )
            plt.plot( bins, mmean, label = "mean" )

        plt.plot( bins, transform( counter.mObservedCounts[label], counter.mObservedCounts.mOutOfBounds[label] ), label = "observed" ) 
        plt.xlim( options.xrange )
        plt.legend()
        plt.title( counter.mName )

        plt.xlabel( "distance from gene / bp" )
        plt.ylabel( "frequency" )

        fig.suptitle( str(label) )

        if options.logscale:
            if "x" in options.logscale:
                plt.gca().set_xscale('log')
            if "y" in options.logscale:
                plt.gca().set_yscale('log')

        if options.hardcopy:
            plt.savefig( os.path.expanduser(options.hardcopy % label) )   

def findMedian( dist ):
    """find median in cumulative and normalized distribution."""
    x = 0
    while dist[x] < 0.5: x+=1
    return x

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: annotator_distance.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-a", "--filename-annotations", dest="filename_annotations", type="string",
                      help="filename mapping gene ids to annotations (a tab-separated table with two-columns) [default=%default]." )

    parser.add_option("-r", "--resolution", dest="resolution", type="int",
                      help="resolution of count vector [default=%default]." )

    parser.add_option("-b", "--num-bins", dest="num_bins", type="int",
                      help="number of bins in count vector [default=%default]." )

    parser.add_option("-i", "--num-samples", dest="num_samples", type="int",
                      help="sample size to compute [default=%default]." )
    
    parser.add_option("-w", "--workspace", dest="filename_workspace", type="string",
                      help="filename with workspace information [default=%default]." )
    
    parser.add_option( "--workspace-builder", dest="workspace_builder", type="choice",
                       choices=("gff", "gtf-intergenic", "gtf-intronic", "gtf-genic" ),
                       help="given a gff/gtf file build a workspace [default=%default]." )

    parser.add_option( "--workspace-labels", dest="workspace_labels", type="choice",
                       choices=("none", "direction", "annotation"),
                       help="labels to use for the workspace workspace [default=%default]." )

    parser.add_option( "--sampler", dest="sampler", type="choice",
                       choices=("permutation", "gaps"),
                       help="sampler to use. The sampler determines the null model of how segments are distributed in the workspace  [default=%default]" )
    
    parser.add_option( "--counter", dest="counters", type="choice", action="append",
                       choices=("transcription", "closest-distance", "all-distances"),
                       help="counter to use. The counter computes the quantity of interest [default=%default]" )

    parser.add_option( "--analysis", dest="analysis", type="choice", action="append",
                       choices=("proximity", "area-under-curve"),
                       help="analysis to perform [default=%default]" )

    parser.add_option( "--transform-counts", dest="transform_counts", type="choice",
                       choices=("raw", "cumulative" ),
                       help="cumulate counts [default=%default]." )

    parser.add_option("-s", "--segments", dest="filename_segments", type="string",
                      help="filename with segment information [default=%default]." )

    parser.add_option( "--xrange", dest="xrange", type="string",
                       help="xrange to plot [default=%default]")

    parser.add_option("-o", "--logscale", dest="logscale", type="string",
                      help="use logscale on x, y or xy [default=%default]"  )

    parser.add_option("-p", "--plot", dest="plot", action="store_true",
                      help="output plots [default=%default]"  )

    parser.add_option(  "--hardcopy", dest="hardcopy", type="string",
                      help="output hardcopies to file [default=%default]"  )

    parser.add_option( "--no-fdr", dest="do_fdr", action="store_false",
                      help="do not compute FDR rates [default=%default]"  )

    parser.add_option( "--segments-format", dest="segments_format", type="choice",
                       choices=("gtf", "bed" ),
                       help="format of segments file [default=%default]." )

    parser.add_option( "--truncate", dest="truncate", action="store_true",
                       help="truncate segments extending beyond a workspace [default=%default]"  )

    parser.add_option( "--remove-overhangs", dest="remove_overhangs", action="store_true",
                       help="remove segments extending beyond a workspace[default=%default]"  )

    parser.add_option( "--keep-ambiguous", dest="keep_ambiguous", action="store_true",
                       help="keep segments extending to more than one workspace [default=%default]"  )
    
    parser.set_defaults(
        filename_annotations = None,
        filename_workspace = "workspace.gff",
        filename_segments = "FastDown.gtf",
        filename_annotations_gtf = "../data/tg1_territories.gff",
        workspace_builder = "gff",
        workspace_labels = "none",
        sampler = "permutation",
        truncate = False,
        num_bins = 10000,
        num_samples = 10,
        resolution = 100,
        plot_samples = False,
        plot_envelope = True,
        counters = [],
        transform_counts = "raw",
        xrange = None,
        plot = False,
        logscale = None,
        output_all = False,
        do_test = False,
        analysis = [],
        do_fdr = True,
        hardcopy = "%s.png",
        segments_format = "gtf",
        remove_overhangs = False,
        )

    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    ###########################################
    ## setup options
    if options.sampler == "permutation":
        sampler = SamplerPermutation
    elif options.sampler == "gaps":
        sampler = SamplerGaps

    if options.xrange: options.xrange = map(float, options.xrange.split(",") )

    if len(options.counters) == 0:
        raise ValueError("please specify at least one counter.")

    if len(options.analysis) == 0:
        raise ValueError("please specify at least one analysis.")

    if options.workspace_labels == "annotation" and not options.filename_annotations:
        raise ValueError("please specify --filename-annotations is --workspace-labels=annotations.")

    ###########################################
    ## read data
    if options.workspace_labels == "annotation":
        def constant_factory(value):
            return itertools.repeat(value).next

        def dicttype():
            return collections.defaultdict( constant_factory( ("unknown",) ) )

        map_id2annotations = IOTools.readMultiMap( open( options.filename_annotations, "r"),
                                                   dtype = dicttype )
    else:
        map_id2annotations = {}

    workspace = readWorkspace( open( options.filename_workspace, "r" ), 
                               options.workspace_builder,
                               options.workspace_labels, 
                               map_id2annotations )

    
    E.info( "read workspace for %i contigs" % (len(workspace)) )

    indexed_workspace = indexIntervals( workspace, with_values = True )
    segments = readSegments( open( options.filename_segments, "r" ), indexed_workspace, 
                             format = options.segments_format,
                             keep_ambiguous = options.keep_ambiguous,
                             truncate = options.truncate,
                             remove_overhangs = options.remove_overhangs )

    nsegments = 0
    for contig, vv in segments.iteritems():
        nsegments += len(vv)
    
    E.info( "read %i segments for %i contigs" % (nsegments, len(workspace)) )
    indexed_segments = indexIntervals( segments, with_values = False)

    if nsegments == 0:
        E.warn( "no segments read - no computation done.")
        E.Stop()
        return
    
    # build labels
    labels = collections.defaultdict( int )
    for contig, vv in workspace.iteritems():
        for start, end, v in vv:
            for l in v[0]: labels[l] += 1
            for l in v[1]: labels[l] += 1

    E.info( "found %i workspace labels" % len(labels ) )
        
    ###########################################
    ## setup counting containers
    counters = []
    for cc in options.counters:
        
        if cc == "transcription":
            counter = CounterTranscription
        elif cc == "closest-distance":
            counter = CounterClosestDistance
        elif cc == "all-distances":
            counter = CounterAllDistances
            
        if nsegments < 256: 
            dtype = numpy.uint8
        elif nsegments < 65536: 
            dtype = numpy.uint16
        elif nsegments < 4294967296:
            dtype = numpy.uint32
        else:
            dtype = numpy.int

        E.debug( "choosen dtype %s" % str(dtype))

        E.info( "samples space is %i bases: %i bins at %i resolution" % \
                    (options.num_bins * options.resolution,
                     options.num_bins, 
                     options.resolution,
                     ) )

        E.info( "allocating counts: %i bytes (%i labels, %i samples, %i bins)" %\
                ( options.num_bins * len(labels) * dtype().itemsize * (options.num_samples + 1),
                  len(labels), 
                  options.num_samples,
                  options.num_bins,
                  ) )

        c = CountingResults( labels )
        c.mObservedCounts = counter( labels, options.num_bins, options.resolution, dtype = dtype )
        
        simulated_counts = []
        for x in range(options.num_samples):
            simulated_counts.append( counter( labels, options.num_bins, options.resolution, dtype = dtype ) )
        c.mSimulatedCounts = simulated_counts
        c.mName = c.mObservedCounts.mName

        counters.append( c )

        E.info( "allocated memory successfully" )
        
    segments_per_workspace = []
    segment_sizes = []
    segments_per_label = collections.defaultdict( int )
    workspaces_per_label = collections.defaultdict( int )

    ############################################
    # get observed and simpulated counts
    nworkspaces, nempty_workspaces, nempty_contigs, nmiddle = 0, 0, 0, 0
    iteration2 = 0
    for contig, vv in workspace.iteritems():
        
        iteration2 += 1
        E.info( "counting %i/%i: %s %i segments" % \
                    ( iteration2, 
                      len(workspace), 
                      contig, 
                      len(vv) ) )

        if len(vv) == 0: continue

        iteration1 = 0
        for work_start, work_end, v in vv:
            
            left_labels, right_labels = v[0], v[1]

            iteration1 += 1

            # ignore empty segments
            if contig not in indexed_segments: 
                nempty_contigs += 1
                continue 

            r = indexed_segments[contig].find( work_start, work_end )
            segments_per_workspace.append( len(r) )

            if not r:
                nempty_workspaces += 1
                continue
            
            # collect segments and stats
            nworkspaces += 1
            observed = [(x.start, x.end) for x in r ]
            observed.sort()
            segments_per_workspace.append( len(observed) )
            segment_sizes.extend( [x[1] - x[0] for x in observed] )

            # collect basic counts
            for label in list(left_labels) + list(right_labels):
                workspaces_per_label[label] += 1
                segments_per_label[label] += len(observed)

            # add observed counts
            for counter in counters:
                counter.mObservedCounts.addCounts( observed, work_start, work_end, left_labels, right_labels )

            # create sampler
            s = sampler( observed, work_start, work_end )

            # add simulated counts
            for iteration in range(options.num_samples):
                simulated = s.sample()
                for counter in counters:
                    counter.mSimulatedCounts[iteration].addCounts( simulated, work_start, work_end, left_labels, right_labels )

    E.info( "counting finished" )
    E.info( "nworkspaces=%i, nmiddle=%i, nempty_workspaces=%i, nempty_contigs=%i" % (nworkspaces, nmiddle, nempty_workspaces, nempty_contigs) )

    ######################################################
    ## transform counts

    if options.transform_counts == "cumulative":
        transform = cumulative_transform
    elif options.transform_counts == "raw":
        transform = normalize_transform

    ####################################################
    ## analysis

    if "proximity" in options.analysis:
        outfile_proximity = E.openOutputFile( "proximity" )
        outfile_proximity.write( "\t".join( ( "label","observed","pvalue",
                                              "expected","CIlower","CIupper","qvalue","segments","workspaces" )) + "\n" )
    else:
        outfile_proximity = None

    if "area-under-curve" in options.analysis:
        outfile_auc = E.openOutputFile( "auc" )
        outfile_auc.write( "label\tobserved\texpected\tCIlower\tCIupper\n" )
    else:
        outfile_auc = None

    ## qvalue: expected false positives at p-value
    ## qvalue = expected false positives / 
    if options.do_fdr:
        E.info( "computing pvalues for fdr" )
        for counter in counters:
            for label in labels:
                E.info( "working on counter:%s label:%s" % (counter, label))

                # collect all P-Values of simulated results to compute FDR
                sim_pvalues = []
                medians = counter.getMedians( label )
                
                for median in medians:
                    pvalue = float(scipy.stats.percentileofscore( medians, median )) / 100.0 
                    sim_pvalues.append( pvalue )
                    
        sim_pvalues.sort()
    else:
        sim_pvalues = []

    ## compute observed p-values
    for counter in counters:
        counter.update()
        
    obs_pvalues = []
    for counter in counters:
        for label in labels:
            obs_pvalues.append( counter.mStats[label].pvalue )
        obs_pvalues.sort()

    ## compute observed p-values
    if options.do_fdr:
        for counter in counters:
            counter.updateFDR( obs_pvalues, sim_pvalues )

    for counter in counters:

        outofbounds_sim, totals_sim = 0, 0
        outofbounds_obs, totals_obs = 0, 0
        for label in labels:
            for sample in range(options.num_samples):
                if counter.mSimulatedCounts[sample].mOutOfBounds[label]:
                    E.debug("out of bounds: sample %i, label %s, counts=%i" % \
                                ( sample, label, counter.mSimulatedCounts[sample].mOutOfBounds[label]) )
                    outofbounds_sim += counter.mSimulatedCounts[sample].mOutOfBounds[label]
                totals_sim += counter.mSimulatedCounts[sample].mTotals[label]

            outofbounds_obs += counter.mObservedCounts.mOutOfBounds[label]
            totals_obs += counter.mObservedCounts.mTotals[label]

        E.info( "out of bounds observations: observed=%i/%i (%5.2f%%), simulations=%i/%i (%5.2f%%)" % \
                (outofbounds_obs, totals_obs,
                 100.0 * outofbounds_obs / totals_obs,
                 outofbounds_sim, totals_sim,
                 100.0 * outofbounds_sim / totals_sim,
                 ) )

        for label in labels:
            
            if outfile_auc:
                mmin, mmax, mmean = counter.getEnvelope( label, transform = normalize_transform )
                obs = normalize_transform(counter.mObservedCounts[label], counter.mObservedCounts.mOutOfBounds[label])

                def block_iterator( a1, a2, a3, num_bins ):
                    x = 0
                    while x < num_bins:
                        while x < num_bins and a1[x] <= a2[x]: x += 1
                        start = x
                        while x < options.num_bins and a1[x] > a2[x]: x += 1
                        end = x
                        total_a1 = a1[start:end].sum()
                        total_a3 = a3[start:end].sum()
                        if total_a1 > total_a3:
                            yield (total_a1 - total_a3, start,end,total_a1, total_a3)

                blocks =  list(block_iterator( obs, mmax, mmean, options.num_bins))

                if options.output_all:
                    for delta, start,end,total_obs,total_mean in blocks:
                        if end - start <= 1: continue
                        outfile_auc.write( "%s\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\n" %\
                                               ( label,
                                                 start * options.resolution, 
                                                 end * options.resolution, 
                                                 (end-start) * options.resolution,
                                                 total_obs, 
                                                 total_mean, 
                                                 delta,
                                                 total_obs / total_mean,
                                                 100.0 * (total_obs / total_mean - 1.0) ) )

                # output best block
                blocks.sort()
                delta, start,end,total_obs,total_mean = blocks[-1]

                outfile_auc.write( "%s\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\n" %\
                                        ( label,
                                          start * options.resolution, 
                                          end * options.resolution, 
                                          (end-start) * options.resolution,
                                          total_obs, 
                                          total_mean, 
                                          delta,
                                          total_obs / total_mean,
                                          100.0 * (total_obs / total_mean - 1.0) ) )

            if outfile_proximity:

                # find error bars at median
                st = counter.mStats[label]
                outfile_proximity.write( "%s\t%i\t%f\t%i\t%i\t%i\t%s\t%i\t%i\n" %\
                                             ( label,
                                               st.observed * options.resolution,
                                               st.pvalue,
                                               st.expected * options.resolution,
                                               st.ci95lower * options.resolution,
                                               st.ci95upper * options.resolution,
                                               IOTools.prettyFloat( st.qvalue ),
                                               segments_per_label[label],
                                               workspaces_per_label[label],
                                               ) )

    if options.plot:
        
        for counter in counters:
            plotCounts( counter, options, transform )

        # plot summary stats
        plt.figure()
        plt.title("distribution of workspace length")
        data = []
        for contig, segs in workspace.iteritems():
            if len(segs) == 0: continue
            data.extend( [ x[1] - x[0] for x in segs ] )

        vals, bins = numpy.histogram( data, bins=numpy.arange(0,max(data),100), new = True)

        t = float(sum(vals))
        plt.plot( bins[:-1], numpy.cumsum( vals ) / t )
        plt.gca().set_xscale('log')
        plt.legend()
        t = float(sum(vals))
        plt.xlabel( "size of workspace")
        plt.ylabel( "cumulative relative frequency")
        if options.hardcopy:
            plt.savefig( os.path.expanduser(options.hardcopy % "workspace_size") )


        plt.figure()
        plt.title("segments per block" )
        vals, bins = numpy.histogram( segments_per_workspace, bins=numpy.arange(0,max(segments_per_workspace),1), new = True)
        plt.plot( bins[:-1], vals )
        plt.xlabel( "segments per block")
        plt.ylabel( "absolute frequency")
        if options.hardcopy:
            plt.savefig( os.path.expanduser(options.hardcopy % "segments_per_block") )

        plt.figure()
        plt.title("workspaces per label" )
        plt.barh( range(0,len(labels) ), [ workspaces_per_label[x] for x in labels ], height=0.5 )
        plt.yticks( range(0,len(labels) ), labels )
        plt.ylabel( "workspaces per label")
        plt.xlabel( "absolute frequency")
        plt.gca().set_xscale('log')

        if options.hardcopy:
            plt.savefig( os.path.expanduser(options.hardcopy % "workspaces_per_label") )

        plt.figure()
        plt.title("segments per label" )
        plt.barh( range(0,len(labels) ), [ segments_per_label[x] for x in labels ], height=0.5 )
        plt.yticks( range(0,len(labels) ), labels )
        plt.ylabel( "segments per label")
        plt.xlabel( "absolute frequency")
        plt.xticks( range(0,len(labels) ), labels )
        if options.hardcopy:
            plt.savefig( os.path.expanduser(options.hardcopy % "segments_per_label") )

        if not options.hardcopy:
            plt.show()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

