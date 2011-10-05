################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
'''
gtf2table.py - annotate genes/transrcipts
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

compute sequence properties for genes of given by a gtf file and output them in 
tabular format.

Usage
-----

Example::

   python gtf2table.py --help

Type::

   python gtf2table.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os, sys, string, re, optparse, math, time, tempfile, subprocess, types, bisect, array, collections
import itertools
import GFF, GTF, Bed, IOTools
import Experiment as E
import IndexedFasta
import Stats
import SequenceProperties
import Genomics
import Intervals

import bx
import bx.bbi.bigwig_file
import bx.intervals.io
import bx.intervals.intersection
import alignlib
import numpy
import IndexedGenome
import pysam

def readIntervalsFromGFF( filename_gff, 
                          source, 
                          feature, 
                          with_values = False, 
                          with_records = False, 
                          fasta = None, 
                          merge_genes = False,
                          format = "gtf",
                          use_strand = False ):
    """read intervals from a file or list.
    """

    assert not (with_values and with_records), "both with_values and with_records are true."

    ninput = 0

    if format in ("gtf", "gff"):
        infile = None
        # read data without value
        if type(filename_gff) == types.StringType:
            E.info(  "loading data from %s for source '%s' and feature '%s'" % (filename_gff, source, feature) )

            infile = IOTools.openFile( filename_gff, "r")        
            if format == "gtf":
                iterator_gff = GTF.iterator(infile)
            elif format == "gff":
                iterator_gff = GTF.iterator(infile)

        elif type(filename_gff) in (types.TupleType, types.ListType):

            E.info( "loading data from cache for source '%s' and feature '%s'" % (source, feature) )

            # from preparsed gff entries
            iterator_gff = filename_gff

        gff_iterator = GFF.iterator_filtered( iterator_gff,
                                              feature = feature,
                                              source = source )

        if format == "gtf":
            e = GTF.readAsIntervals( gff_iterator, 
                                     with_values = with_values, 
                                     with_records = with_records,
                                     merge_genes = merge_genes,
                                     use_strand = use_strand )
        elif format == "gff":
            e = GFF.readAsIntervals( gff_iterator, 
                                     with_values = with_values, 
                                     with_records = with_records,
                                     use_strand = use_strand )

        if infile: infile.close()

    elif format == "bed":
        if merge_genes: raise ValueError("can not merge genes from bed format" )
        if use_strand: raise NotImplementedError( "stranded comparison not implemented for bed format")
        iterator = Bed.iterator( IOTools.openFile(filename_gff, "r") )
        e = collections.defaultdict( list )
        if with_values:
            for bed in iterator:
                ninput += 1
                e[bed.contig].append( (bed.start,bed.end,bed.mFields[0]) )
        elif with_records:
            for bed in iterator:
                ninput += 1
                bed.gene_id = bed.mFields[0]
                bed.transcript_id = bed.gene_id
                e[bed.contig].append( (bed.start,bed.end,bed) )
        else:
            for bed in iterator:
                ninput += 1
                e[bed.contig].append( (bed.start,bed.end) )
        E.info("read intervals for %i contigs from %s: %i intervals" % (len(e), filename_gff, ninput) )

    else:
        raise ValueError("unknown format %s" % format )

    # translate names of contigs
    if fasta:
        if use_strand:
            for contig,strand in e.keys():
                if  contig in fasta:
                    x = e[contig]
                    del e[contig,strand]
                    e[fasta.getToken(contig),strand] = x
        else:
            for contig in e.keys():
                if  contig in fasta:
                    x = e[contig]
                    del e[contig]
                    e[fasta.getToken(contig)] = x

    return e

class Counter:
    """
    This class does not remove small exons/introns,
    so beware ENSEMBL that stores frameshifts as mini-exons.
    """

    header = [ "contig", "strand" ]

    mMinIntronSize = 10

    def __init__(self, fasta = None, section = None, options = None, prefix = None):
        self.mFasta = fasta
        self.mSection = section
        self.mOptions = options
        
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
        if self.mSection:
            return "\t".join( ["%s_%s" % (self.mSection, x) for x in self.header] )
        else:
            return "\t".join( self.header )

    def count(self):
        self.contig = self.getContig()
        self.strand = self.getStrand()

    def getContig(self):
        return self.mGFFs[0].contig

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
            s.append(self.mFasta.getSequence( contig, strand, start, end ))
            
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
        if self.mSection == "exons":
            return self.getExons()
        elif self.mSection == "introns":
            return self.getIntrons()
        else:
            return self.getExons()

class CounterIntronsExons(Counter):
    """count number of introns and exons.
    """

    header = ( "ntranscripts", "nexons", "nintrons", )
    def __init__(self, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs )

    def count(self):
        segments = self.getSegments()
        self.mNTranscripts = len( set( [x.transcript_id for x in self.mGFFs] ) )
        self.mNExons = len(segments)
        self.mNIntrons = len(segments) - 1

    def __str__(self):
        return "\t".join( (str(self.mNTranscripts), str(self.mNSegments), str(self.mNIntrons)) ) 

class CounterPosition(Counter):
    """output the position of the transcript."""
    header = ( "contig", "strand", "start", "end" )

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):

        segments = self.getSegments()
        self.contig = self.getContig()
        self.strand = self.getStrand()

        if len(segments) > 0:
            self.start = min( [ x[0] for x in segments ] )
            self.end = max( [ x[1] for x in segments ] )
        else:
            self.start, self.end = "na", "na"

    def __str__(self):
        return "\t".join( [ str(x) for x in (self.contig, self.strand, self.start, self.end) ] )

##----------------------------------------------------------------
class CounterLengths(Counter):
    header = Stats.Summary().getHeaders()

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):
        segments = self.getSegments()
        self.mResult = Stats.Summary( [x[1] - x[0] for x in segments ], mode="int", allow_empty = True )

    def __str__(self):
        return str(self.mResult)

##----------------------------------------------------------------
class CounterSpliceSites(Counter):

    mIntronTypes = ( ( "U2-GT/AG", "GT", "AG"),
                     ( "U2-nc-GC/AG", "GC", "AG"),
                     ( "U12-AT/AC", "AT", "AC") )

    mNames = [x[0] for x in mIntronTypes ] + ["unknown"]
    header = ["%s" % x for x in mNames ]

    mCheckBothStrands = True

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):

        self.mCounts = dict( [( x[0], 0) for x in self.mIntronTypes ] )
        self.mCounts["unknown"] = 0

        introns = self.getIntrons()
        contig = self.getContig()
        strand = self.getStrand()
        
        for start, end in introns:
            s = self.mFasta.getSequence( contig, "+", start, end )
            
            if self.mCheckBothStrands:
                r = self.getIntronType( s )
                if r == "unknown":
                    r = self.getIntronType( string.translate( s, string.maketrans("ACGTacgt", "TGCAtgca") )[::-1] )
            else:
                if Genomics.IsNegativeStrand( strand ):
                    s = string.translate( s, string.maketrans("ACGTacgt", "TGCAtgca") )[::-1]
                r = self.getIntronType( s )

            self.mCounts[r] += 1
                
    def __str__(self):
        return "\t".join( map(str, [ self.mCounts[x] for x in self.mNames ] ) )

    def getIntronType( self, sequence ):
        """return intron type for an intronic sequence."""

        for name, prime5, prime3 in self.mIntronTypes:
            if sequence[:len(prime5)].upper() == prime5 and \
                    sequence[-len(prime3):].upper() == prime3:
                return name
                
        else:
            return "unknown"

##-----------------------------------------------------------------------------------
class CounterCompositionNucleotides(Counter):
    header = SequenceProperties.SequencePropertiesNA().getHeaders() 

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        s = self.getSequence( ee )
        self.mResult = SequenceProperties.SequencePropertiesNA()
        self.mResult.loadSequence( s )

    def __str__(self):
        return str(self.mResult)

##-----------------------------------------------------------------------------------
class CounterCompositionCpG(Counter):
    header = SequenceProperties.SequencePropertiesCpg().getHeaders() 

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        s = self.getSequence( ee )
        self.mResult = SequenceProperties.SequencePropertiesCpg()
        self.mResult.loadSequence( s )

    def __str__(self):
        return str(self.mResult)

##-----------------------------------------------------------------------------------
class CounterOverlap(Counter):
    """count overlap with segments in another file.

    nover1 and nover2 count "exons".
    """

    headerTemplate = ( "nover1", "nover2", "nover", "pover1", "pover2" )

    ## do not save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = False

    mIsGTF = False

    mUseStrand = False

    def __init__(self, filename_gff, source, feature, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs )

        if feature and source:
            self.header = [ "%s:%s:%s" % (x, source, feature) for x in self.headerTemplate ]
        elif feature:
            self.header = [ "%s:%s" % (x, feature) for x in self.headerTemplate ]
        elif source:
            self.header = [ "%s:%s" % (x, source) for x in self.headerTemplate ]
        else:
            self.header = self.headerTemplate

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )
        
        e = readIntervalsFromGFF( filename_gff[0], source, feature, 
                                  self.mWithValues, self.mWithRecords, 
                                  self.mFasta, 
                                  format = self.mOptions.filename_format,
                                  use_strand = self.mUseStrand )

        # convert intervals to intersectors
        for key in e.keys():
            intersector = bx.intervals.intersection.Intersecter()
            if self.mWithValues or self.mWithRecords:
                for start, end, value in e[key]:
                    intersector.add_interval( bx.intervals.Interval(start,end,value=value) )
            else:
                for start, end in e[key]:
                    intersector.add_interval( bx.intervals.Interval(start,end) )

            e[key] = intersector

        self.mIntersectors = e

        E.info( "loading data finished" )

    def count(self):

        # collect overlapping segments
        segments = []
        contig = self.getContig()
        if self.mFasta: 
            contig = self.mFasta.getToken( contig)
            
        n = 0
        segments = self.getSegments()
        intervals = []
        if contig in self.mIntersectors:
            for start, end in segments:
                r = self.mIntersectors[contig].find( start, end )
                intervals += [ (x.start, x.end) for x in r ]
                if r: n += 1

        self.mNOverlap1 = n
        intervals = list(set(intervals))
        self.mNOverlap2 = len(intervals)

        intervals = Intervals.combineAtDistance( intervals,
                                                 self.mMinIntronSize )

        if n and len(intervals):
            self.mNOverlap = Intervals.calculateOverlap( segments, intervals )

            self.mPOverlap1 = 100.0 * self.mNOverlap / sum( [ end - start for start, end in segments] ) 
            self.mPOverlap2 = 100.0 * self.mNOverlap / sum( [ end - start for start, end in intervals ] )

            self.mONOverlap = "%i" % self.mNOverlap
            self.mOPOverlap1 = "%5.2f" % self.mPOverlap1
            self.mOPOverlap2 = "%5.2f" % self.mPOverlap2
        else:
            self.mONOverlap = 0
            self.mOPOverlap1 = 0
            self.mOPOverlap2 = 0
            self.mPOverlap1 = 0
            self.mPOverlap2 = 0
            self.mNOverlap = 0

    def __str__(self):
        return "\t".join( map(str, (self.mNOverlap1, 
                                    self.mNOverlap2, 
                                    self.mONOverlap, 
                                    self.mOPOverlap1,
                                    self.mOPOverlap2 ) ) )

##-----------------------------------------------------------------------------------
class CounterOverlapStranded(CounterOverlap):
    """count overlap with segments in another file.

    nover1 and nover2 count "exons".

    The overlap is stranded. 

    Negative values of overlap correspond antisense overlap with a feature.

    If there is both sense and antisense overlap, a warning is raised.
    """

    mUseStrand = True

    def __init__(self, *args, **kwargs ):
        CounterOverlap.__init__(self, *args, **kwargs )

    def count(self):

        # collect overlapping segments
        segments = self.getSegments()
        contig, strand = self.getContig(), self.getStrand()
        if self.mFasta: contig = self.mFasta.getToken( contig)
            
        def count( contig, strand ):
            n, intervals = 0, []
            if (contig,strand) in self.mIntersectors:
                for start, end in segments:
                    r = self.mIntersectors[contig,strand].find( start, end )
                    intervals.extend( [ (x.start, x.end) for x in r ] )
                    if r: n += 1
            return n, intervals

        sense = count( contig, strand )
        if strand == "+": strand = "-"
        else: strand = "+"
        antisense = count( contig, strand )

        if sense[0] and antisense[0]:
            E.warn( "%s overlapping both sense and antisense features" % self.getGeneId())
            self.counter.mixed_sense += 1

        is_sense = sense[0] > 0

        if  is_sense: n, intervals = sense
        else: n, intervals = antisense
            
        self.mNOverlap1 = n
        intervals = list(set(intervals))
        self.mNOverlap2 = len(intervals)

        intervals = Intervals.combineAtDistance( intervals,
                                                 self.mMinIntronSize )

        if n and len(intervals):
            self.mNOverlap = Intervals.calculateOverlap( segments, intervals )
            if not is_sense: self.mNOverlap = -self.mNOverlap

            self.mPOverlap1 = 100.0 * self.mNOverlap / sum( [ end - start for start, end in segments] ) 
            self.mPOverlap2 = 100.0 * self.mNOverlap / sum( [ end - start for start, end in intervals ] )

            self.mONOverlap = "%i" % self.mNOverlap
            self.mOPOverlap1 = "%5.2f" % self.mPOverlap1
            self.mOPOverlap2 = "%5.2f" % self.mPOverlap2
        else:
            self.mONOverlap = 0
            self.mOPOverlap1 = 0
            self.mOPOverlap2 = 0
            self.mPOverlap1 = 0
            self.mPOverlap2 = 0
            self.mNOverlap = 0

    def __str__(self):
        return "\t".join( map(str, (self.mNOverlap1, 
                                    self.mNOverlap2, 
                                    self.mONOverlap, 
                                    self.mOPOverlap1,
                                    self.mOPOverlap2 ) ) )

##-----------------------------------------------------------------------------------
class CounterOverlapTranscripts(CounterOverlap):
    """count overlap with segments in another file.

    Nover1 and nover2 now count "transcripts".
    """

    headerTemplate = ( "ngenes", "ntranscripts", "nexons", "nbases", "pover1", "pover2" )

    ## save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = True

    ## input is gtf
    mIsGTF = True

    def __init__(self, *args, **kwargs):
        CounterOverlap.__init__(self, *args, **kwargs )

    def count(self):

        # collect overlapping segments
        segments = []
        contig = self.getContig()
        if self.mFasta: 
            contig = self.mFasta.getToken( contig)
            
        n = 0
        segments = self.getSegments()
        intervals = []
        gene_ids, transcript_ids = set(), set()
        if contig in self.mIntersectors:
            for start, end in segments:
                r = self.mIntersectors[contig].find( start, end )
                intervals += [ (x.start, x.end) for x in r ]
                for x in r: gene_ids.add( x.value.gene_id )
                for x in r: transcript_ids.add( x.value.transcript_id )
                    
                if r: n += 1

        self.mNGenes = len(gene_ids)
        self.mNTranscripts = len(transcript_ids)
        intervals = list(set(intervals))
        self.mNExons = len(intervals)

        intervals = Intervals.combineAtDistance( intervals,
                                                 self.mMinIntronSize )

        if n and len(intervals):

            self.mNBases = Intervals.calculateOverlap( segments, intervals )
            self.mPOverlap1 = 100.0 * self.mNBases / sum( [ end - start for start, end in segments] ) 
            self.mPOverlap2 = 100.0 * self.mNBases / sum( [ end - start for start, end in intervals ] )
            self.mOPOverlap1 = "%5.2f" % self.mPOverlap1
            self.mOPOverlap2 = "%5.2f" % self.mPOverlap2

        else:
            self.mOPOverlap1 = 0
            self.mOPOverlap2 = 0
            self.mNBases = 0

    def __str__(self):
        return "\t".join( map(str, (
                    self.mNGenes,
                    self.mNTranscripts,
                    self.mNExons, 
                    self.mNBases,
                    self.mOPOverlap1,
                    self.mOPOverlap2 ) ) )

##-----------------------------------------------------------------------------------
class CounterCoverage(CounterOverlap):
    """compute base coverage with segments in another file.
    The values are output in 5' to 3' order for each base.
    """

    headerTemplate = ["cov_%s" % x for x in Stats.Summary().getHeaders()+ ("covered", "values",) ]

    ## do not save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = False

    def __init__(self, *args, **kwargs):
        CounterOverlap.__init__(self, *args, **kwargs )

    def count(self):

        # collect overlapping segments
        segments = []
        contig = self.getContig()
        if self.mFasta: 
            contig = self.mFasta.getToken( contig)
            
        n = 0
        segments = self.getSegments()
        segments.sort()

        map_genome2transcript = alignlib.makeAlignmentBlocks()

        x = 0
        for start, end in segments:
            map_genome2transcript.addDiagonal( start,
                                               end, 
                                               x - start)
            x += end - start

        l = sum( [ x[1] - x[0] for x in segments ] )

        counts = numpy.zeros( l, numpy.int )

        intervals = set()
        if contig in self.mIntersectors:
            for start, end in segments:
                r = self.mIntersectors[contig].find( start, end )
                for x in r: intervals.add( (x.start, x.end) )
                
        for start, end in intervals:
            for x in range(start,end):
                y = map_genome2transcript.mapRowToCol(x)
                if y >= 0: counts[ y ] += 1
        
        if Genomics.IsNegativeStrand( self.getStrand() ):
            self.mCounts = counts[::-1].copy()
        else:
            self.mCounts = counts

    def __str__(self):
        s = Stats.Summary( self.mCounts )

        if max(self.mCounts) > 0:
            values = ";".join( map(str,self.mCounts) )
            ncovered = sum( [1 for x in self.mCounts if x > 0 ] )
        else:
            values = "na"
            ncovered = 0

        return "\t".join( (str(s), str(ncovered), values) )

##-----------------------------------------------------------------------------------
class Classifier(Counter):
    """classify transcripts based on a reference annotation.

    This assumes the input is a genome annotation derived from an ENSEMBL gtf file
    created with gff2gtf.py.

    A transcript is classified as (threshold = self.mThresholdMinCoverage)

    known:     overlaps exons of a gene
    unknown:   overlaps any other region but exons
    ambiguous: can't say

    Type of a gene (only applies to known genes)
    pc:        protein coding
    utr:       is a utr transcript (not overlapping the coding part of a gene)
    pseudo:    pseudogene
    npc:       non of the above

    Location of an unkown transcript:
    intronic:     entirely intronic
    associated:   in flank of a protein coding gene, but not in exons or introns
                  (note that this changed from previous releases, where flank referred 
                  to any feature).
    """

    header = [ "is_known", "is_unknown", "is_ambiguous", 
               "is_pc", "is_pseudo", "is_npc", "is_utr", 
               "is_intronic", "is_assoc", "is_intergenic" ] 

    # features to use for classification
    features = ( "CDS", "UTR", "UTR3", "UTR5", "exon", "intronic", "intergenic", "flank", "3flank", "5flank", "telomeric" )
    
    # sources to use for classification
    sources = ("", "protein_coding", "pseudogene", )

    # minimum coverage of a transcript to assign it to a class
    mThresholdMinCoverage = 95

    # full coverage of a transcript to assign it to a class
    mThresholdFullCoverage = 99

    # some coverage of a transcript to assign it to a class
    mThresholdSomeCoverage = 10

    # do not use strand
    mUseStrand = False

    def __init__(self, filename_gff, *args, **kwargs ):

        Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )

        E.info( "loading data from %s" % (filename_gff[0]) )
            
        gffs = []
        infile = IOTools.openFile( filename_gff[0], "r")     
        for g in GTF.iterator(infile):
            gffs.append( g )

        E.info( "loaded data from %s" % (filename_gff[0]) )

        self.mCounters = {}
        self.mKeys = []

        if self.mUseStrand:
            counter = CounterOverlapStranded
        else:
            counter = CounterOverlap

        for source in self.sources:             
            for feature in self.features:
                key = "%s:%s" % (source, feature )
                self.mKeys.append( key )
                self.mCounters[key] = counter( [gffs], 
                                               source = source, 
                                               feature = feature, 
                                               fasta = self.mFasta,
                                               options = self.mOptions )
                
    def count(self):
        
        for key in self.mKeys:
            self.mCounters[key](self.mGFFs)

        def s_min( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdMinCoverage

        def s_excl( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) < (100 - self.mThresholdMinCoverage)

        def s_full( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdFullCoverage

        def s_some( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdSomeCoverage

        # classify wether it is know or unknown
        self.mIsKnown = s_min( ":exon", ":CDS", ":UTR", ":UTR3", ":UTR5" ) 
        self.mIsUnknown = s_excl( ":exon", ":CDS", ":UTR", ":UTR3", ":UTR5" )
        self.mIsAmbiguous = not( self.mIsUnknown or self.mIsKnown )

        # check type of gene:
        self.mIsNPC, self.mIsPC, self.mIsPseudo, self.mIsUTR = False, False, False, False
        if self.mIsKnown:
            self.mIsUTR = s_min( "protein_coding:UTR", "protein_coding:UTR3", "protein_coding:UTR5" )
            # for assigning as protein coding, also include the known UTR
            # do not include intronic sequence, as this conflicts with intronic
            # pseudo genes, which then leads to double assignments
            if not self.mIsUTR:
                self.mIsPC = s_min( "protein_coding:CDS", "protein_coding:UTR", "protein_coding:UTR3", "protein_coding:UTR5" )
                if not self.mIsPC:
                    self.mIsPseudo = s_min( "pseudogene:exon" ) 
            self.mIsNPC = not (self.mIsPC or self.mIsPseudo or self.mIsUTR )

        # classify location of unknown transcripts
        self.mIsIntronic, self.mIsAssociated, self.mIsIntergenic = False, False, False
        if self.mIsUnknown:
            self.mIsAssociated = s_some( "protein_coding:5flank", "protein_coding:3flank", "protein_coding:flank") 
            if not self.mIsAssociated:
                self.mIsIntronic = s_full( ":intronic" )
                self.mIsIntergenic = s_full( ":intergenic", ":telomeric" )

    def __str__(self):

        def to( v ):
            if v: return "1" 
            else: return "0"

        h = [ to(x) for x in (self.mIsKnown, self.mIsUnknown, self.mIsAmbiguous,
                             self.mIsPC, self.mIsPseudo, self.mIsNPC, self.mIsUTR,
                             self.mIsIntronic, self.mIsAssociated, self.mIsIntergenic ) ]

        for key in self.mKeys:
            h.append( str(self.mCounters[key]) )
        return "\t".join( h )

    def getHeader(self):
        h = [Counter.getHeader( self ) ]

        for key in self.mKeys:
            h.append( self.mCounters[key].getHeader() )

        return "\t".join( h )

##-----------------------------------------------------------------------------------
class ClassifierChIPSeq(Classifier):
    """classify ChIPSeq intervals based on a reference annotation.

    This assumes the input is a genome annotation derived from an ENSEMBL gtf file
    created with gff2gtf.py.

    In contrast to transcripts, the intervals are fuzzy. Hence the classification
    is based on a mixture of full/partial overlap.

    An interval is classified as:

    cds
       mostly part of a CDS. These would be intervals fully within a CDS exon.
    utr
       mostly part of UTR. These are intervals fully within the UTR of a gene.
    intergenic
       mostly intergenic. These are intervals fully within the intergenic region
       and more than 1kb from the closest exon.
    upstream
       not any of the above and partly upstream of a gene. These are intervals 
       that might overlap part of the UTR or the 1kb segment before to the 5'-terminal 
       exon of a gene.
    downstream
       not any of the abore and partly downstream of a gene. These are intervals 
       that might overlap part of the UTR or the 1kb segment after to the 3'-terminal 
       exon of a gene.
    intronic
       not any of the above and partly intronic. Note that these could also include
       promotors of short alternative transcripts that skip one or more of the first
       exons.
    ambiguous
       none of the above
    """

    header = [ "is_cds", "is_utr", "is_upstream", "is_downstream", "is_intronic", "is_intergenic", "is_flank", "is_ambiguous" ]

    # sources to use for classification
    sources = ("", ) # "protein_coding", "pseudogene", )

    # minimum coverage of a transcript to assign it to a class
    mThresholdMinCoverage = 95

    # full coverage of a transcript to assign it to a class
    mThresholdFullCoverage = 99

    # some coverage of a transcript to assign it to a class
    mThresholdSomeCoverage = 10

    def count(self):
        
        for key in self.mKeys:
            self.mCounters[key](self.mGFFs)

        def s_min( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdMinCoverage

        def s_excl( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) < (100 - self.mThresholdMinCoverage)

        def s_full( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdFullCoverage

        def s_some( *args ):
            return sum( [ abs(self.mCounters[x].mPOverlap1) for x in args ] ) >= self.mThresholdSomeCoverage

        self.mIsCDS, self.mIsUTR, self.mIsIntergenic = False, False, False
        self.mIsUpStream, self.mIsDownStream, self.mIsIntronic = False, False, False
        self.mIsFlank, self.mIsAmbiguous = False, False

        self.mIsCDS = s_full( ":CDS" )        
        self.mIsUTR = s_full( ":UTR", ":UTR3", ":UTR5" ) 
        self.mIsIntergenic = s_full( ":intergenic", ":telomeric" )

        if not(self.mIsCDS or self.mIsUTR or self.mIsIntergenic):
            self.mIsUpStream = s_some( ":5flank", ":UTR5" )
            if not self.mIsUpStream: 
                self.mIsDownStream = s_some( ":3flank", ":UTR3" )
                if not self.mIsDownStream:
                    self.mIsIntronic = s_some( ":intronic" )
                    if not self.mIsIntronic:
                        self.mIsFlank = s_some( ":flank" )

        self.mIsAmbiguous = not( self.mIsUTR or \
                                     self.mIsIntergenic or self.mIsIntronic or self.mIsCDS or \
                                     self.mIsUpStream or self.mIsDownStream or self.mIsFlank)

    def __str__(self):

        def to( v ):
            if v: return "1" 
            else: return "0"

        h = [ to(x) for x in (self.mIsCDS, 
                              self.mIsUTR, 
                              self.mIsUpStream,
                              self.mIsDownStream,
                              self.mIsIntronic,
                              self.mIsIntergenic,
                              self.mIsFlank,
                              self.mIsAmbiguous,
                              ) ]

        for key in self.mKeys:
            h.append( str(self.mCounters[key]) )
        return "\t".join( h )


##-----------------------------------------------------------------------------------
class ClassifierRNASeq(Counter):
    """classify RNASeq transcripts based on a reference annotation.

    Transcripts are classified by checking overlap with all known
    transcripts.

    If multiple transcripts overlap, select the one that is best matching.

    +--------------------+------------------------------------------------------------+
    |*Labels*            |*Contents*                                                  |
    +--------------------+------------------------------------------------------------+
    |complete            |Intron structure match - all exons are present, though first|
    |                    |and last exon might be different.                           |
    +--------------------+------------------------------------------------------------+
    |extended-fragment   |At least one intron boundary shared. ``predicted``          |
    |                    |transcript extends ``known`` transcript, but at the same is |
    |                    |incomplete.                                                 |
    +--------------------+------------------------------------------------------------+
    |extension           |At least one intron boundary shared. ``predicted``          |
    |                    |transcript extends ``known`` transcript.                    |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |exon-fragment       |Transcript overlap a single exon of a multi-exonic          |
    |                    |transcript.                                                 |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |fragment            |At least one intron boundary shared. ``predicted``          |
    |                    |transcript is shorter than ``known`` transcript.            |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |alternative         |At least one intron boundary shared. ``predicted``          |
    |                    |transcript has additional/missing exons/introns.            |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |unknown             |``predicted`` and ``known`` transcript overlap, but do not  |
    |                    |fall in any of the above categories.                        |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+

    Furthermore, ``predicted`` transcripts that do not overlap the exons of a ``known`` 
    transcript, but only introns, are classed as ``intronic``. All other transcripts
    that overlap neither introns nor exons of a ``known`` transcript are labeled 
    ``intergenic``.

    If the ``known`` transcript is protein coding, the ``predicted`` transcript is further 
    checked if it overlaps with the ``known`` UTR only. These transcripts are labelled ``utr5``
    and ``utr3``.

    Additionaly, the strandedness of the overlap is recorded as well (senes and antisense).

    To decide, which transcript is the closest match, the resultant class are sorted as above 
    in a priority list, where sense orientation has higher priority than anti-sense orientation.
    For example, a ``predicted`` transcript will be rather labeled as ``sense intronic`` rather
    than ``antisense fragment``.

    """

    header = [ "noverlap_transcripts", 
               "noverlap_genes",
               "match_transcript_id", "match_gene_id", "source", "class", "sense" ]

    # number of residues that are permitted for negligible overlap
    tolerance = 10

    # size of flanking region
    flank = 5000

    # priority of classifications to select best match
    priority = ( ( 's', "complete", ),
                 ( 's', "exon-fragment", ),
                 ( 's', "fragment", ),
                 ( 's', "extended-fragment", ),
                 ( 's', "extension", ),
                 ( 's', "alternative", ),
                 ( 's', "unknown", ),
                 ( 's', "intronic", ),
                 ( 's', "utr5", ),
                 ( 's', "utr3", ),
                 ( 's', "flank5", ),
                 ( 's', "flank3", ),
                 ( 'n', "complete", ),
                 ( 'n', "fragment", ),
                 ( 'n', "exon-fragment", ),
                 ( 'n', "extended-fragment", ),
                 ( 'n', "extension", ),
                 ( 'n', "alternative", ),
                 ( 'n', "unknown", ),
                 ( 'n', "intronic", ),
                 ( 'n', "utr5", ),
                 ( 'n', "utr3", ),
                 ( 'n', "flank5", ),
                 ( 'n', "flank3", ),
                 ( 'a', "complete", ),
                 ( 'a', "exon-fragment", ),
                 ( 'a', "fragment", ),
                 ( 'a', "extended-fragment", ),
                 ( 'a', "extension", ),
                 ( 'a', "alternative", ),
                 ( 'a', "unknown", ),
                 ( 'a', "intronic", ),
                 ( 'a', "utr5", ),
                 ( 'a', "utr3", ),
                 ( 'a', "flank5", ),
                 ( 'a', "flank3", ),
                 ( 'n', "flank", ),
                 ( 'n', "intergenic" ) )

    def __init__(self, filename_gff, *args, **kwargs ):

        Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )

        E.info( "loading data from %s" % (filename_gff[0]) )

        map_transcript2gene = {}
        transcripts = {}
        transcript_intervals = IndexedGenome.Quicksect()

        f = IOTools.openFile( filename_gff[0]) 

        for t in GTF.transcript_iterator(GTF.iterator( f )):
            t.sort( key = lambda x: x.start )
            transcript_id, gene_id = t[0].transcript_id, t[0].gene_id
            map_transcript2gene[transcript_id] = gene_id
            transcripts[transcript_id] = t
            transcript_intervals.add( t[0].contig, t[0].start, t[-1].end, transcript_id )

        f.close()

        E.info( "loaded data from %s" % (filename_gff[0]) )

        self.transcripts = transcripts
        self.transcript_intervals = transcript_intervals
        self.map_transcript2gene = map_transcript2gene
        self.map_gene2transcripts = dict( [ (y,x) for x,y in map_transcript2gene.iteritems() ] )

        self.mapClass2Priority = dict( [(y,x) for x,y in enumerate(self.priority) ] )

    def get_sense( self, strand, transcript_strand):
        '''encode sense as (s)ense, (a)ntisense and (n)ot available'''
        
        if strand == "." or transcript_strand == ".": sense = "n"
        elif strand == transcript_strand: sense = "s"
        else: sense = "a"
        return sense


    def classify_nonoverlap( self, exons ):
        '''classify_nonoverlapping transcripts.'''
        
        contig = self.getContig()
        strand = self.getStrand()
        start, end = exons[0][0], exons[-1][1]

        # get closest gene upstream and downstream
        before = self.transcript_intervals.before( contig, start, end, num_intervals = 1, max_dist = self.flank )
        after = self.transcript_intervals.after( contig, start, end, num_intervals = 1, max_dist = self.flank )

        # convert to id, distance, is_before
        associated = []
        if before:
            d = before[0]
            associated.append( (d[2], start - d[1], True ) )
            
        if after:
            d = after[0]
            associated.append( (d[2], d[0] - end, False ) )
                
        results = []

        for transcript_id, distance, is_before in associated:
            transcript_strand = self.transcripts[transcript_id][0].strand
            sense = self.get_sense( strand, transcript_strand )

            if transcript_strand == "+":
                if is_before: cls = "flank3"            
                else: cls = "flank5"
            elif transcript_strand == "-":
                if is_before: cls = "flank5"
                else: cls = "flank3"
            else:
                cls= "flank"
                
            results.append( (cls, sense, transcript_id ) )

        if len(results) == 0:
            results.append( ("intergenic", "n", "" ) )

        return results

    def classify_overlap( self, exons, transcript_id ):
        '''classify overlapping transcripts.'''

        introns = Intervals.complement( exons )
        strand = self.getStrand()
        lexons = Intervals.getLength( exons )

        start, end = exons[0][0], exons[-1][1]

        gtfs = self.transcripts[transcript_id]
        transcript_strand = gtfs[0].strand
        transcript_exons = [ (x.start, x.end) for x in gtfs if x.feature == "exon" ]
        transcript_start, transcript_end = transcript_exons[0][0], transcript_exons[-1][1]
        transcript_lexons = Intervals.getLength( transcript_exons )
        transcript_introns = Intervals.complement( transcript_exons )
        transcript_lintrons = Intervals.getLength( transcript_introns )

        tolerance = self.tolerance

        # check if transcript purely intronic
        overlap_exons = Intervals.calculateOverlap( exons, transcript_exons )

        if overlap_exons == 0: is_intronic = True
        else: is_intronic = False

        # determine sense-ness
        # if no strand is given, set to sense 
        if strand == "." or transcript_strand == ".": sense = True
        else: sense = strand == transcript_strand

        cls = "unclassified"

        if is_intronic:
            cls = "intronic"
        else:
            # count (exactly) shared introns
            shared_introns = [ x for x in introns if x in transcript_introns ]

            # count shared boundaries
            boundaries = sorted([ (x[0] - tolerance, x[0] + tolerance) for x in introns ] + 
                                [ (x[1] - tolerance, x[1] + tolerance) for x in introns ])
            transcript_boundaries = sorted([ (x[0], x[0] + 1) for x in transcript_introns ] + 
                                           [ (x[1], x[1] + 1) for x in transcript_introns ])
            shared_boundaries = Intervals.intersect( boundaries, transcript_boundaries )

            # if there are no introns, matched_structure will be True
            if len(exons) == 1 and len(transcript_exons) == 1:
                matched_structure = approx_structure = True
            elif len(exons) == 1:
                matched_structure = False
                approx_structure = True
            else:
                matched_structure = len(shared_introns) == len(introns)
                approx_structure = len(shared_boundaries) > 0
            
            if matched_structure and abs(overlap_exons - transcript_lexons) < tolerance:
                cls = "complete"
            elif approx_structure and len(exons) > 1 and abs(overlap_exons - transcript_lexons) < tolerance:
                cls = "complete"
            elif approx_structure and (start < transcript_start or end > transcript_end):
                if lexons < transcript_lexons:
                    cls = "extended-fragment"
                else:
                    cls = "extension"
            elif approx_structure and lexons < transcript_lexons:
                if len(exons) == 1 and len(transcript_exons) > 1:
                    cls = "exon-fragment"
                else:
                    cls = "fragment"
            elif approx_structure:
                cls = "alternative"
            else:
                cls = "unknown"

        transcript_cds = [ (x.start, x.end) for x in gtfs if x.feature == "CDS" ]
        transcript_lcds = Intervals.getLength( transcript_cds )
        # intersect with CDS
        if len(transcript_cds) > 0:
            overlap_cds = Intervals.calculateOverlap( transcript_cds, exons )                

            # if no overlap with CDS, check for UTR ovelap
            if overlap_cds < tolerance:
                utrs = Intervals.truncate( transcript_exons, transcript_cds )
                cds_start,cds_end = transcript_cds[0][0], transcript_cds[-1][1]
                utr5 = [ x for x in utrs if x[1] <= cds_start ]
                utr3 = [ x for x in utrs if x[0] >= cds_end ]
                if strand == "-": utr5, utr3= utr3, utr5
                overlap_utr5 = Intervals.calculateOverlap( utr5, exons )
                overlap_utr3 = Intervals.calculateOverlap( utr3, exons )
                if overlap_utr5 >= lexons - tolerance:
                    cls = "utr5"
                elif overlap_utr3 >= lexons - tolerance:
                    cls = "utr3"
            
        return cls, self.get_sense( strand, transcript_strand )

    def count(self):

        contig = self.getContig()
        segments = self.getSegments()
        introns = self.getIntrons()

        try:
            overlaps = list(self.transcript_intervals.get( contig, segments[0][0], segments[-1][1] ))
        except KeyError, msg:
            E.warn( "failed lookup of interval %s:%i-%i: '%s'" % (contig, segments[0][0], segments[-1][1], msg) )
            self.skip = True
            return

        noverlap_transcripts = len(overlaps)
        noverlap_genes = len( set( [self.map_transcript2gene[transcript_id] for start,end,transcript_id in overlaps] ) )

        results = []
        
        if len(overlaps) == 0:
            for cls, sense, transcript_id in self.classify_nonoverlap( segments ):
                if transcript_id != "":
                    source = self.transcripts[transcript_id][0].source
                    gene_id = self.map_transcript2gene[transcript_id]
                else:
                    source, gene_id = "", ""
                results.append( (self.mapClass2Priority[(sense,cls)], 
                                 ( 0, 0, transcript_id, gene_id, source, cls, sense ) ) )
        else:
            for start, end, transcript_id in overlaps:
                cls, sense = self.classify_overlap(segments, transcript_id) 
                source = self.transcripts[transcript_id][0].source
                gene_id = self.map_transcript2gene[transcript_id]
                results.append( (self.mapClass2Priority[(sense,cls)], 
                                 ( noverlap_transcripts, noverlap_genes, transcript_id, gene_id, source, cls, sense ) ) )
        
        results.sort()
        self.result = results[0][1]

    def __str__(self):
        return "\t".join( map(str, self.result) )

##-----------------------------------------------------------------------------------
class CounterOverrun(Counter):
    """count intron overrun. 
    
    This code is similar to counter, but ignores overlap with 
    introns at the 5' or 3' end of a transcript, i.e, only
    introns that are within the transcript range count.

    
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
       RRRRRRRRRR       RRRRRRRRRRRRRR
    ^ignore          ^count

    The code does not treat transcript models joining different genes
    specially.
    """

    header = ( "nover_exonic", "nover_intronic", "nover_external", "length" )

    ## Do not save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = False

    def __init__(self, filename_gff, *args, **kwargs):

        Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )

        source, feature = None, "CDS"
        
        e = readIntervalsFromGFF( filename_gff[0], 
                                  source,
                                  feature, 
                                  with_values = self.mWithValues, 
                                  with_records = self.mWithRecords, 
                                  fasta = self.mFasta,
                                  format = self.mOptions.filename_format )

        # convert intervals to intersectors
        for contig in e.keys():
            intersector = bx.intervals.intersection.Intersecter()
            for start, end in e[contig]:
                intersector.add_interval( bx.intervals.Interval(start,end) )

            e[contig] = intersector

        self.mIntersectors = e
        
        E.info( "loading data finished" )

    def count(self):

        # collect overlapping segments
        segments = []
        contig = self.getContig()
        if self.mFasta: 
            contig = self.mFasta.getToken( contig)
            
        n = 0
        segments = self.getSegments()
        intervals = []
        if contig in self.mIntersectors:
            for start, end in segments:
                r = self.mIntersectors[contig].find( start, end )
                intervals += [ (x.start, x.end) for x in r ]
                if r: n += 1

        intervals = Intervals.combineAtDistance( intervals,
                                                 self.mMinIntronSize )

        self.mNTotalLength = sum( [ x[1] - x[0] for x in segments ] )

        if n and len(intervals):

            # build list with internal introns
            internal_introns = Intervals.complement( intervals )
 
            # Truncate terminal CDS
            #start, end = min( [x[0] for x in segments] ), max( [ x[1] for x in segments ] )
            #intervals = [ ( max( start, x[0]), min( end,x[1]) ) for x in intervals ]

            # remove exons those not overlapping exons and
            # truncate terminal exons in transcripts 
            #start, end = min( [x[0] for x in intervals] ), max( [ x[1] for x in intervals ] )
            #print start, end, segments
            #segments = [ ( max( start, x[0]), min( end,x[1]) ) for x in segments if not (x[0] < end or x[1] > start) ]

            self.mNOverlapExonic = Intervals.calculateOverlap( segments, intervals )
            self.mNOverlapIntronic = Intervals.calculateOverlap( segments, internal_introns )
            self.mNOverlapExternal = sum( [ x[1] - x[0] for x in segments ] ) - self.mNOverlapExonic - self.mNOverlapIntronic

            assert self.mNOverlapExternal >= 0, \
                "intronic=%i, ovl=%i, intervals=%s - segments=%s" % (self.mNOverlapIntronic, 
                                                                     self.mNOverlapExonic,
                                                                     intervals, segments)
        else:
            self.mNOverlapExternal = 0
            self.mNOverlapExonic = 0
            self.mNOverlapIntronic = 0

    def __str__(self):
        return "\t".join( map(str, (self.mNOverlapExonic, 
                                    self.mNOverlapIntronic, 
                                    self.mNOverlapExternal,
                                    self.mNTotalLength,
                                    ) ) )

##-----------------------------------------------------------------------------------
class CounterDistance(Counter):
    """counter for computing the distance to features.

    The columns output are:
    distance: distance to closest feature in any direction
    id: id of closest feature *(
    dist5: distance to closest feature in 5' direction
    strand5: strand of feature in 5' direction
    id5: id of feature in 5' direction
    dist3: distance to closest feature in 3' direction
    strand3: strand of feature in 3' direction
    id3: id of feature in 5' direction

    This counter outputs the gene_id of the closest feature.

    The distance between transcripts and features that overlap
    is set to 0.
    """

    headerTemplate = ( "distance", "id", "dist5", "strand5", "id5", "dist3", "strand3", "id3" )

    ## do not save value for intervals
    mWithValues = False

    ## save records for intervals in order to get strand information
    mWithRecords = True

    ## also compute overlap
    mWithOverlap = False

    def __init__(self, filename_gff, source, feature, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs )

        if feature and source:
            self.header = [ "%s:%s:%s" % (x, source, feature) for x in self.headerTemplate ]
        elif feature:
            self.header = [ "%s:%s" % (x, feature) for x in self.headerTemplate ]
        elif source:
            self.header = [ "%s:%s" % (x, source) for x in self.headerTemplate ]
        else:
            self.header = self.headerTemplate

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )

        e = self.readIntervals( filename_gff[0], source, feature )

        # collect start and end points and points of intersection
        self.startPoints, self.startValues = {}, {}
        self.endPoints, self.endValues = {}, {}
        self.mIntervals = IndexedGenome.IndexedGenome()

        for contig, values in e.items():
            for start, end, value in values:
                self.mIntervals.add( contig, start, end, value )

            values.sort( key = lambda x: x[0] )
            self.startPoints[contig] = [ x[0] for x in values ]
            self.startValues[contig] = [ x[2] for x in values ]
            values.sort( key = lambda x: x[1] )
            self.endPoints[contig] = [ x[1] for x in values ]
            self.endValues[contig] = [ x[2] for x in values ]

        E.info( "loading data finished" )

    def readIntervals( self, filename_gff, source, feature ):

        return readIntervalsFromGFF( filename_gff, source, feature, 
                                     self.mWithValues, self.mWithRecords, self.mFasta,
                                     format = self.mOptions.filename_format )

    def count(self):
        """find closest feature in 5' and 3' direction."""

        # collect overlapping segments
        segments = self.getSegments()
        contig = self.getContig()
        if self.mFasta: contig = self.mFasta.getToken( contig)

        start = min( x[0] for x in segments )
        end = max( x[1] for x in segments )

        self.mData = None
        self.mDistance5 = "na"
        self.strand5 = "na"
        self.mData5 = None
        self.mDistance3 = "na"
        self.strand3 = "na"
        self.mData3 = None
        self.mIsNegativeStrand = Genomics.IsNegativeStrand( self.getStrand() )
        self.mDistance = "na"
        has5, has3 = False, False

        self.mOverlaps = []

        try:
            self.mOverlaps = list(self.mIntervals.get(contig, start, end))
        except KeyError:
            # ignore unknown contigs
            return

        if self.mOverlaps:
            # stop getting distance if there is overlap
            self.mDistance = 0
        else:
            if contig in self.startValues:
                entry_before = bisect.bisect( self.endPoints[contig], start ) - 1
                entry_after = bisect.bisect( self.startPoints[contig], end )

                if entry_before >= 0:
                    e = self.endValues[contig][entry_before]
                    self.mDistance5 = start - e.end
                    self.strand5 = e.strand
                    self.mData5 = e
                    has5 = True

                if entry_after < len( self.startPoints[contig] ):
                    e = self.startValues[contig][entry_after]
                    self.mDistance3 = e.start - end
                    self.strand3 = e.strand
                    self.mData3 = e 
                    has3 = True

            if has5 and has3:
                if self.mDistance5 < self.mDistance3:
                    self.mDistance, self.mData = self.mDistance5, self.mData5
                else:
                    self.mDistance, self.mData = self.mDistance3, self.mData3
            elif has5:
                self.mDistance =self.mDistance5 
                self.mData = self.mData5
            elif has3:
                self.mDistance =self.mDistance3
                self.mData = self.mData3

            if self.mIsNegativeStrand:
                self.mDistance5, self.strand5, self.mData5, self.mDistance3, self.strand3, self.mData3 = \
                    self.mDistance3, self.strand3, self.mData3, self.mDistance5, self.strand5, self.mData5

    def __str__(self):
        
        def toId( x ):
            if x: return x.gene_id
            else: return "na"

        return "\t".join( ( str(self.mDistance),
                            toId( self.mData),
                            str(self.mDistance5),
                            self.strand5,
                            toId( self.mData5),
                            str(self.mDistance3),
                            self.strand3,
                            toId( self.mData3) ))

##-----------------------------------------------------------------------------------
class CounterDistanceGenes(CounterDistance):
    """counter for computing the distance to genes.

    There are two distances of interest:
    1. distances towards the closest gene in 5' and 3' direction
       of the transcript
    2. distance towards the closest 5' and 3' of a gene

    The columns output are:
    closest_id: id of closest feature
    closest_dist: distance to closest feature
    closest_strand: strand of closest feature
    dist5: distance to closest feature in 5' direction
    strand5: strand of feature in 5' direction
    id5: gene id of closest gene in 5' direction
    dist3: distance to closest feature in 3' direction
    strand3: strand of feature in 3' direction
    id3: gene id of closest gene in 3' direction
    min5: minimum distance to 5' end of a gene
    min3: minimum distance to 3' end of a gene
    amin5: absolute minimum distance to 5' end of a gene (set to NULL if 3' end of a gene is closer)
    amin3: absolute minimum distance to 3' end of a gene (set to NULL if 5' end of a gene is closer)

    In order to check if gene is closer to 5' or 3' end of a gene, use::

       CASE WHEN amin3 > 0 THEN '3' ELSE '5' END

    TODO: check if all is correct.
    """

    headerTemplate = ( "closest_id", "closest_dist", "closest_strand", "id5", "dist5", "strand5", "id3", "dist3", "strand3", "min5", "min3", "amin5", "amin3" )

    ## do not save value for intervals
    mWithValues = False

    ## save records for intervals in order to get strand information
    mWithRecords = True

    def __init__(self, *args, **kwargs):
        CounterDistance.__init__(self, *args, **kwargs )

    def readIntervals( self, filename_gff, source, feature ):

        return readIntervalsFromGFF( filename_gff, source, feature, 
                                     self.mWithValues, self.mWithRecords, self.mFasta,
                                     merge_genes = True,
                                     format = self.mOptions.filename_format )

    def __str__(self):

        id5, id3, d3, d5, m3, m5 = "na", "na", "na" , "na", "na", "na"

        if self.mData5: 
            id5 = self.mData5.gene_id
        if self.mData3: 
            id3 = self.mData3.gene_id

        ## get distances to closest feature in 3' and 5' directions
        d5, d3 = [], []
        dist5, dist3, strand5, strand3 = self.mDistance5, self.mDistance3, self.strand5, self.strand3

        if self.mIsNegativeStrand:
            dist5, dist3 = dist3, dist5
            strand5, strand3 = strand3, strand5

        if not Genomics.IsNegativeStrand(strand5):
            d3.append( dist5 )
        else:
            d5.append( dist5 )

        if not Genomics.IsNegativeStrand(strand3):
            d5.append( dist3 )
        else:
            d3.append( dist3 )

        # get minimum distance to 3' or 5' end of a gene
        if d3: any3 = d3 = min(d3) 
        else: any3 = d3 = "na"
        if d5: any5 = d5 = min(d5)
        else: any5 = d5 = "na"
            
        # record the closest distance to a gene in any direction
        if any3 != "na" and any5 != "na":
            if any3 < any5: any5 = "na"
            elif any5 < any3: any3 = "na"

        # record the shortest distance
        if self.mDistance3 == None or self.mDistance5 == None:
            if self.mDistance3 == None:
                closest_id, closest_dist, closest_strand = self.mData5.gene_id, self.mDistance5, self.strand5
            elif self.mDistance5 == None:
                closest_id, closest_dist, closest_strand = self.mData3.gene_id, self.mDistance3, self.strand3
            else:
                closest_id, closest_dist, closest_strand = "na", "na", "na"
        elif self.mDistance3 < self.mDistance5:
            closest_id, closest_dist, closest_strand = self.mData3.gene_id, self.mDistance3, self.strand3
        elif self.mDistance5 < self.mDistance3:
            closest_id, closest_dist, closest_strand = self.mData5.gene_id, self.mDistance5, self.strand5
        else:
            closest_id, closest_dist, closest_strand = "na", "na", "na"

        return "\t".join( ( 
                closest_id,
                str(closest_dist),
                str(closest_strand),
                id5,
                str(self.mDistance5),
                self.strand5,
                id3,
                str(self.mDistance3),
                self.strand3,
                str(d5),
                str(d3),
                str(any5),
                str(any3)) )

##-----------------------------------------------------------------------------------
class CounterDistanceTranscriptionStartSites(CounterDistance):
    """counter for computing the distance to transcription start sites.

    There are two distances of interest:
    1. distances towards the closest tss.
    2. distances towards the closest promoter 5' or 3' end
       
    The columns output are:
    closest_id: id of closest feature
    closest_dist: distance to closest feature
    closest_strand: strand of closest feature
    id5: gene_id of 5' tss 
    dist5: minimum distance to 5' end of a tss
    id3: gene_id of 3' tss
    dist3: minimum distance to 3' end of a tss

    If the interval overlaps with a tss, the distances are set
    towards the closest end and the flag ``is_overlap`` is set.

    In order to check if gene is closer to 5' or 3' end of a gene, use::

       CASE WHEN amin3 > 0 THEN '3' ELSE '5' END

    TODO: check if all is correct.
    """


    headerTemplate = ( "closest_id", "closest_dist", "closest_strand", "id5", "dist5", "id3", "dist3", "is_overlap" )

    ## do not save value for intervals
    mWithValues = False

    ## save records for intervals in order to get strand information
    mWithRecords = True

    ## also compute overlap
    mWithOverlap = True

    def __init__(self, *args, **kwargs):
        CounterDistance.__init__(self, *args, **kwargs )

    def readIntervals( self, filename_gff, source, feature ):

        return readIntervalsFromGFF( filename_gff, source, feature, 
                                     self.mWithValues, self.mWithRecords, self.mFasta,
                                     merge_genes = False,
                                     format = self.mOptions.filename_format )
    def __str__(self):

        closest_id, closest_dist, closest_strand = ["na"] * 3
        d5_id, d5_dist, d3_id, d3_dist = ["na"] * 4

        overlaps = self.mOverlaps
        if overlaps:
            gene_ids = set([ x[2].gene_id for x in overlaps ])
            is_overlap = str(len(gene_ids))
            closest_id, closest_dist, closest_strand = ",".join( gene_ids), 0, "."            
        else:
            is_overlap = "0"

            ## get distances to closest features depending which
            ## direction they are pointing.
            ## The distance is modified, such that always the
            ## distance towards the 5' end is recorded.
            d5, d3 = [], []
            # 5'-----[ this ] ----- 3'

            if self.mData5:
                if Genomics.IsPositiveStrand(self.strand5):
                    # 5'->3'-----[this] 
                    # feature points in same direction
                    d3.append( (self.mDistance5 + (self.mData5.end - self.mData5.start), self.mData5.gene_id, self.strand5) )
                else:
                    # 3'<-5'-----[this] 
                    # feature points in different direction
                    d5.append( (self.mDistance5, self.mData5.gene_id, self.strand5) )

            if self.mData3:
                if Genomics.IsPositiveStrand(self.strand3):
                    # [this]-----5'->3'
                    d5.append( (self.mDistance3, self.mData3.gene_id, self.strand3) )
                else:
                    # [this]-----3'<-5'
                    d3.append( (self.mDistance3 + (self.mData3.end - self.mData3.start), self.mData3.gene_id, self.strand3) )

            # get minimum distance to 3' or 5' end of a gene
            if d3: d3 = min(d3) 
            else: d3 = []
            if d5: d5 = min(d5)
            else: d5 = []

            if d3 and d5:
                closest_dist, closest_id, closest_strand = min( d3, d5 )
            elif d3:
                closest_dist, closest_id, closest_strand = d3
            elif d5:
                closest_dist, closest_id, closest_strand = d5

            if d3: d3_dist, d3_id = d3[:2]
            if d5: d5_dist, d5_id = d5[:2]

        return "\t".join( ( 
                closest_id,
                str(closest_dist),
                str(closest_strand),
                d5_id,
                str(d5_dist),
                d3_id,
                str(d3_dist),
                is_overlap) )    

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterSpliceSiteComparison(CounterOverlap):
    """compare splice sites to reference
    """

    headerTemplate = [ "splice_%s" % x for x in ( "total", "found", "missed", "perfect", "partial", "incomplete", "exon_skipping" ) ]

    def __init__(self, *args, **kwargs):
        CounterOverlap.__init__(self, *args, **kwargs )

    def count(self):

        # collect overlapping segments with introns
        segments = self.getIntrons()
        contig = self.getContig()
        if self.mFasta: contig = self.mFasta.getToken( contig)

        self.mNMissed, self.mNFound = 0, 0
        self.mNPerfectMatch, self.mNPartialMatch = 0, 0
        self.mNIncompleteMatch, self.mNExonSkippingMatch = 0, 0
        for start, end in segments:
            if contig in self.mIntersectors:
                s = self.mIntersectors[contig].find( start, end )
                if s:
                    self.mNFound += 1
                    if len(s) == 1:
                        xstart, xend = s[0].start, s[0].end
                        if xstart == start and xend == end:
                            self.mNPerfectMatch += 1
                        elif xstart == start or xend == end:
                            self.mNPartialMatch += 1
                        else:
                            self.mNIncompleteMatch += 1
                    else:
                        self.mNExonSkippingMatch += 1
                else:
                    self.mNMissed += 1

        self.mNTotal = self.mNMissed + self.mNFound

    def __str__(self):
        return "\t".join( [str(x) for x in (self.mNTotal,
                                            self.mNFound,
                                            self.mNMissed,
                                            self.mNPerfectMatch,
                                            self.mNPartialMatch,
                                            self.mNIncompleteMatch,
                                            self.mNExonSkippingMatch ) ] )

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterProximity(CounterOverlap):
    """extract ranges in proximity to feature. 

    Outputs arrays and stats of the lengths and values 
    of the features in proximity. 

    This counter is useful for normalizing 
    a feature with features in the neighbourhood .

    Proximal ranges are within distance mMaxDistance and
    may overlap the range in question.
    """

    headerTemplate = Stats.Summary().getHeaders() + ( "length", "lengths", "values", ) 

    ## save value for intervals
    mWithValues = True

    def __init__(self, *args, **kwargs):
        CounterOverlap.__init__(self, *args, **kwargs )

        self.mProximalDistance = self.mOptions.proximal_distance

    def count(self):

        # collect overlapping segments
        segments = self.getSegments()
        contig = self.getContig()
        if self.mFasta: 
            contig = self.mFasta.getToken( contig)

        n = 0
        start = min( x[0] for x in segments )
        end = max( x[1] for x in segments )

        start = max( start - self.mProximalDistance, 0 )
        end += self.mProximalDistance

        if contig in self.mIntersectors:
            self.mSegments = self.mIntersectors[contig].find( start, end )
        else:
            self.mSegments = []

    def __str__(self):
        s = Stats.Summary( [ x.value for x in self.mSegments ] )
        values = ";".join( [ "%s-%s:%s" % (x.start, x.end, x.value) for x in self.mSegments ] )
        length = sum( [ x.end - x.start for x in self.mSegments ] )
        lengths = ";".join( [ str(x.end-x.start) for x in self.mSegments ] )
        return "\t".join( (str(s),str(length),lengths,values ) )

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterProximityExclusive(CounterProximity):
    """extract ranges in proximity with ranges.

    Proximal ranges are within distance mMaxDistance and
    may not overlap the range in question. Nested matches are 
    allowed.

    Use this class to normalize rates with surrounding ancestral 
    repeats.
    """

    def __init__(self, *args, **kwargs):
        CounterProximity.__init__(self, *args, **kwargs )

    def count(self):

        # get all segments within range
        CounterProximity.count(self)

        # remove all segments overlapping with any of these segments
        # quick and dirty exhaustive search against exons (should be only few)
        segments = self.getSegments()
        new = []
        for ss in self.mSegments:
            for s, e in segments:
                if min(ss.end,e)-max(ss.start,s) > 0:
                    break
            else:
                new.append(ss)

        self.mSegments = new

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterProximityLengthMatched(CounterProximity):
    """extract ranges in proximity with ranges and also
    do a length match.

    Segments are declared equal in length if they are within
    10% of the original segments length.

    Proximal ranges are within distance mMaxDistance and
    may not overlap the range in question. Nested matches are
    allowed.
    """

    mSizeDifference = 0.1

    def __init__(self, *args, **kwargs):
        CounterProximity.__init__(self, *args, **kwargs )

    def count(self):

        # get all segments within range
        CounterProximity.count(self)

        # remove all segments overlapping with any of these segments
        # quick and dirty exhaustive search against exons (should be only few)
        segments = self.getSegments()
        total_length = sum( x[1] - x[0] for x in segments )
        min_length = int(math.floor((1.0 - self.mSizeDifference) * total_length))
        max_length = int(math.ceil((1.0 + self.mSizeDifference) * total_length))
        new = []
        for ss in self.mSegments:
            l = ss.end - ss.start
            if l < min_length or l > max_length: 
                continue
            for s, e in segments:
                if min(ss.end,e)-max(ss.start,s) > 0:
                    break
            else:
                new.append(ss)

        self.mSegments = new


##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterNeighbours(CounterProximity):
    """extract features in second stream that are in proximity 
    to feature. 

    This class is similar to CounterProxmitiy, but outputs
    the identifiers of the second stream and the distances
    to the element (0 for overlap, negative distances for
    downstream and positive distances for upstream).

    Proximal ranges are within distance mMaxDistance and
    may overlap the range in question.
    """

    headerTemplate = ( "ids", "dists" ) 

    ## save value for intervals
    mWithValues = False
    mWithRecords = True
    mIsGTF = True 

    def count(self):

        # get all segments within range
        CounterProximity.count(self)

        # compute distance to segments
        segments = self.getSegments()
        start = min( x[0] for x in segments )
        end = max( x[1] for x in segments )

        self.mDistances = []

        for seg in self.mSegments:
            if end < seg.start:
                d = end - seg.start
            elif seg.end < start:
                d = start - seg.end
            else:
                d = 0

            if Genomics.IsNegativeStrand( seg.value.strand ):
                d = -d

            self.mDistances.append( d )

    def __str__(self):

        ids = ";".join( [x.value.gene_id for x in self.mSegments] )
        dists = ";".join( map(str,self.mDistances))
        return "\t".join( (ids,dists) )

##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------
class CounterTerritories(CounterOverlap):
    """compute overlap of a feature with territories.

    status flag:
    "U": unique match - the transcript overlaps a single territory
    "A": ambiguous match - the transcript overlaps several territories
    "0": no match: the transcript overlaps no territory
    """

    headerTemplate = ("status", "nterritories", "territories", )

    ## save pointer for intervals
    mWithRecords = True

    def __init__(self, *args, **kwargs):
        CounterOverlap.__init__(self, *args, **kwargs )

    def count(self):

        # collect overlapping segments
        segments = self.getSegments()
        contig = self.getContig()
        if self.mFasta: 
            try:
                contig = self.mFasta.getToken( contig)
            except KeyError:
                warn( "contig %s not found" % contig )

        n = 0
        start = min( x[0] for x in segments )
        end = max( x[1] for x in segments )

        if contig in self.mIntersectors:
            self.mSegments = self.mIntersectors[contig].find( start, end )
        else:
            self.mSegments = []
            
        if self.mSegments:
            if len(self.mSegments) == 1:
                self.mStatus = "U"
            else:
                self.mStatus = "A"
        else:
            self.mStatus = "0"

        self.mNAssignments = len(self.mSegments)
        self.mAssignments = ":".join( [x.value["gene_id"] for x in self.mSegments] )

    def __str__(self):
        return "\t".join( (self.mStatus, 
                           "%i" % self.mNAssignments,
                           self.mAssignments) )

##-----------------------------------------------------------------------------------
class CounterQuality(Counter):

    header = (Stats.Summary().getHeaders() + ("values",) )
    
    # discard segments with size > mMaxLength in order
    # to avoid out-of-memory
    mMaxLength = 100000

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        self.mResult = self.getQuality( ee )

    def getQuality(self, segments):
        """get sequence from a set of segments."""

        contig = self.getContig()
        # quality is always '+'
        strand = "+"
        
        s = None
        for start,end in segments:
            if end - start > self.mMaxLength: return []
            if not s:
                s = self.mFasta.getSequence( contig, strand, start, end )
            else:
                s.extend(self.mFasta.getSequence( contig, strand, start, end ))
        return s

    def __str__(self):
        s = Stats.Summary( self.mResult, mode = "int" )
        values = ";".join( [ str(x) for x in self.mResult ] )
        return "\t".join( (str(s),values ) )

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

        length = sum( [x[1] - x[0] for x in segments ] )
        counts_sense = numpy.zeros( length )
        counts_antisense = numpy.zeros( length )
        reads_sense, reads_antisense = set(), set()
        contig = self.getContig()
        if self.getStrand() == "+":
            is_reverse = False
        else:
            is_reverse = True

        def __add( counts, positions, offset ):
            for p in positions:
                pos = p - offset
                if 0 <= pos < length: counts[ pos ] += 1
            

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
                        __add( counts_sense, positions, offset )
                        reads_sense.add( read.qname )
                    else:
                        __add( counts_antisense, positions, offset )
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
    '''
    
    header = ( [ "%s_%s" % (x,y) for x,y in itertools.product( ( "sense", "antisense", "anysense"),
                                                                ( "unique_counts", "all_counts") ) ] )
    
    def __init__(self, bamfiles, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bamfiles: raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles

    def count(self):
        segments = self.getSegments()
        contig = self.getContig()

        nsense_unique_counts, nsense_all_counts = 0, 0
        nantisense_unique_counts, nantisense_all_counts = 0, 0
        nanysense_unique_counts, nanysense_all_counts = 0, 0

        if self.getStrand() == "+":
            is_reverse = False
        else:
            is_reverse = True

        for start, end in segments:
            for samfile in self.mBamFiles:
                last_any_pos = None
                last_sense_pos = None
                last_anti_pos = None
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

##-----------------------------------------------------------------------------------
class CounterReadExtension(Counter):
    '''compute read distribution from 3' to 5' end.

    Requires bam files to compute that coverage. Multiple bam
    files can be supplied, these will be summed up.

    Counts are separated into sense, antisense and any sense.

    This method requires a gff-file describing each gene's territory
    to avoid miscounting reads from genes in close proximity.



    Returns the coverage distribution in the territory, the median
    distance and the cumulative distribution every 1kb starting
    from the gene's end.
    '''

    # look at 15kb either way
    max_territory_size = 15000

    # distance increment
    increment = 100

    # minimum mapping quality
    # ignore those with 0 quality
    min_quality = 1

    def __init__(self, bamfiles, filename_gff, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 2:
            raise ValueError("expected two gff files: territories and UTRs" )

        self.labels = ("upstream", "downstream", "firstexon", "lastexon", "utr5", "utr3" )
        self.directions =  ("sense", "antisense", "anysense" )

        self.header = \
            tuple( 
            [ "%s_%s" % (x,y) for x,y in itertools.product( self.labels, ("length", "start", "end")) ] +\
                [ "%s_%s" % (x,y) for x,y in itertools.product(             
                        [ "%s_%s" % (xx, yy) for xx, yy in itertools.product(self.labels, self.directions)],
                        ("pcovered", ) + Stats.Summary().getHeaders() ) ] )

            
        self.outfiles = IOTools.FilePool( options.output_filename_pattern % "readextension_%s" )

        # -1 is the terminal exon
        for x,y in itertools.product( self.labels[:2], self.directions):
            self.outfiles.write( "%s_%s" % (x,y), 
                                 "%s\n" % "\t".join( 
                    ("gene_id", "length", "utr", "exon" ) + tuple(map(str, range( 0, 
                                                                                  self.max_territory_size, 
                                                                                  self.increment ) ) ) ) )

        if not bamfiles: raise ValueError("supply --bam-file options for readcoverage")
        self.mBamFiles = bamfiles

        filename_territories_gff, filename_utrs_gff = filename_gff

        # read territories
        self.territories = {}
        for gtf in GTF.iterator( IOTools.openFile( filename_territories_gff ) ):
            if gtf.gene_id in self.territories:
                raise ValueError( "need territories - multiple entries for gene %s" % gtf.gene_id)
            
            self.territories[gtf.gene_id] = (gtf.contig, gtf.start, gtf.end)
            
        # read known UTRs
        self.UTRs = collections.defaultdict( list )
        for gtf in GTF.iterator( IOTools.openFile( filename_utrs_gff ) ):
            if gtf.feature in ("UTR5", "UTR3", "UTR"):
                self.UTRs[gtf.gene_id].append( (gtf.contig, gtf.start, gtf.end) )

    def count(self):
        
        segments = self.getSegments()
        first_exon_start, last_exon_end = segments[0][0], segments[-1][1]
        first_exon_end, last_exon_start = segments[0][1], segments[-1][0]
        gene_id = self.mGFFs[0].gene_id
        self.strand = self.getStrand()
        self.gene_id = gene_id
        contig = self.getContig()
        if self.strand == "+":
            is_reverse = False
        else:
            is_reverse = True

        try:
            territory_contig, territory_start, territory_end = self.territories[gene_id]
        except KeyError:
            self.skip = True
            return

        # sanity check - is gene within territory?
        assert first_exon_start >= territory_start
        assert last_exon_end <= territory_end
        assert contig == territory_contig

        # truncate territory
        territory_start = max( territory_start, first_exon_start - self.max_territory_size )
        territory_end = min( territory_end, last_exon_end + self.max_territory_size )

        min_quality = self.min_quality

        # get UTRs - reorient them as upstream/downstream
        utrs = [None, None]
        if gene_id in self.UTRs:
            for utr_contig, utr_start, utr_end in self.UTRs[gene_id]:
                assert utr_contig == contig
                if utr_start == last_exon_end:
                    # downstream
                    utrs[1] = (utr_start, utr_end)
                elif utr_end == first_exon_start:
                    # upstream
                    utrs[0] = (utr_start, utr_end)
                else:
                    # ignore - "internal" UTRs
                    pass
#                    raise ValueError( "UTR mismatch for %s:%i-%i: utrs=%s" % \
#                                          (gene_id, 
#                                           first_exon_start, last_exon_end,
#                                           str(self.UTRs[gene_id]) ) )

        #############################################
        # these are pairs of before/after on
        # positive strand coordinates - will be re-oriented
        # later as downstream/upstream
        regions = [ (territory_start, first_exon_start),
                    (last_exon_end, territory_end),
                    (first_exon_start, first_exon_end),
                    (last_exon_start, last_exon_end ) ] + utrs
        
        counts = []

        for region in regions:
            if region:
                start, end = region
            else:
                # create dummy counts vector
                start, end = 0, 1
            
            counts.append( [numpy.zeros( end-start, numpy.int ),
                            numpy.zeros( end-start, numpy.int )] )

        def __add( counts, positions, offset ):
            l = len(counts)
            for p in positions:
                pos = p - offset
                if 0 <= pos < l: counts[ pos ] += 1

        def _update( counts, start, end ):
            counts_sense, counts_antisense = counts
            
            for samfile in self.mBamFiles:
                for read in samfile.fetch( contig, start, end ):
                    # ignore low quality mapped reads
                    if read.mapq < min_quality: continue
                    # only count positions actually overlapping
                    positions = read.positions
                    if not positions: continue
                    if is_reverse == read.is_reverse:
                        __add( counts_sense, positions, start )
                    else:
                        __add( counts_antisense, positions, start )

        for region, cc in zip( regions, counts):
            if region:
                start, end = region
                _update( cc, start, end )

        # invert "upstream" counts so that they are counting from TSS
        for x in xrange( 0, len(regions), 2):
            for y in range( len(counts[x]) ):
                counts[x][y] = counts[x][y][::-1]

        # switch upstream (x) /downstream (x+1)
        if is_reverse:
            for x in xrange( 0, len(regions), 2):
                counts[x], counts[x+1] = counts[x+1], counts[x]
                regions[x], regions[x+1] = regions[x+1], regions[x]

        # add anysense_counts
        for x in range(len(counts)):
            counts[x] = counts[x] + [ counts[x][0] + counts[x][1] ]
            
        self.counts = counts
        self.regions = regions

    def __str__(self):

        if self.skip: return "\t".join( ["na"] * len(self.header))

        r = []

        for label, region in zip( self.labels, self.regions):
            if region:
                start, end = region
                length = end - start
                r.append( "%i" % length )
                r.append( "%i" % start )
                r.append( "%i" % end )
            else:
                r.extend( ["na"] * 3 )

        # max number of distances, +1 for exon
        max_d = len(range(0,self.max_territory_size, self.increment) ) + 1

        # compute distributions for regions 
        for label, region, counts, exon_counts, utr_region in \
                zip( self.labels[:2], 
                     self.regions[:2], 
                     self.counts[:2], 
                     self.counts[2:4],
                     self.regions[4:6],
                     ):
            start, end = region
            length = end - start
            if utr_region:
                utr_start, utr_end = utr_region
                utr_extension = str(utr_end - utr_start)
            else:
                utr_extension = ""

            # output distributions (only for UTR counts)
            for direction, cc, ec in zip( self.directions, counts, exon_counts ):
                if len(ec) == 0: d = [0]
                else: d = [ec.max()]
                for x in range(0, length, self.increment ):
                    d.append( cc[x:min( length, x+self.increment)].max() )
                d.extend( [""] * (max_d - len(d)) )
                self.outfiles.write( "%s_%s" % (label, direction),
                                     "%s\t%i\t%s\t%s\n" % (self.gene_id,
                                                           length, 
                                                           utr_extension,
                                                           "\t".join( map(str, d ) ) ) ) 

        # output coverage stats
        for label, region, counts in zip( self.labels, self.regions, self.counts ):
            if not region:
                r.extend( ["na"] * (1+len(Stats.Summary().getHeaders()) ) )
                continue

            start, end = region
            length = end - start

            # output stats
            for direction, c in zip( self.directions, counts ):
                # compress
                c = c[c > 0]
                # pcovered
                if length == 0: r.append( "na" )
                else: r.append( "%5.2f" % (100.0 * len(c) / length) )

                # coverage stats
                r.append( str( Stats.Summary( c, mode = "int" ) ) )

        return "\t".join( r )

##-----------------------------------------------------------------------------------
class CounterBigwigCounts(Counter):
    '''obtain bigwig values and return summary stats.

    Requires a bigwig files to compute.
    '''
    
    header = ("length", "pcovered",) + Stats.Summary().getHeaders()
    
    # discard segments with size > mMaxLength in order
    # to avoid out-of-memory
    mMaxLength = 100000

    def __init__(self, bigwig_file, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs )
        if not bigwig_file: raise ValueError("supply --bigwig-file options for bigwig")
        self.mBigwigFile = bigwig_file

    def count(self):
        segments = self.getSegments()

        length = sum( [x[1] - x[0] for x in segments ] )
        nreads = 0
        contig = self.getContig()

        t, valid = None, None
        l = 0
        for start, end in segments:
            d = self.mBigwigFile.summarize( contig, start, end, end - start)
            if t != None: 
                t = numpy.append( t, d.sum_data)
                valid = numpy.append( valid, d.valid_count )
            else: t, valid = d.sum_data, d.valid_count

            l += end - start

        self.mTotalLength = length
        valid = valid[valid > 0]

        self.mCovered = len(valid)
        self.mResult = t

    def __str__(self):
        s = Stats.Summary( self.mResult, mode = "int" )
        return "\t".join( (str(self.mTotalLength),
                           "%5.2f" % (100.0 * self.mCovered / self.mTotalLength),
                           str(s),))


##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gtf2table.py 2888 2010-04-07 08:48:36Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-q", "--quality-file", dest="quality_file", type="string",
                      help="filename with genomic base quality information [default=%default]."  )

    parser.add_option("-b", "--bam-file", dest="bam_files", type="string",
                      help="filename with read mapping information. Multiple files can be submitted in a comma-separated list [default=%default]."  )

    parser.add_option("-i", "--bigwig-file", dest="bigwig_file", type="string",
                      help="filename with bigwig information [default=%default]."  )

    parser.add_option("-f", "--filename-gff", dest="filename_gff", type="string", action="append",
                      help="filename with extra gff files. The order is important [default=%default]."  )

    parser.add_option( "--filename-format", dest="filename_format", type="choice",
                       choices=("bed", "gff", "gtf" ),
                       help="format of secondary stream [default=%default]."  )

    parser.add_option( "--gff-source", dest="gff_sources", type="string", action="append",
                      help="restrict to this source in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option( "--gff-feature", dest="gff_features", type="string", action="append",
                      help="restrict to this feature in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option("-r", "--reporter", dest="reporter", type="choice",
                      choices=("genes", "transcripts" ),
                      help="report results for 'genes' or 'transcripts' [default=%default]."  )

    parser.add_option("-s", "--section", dest="sections", type="choice", action="append",
                      choices=("exons", "introns" ),
                      help="select range on which counters will operate [default=%default]."  )

    parser.add_option("-c", "--counter", dest="counters", type="choice", action="append",
                      choices=("length", "splice", "composition-na", "composition-cpg", 
                               "overlap", 
                               "classifier", 
                               "classifier-chipseq",
                               "classifier-rnaseq",
                               "overlap-stranded",
                               "overlap-transcripts",
                               "read-coverage", 
                               "read-extension", 
                               "read-counts",
                               "bigwig-counts",
                               'neighbours',
                               "proximity", "proximity-exclusive", "proximity-lengthmatched",
                               "position", "territories", "splice-comparison", 
                               "distance", "distance-genes", "distance-tss",
                               "coverage", "quality", "overrun" ),
                      help="select counters to apply [default=%default]."  )

    parser.add_option( "--add-gtf-source", dest="add_gtf_source", action="store_true",
                      help="add gtf field of source to output [default=%default]."  )

    parser.add_option( "--proximal-distance", dest="proximal_distance", type="int",
                      help="distance to be considered proximal to an interval [default=%default]."  )

    parser.add_option( "--prefix", dest="prefixes", type="string", action="append",
                      help="add prefix to column headers - prefixes are used in the same order as the counters [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        reporter = "genes",
        with_values = True,
        sections = [],
        counters = [],
        filename_gff = [],
        filename_format = "gtf",
        gff_features = [],
        gff_sources = [],
        add_gtf_source = False,
        proximal_distance = 10000, 
        bam_files = None,
        prefixes = []
        )

    (options, args) = E.Start( parser, add_output_options = True)

    if options.prefixes:
        if len(options.prefixes) != len(options.counters):
            raise ValueError("if any prefix is given, the number of prefixes must be the same as the number of counters" )

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.quality_file:
        quality = IndexedFasta.IndexedFasta( options.quality_file )
        quality.setTranslator( IndexedFasta.TranslatorBytes() )
    else:
        quality = None

    if options.bam_files:
        bam_files = []
        for bamfile in options.bam_files.split(","):
            bam_files.append( pysam.Samfile(bamfile, "rb" ) )
    else:
        bam_files = None

    if options.bigwig_file:
        bigwig_file = bx.bbi.bigwig_file.BigWigFile( open(options.bigwig_file))
    else:
        bigwig_file = None

    counters = []

    if not options.sections: 
        E.info( "counters will use the default section (exons)" )
        options.sections.append( None)

    if not options.gff_sources: options.gff_sources.append( None )
    if not options.gff_features: options.gff_features.append( None )

        

    cc = E.Counter()

    for n, c in enumerate(options.counters):
        if options.prefixes:
            prefix = options.prefixes[n]
        else:
            prefix = None

        if c == "position":
            for section in options.sections:
                counters.append( CounterPosition( section = section, options = options, prefix = prefix ) )
        elif c == "length":
            for section in options.sections:
                counters.append( CounterLengths( section = section, options = options, prefix = prefix ) )
        elif c == "splice":
            counters.append( CounterSpliceSites( fasta=fasta, prefix = prefix ) )
        elif c == "quality":
            counters.append( CounterQuality( fasta=quality, prefix = prefix ) )
        elif c == "overrun":
            counters.append( CounterOverrun( filename_gff = options.filename_gff,
                                             options = options, prefix = prefix ) )
        elif c == "read-coverage":
            counters.append( CounterReadCoverage( bam_files,
                                                  options = options, prefix = prefix ) )
        elif c == "read-extension":
            counters.append( CounterReadExtension( bam_files,
                                                   filename_gff = options.filename_gff,
                                                   options = options, prefix = prefix ) )
        elif c == "read-counts":
            counters.append( CounterReadCounts( bam_files,
                                                options = options, prefix = prefix ) )
        elif c == "bigwig-counts":
            counters.append( CounterBigwigCounts( bigwig_file,
                                                  options = options, prefix = prefix ) )

        elif c == "splice-comparison":
            counters.append( CounterSpliceSiteComparison( fasta=fasta, 
                                                          filename_gff = options.filename_gff,
                                                          feature=None, 
                                                          source=None, 
                                                          options=options, prefix = prefix ) )
        elif c == "composition-na":
            for section in options.sections:
                counters.append( CounterCompositionNucleotides( fasta=fasta,
                                                                section = section,
                                                                options = options, prefix = prefix ) )
        elif c == "composition-cpg":
            for section in options.sections:
                counters.append( CounterCompositionCpG( fasta=fasta,
                                                                section = section,
                                                                options = options, prefix = prefix ) )

        elif c in ( "overlap", 
                    "overlap-stranded",
                    "overlap-transcripts", 
                    "proximity", "proximity-exclusive", "proximity-lengthmatched",
                    "neighbours",
                    "territories", 
                    "distance", "distance-genes", "distance-tss",
                    "coverage" ):
            if c == "overlap":
                template = CounterOverlap
            if c == "overlap-stranded":
                template = CounterOverlapStranded
            elif c == "overlap-transcripts":
                template = CounterOverlapTranscripts
            elif c == "proximity":
                template = CounterProximity
            elif c == "neighbours":
                template = CounterNeighbours
            elif c == "proximity-exclusive":
                template = CounterProximityExclusive
            elif c == "proximity-lengthmatched":
                template = CounterProximityLengthMatched
            elif c == "territories":
                template = CounterTerritories
            elif c == "distance":
                template = CounterDistance
            elif c == "distance-genes":
                template = CounterDistanceGenes
            elif c == "distance-tss":
                template = CounterDistanceTranscriptionStartSites
            elif c == "coverage":
                template = CounterCoverage

            for section in options.sections:
                for source in options.gff_sources:
                    for feature in options.gff_features:
                        counters.append( template( filename_gff = options.filename_gff,
                                                   feature = feature,
                                                   source = source,
                                                   fasta=fasta,
                                                   section = section,
                                                   options = options, prefix = prefix) )

        elif c == "classifier":
            counters.append( Classifier( filename_gff = options.filename_gff,
                                         fasta = fasta,
                                         options = options, prefix = prefix) )
        elif c == "classifier-chipseq":
            counters.append( ClassifierChIPSeq( filename_gff = options.filename_gff,
                                                fasta = fasta,
                                                options = options, prefix = prefix) )

        elif c == "classifier-rnaseq":
            counters.append( ClassifierRNASeq( filename_gff = options.filename_gff,
                                               fasta = fasta,
                                               options = options, prefix = prefix) )

    if options.reporter == "genes":
        iterator = GTF.flat_gene_iterator
        header = ["gene_id"]
        fheader = lambda x: [x[0].gene_id]

    elif options.reporter == "transcripts":
        iterator = GTF.transcript_iterator
        header = ["transcript_id"]
        fheader = lambda x: [x[0].transcript_id]

    if options.add_gtf_source:
        header.append( "source" )
        ffields = lambda x: [x[0].source]
    else:
        ffields = lambda x: []

    options.stdout.write( "\t".join( 
            header + [ x.getHeader() for x in counters] ) + "\n" )

    for gffs in iterator( GTF.iterator(options.stdin) ):
        cc.input += 1
        for counter in counters: counter(gffs)

        skip = len( [x for x in counters if x.skip] ) == len(counters)
        if skip:
            cc.skipped += 1
            continue

        options.stdout.write( "\t".join( \
                fheader(gffs) +\
                ffields(gffs) +\
                [str( counter ) for counter in counters] ) + "\n" )
        
        cc.output += 1

    E.info("%s" %str(cc))
    for counter in counters:
        E.info( "%s\t%s" % (repr(counter), str(counter.counter) ) )
    E.Stop()
