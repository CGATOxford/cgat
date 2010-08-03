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
import GFF, GTF, Bed
import Experiment as E
import IndexedFasta
import Stats
import SequenceProperties
import Genomics
import Intervals

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
                          format = "gtf" ):
    """read intervals from a file or list.
    """

    assert not (with_values and with_records), "both with_values and with_records are true."

    if format in ("gtf", "gff"):
        infile = None
        # read data without value
        if type(filename_gff) == types.StringType:
            E.info(  "loading data from %s for source '%s' and feature '%s'" % (filename_gff, source, feature) )

            infile = open( filename_gff, "r")        
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
                                     merge_genes = merge_genes )
        elif format == "gff":
            e = GFF.readAsIntervals( gff_iterator, 
                                     with_values = with_values, 
                                     with_records = with_records )

        if infile: infile.close()

    elif format == "bed":
        if merge_genes: raise ValueError("can not merge genes from bed format" )
        iterator = Bed.iterator( open(filename_gff, "r") )
        e = collections.defaultdict( list )
        if with_values:
            for bed in iterator:
                e[bed.contig].append( (bed.start,bed.end,bed.mFields[0]) )
        elif with_records:
            for bed in iterator:
                bed.gene_id = bed.mFields[0]
                bed.transcript_id = bed.gene_id
                e[bed.contig].append( (bed.start,bed.end,bed) )
        else:
            for bed in iterator:
                e[bed.contig].append( (bed.start,bed.end) )

    else:
        raise ValueError("unknown format %s" % format )

    # translate names of contigs
    if fasta:
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

    mHeader = [ "contig", "strand" ]

    mMinIntronSize = 10

    def __init__(self, fasta = None, section = None, options = None):
        self.mFasta = fasta
        self.mSection = section
        self.mOptions = options

    def __call__(self, gffs):
        self.mGFFs = gffs
        self.count()

    def __str__(self):
        return "\t".join( (self.contig, self.strand) )

    def getHeader(self):
        if self.mSection:
            return "\t".join( ["%s_%s" % (self.mSection, x) for x in self.mHeader] )
        else:
            return "\t".join( self.mHeader )

    def count(self):
        self.contig = self.getContig()
        self.strand = self.getStrand()

    def getContig(self):
        return self.mGFFs[0].contig

    def getStrand(self):
        return self.mGFFs[0].strand
    
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

    mHeader = ( "ntranscripts", "nexons", "nintrons", )
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
    mHeader = ( "contig", "strand", "start", "end" )

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
    mHeader = Stats.Summary().getHeaders()

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
    mHeader = ["%s" % x for x in mNames ]

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
    mHeader = SequenceProperties.SequencePropertiesNA().getHeaders() 

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
class CounterOverlap(Counter):
    """count overlap with segments in another file.

    Nover1 and nover2 count "exons".
    """

    mHeaderTemplate = ( "nover1", "nover2", "nover", "pover1", "pover2" )

    ## do not save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = False

    mIsGTF = False

    def __init__(self, filename_gff, source, feature, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs )

        if feature and source:
            self.mHeader = [ "%s:%s:%s" % (x, source, feature) for x in self.mHeaderTemplate ]
        elif feature:
            self.mHeader = [ "%s:%s" % (x, feature) for x in self.mHeaderTemplate ]
        elif source:
            self.mHeader = [ "%s:%s" % (x, source) for x in self.mHeaderTemplate ]
        else:
            self.mHeader = self.mHeaderTemplate

        e = readIntervalsFromGFF( filename_gff, source, feature, 
                                  self.mWithValues, self.mWithRecords, 
                                  self.mFasta, 
                                  format = self.mOptions.filename_format )

        # convert intervals to intersectors
        for contig in e.keys():
            intersector = bx.intervals.intersection.Intersecter()
            if self.mWithValues or self.mWithRecords:
                for start, end, value in e[contig]:
                    intersector.add_interval( bx.intervals.Interval(start,end,value=value) )
            else:
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
class CounterOverlapTranscripts(CounterOverlap):
    """count overlap with segments in another file.

    Nover1 and nover2 now count "transcripts".
    """

    mHeaderTemplate = ( "ngenes", "ntranscripts", "nexons", "nbases", "pover1", "pover2" )

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

    mHeaderTemplate = ["cov_%s" % x for x in Stats.Summary().getHeaders()+ ("covered", "values",) ]

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

    mHeader = [ "is_known", "is_unknown", "is_ambiguous", "is_pc", "is_pseudo", "is_npc", "is_utr", "is_intronic", "is_assoc", "is_intergenic" ] 

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

    def __init__(self, filename_gff, *args, **kwargs ):

        Counter.__init__(self, *args, **kwargs )

        E.info( "loading data from %s" % (filename_gff) )
            
        gffs = []
        infile = open( filename_gff, "r")     
        for g in GTF.iterator(infile):
            gffs.append( g )

        E.info( "loaded data from %s" % (filename_gff) )

        self.mCounters = {}
        self.mKeys = []
        for source in self.sources:             
            for feature in self.features:
                key = "%s:%s" % (source, feature )
                self.mKeys.append( key )
                self.mCounters[key] = CounterOverlap( gffs, 
                                                      source = source, 
                                                      feature = feature, 
                                                      fasta = self.mFasta,
                                                      options = self.mOptions )
                
    def count(self):
        
        for key in self.mKeys:
            self.mCounters[key](self.mGFFs)

        def s_min( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdMinCoverage

        def s_excl( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) < (100 - self.mThresholdMinCoverage)

        def s_full( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdFullCoverage

        def s_some( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdSomeCoverage

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

    mHeader = [ "is_cds", "is_utr", "is_upstream", "is_downstream", "is_intronic", "is_intergenic", "is_flank", "is_ambiguous" ]

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
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdMinCoverage

        def s_excl( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) < (100 - self.mThresholdMinCoverage)

        def s_full( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdFullCoverage

        def s_some( *args ):
            return sum( [ self.mCounters[x].mPOverlap1 for x in args ] ) >= self.mThresholdSomeCoverage

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

    mHeader = ( "nover_exonic", "nover_intronic", "nover_external", "length" )

    ## Do not save value for intervals
    mWithValues = False

    ## do not save records for intervals
    mWithRecords = False

    def __init__(self, filename_gff, *args, **kwargs):

        Counter.__init__(self, *args, **kwargs )

        source, feature = None, "CDS"
        
        e = readIntervalsFromGFF( filename_gff, 
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

    mHeaderTemplate = ( "distance", "id", "dist5", "strand5", "id5", "dist3", "strand3", "id3" )

    ## do not save value for intervals
    mWithValues = False

    ## save records for intervals in order to get strand information
    mWithRecords = True

    ## also compute overlap
    mWithOverlap = False

    def __init__(self, filename_gff, source, feature, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs )

        if feature and source:
            self.mHeader = [ "%s:%s:%s" % (x, source, feature) for x in self.mHeaderTemplate ]
        elif feature:
            self.mHeader = [ "%s:%s" % (x, feature) for x in self.mHeaderTemplate ]
        elif source:
            self.mHeader = [ "%s:%s" % (x, source) for x in self.mHeaderTemplate ]
        else:
            self.mHeader = self.mHeaderTemplate

        e = self.readIntervals( filename_gff, source, feature )

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

    mHeaderTemplate = ( "closest_id", "closest_dist", "closest_strand", "id5", "dist5", "strand5", "id3", "dist3", "strand3", "min5", "min3", "amin5", "amin3" )

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


    mHeaderTemplate = ( "closest_id", "closest_dist", "closest_strand", "id5", "dist5", "id3", "dist3", "is_overlap" )

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

    mHeaderTemplate = [ "splice_%s" % x for x in ( "total", "found", "missed", "perfect", "partial", "incomplete", "exon_skipping" ) ]

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

    mHeaderTemplate = Stats.Summary().getHeaders() + ( "length", "lengths", "values", ) 

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

    mHeaderTemplate = ( "ids", "dists" ) 

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

    mHeaderTemplate = ("status", "nterritories", "territories", )

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

    mHeader = (Stats.Summary().getHeaders() + ("values",) )
    
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
    '''
    
    mHeader = ("length", "pcovered",) + Stats.Summary().getHeaders()
    
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
        counts = numpy.zeros( length )
        nreads = 0
        contig = self.getContig()

        l = 0
        for start, end in segments:
            if end - start > self.mMaxLength: return []
            offset = start - l
            for samfile in self.mBamFiles:
                for read in samfile.fetch( contig, start, end ):
                    nreads += 1
                    rstart = max( 0, read.pos - offset )
                    rend = min( length, read.pos - offset + read.rlen ) 
                    counts[ rstart:rend ] += 1
            l += end - start

        self.mTotalLength = length
        counts = counts[ counts > 0]
        self.mResult = counts
        self.mCovered = len(counts)

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

    parser.add_option("-f", "--filename-gff", dest="filename_gff", type="string",
                      help="filename with extra gff file (for counter: overlap) [default=%default]."  )

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
                      choices=("length", "splice", "composition-na", "overlap", 
                               "classifier", "classifier-chipseq",
                               "overlap-transcripts",
                               "read-coverage",
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

    parser.set_defaults(
        genome_file = None,
        reporter = "genes",
        with_values = True,
        sections = [],
        counters = [],
        filename_gff = None,
        filename_format = "gtf",
        gff_features = [],
        gff_sources = [],
        add_gtf_source = False,
        proximal_distance = 10000, 
        bam_files = None,
        )

    (options, args) = E.Start( parser )

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
        bamfiles = []
        for bamfile in options.bam_files.split(","):
            bamfiles.append( pysam.Samfile(bamfile, "rb" ) )
    else:
        bamfiles = None

    counters = []

    if not options.sections: 
        E.info( "counters will use the default section (exons)" )
        options.sections.append( None)

    if not options.gff_sources: options.gff_sources.append( None )
    if not options.gff_features: options.gff_features.append( None )

    for c in options.counters:
        if c == "position":
            for section in options.sections:
                counters.append( CounterPosition( section = section, options = options ) )
        elif c == "length":
            for section in options.sections:
                counters.append( CounterLengths( section = section, options = options ) )
        elif c == "splice":
            counters.append( CounterSpliceSites( fasta=fasta ) )
        elif c == "quality":
            counters.append( CounterQuality( fasta=quality ) )
        elif c == "overrun":
            counters.append( CounterOverrun( filename_gff = options.filename_gff,
                                             options = options ) )
        elif c == "read-coverage":
            counters.append( CounterReadCoverage( bamfiles,
                                                  options = options ) )

        elif c == "splice-comparison":
            counters.append( CounterSpliceSiteComparison( fasta=fasta, 
                                                          filename_gff = options.filename_gff,
                                                          feature=None, 
                                                          source=None, 
                                                          options=options ) )
        elif c == "composition-na":
            for section in options.sections:
                counters.append( CounterCompositionNucleotides( fasta=fasta,
                                                                section = section,
                                                                options = options ) )

        elif c in ( "overlap", "overlap-transcripts", 
                    "proximity", "proximity-exclusive", "proximity-lengthmatched",
                    "neighbours",
                    "territories", 
                    "distance", "distance-genes", "distance-tss",
                    "coverage" ):
            if c == "overlap":
                template = CounterOverlap
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
                                                   options = options) )

        elif c == "classifier":
            counters.append( Classifier( filename_gff = options.filename_gff,
                                         fasta = fasta,
                                         options = options) )
        elif c == "classifier-chipseq":
            counters.append( ClassifierChIPSeq( filename_gff = options.filename_gff,
                                                fasta = fasta,
                                                options = options) )

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
        for counter in counters: counter(gffs)
        options.stdout.write( "\t".join( \
                fheader(gffs) +\
                ffields(gffs) +\
                [str( counter ) for counter in counters] ) + "\n" )

    E.Stop()
