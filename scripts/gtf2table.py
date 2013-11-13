'''
gtf2table.py - annotate genes/transcripts
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Genesets GTF Annotation

Purpose
-------

compute properties of genes or transcripts given in a :term:`gtf` 
formatted file and output them in tabular format. 

Quantities can be either computed per gene (all exons across all transcripts)
or per transcript. The input needs to be sorted accordingly.

The following methods (see option ``--counter``) are available:

bigwig-counts 
   collect density values from bigwig file and output
   summary statistics and percentage of bases covered (``pcovered``)
   by value in bigwig file.
   
binding-pattern 
   given a list of intervals, determine the binding
   pattern within and surrounding the gene. For each gene, intervals
   overlapping the CDS, introns, UTRs and the flank are collected and
   recorded. The binding is summarized with a binding pattern, a
   binary pattern indicating overlap/no overlap with 5' flank, 5' UTR,
   CDS, Introns, 3' UTR, 3' flank.

classifier
   classify transcripts according to genomic annotation 
   Requires a :term:`gff` file with genomic annotations
   (see :doc:`gtf2gff`.)

classifier-rnaseq
   classify rnaseq transcripts with respect to a reference
   geneset. Requires a :term:`gtf` file with a reference
   gene set.

classifier-polii classify 
   according to PolII transcripts. A
   gene/transcript is transcribed, if it is covered by large PolII
   intervals over 80% of its length. A gene/transript is primed if its
   promotor/UTR is covered by 50% of its length, while the rest of the
   gene body isn't.

composition-na
   output nucleotide composition of gene

composition-cgp
   output cpg composition of gene

coverage
   compute nucleotide coverage of input with segments in another file.
   The values are output in 5' to 3' order for each nucleotide. Requires
   a second :term:`gff` formatted file with features to cover.

distance
   compute distance of genes to features in a second file. Requires
   a second :term:`gff` formatted file with transcripts. The strand
   information of the features is ignored.

distance-genes
   compute distance of genes to genes in a second file. Requires a 
   second :term:`gtf` formatted file with genes. The counter distinguishes
   a variety of cases (closest upstream/downstream).

distance-tss
   compute distance of genes to transcription start sites. Requires a
   second :term:`gtf` formatted file with genes.

length
   output exon length summary of gene.

neighbours
    output features in second stream that are in proximity to genes
    in input. Requires a :term:`gtf` formatted file with genes.

overlap
    compute overlap of genes in input with features in second stream.
    Requires a :term:`gff` formatted file with gene territories.

overlap-stranded
    count overlap with genomic features in second another file. Outputs
    the number of overlapping exons. Records the direction of overlap
    (sense/antisense). Requires a :term:`gff` formatted file with 
    features.

overlap-transcripts
    count overlap of genes with transcripts in another set.
    Requires a :term:`gtf` formatted file.

overrun
   output intron overrun, exons in the input gene set extending
   into the introns of a reference gene set. Requries a :term:`gtf`
   formatted file with a reference gene set.

position
   output genomic coordinates of gene

proximity
   report summary stats (lengths,values) of features in proximity to 
   genes input gene set. Requires a :term:`gff` formatted file with
   genomic features.

proximity-exclusive
   as proximity, but exclude any ranges overlapping the gene set.

proximity-lengthmatched
  as proximity-exclusive, but length-match features with genes.

quality
   output base-quality information summary of gene. Needs quality scores.

read-coverage
   output read coverage summary statistics of gene

read-extension

read-counts
   count number of reads overlapping a gene or transcript. Counts
   uniquely by read name and counts duplicate and non-duplicate (reads)
   separately.

splice
   output splicing summary of gene

splice-comparison

territories

Usage
-----

Example::

   python gtf2table.py --counter=length < geneset.gtf > geneset.tsv

Type::

   python gtf2table.py --help

for command line help.

Command line options
--------------------

'''

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

try:
    import alignlib_lite
except ImportError:
    pass

import numpy
import CGAT.IndexedGenome as IndexedGenome
import pysam

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _gtf2table
except ImportError:
    import CGAT._gtf2table as _gtf2table

def readIntervalsFromGFF( filename_gff, source, feature, 
			  with_values = False, with_records = False, fasta = None, 
			  merge_genes = False, format = "gtf", use_strand = False ):
    """read intervals from a file or list.
    """

    assert not (with_values and with_records), "both with_values and with_records are true."

    ninput = 0

    if format == None:
        if type(filename_gff) == types.StringType:
            fn = filename_gff
            if fn.endswith( ".gtf" ) or fn.endswith( ".gtf.gz"):
                format = "gtf"
            elif fn.endswith( ".gff" ) or fn.endswith( ".gff.gz"):
                format = "gff"
            elif fn.endswith( ".bed" ) or fn.endswith( ".bed.gz"):
                format = "bed"
        else:
            format = "gff"

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

        gff_iterator = GTF.iterator_filtered( iterator_gff,
                                              feature = feature,
                                              source = source )

        if format == "gtf":
            e = GTF.readAsIntervals( gff_iterator, 
                                     with_values = with_values, 
                                     with_records = with_records,
                                     merge_genes = merge_genes,
                                     use_strand = use_strand )
        elif format == "gff":
            e = GTF.readAsIntervals( gff_iterator, 
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
                e[bed.contig].append( (bed.start,bed.end,bed.fields[0]) )
        elif with_records:
            for bed in iterator:
                ninput += 1
                bed.gene_id = bed.fields[0]
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

class CounterIntronsExons(_gtf2table.Counter):
    """count number of introns and exons.
    """

    header = ( "ntranscripts", "nexons", "nintrons", )
    def __init__(self, *args, **kwargs):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):
        segments = self.getSegments()
        self.mNTranscripts = len( set( [x.transcript_id for x in self.mGFFs] ) )
        self.mNExons = len(segments)
        self.mNIntrons = len(segments) - 1

    def __str__(self):
        return "\t".join( (str(self.mNTranscripts), str(self.mNSegments), str(self.mNIntrons)) ) 

class CounterPosition(_gtf2table.Counter):
    """output the position of the transcript."""
    header = ( "contig", "strand", "start", "end" )

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

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
class CounterLengths(_gtf2table.Counter):
    header = Stats.Summary().getHeaders()

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):
        segments = self.getSegments()
        self.result = Stats.Summary( [x[1] - x[0] for x in segments ], mode="int", allow_empty = True )

    def __str__(self):
        return str(self.result)

##----------------------------------------------------------------
class CounterSpliceSites(_gtf2table.Counter):

    mIntronTypes = ( ( "U2-GT/AG", "GT", "AG"),
                     ( "U2-nc-GC/AG", "GC", "AG"),
                     ( "U12-AT/AC", "AT", "AC") )

    mNames = [x[0] for x in mIntronTypes ] + ["unknown"]
    header = ["%s" % x for x in mNames ]

    mCheckBothStrands = True

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):

        self.mCounts = dict( [( x[0], 0) for x in self.mIntronTypes ] )
        self.mCounts["unknown"] = 0

        introns = self.getIntrons()
        contig = self.getContig()
        strand = self.getStrand()
        
        for start, end in introns:
            s = self.fasta.getSequence( contig, "+", start, end )
            
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
class CounterCompositionNucleotides(_gtf2table.Counter):
    header = SequenceProperties.SequencePropertiesNA().getHeaders() 

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        s = self.getSequence( ee )
        self.result = SequenceProperties.SequencePropertiesNA()
        self.result.loadSequence( s )

    def __str__(self):
        return str(self.result)

##-----------------------------------------------------------------------------------
class CounterCompositionCpG(_gtf2table.Counter):
    '''compute CpG frequencies as well as nucleotide frequencies.

    Note that CpG density is calculated across the merged exons
    of a transcript. Thus, there might be difference between the CpG 
    on a genomic level and on the transrcipt level depending on how
    many genomic CpG are lost across an intron-exon boundary or how
    many transcript CpG are created by exon fusion.
    '''

    header = SequenceProperties.SequencePropertiesCpg().getHeaders() 

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        s = self.getSequence( ee )
        self.result = SequenceProperties.SequencePropertiesCpg()
        self.result.loadSequence( s )

    def __str__(self):
        return str(self.result)

##-----------------------------------------------------------------------------------
class CounterOverlap(_gtf2table.Counter):
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
        _gtf2table.Counter.__init__(self, *args, **kwargs )

        if feature and source:
            self.header = [ "%s:%s:%s" % (x, source, feature) for x in self.headerTemplate ]
        elif feature:
            self.header = [ "%s:%s" % (x, feature) for x in self.headerTemplate ]
        elif source:
            self.header = [ "%s:%s" % (x, source) for x in self.headerTemplate ]
        else:
            self.header = self.headerTemplate

        if len(filename_gff) != 1:
            raise ValueError("expected one gff file" )
        
        e = readIntervalsFromGFF( filename_gff[0], 
				  source, 
				  feature, 
                                  self.mWithValues, 
				  self.mWithRecords, 
                                  self.fasta, 
                                  format = self.options.filename_format,
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
        if self.fasta: 
            contig = self.fasta.getToken( contig)
            
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
        if self.fasta: contig = self.fasta.getToken( contig)
            
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

    headerTemplate = ( "ngenes", "ntranscripts", 
		       "nexons", "nbases", 
		       "pover1", "pover2" )

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
        if self.fasta: 
            contig = self.fasta.getToken( contig)
            
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
    """compute nucleotide coverage of input with segments in another file.
    The values are output in 5' to 3' order for each nucleotide.
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
        if self.fasta: 
            contig = self.fasta.getToken( contig)
            
        n = 0
        segments = self.getSegments()
        segments.sort()

        map_genome2transcript = alignlib_lite.py_makeAlignmentBlocks()

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
class Classifier(_gtf2table.Counter):
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

        _gtf2table.Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 1:
            raise ValueError("expected one gff file, received %i" % len(filename_gff) )

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
                                               fasta = self.fasta,
                                               options = self.options )
                
    def count(self):
        
        for key in self.mKeys:
            self.mCounters[key].update(self.mGFFs)

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
        h = [_gtf2table.Counter.getHeader( self ) ]

        for key in self.mKeys:
            h.append( self.mCounters[key].getHeader() )

        return "\t".join( h )


##-----------------------------------------------------------------------------------
class ClassifierRNASeq(_gtf2table.Counter):
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
                 ( 's', "utr5", ),
                 ( 's', "utr3", ),
                 ( 's', "intronic", ),
		 ( 's', "flank5", ),
                 ( 's', "flank3", ),
                 ( 'n', "complete", ),
                 ( 'n', "fragment", ),
                 ( 'n', "exon-fragment", ),
                 ( 'n', "extended-fragment", ),
                 ( 'n', "extension", ),
                 ( 'n', "alternative", ),
                 ( 'n', "unknown", ),
                 ( 'n', "utr5", ),
                 ( 'n', "utr3", ),
		 ( 'n', "intronic", ),
                 ( 'n', "flank5", ),
                 ( 'n', "flank3", ),
                 ( 'a', "complete", ),
                 ( 'a', "exon-fragment", ),
                 ( 'a', "fragment", ),
                 ( 'a', "extended-fragment", ),
                 ( 'a', "extension", ),
                 ( 'a', "alternative", ),
                 ( 'a', "unknown", ),
		 ( 'a', "utr5", ),
                 ( 'a', "utr3", ),
		 ( 'a', "intronic", ),
		 ( 'a', "flank5", ),
                 ( 'a', "flank3", ),
                 ( 'n', "flank", ),
                 ( 'n', "intergenic" ) )

    def __init__(self, filename_gff, *args, **kwargs ):

        _gtf2table.Counter.__init__(self, *args, **kwargs )

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
class ClassifierRNASeqSplicing(_gtf2table.Counter):
    """This is IMSs new style transcript classifier. It aims to give classifications
    that make more sense to biologists involved in splicing by using familiar catagories.

    classify RNASeq transcripts based on a reference annotation.

    Transcripts are classified by checking overlap with all known
    transcripts.

    If multiple transcripts overlap, select the one that is best matching.

    +--------------------+------------------------------------------------------------+
    |*Labels*            |*Contents*                                                  |
    +--------------------+------------------------------------------------------------+
    |complete            |Intron chains in transcripts are identical                  |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |retained-intron     |All introns in test transcript are present in an existing   |
    |                    |transcript , but not all introns between the start and end  |
    |                    |co-ordinates of the test transcript are included            |
    +--------------------+------------------------------------------------------------+
    |extended            |Intron chain in the known transcripts is the same as the    |
    |                    |predicted transcirpt between the strat and end of the       |
    |                    |predicted transcript, but the first and/or last exon extends|
    |                    |beyond the boundaries of the first/last included exon of the|
    |                    |known transcirpt                                            |
    +--------------------+------------------------------------------------------------+
    |alternate-5prime    |All of the exon boundaries of the known transcirpt between  |
    |alternate-3prime    |the start and end of the predicted transcirpt are shared    |
    |novel-exon          |with the predicted transcript but the predicted transcript  |
    |                    |contains extra boundaries                                   |
    +--------------------+------------------------------------------------------------+
    |skipped-exon        |All of the exon boundaries in the predicted transcript are  |
    |                    |present in the known transcript, but the known transcript   |
    |                    |contains boundaires that the predicted transcript does not  |
    +--------------------+------------------------------------------------------------+
    |alternative         |Known and predicted transcripts share at least one exon     |
    |                    |boundary                                                    |
    |                    |                                                            |
    +--------------------+------------------------------------------------------------+
    |unknown             |``predicted`` and ``known`` transcript overlap, but do not  |
	|                    |fall in any of the above categories.                        |
	|                    |                                                            |
	+--------------------+------------------------------------------------------------+

    Any of these forms can also exist as a -fragment. This means that the comparison
	known and the predicted transcript did not make the criteria for a class, but the 
	comparison of the predicted and the part of the known transcirpt bewteen the start and
	end of the known transcript do. 

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
    tolerance = 1

    # size of flanking region
    flank = 5000

    # priority of classifications to select best match
    priority = ( ( 's', "complete", ),
		 ( 's', "fragment", ),
		 ( 's', "extended", ),
		 ( 's', "extended-fragment", ),
		 ( 's', "retained-intron", ),
		 ( 's', "retained-intron-fragment", ),
		 ( 's', "exon-boundary-change", ),
		 ( 's', "alternate-5prime", ),
		 ( 's', "alternate-3prime", ),
		 ( 's', "novel-intron",),
		 ( 's', "skipped-exon" ,),
		 ( 's', "novel-exon", ),
		 ( 's', "alternate-exon", ),
		 ( 's', "exon-boundary-change-fragment", ),
		 ( 's', "alternate-5prime-fragment", ),
		 ( 's', "alternate-3prime-fragment", ),
		 ( 's', "novel-intron-fragment",),
		 ( 's', "skipped-exon-fragment", ),
		 ( 's', "novel-exon-fragment", ),
		 ( 's', "alternate-exon-fragment", ),
		 ( 's', "alternative", ),
		 ( 's', "unknown", ),
		 ( 's', "utr5", ),
		 ( 's', "utr3", ),
		 ( 's', "intronic", ),
		 ( 's', "flank5", ),
		 ( 's', "flank3", ),
		 ( 's', "intergenic", ),
		 ( 'n', "complete", ),
		 ( 'n', "fragment", ),
		 ( 'n', "extended", ),
		 ( 'n', "extended-fragment", ),
		 ( 'n', "retained-intron", ),
		 ( 'n', "retained-intron-fragment", ),
		 ( 'n', "exon-boundary-change", ),
		 ( 'n', "alternate-5prime", ),
		 ( 'n', "alternate-3prime", ),
		 ( 'n', "novel-intron",),
		 ( 'n', "skipped-exon" ,),
		 ( 'n', "novel-exon", ),
		 ( 'n', "alternate-exon", ),
		 ( 'n', "exon-boundary-change-fragment", ),
		 ( 'n', "alternate-5prime-fragment", ),
		 ( 'n', "alternate-3prime-fragment", ),
		 ( 'n', "novel-intron-fragment",),
		 ( 'n', "skipped-exon-fragment", ),
		 ( 'n', "novel-exon-fragment", ),
		 ( 'n', "alternate-exon-fragment", ),
		 ( 'n', "alternative", ),
		 ( 'n', "unknown", ),
		 ( 'n', "utr5", ),
		 ( 'n', "utr3", ),
		 ( 'n', "intronic", ),
		 ( 'n', "flank5", ),
		 ( 'n', "flank3", ),
		 ( 'n', "flank", ),
		 ( 'n', "intergenic", ),
		 ( 'a', "complete", ),
		 ( 'a', "fragment", ),
		 ( 'a', "extended", ),
		 ( 'a', "extended-fragment", ),
		 ( 'a', "retained-intron", ),
		 ( 'a', "retained-intron-fragment", ),
		 ( 'a', "exon-boundary-change", ),
		 ( 'a', "alternate-5prime", ),
		 ( 'a', "alternate-3prime", ),
		 ( 'a', "novel-intron",),
		 ( 'a', "skipped-exon" ,),
		 ( 'a', "novel-exon", ),
		 ( 'a', "alternate-exon", ),
		 ( 'a', "exon-boundary-change-fragment", ),
		 ( 'a', "alternate-5prime-fragment", ),
		 ( 'a', "alternate-3prime-fragment", ),
		 ( 'a', "novel-intron-fragment",),
		 ( 'a', "skipped-exon-fragment", ),
		 ( 'a', "novel-exon-fragment", ),
		 ( 'a', "alternate-exon-fragment", ),
		 ( 'a', "alternative", ),
		 ( 'a', "unknown", ),
		 ( 'a', "utr5", ),
		 ( 'a', "utr3", ),
		 ( 'a', "intronic", ),
		 ( 'a', "flank5", ),
		 ( 'a', "flank3", ),
		 ( 'a', "intergenic", )	 )

    def __init__(self, filename_gff, *args, **kwargs ):

        _gtf2table.Counter.__init__(self, *args, **kwargs )

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
            t = [x for x in t if x.feature == "exon"]
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
                approx_structure = False
            else:
				# IMS structure is an exact match if and only if all intron boundaries in both exons and 
				# transcripts are shared. This allows extensions of the 3' or 5' UTRs.
                matched_structure = len(shared_introns) == len(introns)
                approx_structure = len(shared_boundaries) > 0
             

            included_exons = [x for x in transcript_exons if x[1] > start and x[0] < end ]
            included_transcript_introns = [x for x in transcript_introns if x[0] > start and x[1] <= end]
            included_boundaries = sorted([x for x in transcript_boundaries if x[1] >= exons[0][1]  and x[0] <= exons[-1][0]  ])
            shared_included_boundaries = Intervals.intersect( boundaries, included_boundaries )
            
			# If there is a matched structure, i.e. all of the introns in the gene model are in an existing gene
			# model, then either we have a fragment of the complete transcript or some sort of retain intron 
			# (or both)
            #if matched_structure and abs(overlap_exons - transcript_lexons) < tolerance:
            if matched_structure and len(shared_introns) == len(transcript_introns):
                cls = "complete"
            elif matched_structure:
                shared_included = [x for x in introns if x in included_transcript_introns]

                cls = []

                if not len(shared_included) == len(included_transcript_introns):   
                    cls.append("retained-intron")
                elif start < included_exons[0][0] or end > included_exons[-1][1]:
                    cls.append("extended")
                    
                if (Intervals.calculateOverlap([boundaries[0]], transcript_boundaries[1:]) > 0 or
                    Intervals.calculateOverlap([boundaries[-1]], transcript_boundaries[:-2]) > 0):
                
                    cls.append("fragment")
                
                cls = "-".join(cls)         
                
            #elif approx_structure and len(exons) > 1 and abs(overlap_exons - transcript_lexons) < tolerance:
            #    cls = "complete"

        
            elif approx_structure and (len(shared_included_boundaries) == len(included_boundaries) and
                  len(shared_included_boundaries) < len(boundaries)):
                #there are boundaries present in the test model that are not present in the transcript 
                #therefore there is a novel exon. Could be fragement or complete.
                
                #print exons[0]
		#print transcript_exons
		#print Intervals.calculateOverlap([exons[0]],transcript_exons)

                if (Intervals.calculateOverlap([exons[0]],transcript_exons) == 0):
                    #print "hello"
                    cls = "alternate-5prime"
                        
                    if len(Intervals.intersect([boundaries[-1]], transcript_boundaries[0:-1])) > 0 :
                    # is a fragment if the final boundary is in the transcript but isn't its final one.
                        cls = cls + "-fragment"
                    if not len(boundaries) == len(shared_included_boundaries) + 1:
                        cls = "alternative"
                elif (Intervals.calculateOverlap([exons[-1]],transcript_exons) == 0):
                    cls = "alternate-3prime"
                    if len(Intervals.intersect([boundaries[-1]], transcript_boundaries[2:])) > 0:
                        #is a fragment if the first boundary is in the transcript but is the first one.
                        cls = cls + "-fragment"
                    if not len(boundaries) == len(shared_included_boundaries) + 1:
                        cls = "alternative"
                elif (len(shared_boundaries) == len(transcript_boundaries) and 
		      len(shared_boundaries) < len (boundaries) ):
                    
                    novel_exons = [exon for exon in exons if Intervals.calculateOverlap([exon],transcript_exons) ==0]
                    if len(novel_exons) == 1:
                        cls = "novel-exon"
                    elif len(novel_exons) == 0:
                        cls = "novel-intron"
                    else:
                        cls = "alternative"
                else:
                    novel_exons = [exon for exon in exons if Intervals.calculateOverlap([exon],transcript_exons) ==0]
                    if len(novel_exons) == 1:
                        cls = "novel-exon-fragment"
                    elif len(novel_exons) == 0:
                        cls = "novel-intron-fragment"
                    else:
                        cls = "alternative"
                
            elif approx_structure and (len(shared_included_boundaries) == len (boundaries) and 
                  len(shared_included_boundaries) < len (included_boundaries)):
                #there are boundaries in the transcript that are not in the test model. There for an exon has been 
                #skipped. Could be fragment or complete.
                if (Intervals.calculateOverlap([boundaries[0]],[transcript_boundaries[0]]) > 0 and
                    Intervals.calculateOverlap([boundaries[-1]],[transcript_boundaries[-1]]) > 0):
                    cls = "skipped-exon"
                else:
                    cls = "skipped-exon-fragment"
                
            elif approx_structure and (len(Intervals.combine(exons + included_exons)) == len(exons) and
                  len(Intervals.combine(exons + included_exons)) == len(included_exons)):
                #same exons exist in each (where two exons are the "same" exon if they overlap with each other, but
                #only each other).
                cls = "exon-boundary-change"
                if not (len(Intervals.combine(exons + transcript_exons)) == len (exons) and
                        len(Intervals.combine(exons + transcript_exons)) == len (transcript_exons)):
                    cls = cls +"-fragment"

            elif approx_structure and (len(Intervals.combine(exons + included_exons)) == len(exons) +1 and
                                       len(Intervals.combine(exons + included_exons)) == len(included_exons)+1):
                #there is exactly one exon in test model that is not in transcript and vice versa
                cls = "alternate-"
                if (Intervals.calculateOverlap([exons[0]], transcript_exons) == 0):
                    if len(boundaries) == len(shared_boundaries) + 1:
                        cls = cls + "5prime"
                    else:
                        cls = "alternative"

                elif (Intervals.calculateOverlap([exons[-1]], transcript_exons) == 0):
                    if len(boundaries) == len(shared_boundaries) +1:
                        cls = cls + "3prime"
                    else:
                        cls = "alternative"
                elif len(boundaries) == len(shared_boundaries) +2:
                    cls = cls + "exon"
                else:
                    cls = "alternative"
                    

                if len(Intervals.intersect([boundaries[-1]], transcript_boundaries[0:-1])) > 0 :
                    # is a fragment if the final boundary is in the transcript but isn't its final one.
                    cls = cls + "-fragment"
                elif len(Intervals.intersect([boundaries[-1]], transcript_boundaries[2:])) > 0:
                    #is a fragment if the first boundary is in the transcript but is the first one.
                    cls = cls + "-fragment"
                if cls == "alternative-fragment":
                    cls = "alternative"

            elif approx_structure:
                cls = "alternative"

            elif overlap_exons == lexons:
                cls = "fragment"
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
                if overlap_utr5 > 0:
                    cls = "utr5"
                elif overlap_utr3 > 0:
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
class ClassifierIntervals(CounterOverlap):
    """classify transcripts based on a list of intervals.
    """

    header = [ "is_known", "is_unknown", "is_ambiguous", 
               "is_pc", "is_pseudo", "is_npc", "is_utr", 
               "is_intronic", "is_assoc", "is_intergenic" ] 

    # do not use strand
    mUseStrand = False

    def __init__(self, *args, **kwargs ):

        CounterOverlap.__init__(self, *args, **kwargs )

##-----------------------------------------------------------------------------------
class ClassifierPolII(ClassifierIntervals):
    """classify transcripts based on PolII annotation.
    
    The input is a list of PolII intervals.
    
    An interval is classified as:
    
    is_transcribed: 
        the longest interval covers *threshold_min_coverage* of the 
        gene body.

    """

    headerTemplate = [ "is_transcribed", "is_primed", "noverlap", 
                       "largest_size", "largest_overlap", "largest_coverage",
                       "promotor_overlap", "promotor_coverage" ]

    # minimum coverage of a transcript to accept
    threshold_min_coverage_transcript = 80

    # minimum coverage of a promotor to be accepted
    threshold_min_coverage_promotor = 50

    # size of promotor - 1kb upstream of tss
    promotor_upstream = 200
    
    # size of promotor - 1kb downstream of tss
    promotor_downstream = 200

    def count(self):

        self.is_transcribed = False
        self.is_primed = False
        self.noverlap = 0
        self.largest = 0
        self.largest_interval = 0
        self.largest_overlap = 0
        self.largest_coverage = 0
        self.promotor_overlap = 0
        self.promotor_coverage = 0

        contig = self.getContig()

        if contig not in self.mIntersectors: return
        
        # count overlap over full gene body
        segments = self.getExons()
        start, end = segments[0][0], segments[-1][1]
        
        promotor_max = max( self.promotor_upstream, self.promotor_downstream )

        r = list(self.mIntersectors[contig].find( start - promotor_max, end + promotor_max ))
        if len(r) == 0: return
        self.noverlap = len(r)

        r.sort( key = lambda x: x.end-x.start )

        largest_start, largest_end = r[-1].start, r[-1].end
        self.largest_interval = largest_end - largest_start
        self.largest_overlap = min( largest_end, end) - max(largest_start, start)
        self.largest_coverage = 100.0 * float( self.largest_overlap ) / (end - start)

        strand = self.getStrand()
        if strand == "+":
            promotor_start, promotor_end = start - self.promotor_upstream, start + self.promotor_downstream
        else:
            promotor_start, promotor_end = end - self.promotor_downstream, end + self.promotor_upstream
            
        self.promotor_overlap = Intervals.calculateOverlap( [ (x.start, x.end) for x in r ],
                                                            [ (promotor_start, promotor_end), ] )

        self.promotor_coverage = 100.0 * float( self.promotor_overlap ) / (promotor_end - promotor_start)

        if self.largest_coverage >= self.threshold_min_coverage_transcript:
            self.is_transcribed = True
        elif self.promotor_coverage >= self.threshold_min_coverage_promotor:
            self.is_primed = True

    def __str__(self):

        def to( v ):
            if v: return "1" 
            else: return "0"

        h = [ to(x) for x in (self.is_transcribed, 
                              self.is_primed ) ]
        h.extend( [ "%i" % self.noverlap,
                    "%i" % self.largest_interval,
                    "%i" % self.largest_overlap,
                    "%5.2f" % self.largest_coverage,
                    "%i" % self.promotor_overlap,
                    "%5.2f" % self.promotor_coverage ] )
        
        return "\t".join( h )

##-----------------------------------------------------------------------------------
class CounterBindingPattern(CounterOverlap):
    """compute overlaps between gene and various tracks given
    by one or more gff files.

    pattern
       a binding pattern.
    overlap
       number of intervals overlapping the genic region
    
    Single exon genes have no intron.
    Genes with one intron have only a first intron.
    Genes with two introns have a first and a middle intron.

    This method also computes the expected number of times 
    an interval would overlap with a certain feature. The 
    probability is only exact if it is assumed that the 
    intervals are point features.
    """

    headerTemplate = [ "pattern", "overlap" ]


    # do not use strand
    mUseStrand = False

    mWithValues = False
    
    mWithRecords = False

    # size of flank, make sure the headers are updated above.
    flank = 10000

    # number of bins in the flank
    flank_bins = 10

    def __init__(self, *args, **kwargs ):
        CounterOverlap.__init__(self, *args, **kwargs )

        increment = self.flank // self.flank_bins
	
        self.headerTemplate.extend( 
            [ "%s_%s" % (x,y) for x,y in itertools.product( 
                    ("cds", "first_exon", "exon", "utr5", "utr3", "first_intron", "middle_intron", "last_intron", "intron" ) +\
                        tuple( ["flank5_%05i" % x for x in range(0, self.flank, increment) ] ) +\
                        tuple( ["flank3_%05i" % x for x in range(0, self.flank, increment) ] ),
                    ("overlap", "poverlap" ) ) ] )

    def count( self ):

        self.overlap_intron = 0
        self.overlap_exon = 0
        self.overlap_utr5 = 0
        self.overlap_utr3 = 0
        self.overlap_cds = 0

        self.overlap_first_intron = 0
        self.overlap_last_intron = 0
        self.overlap_middle_intron = 0
        self.overlap_first_exon = 0
        self.overlap_flank5 = [0] * self.flank_bins
        self.overlap_flank3 = [0] * self.flank_bins

        self.poverlap_intron = 0
        self.poverlap_exon = 0
        self.poverlap_utr5 = 0
        self.poverlap_utr3 = 0
        self.poverlap_cds = 0

        self.poverlap_first_intron = 0
        self.poverlap_last_intron = 0
        self.poverlap_middle_intron = 0
        self.poverlap_first_exon = 0
        self.poverlap_flank5 = [0] * self.flank_bins
        self.poverlap_flank3 = [0] * self.flank_bins
        self.overlap = 0

        # pattern takes only first bin for flank5 and flank3
        self.pattern = "0" * 6

        contig = self.getContig()
        strand = self.getStrand()
        if contig not in self.mIntersectors: return

        ######################################
        ## build sub-intervals to count
        utr5, utr3 = self.getUTRs()
        cds = self.getCDS()
        exons = self.getExons()
        introns = self.getIntrons()
        start, end = exons[0][0], exons[-1][1]

        ######################################
        ## get overlapping intervals
        extended_start, extended_end = start - self.flank, end + self.flank
        overlaps = list(self.mIntersectors[contig].find( extended_start, extended_end ))
        if len(overlaps) == 0 : return 
        intervals = [(x.start, x.end) for x in overlaps ]

        self.overlap = len(intervals )

        ######################################
        ## build special sets
        first_intron, first_exon = [], []
        last_intron, middle_intron = [], []
        flank_increment = self.flank // self.flank_bins

        # flank is ordered such that indices move away from the tss or tes
        if strand == "+":
            flank5 = [ [(x-flank_increment, x)] for x in range(start, start - self.flank, -flank_increment) ]
            flank3 = [ [(x,x+flank_increment)] for x in range(end, end + self.flank, flank_increment) ]
            if introns: 
                first_intron = [introns[0]]
                if len(introns) > 1:
                    last_intron = [introns[-1]]
                if len(introns) > 2:
                    middle_intron = [introns[len(introns) // 2 ] ] 
            first_exon = [exons[0]]
            del exons[0]
        else:
            flank3 = [ [(x-flank_increment, x)] for x in range(start, start - self.flank, -flank_increment) ]
            flank5 = [ [(x,x+flank_increment)] for x in range(end, end + self.flank, flank_increment) ]
            if introns: 
                first_intron = [introns[-1]]
                if len(introns) > 1:
                    last_intron = [introns[0]]
                if len(introns) > 2:
                    middle_intron = [introns[len(introns) // 2 ] ] 

            first_exon = [exons[-1]]
            del exons[-1]

        #######################################
        ## calculate overlap
        self.overlap_intron = Intervals.calculateOverlap( introns, intervals )
        self.overlap_exon = Intervals.calculateOverlap( exons, intervals )
        self.overlap_utr5 = Intervals.calculateOverlap( utr5, intervals )
        self.overlap_utr3 = Intervals.calculateOverlap( utr3, intervals )
        self.overlap_cds = Intervals.calculateOverlap( cds, intervals )
        self.overlap_first_exon = Intervals.calculateOverlap( first_exon, intervals )
        self.overlap_first_intron = Intervals.calculateOverlap( first_intron, intervals )
        self.overlap_last_intron = Intervals.calculateOverlap( last_intron, intervals )
        self.overlap_middle_intron = Intervals.calculateOverlap( middle_intron, intervals )

        for x, i in enumerate( flank5):
            self.overlap_flank5[x] = Intervals.calculateOverlap( i, intervals )
        for x, i in enumerate( flank3):
            self.overlap_flank3[x] = Intervals.calculateOverlap( i, intervals )
        
        #######################################
        ## calculate percent overlap
        pp = IOTools.prettyPercent
        self.poverlap_intron = pp( self.overlap_intron, Intervals.getLength( introns ), na = 0)
        self.poverlap_exon = pp(self.overlap_exon, Intervals.getLength( exons ), na = 0)
        self.poverlap_cds = pp(self.overlap_cds, Intervals.getLength( cds ), na = 0)
        self.poverlap_utr5 = pp( self.overlap_utr5, Intervals.getLength( utr5 ), na = 0)
        self.poverlap_utr3 = pp( self.overlap_utr3, Intervals.getLength( utr3 ), na = 0)
        self.poverlap_first_exon = pp( self.overlap_first_exon, Intervals.getLength( first_exon ), na = 0)
        self.poverlap_first_intron = pp( self.overlap_first_intron, Intervals.getLength( first_intron ), na = 0)
        self.poverlap_middle_intron = pp( self.overlap_middle_intron, Intervals.getLength( middle_intron ), na = 0)
        self.poverlap_last_intron = pp( self.overlap_last_intron, Intervals.getLength( last_intron ), na = 0)

        for x, v in enumerate( self.overlap_flank5 ):
            self.poverlap_flank5[x] = pp( v, Intervals.getLength( flank5[x] ), na = 0)
            
        for x, v in enumerate( self.overlap_flank3 ):
            self.poverlap_flank3[x] = pp( v, Intervals.getLength( flank3[x] ), na = 0)

        #######################################
        ## build binding pattern
        pattern = []
        for x in ( self.overlap_flank5[0], self.overlap_utr5, self.overlap_cds, 
                   self.overlap_intron, self.overlap_utr3, self.overlap_flank3[0]):
            if x: pattern.append("1")
            else: pattern.append("0")

        self.pattern = "".join(pattern)

    def __str__(self):
        
        data = [ self.pattern,
                 self.overlap,
                 self.overlap_cds,
                 self.poverlap_cds,
                 self.overlap_first_exon,
                 self.poverlap_first_exon,
                 self.overlap_exon,
                 self.poverlap_exon,
                 self.overlap_utr5,
                 self.poverlap_utr5,
                 self.overlap_utr3,
                 self.poverlap_utr3,
                 self.overlap_first_intron,
                 self.poverlap_first_intron,
                 self.overlap_middle_intron,
                 self.poverlap_middle_intron,
                 self.overlap_last_intron,
                 self.poverlap_last_intron,
                 self.overlap_intron,
                 self.poverlap_intron,
		 ]

        for x in range(len(self.overlap_flank5)):
            data.append( self.overlap_flank5[x] )
            data.append( self.poverlap_flank5[x] )

        for x in range(len(self.overlap_flank3)):
            data.append( self.overlap_flank3[x] )
            data.append( self.poverlap_flank3[x] )

        return "\t".join( map(str, data ) )

##-----------------------------------------------------------------------------------
class CounterOverrun(_gtf2table.Counter):
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

        _gtf2table.Counter.__init__(self, *args, **kwargs )

        if len(filename_gff) != 1:
            raise ValueError("expected only one gff file" )

        source, feature = None, "CDS"
        
        e = readIntervalsFromGFF( filename_gff[0], 
                                  source,
                                  feature, 
                                  with_values = self.mWithValues, 
                                  with_records = self.mWithRecords, 
                                  fasta = self.fasta,
                                  format = self.options.filename_format )

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
        if self.fasta: 
            contig = self.fasta.getToken( contig)
            
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
class CounterDistance(_gtf2table.Counter):
    """counter for computing the distance to features on either side 
    of a feature.

    The strand of a feature is taken into account.

    The columns output are:

    distance
       distance to closest feature in any direction
    id
       id of closest feature 
    dist5
       distance to closest feature in 5' direction
    strand5
       strand of feature in 5' direction
    id5
       id of feature in 5' direction
    dist3
       distance to closest feature in 3' direction
    strand3
       strand of feature in 3' direction
    id3
       id of feature in 5' direction

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
        _gtf2table.Counter.__init__(self, *args, **kwargs )

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
                                     self.mWithValues, self.mWithRecords, self.fasta,
                                     format = self.options.filename_format )

    def count(self):
        """find closest feature in 5' and 3' direction."""

        # collect overlapping segments
        segments = self.getSegments()
        contig = self.getContig()
        if self.fasta: contig = self.fasta.getToken( contig)

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
            E.warn( "unknown contig %s" % contig )
            return

        if self.mOverlaps:
            # stop getting distance if there is overlap
            self.mDistance = 0
            self.mData = [x[2] for x in self.mOverlaps]
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

    If the transcript overlaps a gene, the distance is set to 0. If it overlaps multiple genes,
    only one (the first) is output.
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
                                     self.mWithValues, self.mWithRecords, self.fasta,
                                     merge_genes = True,
                                     format = self.options.filename_format )

    def __str__(self):

        id5, id3, d3, d5, m3, m5 = "na", "na", "na" , "na", "na", "na"
        any3, any5 = "", ""

        if self.mOverlaps:
            closest_id = self.mData[0].gene_id
            closest_dist = 0
            closest_strand = "."
        else:
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
    """compute the distance to transcription start sites.

    Assumes that the intervals supplied are transcription start sites.

    An interval is associated with various tss:

       * the closest tss overall (prefix = closest)
       * the closest tss for which this interval is upstream of (prefix=upstream)
       * the closest tss for which this interval is downstream of (prefix=downstream)

    Note that if the tss is upstream to two genes (on either side), then
    there will no downstream entry (and vice versa).

    For each of these tss', the following features are output:
       * id: the transcript identifier of this tss
       * dist: the distance to the tss
       * strand: the strand of the tss

    If the interval overlaps with a tss, the flag ``is_overlap`` is set and
    the distance is set to 0.

    For closest_tss_dist, the distances are positive if the interval is downstream of
    the tss. It is negative for upstream distances.
    """


    headerTemplate = ( "closest_id", "closest_dist",                        
                       "closest_strand", 
                       "upstream_id", "upstream_dist",
                       "downstream_id", "downstream_dist",
                       "is_overlap" )

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
                                     self.mWithValues, self.mWithRecords, self.fasta,
                                     merge_genes = False,
                                     format = self.options.filename_format )
    def __str__(self):

        closest_id, closest_dist, closest_strand = ["na"] * 3
        upstream_id, upstream_dist, downstream_id, downstream_dist = ["na"] * 4

        overlaps = self.mOverlaps
        if overlaps:
            gene_ids = set([ x[2].gene_id for x in overlaps ])
            is_overlap = str(len(gene_ids))
            closest_id, closest_dist, closest_strand = ",".join( gene_ids), 0, "."            
            closest_id, closest_dist, closest_strand = list(gene_ids)[0], 0, "."
        else:
            is_overlap = "0"

            ## get distances to closest features depending which
            ## direction they are pointing.

            ## distance is always distance to 5' end

            upstream, downstream = [], []
            # 5'-----[ this ] ----- 3'

            # feature closest in 5' direction of this feature
            if self.mData5:
                if Genomics.IsPositiveStrand(self.strand5):
                    # 5'->3'-----[this] 
                    # feature points in same direction
                    downstream.append( (self.mDistance5 + (self.mData5.end - self.mData5.start), 
                                        self.mData5.gene_id, self.strand5) )
                else:
                    # 3'<-5'-----[this] 
                    # feature points in different direction
                    upstream.append( (self.mDistance5, self.mData5.gene_id, self.strand5) )

            if self.mData3:
                if Genomics.IsPositiveStrand(self.strand3):
                    # [this]-----5'->3'
                    upstream.append( (self.mDistance3, self.mData3.gene_id, self.strand3) )
                else:
                    # [this]-----3'<-5'
                    downstream.append( (self.mDistance3 + (self.mData3.end - self.mData3.start), self.mData3.gene_id, self.strand3) )

            # get minimum distance to 3' or 5' end of a gene
            if downstream: downstream = min(downstream) 
            else: downstream = []
            if upstream: upstream = min(upstream)
            else: upstream = []

            if downstream and upstream:
                if downstream < upstream:
                    closest_dist, closest_id, closest_strand = downstream
                else:
                    closest_dist, closest_id, closest_strand = upstream
                    closest_dist *= -1
            elif downstream:
                closest_dist, closest_id, closest_strand = downstream
            elif upstream:
                closest_dist, closest_id, closest_strand = upstream

            if downstream: downstream_dist, downstream_id = downstream[:2]
            if upstream: upstream_dist, upstream_id = upstream[:2]

        return "\t".join( ( 
                closest_id,
                str(closest_dist),
                str(closest_strand),
                upstream_id,
                str(upstream_dist),
                downstream_id,
                str(downstream_dist),
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
        if self.fasta: contig = self.fasta.getToken( contig)

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

        self.mProximalDistance = self.options.proximal_distance

    def count(self):

        # collect overlapping segments
        segments = self.getSegments()
        contig = self.getContig()
        if self.fasta: 
            contig = self.fasta.getToken( contig)

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
        if self.fasta: 
            try:
                contig = self.fasta.getToken( contig)
            except KeyError:
                E.warn( "contig %s not found" % contig )

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
class CounterQuality(_gtf2table.Counter):

    header = (Stats.Summary().getHeaders() + ("values",) )
    
    # discard segments with size > mMaxLength in order
    # to avoid out-of-memory
    mMaxLength = 100000

    def __init__(self, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )

    def count(self):
        ee = self.getSegments()
        self.result = self.getQuality( ee )

    def getQuality(self, segments):
        """get sequence from a set of segments."""

        contig = self.getContig()
        # quality is always '+'
        strand = "+"
        
        s = None
        for start,end in segments:
            if end - start > self.mMaxLength: return []
            if not s:
                s = self.fasta.getSequence( contig, strand, start, end )
            else:
                s.extend(self.fasta.getSequence( contig, strand, start, end ))
        return s

    def __str__(self):
        s = Stats.Summary( self.result, mode = "int" )
        values = ";".join( [ str(x) for x in self.result ] )
        return "\t".join( (str(s),values ) )


##-----------------------------------------------------------------------------------
class CounterReadExtension(_gtf2table.Counter):
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

        _gtf2table.Counter.__init__(self, *args, **kwargs )

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
            
        self.outfiles = IOTools.FilePool( self.options.output_filename_pattern % "readextension_%s" )

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
        assert last_exon_end <= territory_end, "assertion failed for %s: exon %i > territory %i" % (self.gene_id,
                                                                                                    last_exon_end,
                                                                                                    territory_end)
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
class CounterBigwigCounts(_gtf2table.Counter):
    '''obtain bigwig values and return summary stats.

    Requires a bigwig files to compute.
    '''
    
    header = ("length", "pcovered",) + Stats.Summary().getHeaders()
    
    # discard segments with size > mMaxLength in order
    # to avoid out-of-memory
    mMaxLength = 100000

    def __init__(self, bigwig_file, *args, **kwargs ):
        _gtf2table.Counter.__init__(self, *args, **kwargs )
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
        self.result = t

    def __str__(self):
        s = Stats.Summary( self.result, mode = "int" )
        return "\t".join( (str(self.mTotalLength),
                           "%5.2f" % (100.0 * self.mCovered / self.mTotalLength),
                           str(s),))

##------------------------------------------------------------
def main( argv = None ):

    parser = E.OptionParser( version = "%prog version: $Id: gtf2table.py 2888 2010-04-07 08:48:36Z andreas $", 
			     usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-q", "--quality-file", dest="quality_file", type="string",
                      help="filename with genomic base quality information [default=%default]."  )

    parser.add_option("-b", "--bam-file", dest="bam_files", type="string", metavar="bam",
                      help="filename with read mapping information. Multiple files can be submitted in a comma-separated list [default=%default]."  )

    parser.add_option("-i", "--bigwig-file", dest="bigwig_file", type="string", metavar="bigwig",
                      help="filename with bigwig information [default=%default]."  )

    parser.add_option("-f", "--filename-gff", dest="filename_gff", type="string", action="append", metavar='bed',
                      help="filename with extra gff files. The order is important [default=%default]."  )

    parser.add_option( "--filename-format", dest="filename_format", type="choice",
                       choices=("bed", "gff", "gtf" ),
                       help="format of secondary stream [default=%default]."  )

    parser.add_option( "--gff-source", dest="gff_sources", type="string", action="append",
                      help="restrict input to this 'source' in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option( "--gff-feature", dest="gff_features", type="string", action="append",
                      help="restrict input to this 'feature' in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option("-r", "--reporter", dest="reporter", type="choice",
                      choices=("genes", "transcripts" ),
                      help="report results for 'genes' or 'transcripts' [default=%default]."  )

    parser.add_option("-s", "--section", dest="sections", type="choice", action="append",
                      choices= ("exons", "introns" ), 
                      help="select range on which counters will operate [default=%default]."  )

    parser.add_option("-c", "--counter", dest="counters", type="choice", action="append",
                      choices=(	"bigwig-counts",
				"binding-pattern",
				"classifier", 
				"classifier-rnaseq",
				"classifier-rnaseq-splicing",
				"classifier-polii",
				"composition-na", 
				"composition-cpg", 
				"coverage", 
				"distance", 
				"distance-genes", 
				"distance-tss",
				"length", 
				'neighbours',
				"overlap", 
				"overlap-stranded",
				"overlap-transcripts",
				"overrun",
				"position", 
				"proximity", 
				"proximity-exclusive", 
				"proximity-lengthmatched",
				"quality",
				"read-coverage", 
				"read-extension", 
				"read-counts",
				"splice", 
				"splice-comparison", 
				"territories"),
		      help="select counters to apply to input [default=%default]."  )

    parser.add_option( "--add-gtf-source", dest="add_gtf_source", action="store_true",
                      help="add gtf field of source to output [default=%default]."  )

    parser.add_option( "--proximal-distance", dest="proximal_distance", type="int",
                      help="distance to be considered proximal to an interval [default=%default]."  )

    parser.add_option( "--weight-multi-mapping", dest="weight_multi_mapping", action="store_true",
		       help="weight multi-mapping reads is read-counts counter. Requires "
		       "the NH flag to be set by the mapper [default=%default]."  )

    parser.add_option( "--prefix", dest="prefixes", type="string", action="append",
                      help="add prefix to column headers - prefixes are used in the same order as the counters [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        reporter = "genes",
        with_values = True,
        sections = [],
        counters = [],
        filename_gff = [],
        filename_format = None,
        gff_features = [],
        gff_sources = [],
        add_gtf_source = False,
        proximal_distance = 10000, 
        bam_files = None,
	weight_multi_mapping = False,
        prefixes = []
        )

    if not argv: argv = sys.argv

    (options, args) = E.Start( parser, add_output_options = True, argv = argv )

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
            counters.append( _gtf2table.CounterReadCoverage( bam_files,
                                                  options = options, prefix = prefix ) )
        elif c == "read-extension":
            counters.append( CounterReadExtension( bam_files,
                                                   filename_gff = options.filename_gff,
                                                   options = options, 
                                                   prefix = prefix ) )
        elif c == "read-counts":
            counters.append( _gtf2table.CounterReadCounts( bam_files,
							   weight_multi_mapping = options.weight_multi_mapping,
							   options = options, 
							   prefix = prefix ) )
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
                    "distance", 
                    "distance-genes", 
                    "distance-tss",
                    "binding-pattern",
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
            elif c == "binding-pattern":
                template = CounterBindingPattern

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

        elif c == "classifier-rnaseq":
            counters.append( ClassifierRNASeq( filename_gff = options.filename_gff,
                                               fasta = fasta,
                                               options = options, prefix = prefix) )
        elif c == "classifier-rnaseq-splicing":
            counters.append( ClassifierRNASeqSplicing (filename_gff = options.filename_gff,
						       fasta = fasta,
						       options = options, prefix = prefix) )
        elif c == "classifier-polii":
            counters.append( ClassifierPolII( filename_gff = options.filename_gff,
                                              feature = None,
                                              source = None,
                                              fasta = fasta,
                                              options = options, prefix = prefix) )
        elif c == "binding-pattern":
            counters.append( CounterBindingPattern( filename_gff = options.filename_gff,
                                                    feature = None,
                                                    source = None,
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

	for counter in counters: counter.update(gffs)
	    
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

if __name__ == "__main__":
	sys.exit( main( sys.argv) )



