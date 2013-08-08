####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: psl2assembly.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####
'''
psl2assembly.py - assemble (long) reads on genome
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Assemble reads on a reference genome. Overlapping reads are merged.

This script is useful for building transcript models from
long aligned reads (for example, 454 reads or ESTs).

The script is aware of splice motifs and modifies alignments 
around intron boundaries to conform to splice site signals.

Input is a psl formatted file sorted by target contig and target start.

Output is a table linking reads to transcripts. Optional
output files are:

.transcripts: a psl formatted output file containing transcripts

.pileup: fasta file with pileup alignments of reads on transcripts

.regions: a tab separated file of transcript coordinates

.coverage: coverage statistics for each transcript

.sbjct_coverage: coverage statistics of the target genome

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----
'''

import sys
import re
import string
import optparse
import time
import os
import tempfile
import shutil
import collections
import math

import scipy
import scipy.stats
import numpy

import CGAT.Intervals as Intervals
import CGAT.Experiment as E
import CGAT.Histogram as Histogram
import CGAT.Blat as Blat
import CGAT.IndexedFasta as IndexedFasta
import alignlib
import CGAT.Mali as Mali
import CGAT.Stats as Stats
import CGAT.Genomics as Genomics



class Processor:

    def __init__(self, genome_fasta, genome_queries, options ):
        self.mGenomeFasta = genome_fasta

        self.mQueriesFasta = queries_fasta

        self.options = options

        # map of id to coverage information
        self.mMapId2Coverage = None

    def setMapId2Coverage( self, d ):
        """set mapping from ids to coverage information from previous run."""
        self.mMapId2Coverage = d

    def getReadCounts( self, matches, lgenome ):
        """return three arrays with the number of reads per position

        The three arrays are:
        1. counts for exons
        2. counts for introns
        3. counts for termini (intron/exon boundaries as well)
        ."""

        is_exon = numpy.zeros( lgenome, numpy.int )
        is_intron = numpy.zeros( lgenome, numpy.int )
        is_terminal = numpy.zeros( lgenome, numpy.int )
        
        ###############################################
        # map matches onto genomic sequence counting
        # introns and exons
        for match in matches:

            if self.options.loglevel >= 4:
                self.options.stdlog.write( "# "+ str( match ) + "\n" )
                self.options.stdlog.flush()

            was_exon = False


            assert lgenome >= match.mMapTarget2Query.getRowTo(), "alignment for match %s is out of bounds (%i > %i)" % \
                (str(match), match.mMapTarget2Query.getRowTo(), lgenome ) 

            # collapse the alignment
            for x in range(match.mMapTarget2Query.getRowFrom(), match.mMapTarget2Query.getRowTo() ):
                if match.mMapTarget2Query.mapRowToCol( x ) >= 0:            
                    is_exon[ x ] += 1
                    if not was_exon: is_terminal[x] += 1
                    was_exon = True
                else:
                    is_intron[ x ] += 1
                    if was_exon: is_terminal[x-1] += 1
                    was_exon = False
                    
            is_terminal[match.mMapTarget2Query.getRowTo()-1] += 1

        return is_exon, is_intron, is_terminal

class Builder(Processor):

    def __init__(self, *args, **kwargs ):
        Processor.__init__( self, *args, **kwargs )

        self.mIdFormat = options.output_format

        if self.options.output_filename_pattern:
            self.mOutFile = open( self.options.output_filename_pattern % self.mName, "w" )
        else:
            self.mOutFile = self.options.stdout

        self.mId = 0

    def __call__( self, id, *args ):

        self.mId = id

        if self.options.loglevel >= 3:
            self.options.stdlog.write("# module %s for %s started\n" % (self.mName, self.mId ) )
            self.options.stdlog.flush()
                                 
        self.mOutputId = self.mIdFormat % self.mId 
        self.process( *args )
        self.mOutFile.flush()

        if self.options.loglevel >= 3:
            self.options.stdlog.write("# module %s for %s finished\n" % (self.mName, self.mId ) )
            self.options.stdlog.flush()

    def printHeader( self ):
        header = self.getHeader()
        if header:
            self.mOutFile.write( header + "\n" )

    def getHeader( self ):
        return None

    def buildPileUp( self, sbjct_id, sbjct_start, sbjct_end, matches, map_genome2transcript = None ):
        """build a pileup multiple alignment of all matches on the genomic sequence.
        
        If map_genome2transcript is given, then the alignment is mapped from genomic coordinates
        (starting a 0) to a different coordinate system. Use this mapper to remove gaps or
        undesired regions from the sequence.
        """

        if self.options.loglevel >= 4:
            self.options.stdlog.write("# pileup started\n" )
            self.options.stdlog.flush()

        mali = alignlib.makeMultAlignment()

        genome = self.mGenomeFasta.getSequence( sbjct_id, "+", sbjct_start, sbjct_end )        

        if map_genome2transcript:
            mapped = alignlib.makeAlignmentBlocks()
            genome = "".join( [ genome[x] for x in range(len(genome)) if map_genome2transcript.mapRowToCol(x) != -1 ] )

        identifiers = ["%s:%s" % (self.mOutputId, "genome" ) ]
        ali = alignlib.makeAlignmentBlocks()
        lgenome = sbjct_end - sbjct_start
        ali.addDiagonal( 0, lgenome, 0 )
        seqs = alignlib.StringVector()
        seqs.append( genome )
        mali.add( ali )

        x = 0

        for match in matches:
            if self.options.loglevel >= 3:
                self.options.stdlog.write( "# "+ str( match ) + "\n" )
            
            query = self.mQueriesFasta.getSequence( match.mQueryId, match.strand )
            
            if self.options.loglevel >= 4:

                a = str(alignlib.AlignmentFormatExplicit( match.mMapTarget2Query,
                                                          alignlib.makeSequence( genome ),
                                                          alignlib.makeSequence( query ) ))

                # self.options.stdlog.write( str( alignlib.AlignmentFormatExplicit( match.mMapTarget2Query, 
#                                                                                  alignlib.makeSequence( genome ),
#                                                                                  alignlib.makeSequence( query ) ) ) + "\n" )
            
            if map_genome2transcript:
                alignlib.combineAlignment( mapped, 
                                           map_genome2transcript,
                                           match.mMapTarget2Query,
                                           alignlib.RR )
                                            
                mali.add( mapped )
            else:
                mali.add( match.mMapTarget2Query )
            
            seqs.append( query )

            x += 1
            identifiers.append( "%s:%s" % (self.mOutputId, match.mQueryId ))
        
        if self.options.loglevel >= 4:
            self.options.stdlog.write("# pileup finished\n" )
            self.options.stdlog.flush()

        return mali, identifiers, seqs 
    
    def finish( self ):
        self.options.stdlog.write( "# %s\tfinished\n" % (self.mName ))

class BuilderRegion(Builder):

    mName = "regions"

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

    def getHeader( self ):
        return "id\tcontig\tstart\tend\tsize\tnmatches\tnpos\tnneg"

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        npositive = sum( Genomics.IsPositiveStrand( x.strand ) for x in matches )
        nnegative = sum( Genomics.IsNegativeStrand( x.strand ) for x in matches )
        self.mOutFile.write( "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\n" % (self.mOutputId, sbjct_id, sbjct_start, sbjct_end, sbjct_end - sbjct_start, len(matches), npositive, nnegative ))
        self.mOutFile.flush()


class BuilderTranscribedLocus(Builder):

    mName = "locus.psl"

    """build a transcript from all aligned matches."""
    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

        if self.options.output_filename_pattern:
            self.mOutFileFasta = open( self.options.output_filename_pattern % self.mName + ".fasta", "w" )
            if self.options.output_pileup:
                self.mOutFilePileup = open( self.options.output_filename_pattern % self.mName + ".pileup", "w" )
            self.mOutFileStrand = open( self.options.output_filename_pattern % self.mName + ".strand", "w" )
        else:
            self.mOutFileFasta = self.options.stdout
            self.mOutFileStrand = self.options.stdout
            if self.options.output_pileup:
                self.mOutFilePileup = self.options.stdout

        self.mOutFileStrand.write( "gene_id\tstrand\tnmatches\tnpos\tnneg\tnintrons\tnintrons_pos\tnintrons_neg\tnintrons_unknown\n" )

    def getHeader( self ):
        return Blat.Match().getHeader()

    def guessStrand( self, matches, intron_counts = None ):
        """guess the strand of a transcript.

        The first guess is based on the directionality of splice site motifs
        given by the dict intron_counts. If that is inconclusive, the directionality
        of the matches is considered.

        Only clear votes count, i.e. if one splice site/read disagrees, the 
        strand is labelled ambiguous "."
        """
        strand = None
        if intron_counts and sum( intron_counts.values() ) > 0:
            nintrons_positive = intron_counts["+"]
            nintrons_negative = intron_counts["-"]
            nintrons_unknown = intron_counts["."]

            if nintrons_positive and nintrons_negative == 0:
                strand = "+"
            elif nintrons_negative and nintrons_positive == 0:
                strand = "-"
            else:
                strand = "."
        else:
            nintrons_positive = 0
            nintrons_negative = 0
            nintrons_unknown = 0

        npositive = sum( Genomics.IsPositiveStrand( x.strand ) for x in matches )
        nnegative = sum( Genomics.IsNegativeStrand( x.strand ) for x in matches )

        if not strand:
            if npositive and nnegative == 0:
                strand = "+"
            elif nnegative and npositive == 0:
                strand = "-"
            else:
                strand = "."

        self.mOutFileStrand.write( "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n" % \
                                        (self.mOutputId, 
                                         strand,
                                         len(matches),
                                         npositive,
                                         nnegative,
                                         nintrons_positive + nintrons_negative + nintrons_unknown,
                                         nintrons_positive,
                                         nintrons_negative,
                                         nintrons_unknown ) )

        return strand

    def getTranscriptMap( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        
        genome = self.mGenomeFasta.getSequence( sbjct_id, "+", sbjct_start, sbjct_end )
        
        residues = set()
        for match in matches:
            if self.options.loglevel >= 3:
                self.options.stdlog.write( "# "+ str( match ) + "\n" )
                self.options.stdlog.flush()

            # collapse the alignment
            for x in range(match.mMapTarget2Query.getRowFrom(), match.mMapTarget2Query.getRowTo() ):
                if match.mMapTarget2Query.mapRowToCol( x ) >= 0:            
                    residues.add( x )
        
        residues = list( residues )
        residues.sort()
        
        map_transcript2genome = alignlib.makeAlignmentBlocks()

        sequence = []

        for x in range(len(residues)):
            sequence.append( genome[residues[x]] )
            map_transcript2genome.addPair( x, residues[x] )
            
        map_transcript2genome.moveAlignment( 0, sbjct_start )

        return map_transcript2genome, "".join(sequence ), self.guessStrand( matches )

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        map_transcript2genome, sequence, strand = self.getTranscriptMap( sbjct_id, sbjct_start, sbjct_end, matches )

        match = Blat.Match()
        match.fromMap( map_transcript2genome )
        ## ambiguous strands are labelled "+" as psl does not know
        if strand == "-":
            match.strand = "-"
        else:
            match.strand = "+"

        match.mSbjctId = sbjct_id
        match.mQueryId = self.mOutputId
        match.mQueryLength = map_transcript2genome.getRowTo()
        match.mSbjctLength = self.mGenomeFasta.getLength( sbjct_id )

        self.mOutFile.write( str( match ) + "\n" )
        self.mOutFile.flush()
        self.mOutFileFasta.write( ">%s\n%s\n" % (self.mOutputId, sequence) )
        self.mOutFileFasta.flush()

        if self.options.output_pileup:
            self.outputPileUp( sbjct_id, sbjct_start, sbjct_end, matches, match, sequence )

    def outputPileUp( self, sbjct_id, sbjct_start, sbjct_end, matches, match, sequence ):
        """output a pileup multiple alignment for debugging purposes.

        Matches: the reads
        match: the transcript match
        """

        mali, identifiers, seqs = self.buildPileUp( sbjct_id, sbjct_start, sbjct_end, matches )
        
        map_genome2query = match.getMapTarget2Query()
        map_genome2query.moveAlignment( -sbjct_start, 0 ) 
        
        mali.add( map_genome2query )
        seqs.append( sequence )

        identifiers.append( "%s:%s" % (self.mOutputId, "transcript" ))

        if options.mali_output_format == "fasta":
            mm = Mali.convertAlignlib2Mali( mali, identifiers, seqs )
            mm.writeToFile( self.mOutFilePileup, format = "fasta" )

class BuilderTranscriptSpliced(BuilderTranscribedLocus):

    mName = "transcripts.psl"

    """build a transcript from all aligned matches.

    This processor includes a heuristics to correct for errors
    in the alignment/sequencing stage. The following errors are
    corrected:

    1. small gaps in the transcript. These are filled with genomic sequence
    2. splice site overrun = transcripts ending close to a splice site sometimes
       overshoot the splice site.

       Check for the existence of a splice site close by that is supported by at least one
       read.
       
    """
    mMaxGapLength = 15
    mFixIntronDistance = 10
    ## intron motifs [ forward strand ]
    mIntronTypes = ( ( "U2-GT/AG", "GT", "AG"),
                     ( "U2-nc-GC/AG", "GC", "AG"),
                     ( "U12-AT/AC", "AT", "AC") )
    mMaxSpliceMotifLength = 2

    def __init__(self, *args, **kwargs ):
        BuilderTranscribedLocus.__init__( self, *args, **kwargs )

        self.mFixedOverruns = 0
        self.mFixedGaps = 0
        self.mFixedCoding = 0
        self.mNIntrons = 0

        self.mIntronCounts = collections.defaultdict( int )

    def getIntronType( self, s ):
        """return intron type for an intronic sequence. 

        Checks both directions and returns the type of the intron and the strand.
        If the intron-type is not known, type is set to "unknown"
        """

        r = string.translate( s, string.maketrans("ACGTacgt", "TGCAtgca") )[::-1]

        for name, prime5, prime3 in self.mIntronTypes:
            if s[:len(prime5)].upper() == prime5 and \
                    s[-len(prime3):].upper() == prime3:
                return name, "+"
            elif r[:len(prime5)].upper() == prime5 and \
                    r[-len(prime3):].upper() == prime3:
                return name, "-"
        else:
            return "unknown", "."

    def getTranscriptMap( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        genome = self.mGenomeFasta.getSequence( sbjct_id, "+", sbjct_start, sbjct_end ).upper()

        lgenome = sbjct_end - sbjct_start
        is_exon, is_intron, is_terminal = self.getReadCounts( matches, lgenome )

        ###############################################
        # fill small gaps
        is_gap = False
        lgap = 0
        for x in xrange( 0, lgenome ):
            if is_exon[x] == 0:
                lgap += 1
            else:
                if lgap:
                    if lgap <= self.mMaxGapLength:
                        if self.options.loglevel >= 3:
                            self.options.stdlog.write( "# fixing small gap of length %i: %i-%i\n" % ( lgap, x-lgap, x) )

                        for y in xrange( x - lgap, x):
                            is_exon[y] = 1
                            is_intron[y] = 0
                        self.mFixedGaps += 1

                lgap = 0

        
        ###############################################
        # deal with intron overrun
        def hasSpliceMotif( sequence ):
            """find a splice motif in seq."""
            s = sequence.upper()
            r = string.translate( s.upper(), string.maketrans("ACGTacgt", "TGCAtgca") )[::-1]

            for name, prime5, prime3 in self.mIntronTypes:
                if s.startswith( prime5 ) or r.startswith( prime5 ) or \
                        s.startswith( prime3 ) or r.startswith( prime3 ):
                    return True
            return False

        def findSpliceStart( start, end ):
            """find a splice site start within start,end with more evidence than the original site."""
            for x in range(end,start,-1):
                if is_intron[x] > is_intron[x-1] and hasSpliceMotif( genome[x:x+self.mMaxSpliceMotifLength] ):
                    return x
            return None

        def findSpliceEnd( start, end ):
            """find a splice site end within start,end with more evidence than the original site."""
            for x in range(start,end):
                if is_intron[x-1] > is_intron[x] and hasSpliceMotif( genome[x-self.mMaxSpliceMotifLength:x] ):
                    return x-1
            return None

        for x in xrange( 0, lgenome ):
            # check for read ends within introns
            if self.options.loglevel >= 5:
                self.options.stdlog.write( "# %i(%i): is_term=%i, is_exon=%i, is_intron=%i, check=%i\n" %\
                                               (x, x+sbjct_start,
                                                is_terminal[x], is_exon[x], is_intron[x],
                                                is_terminal[x] and is_terminal[x] == is_exon[x] and is_intron[x] ) )

            if is_terminal[x] and is_terminal[x] == is_exon[x] and is_intron[x]:
                left = max(0, x-self.mFixIntronDistance )
                right = min( lgenome-1, x+self.mFixIntronDistance )
                if is_exon[left] > 0:
                    if is_exon[left] > is_exon[x]:
                        z = findSpliceStart( left, x )
                        if z != None:
                            if self.options.loglevel >= 3:
                                self.options.stdlog.write( "# correcting intron boundary: %i-%i %s %s %s\n" %\
                                                               ( z, 
                                                                 x + 1, 
                                                                 genome[z-5:z],
                                                                 genome[z:x+1],
                                                                 genome[x+1:x+6]))

                            for y in range( z, x + 1):
                                is_exon[y] = 0
                                is_intron[y] = 1
                            self.mFixedOverruns += 1

                elif is_exon[right] > 0:
                    if is_exon[right] > is_exon[x]:
                        z = findSpliceEnd( x, right )
                        if z != None:
                            if self.options.loglevel >= 3:
                                self.options.stdlog.write( "# correcting intron boundary: %i-%i: %s %s %s\n" %\
                                                           ( x, z + 1, 
                                                             genome[x-5:x],
                                                             genome[x:z+1],
                                                             genome[z+1:z+6]))

                            for y in range( x, z + 1):
                                is_exon[y] = 0
                                is_intron[y] = 1
                            self.mFixedOverruns += 1
                
        map_transcript2genome = alignlib.makeAlignmentBlocks()

        sequence = []

        y = 0
        for x in range(len(is_exon)):
            if is_exon[x]:
                sequence.append( genome[x] )
                map_transcript2genome.addPair( y, x )
                y += 1

        map_transcript2genome.moveAlignment( 0, sbjct_start )

        ## decide upon strand based on introns
        ## and do counting
        strands = { ".": 0, "+" : 0, "-": 0 }
        for start, end in self.getIntrons( is_exon ):
            name, strand = self.getIntronType( genome[start:end] )
            strands[strand] += 1
            if options.loglevel >= 3:
                self.options.stdlog.write( "# intron check: %i-%i: %i-%i: %s %s\n" %\
                                               ( start, end, 
                                                 start+sbjct_start,
                                                 end+sbjct_start,
                                                 name, strand ) )


            self.mIntronCounts[name] += 1
            
        return map_transcript2genome, "".join(sequence), self.guessStrand( matches, intron_counts = strands )

    def getIntrons( self, is_exon ):
        """return a list of intron intervals."""
        introns = []
        l = len(is_exon)
        start, end = 0, l - 1
        while start < l and is_exon[start] == 0: start += 1
        while end >= 0 and is_exon[end] == 0: end -= 1
            
        was_intron = False
        f = 0
        for x in xrange( start, end ):
            if is_exon[x] and was_intron:
                introns.append( (f, x) )
            elif not is_exon[x] and not was_intron:
                f = x
            was_intron = not is_exon[x]

        if was_intron: introns.append( (f,end) )
        return introns

    def finish( self ):

        if self.options.loglevel >= 1:
            self.options.stdlog.write( "# %s\tfinished: fixed_gaps=%i, fixed_overruns=%i\n" % (self.mName, self.mFixedGaps, self.mFixedOverruns ) )
            total = sum( self.mIntronCounts.values() )
            if total:
                self.options.stdlog.write( "# %s\tstats:\n" % self.mName )
                self.options.stdlog.write( "# %s\ttype\tcounts\tfreq\n" % (self.mName) )
                for key,val in sorted(self.mIntronCounts.items()):
                    self.options.stdlog.write( "# %s\t%s\t%i\t%5.2f\n" % (self.mName, key, val, 100.0 * val / total ) )
                self.options.stdlog.write( "# %s\t%s\t%i\t%5.2f\n" % ( self.mName, "total", total, 100.0 * total / total ) )

class BuilderCoverage(Builder):
    """compute residue read coverage over all genomic bases that have at least one aligned
    match.

    If previous coverage information is available, aggregate information is output.
    (updated mean and pooled standard deviation).
    """

    mName = "coverage"

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )
       
    def getHeader( self ):
        ## if using previous coverage values, only mean is available
        if self.mMapId2Coverage:
            return "id\tcontig\tstart\tend\tsize\tnmatches\tncovered\tmean\tstddev" 
        else:
            return "id\tcontig\tstart\tend\tsize\tnmatches\tncovered\t%s" % ("\t".join(Stats.DistributionalParameters().getHeaders()) )

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        residues = collections.defaultdict(int)
        for match in matches:
            if self.options.loglevel >= 4:
                self.options.stdlog.write( "# "+ str( match ) + "\n" )
            
            # collapse the alignment
            # offset = match.mSbjctFrom - sbjct_start
            for x in range(match.mMapTarget2Query.getRowFrom(), match.mMapTarget2Query.getRowTo() ):
                if match.mMapTarget2Query.mapRowToCol( x ) >= 0:            
                    residues[x] += 1

        ncovered = len(residues)

        if self.mMapId2Coverage:
            c, t, nmatches = 0, 0, 0

            vars = []
            for match in matches:
                try:
                    m = self.mMapId2Coverage[match.mQueryId]
                except KeyError:
                    E.warn("no coverage information for %s" % match.mQueryId)
                    continue
                c += m.mMean * m.mNCovered
                t += m.mNCovered
                nmatches += m.mNMatches
                vars.append( (m.mNCovered, m.mStddev * m.mStddev)  )
                
            pooled_variance = Stats.getPooledVariance( vars )
            extra_info = "%6.4f\t%6.4f" % ( (float(c) / t), math.sqrt( pooled_variance) ) 

        else:
            nmatches = len(matches)
            extra_info = str(Stats.DistributionalParameters( residues.values() ))

        self.mOutFile.write( "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\n" % (self.mOutputId, sbjct_id, sbjct_start, sbjct_end, 
                                                                   sbjct_end - sbjct_start, 
                                                                   nmatches,
                                                                   ncovered,
                                                                   extra_info))
        self.mOutFile.flush()



class BuilderSbjctCoverage(Builder):
    """compute residue read coverage over all genomic bases in the sbjct
    irrespective if they match zero or more reads.
    """

    mName = "sbjct_coverage"

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )
       
    def getHeader( self ):
        return "id\tcontig\tstart\tend\tsize\tnmatches\tncovered\t%s" % ("\t".join(Stats.DistributionalParameters().getHeaders()) )

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        sbjct_length = self.mGenomeFasta.getLength( sbjct_id )

        residues = collections.defaultdict(int)

        for match in matches:
            if self.options.loglevel >= 4:
                self.options.stdlog.write( "# "+ str( match ) + "\n" )
            
            # collapse the alignment
            # offset = match.mSbjctFrom - sbjct_start
            for x in range(match.mMapTarget2Query.getRowFrom(), match.mMapTarget2Query.getRowTo() ):
                if match.mMapTarget2Query.mapRowToCol( x ) >= 0:            
                    residues[x] += 1

        self.mOutFile.write( "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\n" % (self.mOutputId, sbjct_id, sbjct_start, sbjct_end, 
                                                                   sbjct_end - sbjct_start, 
                                                                   len(matches),
                                                                   len(residues),
                                                                   str(Stats.DistributionalParameters( residues.values() + [0.0] * (sbjct_length - len(residues)) ) )))
        self.mOutFile.flush()

class BuilderIndels(Builder):
    """indel variation in the genome or the reads.
    """

    mName = "indels"

    # region to examine on either side to call an indel in the genome
    mRegionSize = 5

    # minimum support (= number of reads) on either side of an indel
    mMinSupport = 3

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

        assert self.mQueriesFasta, "%s requires a sequence database with the queries" % self.mName
        assert self.mGenomeFasta, "%s requires a sequence database with the genome" % self.mName
       
        self.mNDeletions = 0
        self.mNInsertions = 0

        if self.options.loglevel >= 1:
            self.options.stdlog.write("# %s\tstarted: mRegionSize=%i, mMinSupport=%i\n" % (self.mName, self.mRegionSize, self.mMinSupport) )

    def getHeader( self ):
        return "id\tcontig\tpos\tposT\tposM\tchar\tnmatches\tnaligned\tngaps\tsupport"

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        """output info on positions that are gap-chars in genome.

        This script will only report single base indels.
        """

        lgenome = sbjct_end - sbjct_start
        is_exon, is_intron, is_terminal = self.getReadCounts( matches, lgenome )
        mali, identifiers, seqs = self.buildPileUp( sbjct_id, sbjct_start, sbjct_end, matches )

        if self.options.loglevel >= 4:
            self.options.stdlog.write("# converting mali started\n" )
            self.options.stdlog.flush()
        
        columns = Mali.convertAlignlib2Mali( mali, identifiers, seqs ).getColumns()

        if self.options.loglevel >= 4:
            self.options.stdlog.write("# converting finished\n" )
            self.options.stdlog.flush()

        assert identifiers[0] == "%s:genome" % self.mOutputId

        nmatches = len(matches)

        min_support = self.mMinSupport

        # c: index in alignment, g: index on genome
        # skip first and last position
        c, g = 1, 1
        for column in columns[1:-1]:

            ngaps = column[1:].count("-")
            # can not use is_exons, as this has only positions on genome
            naligned = nmatches - ngaps
            is_gap = column[0] == "-"
            # deletion in genome: support = reads with inserted base
            if (is_gap and naligned >= min_support) or \
                    (not is_gap and is_exon[g-1] >= min_support and naligned == 0 and is_exon[g+1] >= min_support):

                naligned_left = min( is_exon[ max(0, g-self.mRegionSize):g] )
                naligned_right = min( is_exon[ g+1:min(lgenome, g+self.mRegionSize+1)] ) 
                support = min( (naligned_left, naligned_right) )

                if support >= min_support:
                    if is_gap:
                        self.mNDeletions += 1
                    else:
                        self.mNInsertions += 1

                    self.mOutFile.write( "%s\t%s\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\n" % ( self.mOutputId, 
                                                                                    sbjct_id, 
                                                                                    sbjct_start + g, 
                                                                                    g,
                                                                                    c,
                                                                                    column[0],
                                                                                    nmatches,
                                                                                    naligned,
                                                                                    ngaps,
                                                                                    support ) )
                    self.mOutFile.flush()

            # increment indices: skip over gaps in genome
            if not is_gap: g += 1
            c += 1

    def finish( self ):

        if self.options.loglevel >= 1:
            self.options.stdlog.write( "# %s\tfinished: ndeletions=%i, ninsertions=%i\n" % (self.mName, self.mNDeletions, self.mNInsertions ) )

##-------------------------------------------------------------------------
class BuilderConsensus(Builder):
    """builds a pileup alignment of reads on the transcript sequence 
    and outputs the consensus sequence (without introns).

    This module can be used to build transcripts of cross-species
    mappings. Note that the sequence will no correspond to the
    one output by BuilderTranscriptsIntrons as it does not 
    correct for exon overrun.
    """
    
    mName = "consensus"

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

        assert self.mQueriesFasta, "%s requires a sequence database with the queries" % self.mName
        assert self.mGenomeFasta, "%s requires a sequence database with the genome" % self.mName

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        lgenome = sbjct_end - sbjct_start
        is_exon, is_intron, is_terminal = self.getReadCounts( matches, lgenome )
        
        map_genome2transcript = alignlib.makeAlignmentVector()

        c = 0
        for x in range(lgenome):
            if is_exon[x] > 0:
                map_genome2transcript.addPair( x, c, 0 )
                c+=1
                
        if options.loglevel >= 3:
            options.stdlog.write( "# %s\n" % str(alignlib.AlignmentFormatEmissions( map_genome2transcript)) )
            options.stdlog.flush()

        mali, identifiers, seqs = self.buildPileUp( sbjct_id, sbjct_start, sbjct_end, matches,
                                                    map_genome2transcript = map_genome2transcript)

        if options.loglevel >= 4:
            print str(mali)
            options.stdlog.flush()

        if not identifiers:
            identifiers = [ "%i" % x for x in range( mali.getNumSequences() ) ]

        if options.mali_output_format == "fasta":
            mm = Mali.convertAlignlib2Mali( mali, identifiers, seqs )
            mm.writeToFile( self.mOutFile, format = "fasta" )

##-------------------------------------------------------------------------
class BuilderPileUp(Builder):

    mName = "pileup"

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

        assert self.mQueriesFasta, "%s requires a sequence database with the queries" % self.mName
        assert self.mGenomeFasta, "%s requires a sequence database with the genome" % self.mName

    @E.benchmark
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):

        mali, identifiers, seqs = self.buildPileUp( sbjct_id, sbjct_start, sbjct_end, matches )

        if not identifiers:
            identifiers = [ "%i" % x for x in range( mali.getNumSequences() ) ]

        if options.mali_output_format == "fasta":
            mm = Mali.convertAlignlib2Mali( mali, identifiers, seqs )
            mm.writeToFile( self.mOutFile, format = "fasta" )


##-------------------------------------------------------------------------
class BuilderPolyA(Builder):
    """predict polyA tails.

    The output is a tab-separated table with the following columns

    self.mOutFile.write( "id\ttails\tnremoved\tnstrands\tnmotifs\tninconsistent\n" )

    id
       transcript id
    tails
       number of predicted tails    
    nremoved
       number of tails that were joined with other tails
    nstrands
       number of strands that tails are pointing to
    nmotifs
       number of detected polyA motifs (AATAAA)
    ninconsistent
       number of inconsistent tails

    """    
    mName = "polyA"

    mMotifs = { "AATAAA": ("AATAAA", "TTTATT") }

    # number of residues to search for polyA motifs
    mMotifArea = 50

    # combine polyA starts that are closer than x
    mMaxAggregateArea = 10

    # max residues missing from non polyA end
    mMaxUnaligned = 3
    # min residues in tail
    mMinUnaligned = 10
    # min percent residues that are A/T in tail
    mMinPercent = 70.0

    def __init__(self, *args, **kwargs ):
        Builder.__init__( self, *args, **kwargs )

        assert self.mQueriesFasta, "%s requires a sequence database with the queries" % self.mName
        assert self.mGenomeFasta, "%s requires a sequence database with the genome" % self.mName

        self.mOutFile.write( "id\ttails\tnremoved\tnstrands\tnmotifs\tninconsistent\n" )

        self.mOutFileTails = open( self.options.output_filename_pattern % self.mName + ".tails", "w" )
        self.mOutFileTails.write( "id\tpos\tstrand\tnmotifs\tmotifs\tnids\tids\n" )

    @E.benchmark        
    def process( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        """detect PolyA tails in transcripts.

        1. collect all matches with polyA tail - the remainder is unchanged
        2. declare polyA tail by the shortest! unaligned segment and compute
        coverage for each match appropriately. 
        
        The method checks whether the tails are consistent (always at the same end).
        If not, an AssertionError is thrown
        """

        genome = self.mGenomeFasta.getSequence( sbjct_id, "+", sbjct_start, sbjct_end ).upper()

        tails = collections.defaultdict( list )

        for match in matches:

            if options.loglevel >= 5:
                options.stdlog.write( "#%s\n" % str( match) )

            query_id = match.mQueryId
            missing_start = match.mQueryFrom
            missing_end = match.mQueryLength - match.mQueryTo
            if missing_start < missing_end:
                smaller = missing_start
                larger = missing_end
                start,end = match.mQueryTo,match.mQueryLength
                if match.strand == "+":
                    poly_start = match.mQueryTo - 1
                    direction = "+"
                else:
                    poly_start = match.mQueryLength - match.mQueryTo
                    direction = "-"
            else:
                smaller = missing_end
                larger = missing_start
                start,end = 0, match.mQueryFrom
                if match.strand == "+":
                    poly_start = match.mQueryFrom
                    direction = "-"
                else:
                    poly_start = match.mQueryLength - match.mQueryFrom - 1
                    direction = "+"

            # check if tail is at least polyA_min_aligned and at most polyA_max_unaligned
            # are missing from the other end.
            if not(smaller < self.mMaxUnaligned and larger > self.mMinUnaligned):
                continue

            tail = queries_fasta.getSequence( query_id )[start:end]

            counts = {"A": 0, "T": 0, "N": 0}
            for c in tail.upper(): counts[c] = counts.get(c, 0) + 1
            total = end-start
            pA = 100.0 * (counts["A"] + counts["N"]) / total
            pT = 100.0 * (counts["T"] + counts["N"]) / total
            
            if max(pA,pT) < self.mMinPercent:
                continue

            map_query2genome = match.getMapQuery2Target()
            map_query2genome.moveAlignment( 0, -sbjct_start ) 
            
            mapped = map_query2genome.mapRowToCol( poly_start )

            tails[(mapped,direction)].append( (query_id,len(tail)) )

            if options.loglevel >= 5:
                options.stdlog.write( "# polyA detection: found=%s:%i-%i pA=%5.2f pT=%5.2f pos=%i:%s tail=%i:%s\n" % (query_id, start, end, pA, pT, mapped, direction,len(tail),tail))
                
        if not tails: 
            self.mOutFile.write( "%s\t%i\t%i\t%i\t%i\t%i\n" % \
                                     (self.mOutputId,
                                      len(tails),
                                      0, 0, 0, 0))
            return

        # filter polyA tails
        if options.loglevel >= 3:
            options.stdlog.write( "# putative polyA tails: %s\n" % (str(tails)))

        # aggregate tails
        k = sorted(tails.keys())
        data = tails[k[0]]
        last_mapped, last_strand = k[0]
        new_tails = { (last_mapped, last_strand): data }
        for this_key in k[1:]:
            mapped,strand =this_key

            if strand == last_strand and mapped - last_mapped < self.mMaxAggregateArea:
                data += tails[ this_key ]
            else:
                data = tails[ this_key]
                new_tails[ this_key ] = data

            last_strand, last_mapped = strand, mapped

        nreduced = len(tails) - len(new_tails)
        tails = new_tails

        # check direction
        directions = set([ x[1] for x in tails.keys() ])

        if len(directions) > 1:

            # use number of residues in tail to decide upon strand
            #nforward = sum([ sum([ y[1] for y in x[1]]) for x in tails.items() if x[0][1] == "+"] )
            #nreverse = sum([ sum([ y[1] for y in x[1]]) for x in tails.items() if x[0][1] == "-"] )
            nforward = sum([ len(x[1]) for x in tails.items() if x[0][1] == "+"] )
            nreverse = sum([ len(x[1]) for x in tails.items() if x[0][1] == "-"] )

            ninconsistent = min(nforward, nreverse )

            if options.loglevel >= 1:
                options.stdlog.write("# WARNING: %s: inconsistent polyA tails: forward=%i, reverse=%i from %s\n" % \
                                         (self.mOutputId,
                                          nforward, nreverse, str(tails)) )

            #if nforward >= nreverse:
            #    tails = dict( [ (x[0],x[1]) for x in tails.items() if x[0][1] == "+"] )
            #else:
            #    tails = dict( [ (x[0],x[1]) for x in tails.items() if x[0][1] == "-"] )
            
        else:
            ninconsistent = 0

        nmotifs = 0
        for key, ids in tails.items():

            mapped, direction = key

            distances = []
            for name, motifs in self.mMotifs.items():
                if direction == "+":
                    pos = genome.rfind( motifs[0], mapped-self.mMotifArea,mapped )
                else:
                    pos = genome.find( motifs[1], mapped, mapped+self.mMotifArea )

                if pos >= 0:
                    distances.append( (name, abs(pos - mapped)) )
            
            self.mOutFileTails.write( "%s\t%i\t%s\t%i\t%s\t%s\t%s\n" % \
                                          (self.mOutputId,
                                           mapped,
                                           direction,
                                           len(distances),
                                           ",".join( [ "%s:%i" % x for x in distances ] ),
                                           len(ids),
                                           ",".join( [ "%s:%i" % x for x in ids] )))

            if len(distances): nmotifs += 1

        self.mOutFile.write( "%s\t%i\t%i\t%i\t%i\t%i\n" % \
                                 (self.mOutputId,
                                  len(tails),
                                  nreduced,
                                  len(directions),
                                  nmotifs,
                                  ninconsistent))
        

        self.mOutFileTails.flush()
        self.mOutFile.flush()

##-------------------------------------------------------------------------
class Filter(Processor):

    mName = "Filter"

    def __init__(self, *args, **kwargs ):
        Processor.__init__( self, *args, **kwargs )

        if self.options.output_filename_pattern:
            self.mOutFile = open( self.options.output_filename_pattern % self.mName, "w" )
        else:
            self.mOutFile = self.options.stdout

    def __call__( self, *args ):
        return self.filter( *args )

class FilterImperfect(Filter):
    """remove imperfect matches, if there is no evidence from other matches."""

    mName = "FilterImperfect"

    mMaxNMismatches = 0
    mMaxNMissingResidues = 0

    mMaxQueryNGapsCounts = 0
    mMaxQueryNGapsBases = 0
    mMaxSbjctNGapsCounts = 0
    mMaxSbjctNGapsBases = 0

    def __init__(self, *args, **kwargs ):
        Filter.__init__( self, *args, **kwargs )

    def filter( self, matches ):
        """filter matches."""
        
        if self.options.loglevel >= 2:
            options.stdlog.write( "# started filter %s with %i matches\n" %\
                                      (self.mName, len(matches) ) )
            options.stdlog.flush()

        if len(matches) == 1:
            match = matches[0]
            if match.mQueryLength - (match.mQueryTo - match.mQueryFrom) > self.mMaxNMissingResidues or \
                    match.mNMismatches > self.mMaxNMismatches or \
                    match.mQueryNGapsCounts > self.mMaxQueryNGapsCounts or \
                    match.mQueryNGapsBases > self.mMaxQueryNGapsBases or \
                    match.mSbjctNGapsCounts > self.mMaxSbjctNGapsCounts or \
                    match.mSbjctNGapsBases > self.mMaxSbjctNGapsBases:
                return []
            
        return matches

class FilterIntrons(Filter):
    """remove intronic matches that have little evidence

    An intron needs to be supported by at least two matches and
    the intron boundaries need to coincide perfectly.

    Introns need to be at least mIntronSize residues long.
    """

    # minimum intron size
    mIntronSize = 30

    def __init__(self, *args, **kwargs ):
        Filter.__init__( self, *args, **kwargs )

    def filter( self, matches ):
        """filter matches.
        not fully implemented.
        """
        
        raise "not implemented"
        matches_without_introns, matches_with_introns = [], []

        for match in matches:
            if match.mSbjctNGapsBases < self.mIntronSize:
                matches_without_introns.append( match )
            else:
                # ignore matches with two introns or two introns and a gap
                if match.mSbjctNGapsCounts <= 1:
                    matches_with_introns.append( match )
            
        new_matches_with_introns = []
        for match in matches_with_introns:

            ## find intron
            b = match.mBlockSizes[0] 
            last_q = match.mQueryBlockStarts[0] + b
            last_s = match.mSbjctBlockStarts[0] 
            for q, s, b in zip( match.mQueryBlockStarts[1:], match.mSbjctBlockStarts[1:], match.mBlockSizes[1:]):
                if last_s - s >= self.mIntronSize:
                    is_confirmed == False
                    # positions of exon boundaries on genome in alignment coordinates
                    pos1, pos2 = last_s - 1 - match.mSbjctFrom, s - match.mSbjctFrom
                    for mm in matches_with_introns:
                        if match.mQueryToken == mm.mQueryToken: continue
                        offset = match.mSbjctFrom - mm.mSbjctFrom
                        pp1 = pos1 + offset
                        pp2 = pos2 + offset
                        if pp1 >= 0 and pp2 >= 0 and \
                                match.mMapTarget2Query.mapRowToCol( pp1 ) >= 0 and \
                                match.mMapTarget2Query.mapRowToCol( pp2 ) >= 0 and \
                                match.mMapTarget2Query.mapRowToCol( pp2 ) - match.mMapTarget2Query.mapRowToCol( pp1 ) == 1:
                            is_confirmed = True
                            break
                    if is_confirmed:
                        new_matches_with_introns.append( match )
                        break

                new_matches_with_introns.append( match )
            
        return matches_without_introns + matches_with_introns

class FilterExonExtenders(Filter):
    """remove matches without introns that extend over an intron/exon boundary 
    created by another read.

    Introns need to be at least mMinIntronSize residues long.
    """
    mName = "filter_exon_extenders"

    # minimum intron size
    mMinIntronSize = 30

    # at least this number of residues are extended
    mMinLengthExonExtender = 10

    def __init__(self, *args, **kwargs ):
        Filter.__init__( self, *args, **kwargs )

        self.mOutFile.write( "tName\ttStart\ttEnd\tqName\trName\tintron\tnOverlap\n" )

    @E.benchmark
    def filter( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        """filter matches.
        not fully implemented.
        """

        if self.options.loglevel >= 2:
            options.stdlog.write( "# started filter %s with %i matches\n" %\
                                      (self.mName, len(matches) ) )
            options.stdlog.flush()
        
        matches_without_introns, matches_with_introns = [], []

        for match in matches:
            if match.mSbjctNGapsBases < self.mMinIntronSize:
                matches_without_introns.append( match )
            else:
                matches_with_introns.append( match )
                    
        new_matches_without_introns = []
        for match_without_introns in matches_without_introns:
            for match_with_introns in matches_with_introns:
                d = alignlib.getAlignmentShortestDistance( \
                    match_without_introns.mMapTarget2Query, 
                    match_with_introns.mMapTarget2Query, 
                    alignlib.RR )

                if d > 0: continue

                do_delete = False

                for intron_start, intron_end in match_with_introns.iterator_introns():

                    if intron_end - intron_start < self.mMinIntronSize:
                        continue

                    # check overlap with this match
                    intron_map = alignlib.makeAlignmentBlocks()
                    intron_map.addDiagonal( intron_start - sbjct_start, 
                                            intron_end - sbjct_start, 
                                            0 )

                    d = alignlib.getAlignmentShortestDistance( \
                        match_without_introns.mMapTarget2Query, 
                        intron_map,
                        alignlib.RR )

                    if d > 0: continue
                    intron_sequence = self.mGenomeFasta.getSequence( sbjct_id, "+", intron_start, intron_end )
                    intron_type, prime5, prime3 = Genomics.GetIntronType( intron_sequence, both_strands = True )
                    
                    overlap = alignlib.getAlignmentOverlap( match_without_introns.mMapTarget2Query, 
                                                            intron_map,
                                                            alignlib.RR )

                    if intron_type != "unknown" and overlap > self.mMinLengthExonExtender:
                        
                        if self.options.loglevel >= 1:
                            self.options.stdlog.write( "# deleting read: %s: extends an exon into a %s intron in %s by %i residues\n" %\
                                                           (match_without_introns.mQueryId, intron_type, 
                                                            match_with_introns.mQueryId, overlap ) )
                            
                        do_delete = True
                        self.mOutFile.write( "%s\t%i\t%i\t%s\t%s\t%s\t%i\n" % 
                                             (sbjct_id, 
                                              sbjct_start, 
                                              sbjct_end, 
                                              match_without_introns.mQueryId,
                                              match_with_introns.mQueryId,
                                              intron_type,
                                              overlap) )
                        self.mOutFile.flush()
                        break

                if do_delete: 
                    break
            else:
                new_matches_without_introns.append( match_without_introns )

        return matches_with_introns + new_matches_without_introns

class FilterTranscriptMergers(Filter):
    """remove matches with introns that are not well supported.

    Introns need to be at least mMinIntronSize residues long.
    """
    
    mName = "filter_transcript_mergers"

    # minimum intron size to be deleted
    mMinIntronSize = 300

    # maximum minimum coverage in intron to be deleted
    mMaxCoverageIntron = 1.0

    # minimum avegare coverage in adjacent exons to consider intron for deletion
    mMinCoverageExons = 2.0

    def __init__(self, *args, **kwargs ):
        Filter.__init__( self, *args, **kwargs )

        self.mOutFile.write( "tName\ttStart\ttEnd\tqName\tiStart\tiEnd\tiLen\tcExons\tcIntrons\n" )
                             
    @E.benchmark
    def filter( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        """filter matches.
        not fully implemented.
        """

        if self.options.loglevel >= 2:
            options.stdlog.write( "# started filter %s with %i matches\n" %\
                                      (self.mName, len(matches) ) )
            options.stdlog.flush()

        lgenome = sbjct_end - sbjct_start
        is_exon, is_intron, is_terminal = self.getReadCounts( matches, lgenome )

        new_matches = []

        for match in matches:

            keep = True

            exons = match.iterator_exons()
            last_start,last_end = exons.next()
            last_start -= sbjct_start
            last_end -= sbjct_start
            for start,end in exons:
                start -= sbjct_start
                end -= sbjct_start

                if start - last_end < self.mMinIntronSize: continue

                coverage_exons = float(numpy.sum( is_exon[last_start:last_end] ) + numpy.sum( is_exon[start:end] )) / \
                    ( last_end - last_start + end - start )

                v = [ is_intron[x] for x in filter( lambda x: is_exon[x] == 0, range(last_end,start)) ] 
                if len(v) == 0: continue
                coverage_intron = min ( v )
                
                lintron = start - last_end
                
                if coverage_intron <= self.mMaxCoverageIntron and coverage_exons >= self.mMinCoverageExons:
                    if options.loglevel >= 3:
                        options.stdlog.write( "# read: %s with intron %i-%i (l=%i): mean_cov_exons=%5.2f, min_cov_intron=%5.2f\n" %\
                                                  (match.mQueryId, 
                                                   sbjct_start + last_end, 
                                                   sbjct_start + start, 
                                                   lintron,
                                                   coverage_exons,
                                                   coverage_intron ) )

                    self.mOutFile.write( "%s\t%i\t%i\t%s\t%i\t%i\t%i\t%5.2f\t%5.2f\n" % 
                                         (sbjct_id, sbjct_start, sbjct_end, 
                                          match.mQueryId,
                                          sbjct_start + last_end, 
                                          sbjct_start + start, 
                                          lintron,
                                          coverage_exons,
                                          coverage_intron ) )

                    self.mOutFile.flush()
                    keep = False
                    break

                last_start, last_end = start, end

            if keep: new_matches.append( match )
            
        if options.loglevel >= 3:
            options.stdlog.write( "# %s - %s:%i..%i before=%i after=%i\n" %\
                                      (self.mName,
                                       sbjct_id, 
                                       sbjct_start,
                                       sbjct_end,
                                       len(matches),
                                       len(new_matches)))


        return new_matches

class FilterDuplicates(Filter):
    """remove duplicate reads with identical sequence from the list.
    """
    
    mName = "duplicates"

    # maximum difference in length
    mMinLengthDifference = 5

    mMinDistance = 10
    
    mMinAlignmentOverlap = 10

    mMinProbability = 0.05

    mMaxMatches = 10000

    def __init__(self, *args, **kwargs ):
        Filter.__init__( self, *args, **kwargs )

        self.mOutFile.write( "rep\tmem\tlRep\tlMem\tnoverlap\tdelta\tnreads\tPL\tPS\tP\n" )

        if options.loglevel >= 1:
            options.stdlog.write( "# %s: obtaining length distribution\n" % self.mName )
            options.stdlog.flush()

        ## estimate mean and stddev 
        lengths = self.mQueriesFasta.getLengths()
        summary = Stats.Summary( lengths )
        self.mLengthMean = summary.mMean
        self.mLengthMedian = summary.mMedian
        self.mLengthStd = summary.mSampleStd
                             
        if options.loglevel >= 1:
            options.stdlog.write( "# %s: estimated length distribution: n=%i, mean=%f, std=%f, median=%f\n" %\
                                      (self.mName,
                                       len(lengths),
                                       self.mLengthMean,
                                       self.mLengthStd,
                                       self.mLengthMedian) )
            options.stdlog.flush()

    @E.benchmark
    def filter( self, sbjct_id, sbjct_start, sbjct_end, matches ):
        """remove duplicate reads.

        Duplicate reads are mapped to the same position. The matches are sorted by
        increasing number of mismatches and gaps in query. With homopolymer runs
        in 454 data extra bases are often inserted.
        """

        if self.options.loglevel >= 2:
            options.stdlog.write( "# started filter %s with %i matches\n" %\
                                      (self.mName, len(matches) ) )
            options.stdlog.flush()

        new_matches = []
        keep = True
        eliminated = set()
        matches.sort( key=lambda x: x.mNMismatches + x.mQueryNGapsBases)

        lgenome = sbjct_end - sbjct_start
        is_exon, is_intron, is_terminal = self.getReadCounts( matches, lgenome )

        if len(matches) > self.mMaxMatches:
            return matches

        for x in range(len(matches)):
            if x in eliminated: continue
            mx = matches[x] 
            new_matches.append( mx )
            lx = mx.mQueryLength
            # seqx = self.mQueriesFasta.getSequence( mx.mQueryId )
            for y in range(x+1, len(matches)):
                if y in eliminated: continue
                my = matches[y]
                ly = my.mQueryLength

                # make sure they overlap
                overlap = alignlib.getAlignmentOverlap( \
                    mx.mMapTarget2Query, 
                    my.mMapTarget2Query, 
                    alignlib.RR )

                if overlap < self.mMinAlignmentOverlap: continue

                if mx.mSbjctFrom > my.mSbjctTo or mx.mSbjctTo < my.mSbjctFrom: continue
                if abs(mx.mSbjctFrom - my.mSbjctFrom) > self.mMinDistance or abs(mx.mSbjctTo - my.mSbjctTo) > self.mMinDistance: continue

                # compute 
                delta = abs( lx - ly )
                start = max(mx.mSbjctFrom,my.mSbjctFrom) - sbjct_start
                end = min(mx.mSbjctTo,my.mSbjctTo) - sbjct_start
                nreads = max( is_exon[start:end] )
                # P = 1.0 - ( 1.0 - float(delta+1) * min (lx, ly) / (lx * ly)) ** ((nreads * (nreads -1) ) / 2)
                # The factor 0.5 ensures that there is at least 1 residue in the interval if delta = 0.
                # are the brackets correct?
                PL = scipy.stats.norm.cdf( min(lx, lx + delta) + 0.5, loc = self.mLengthMean, scale = self.mLengthStd ) - \
                    scipy.stats.norm.cdf( max(0, lx - delta - 0.5), loc = self.mLengthMean, scale = self.mLengthStd )
                PS = float(delta+1) * min (lx, ly) / (lx * ly)
                # P = ( 1.0 - min(1.0, float(delta+1) * nreads * ( nreads-1 ) / 2 * min(lx,ly) / (lx * ly) ) ) * plx
                PSP = 1.0 - (1.0 - PS) ** ( nreads * (nreads-1) / 2.0 ) 
                PLP = 1.0 - (1.0 - PL) ** ( nreads * (nreads-1) / 2.0 )

                P = PLP * PSP
                # if abs( mx.mSbjctFrom - my.mSbjctFrom ) + abs(mx.mSbjctTo - my.mSbjctTo ) > self.mMinLengthDifference: continue
                # if abs(mx.mQueryLength - my.mQueryLength) > self.mMinLengthDifference: continue
                # seqy = self.mQueriesFasta.getSequence( my.mQueryId )
                if P <= self.mMinProbability:
                    self.mOutFile.write( "%s\t%s\t%i\t%i\t%i\t%i\t%i\t%f\t%f\t%f\n" % (mx.mQueryId, my.mQueryId, mx.mQueryLength, my.mQueryLength, overlap, delta, nreads, PL, PS, P ))
                    eliminated.add(y)
                    self.mOutFile.flush()

        if options.loglevel >= 3:
            options.stdlog.write( "# %s - %s:%i..%i before=%i after=%i\n" %\
                                      (self.mName,
                                       sbjct_id, 
                                       sbjct_start,
                                       sbjct_end,
                                       len(matches),
                                       len(new_matches)))


        return new_matches

class CoverageInfo: pass

def readMapId2Coverage( infile ):
    """read map of id to coverage information."""
    map_id2coverage = {}
    for line in infile:
        if line.startswith( "id" ): continue
        if line.startswith( "#" ): continue
        data = line[:-1].split("\t")
        x = CoverageInfo()
        id, x.mNMatches, x.mNCovered, x.mMean, x.mStddev = data[0], int(data[5]), int(data[6]), float(data[10]), float(data[12])
        map_id2coverage[id] = x

    return map_id2coverage

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: psl2assembly.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("--input-filename-queries", dest="input_filename_queries", type="string",
                      help="fasta filename with queries [default=%default]."  )

    parser.add_option("--input-filename-coverage", dest="input_filename_coverage", type="string",
                      help="tabular table with coverage information [default=%default]."  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename with histogram information on aggregate coverages [%default].")

    parser.add_option( "--mali-output-format", dest="mali_output_format", type="choice" ,
                       choices=( "fasta", ),
                       help="output format to choose [default=%default].")

    parser.add_option( "--method", dest="methods", type="choice", action="append",
                       choices=( "region", "consensus", "pileup", "transcript", "coverage", "sbjct_coverage", "indels", "polyA", "locus" ),
                       help="methods to apply [%default].")

    parser.add_option( "-z", "--from-zipped", dest="from_zipped", action="store_true",
                       help="input is zipped.")

    parser.add_option( "--output-pileup", dest="output_pileup", action="store_true",
                       help="in the transcript modules, output the pileup alignment with the predicted transcript [default=%default].")

    parser.add_option( "--test", dest="test", type="int",
                       help="test - stop after # rows of parsing[%default]." )

    parser.add_option( "--staggered", dest="staggered_alignments", type="choice",
                       choices=("all", "none", "merge" ),
                       help="how to deal with staggered alignments[%default]." )

    parser.add_option( "--force-merge", dest="force_merge", type="int",
                       help="in case of staggered alignments, force merge if there are # alignments. This avoids costly computation of components in large sets. If 0, do not apply the threshold [%default]." )

    parser.add_option( "--filter", dest="filters", type="choice", action="append",
                       choices=("duplicates", "imperfect", "exon-extenders", "transcript-mergers"),
                       help="filters to apply [%default]." )

    parser.add_option( "--threshold-merge-distance", dest="threshold_merge_distance", type="int",
                       help="distance in nucleotides at which two adjacent reads shall be merged even if they are not overlapping [%default]." )

    parser.add_option( "--threshold-merge-overlap", dest="threshold_merge_overlap", type="int",
                       help="require reads to be overlapping by this amount to be joined into the same transcript [default=%default]." )

    parser.add_option( "--start-at", dest="start_at", type="int",
                       help="start numbering ids at number [default=%default]." )

    parser.set_defaults( input_filename_queries = None,
                         threshold_good_query_coverage = 90.0,
                         threshold_min_pid = 30.0,
                         output_filename_pattern = "%s",
                         mali_output_format = "fasta",
                         from_zipped = False,
                         genome_file = None,
                         test = None,
                         staggered_alignments = "none",
                         methods = [],
                         filters = [],
                         output_format = "%06i",
                         threshold_merge_distance = 1,
                         threshold_merge_overlap = 0,
                         output_pileup = False,
                         start_at = 1,
                         input_filename_coverage = None,
                         force_merge = None,
                         )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) == 1:
        if options.from_zipped or args[0][-3:] == ".gz":
            import gzip
            infile = gzip.open( args[0], "r" )
        else:
            infile = open(args[0], "r" )
    else:
        infile = sys.stdin

    genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
    if options.input_filename_queries:
        queries_fasta =  IndexedFasta.IndexedFasta( options.input_filename_queries )
    else:
        queries_fasta = None

    if options.input_filename_coverage:
        map_id2coverage = readMapId2Coverage( open( options.input_filename_coverage, "r" ) )
        E.info("read coverage information for %i ids" % len(map_id2coverage) )
    else:
        map_id2coverage = None

    ################################################
    ################################################
    ################################################
    ## pick a processor
    ################################################        
    methods = []

    if len(options.methods) == 0:
        raise ValueError("please supply at least one method to apply.")

    for method in options.methods:
        if method == "region":
            methods.append( BuilderRegion( genome_fasta, queries_fasta, options ) )
        elif method == "pileup":
            methods.append( BuilderPileUp( genome_fasta, queries_fasta, options) )
        elif method == "consensus":
            methods.append( BuilderConsensus( genome_fasta, queries_fasta, options) )
        elif method == "locus":
            methods.append( BuilderTranscribedLocus( genome_fasta, queries_fasta, options) )
        elif method == "transcript":
            methods.append( BuilderTranscriptSpliced( genome_fasta, queries_fasta, options) )
        elif method == "coverage":
            methods.append( BuilderCoverage( genome_fasta, queries_fasta, options) )
        elif method == "sbjct_coverage":
            methods.append( BuilderSbjctCoverage( genome_fasta, queries_fasta, options) )
        elif method == "indels":
            methods.append( BuilderIndels( genome_fasta, queries_fasta, options) )
        elif method == "polyA":
            if not queries_fasta:
                raise ValueError( "need --input-filename-queries for polyA tail prediction." )
            methods.append( BuilderPolyA( genome_fasta, queries_fasta, options) )

    for method in methods:
        if map_id2coverage:
            method.setMapId2Coverage( map_id2coverage )
        method.printHeader()

    ################################################
    ################################################
    ################################################
    ## pick filters
    ## pre filters change components, post filters
    ## don't
    ################################################        
    pre_filters = []
    post_filters = []
    for f in options.filters:
        if f == "imperfect":
            post_filters.append( FilterImperfect( genome_fasta, queries_fasta, options ) )
        elif f == "exon-extenders":
            post_filters.append( FilterExonExtenders( genome_fasta, queries_fasta, options ) )
        elif f == "duplicates":
            pre_filters = [ FilterDuplicates( genome_fasta, queries_fasta, options ) ] + pre_filters
        elif f == "transcript-mergers":
            pre_filters.append( FilterTranscriptMergers( genome_fasta, queries_fasta, options ) )

    ################################################
    ################################################
    ################################################
    ## processing of a chunk (matches of same query)
    ################################################        
    ninput, noutput, nskipped, ncomponents = 0, 0, 0, 0

    output_id = options.start_at
    nskipped_components = 0
    n_pre_filtered_reads = 0
    n_post_filtered_reads = 0

    options.stdout.write( "%s\t%s\t%s\n" % ( "id","read", "strand" ))

    def processChunk( sbjct_id, sbjct_start, sbjct_end, matches ):
        global ninput, noutput, nskipped, ncomponents
        global nskipped_components
        global n_pre_filtered_reads
        global n_post_filtered_reads
        global output_id

        if len(matches) == 0:
            raise ValueError("no matches")

        ninput += 1

        if options.loglevel >= 2:
            options.stdlog.write( "## processing region: %s:%i..%i size=%i, matches=%i\n" %\
                                      (sbjct_id, sbjct_start, sbjct_end, sbjct_end - sbjct_start, len(matches)))
            options.stdlog.flush()

        # calculate overlapping components
        addAlignments( matches, shift = -sbjct_start, by_query = False )

        if options.loglevel >= 2:
            options.stdlog.write( "## region: %s:%i..%i: built alignments for %i matches\n" % \
                                      (sbjct_id, sbjct_start, sbjct_end, len(matches)) )

            options.stdlog.flush()


        if options.force_merge and len(matches) >= options.force_merge:
            options.stdlog.write( "## region: %s:%i..%i: forcing merge for %i matches\n" % \
                                      (sbjct_id, sbjct_start, sbjct_end, len(matches)) )
            options.stdlog.flush()
            components = ( range(0,len(matches), ), )
        else:
            components = Blat.getComponents( matches, 
                                            max_distance = options.threshold_merge_distance,
                                            min_overlap = options.threshold_merge_overlap )
            
            if options.loglevel >= 2:
                options.stdlog.write( "## region: %s:%i..%i: computed %i components\n" % \
                                          (sbjct_id, sbjct_start, sbjct_end, len(components)) )

                options.stdlog.flush()

        redo_components = False
        if pre_filters:
            new_matches = []
            for component in components:

                m = [ matches[ x ] for x in component ]

                before = len(m)
                for f in pre_filters:
                    m = f( sbjct_id, sbjct_start, sbjct_end, m )
                    if not m: break
                after = len(m)
                
                n_pre_filtered_reads += before - after
                
                if before - after > 0: 
                    redo_components = True
                
                new_matches += m

            if redo_components:
                matches = new_matches
                if options.force_merge and len(matches) >= options.force_merge:
                    options.stdlog.write( "## region: %s:%i..%i: forcing merge for %i matches\n" % \
                                              (sbjct_id, sbjct_start, sbjct_end, len(matches)) )
                    options.stdlog.flush()
                    components = ( range(0,len(matches), ), )
                else:
                    components = Blat.getComponents( matches, 
                                                     max_distance = options.threshold_merge_distance,
                                                     min_overlap = options.threshold_merge_overlap )

        if len(components) > 1:
            if options.staggered_alignments == "none":
                if options.loglevel >= 1:
                    options.stdlog.write( "## WARNING: more than one component (%i) in region %s:%i..%i (%s) - all are ignored\n" %\
                                              (len(components), sbjct_id, sbjct_start, sbjct_end, 
                                               str( [ len(x) for x in components ] )))
                nskipped += 1
                return
            elif options.staggered_alignments == "all":
                if options.loglevel >= 1:
                    options.stdlog.write( "## WARNING: more than one component (%i) in region %s:%i..%i (%s) - all are output\n" %\
                                              (len(components), sbjct_id, sbjct_start, sbjct_end, 
                                               str( [ len(x) for x in components ] )))
        
            elif options.staggered_alignments == "merge":
                if options.loglevel >= 1:
                    options.stdlog.write( "## WARNING: more than one component (%i) in region %s:%i..%i (%s) - all are merged\n" %\
                                              (len(components), sbjct_id, sbjct_start, sbjct_end, 
                                               str( [ len(x) for x in components ] )))
                    
                components = [ range(len(matches) ), ]

        ## applying filters
        for component in components:
            
            ncomponents += 1
            m = [ matches[ x ] for x in component ]

            before = len(m)
            for f in post_filters:
                m = f( sbjct_id, sbjct_start, sbjct_end, m )
                if not m: break
            after = len(m)

            n_post_filtered_reads += before - after

            if not m:
                if options.loglevel >= 1:
                    options.stdlog.write( "## region %s:%i..%i - component %i - all matches removed by filters\n" %\
                                              (sbjct_id, sbjct_start, sbjct_end, ncomponents ) ) 
                nskipped_components += 1
                continue

            for p in methods:
                p( output_id, sbjct_id, sbjct_start, sbjct_end, m )
                
            for match in m:
                options.stdout.write( "%s\t%s\t%s\n" %\
                                          (options.output_format % output_id,
                                           match.mQueryId,
                                           match.strand ) )
                
            options.stdout.flush()
            output_id += 1
            noutput += 1

    ################################################
    ################################################
    ################################################
    ## main loop: 
    ##
    ## collect overlapping entries on sbjct.
    ## Requirement: the input needs to be sorted
    ## by sbjct_id and sbjct_from
    ################################################        
    nfully_covered = None
    matches = []
    last_sbjct_id = None
    sbjct_start = None
    sbjct_end = None
    ninput_lines = 0

    skip = 0

    iterator = Blat.BlatIterator( sys.stdin )

    processed_contigs = set()

    while 1:

        match = iterator.next()
        
        if match == None: break
        
        ninput_lines += 1

        if options.test and ninput_lines > options.test:
            break
        
        if match.mSbjctId != last_sbjct_id or match.mSbjctFrom >= (sbjct_end + options.threshold_merge_distance):
            if last_sbjct_id:
                processChunk( last_sbjct_id, 
                              sbjct_start, sbjct_end,
                              matches )
            matches = []

            if last_sbjct_id != match.mSbjctId and match.mSbjctId in processed_contigs:
                raise ValueError("input not sorted correctly (contig,start): already encountered %s\n%s" % (match.mSbjctId, str(match)))
            last_sbjct_id = match.mSbjctId
            processed_contigs.add( last_sbjct_id )

            sbjct_start = match.mSbjctFrom
            sbjct_end = match.mSbjctTo
            
        if match.mSbjctFrom < sbjct_start:
            raise ValueError("input not sorted correctly (contig,start): %i < %i\n%s" % (match.mSbjctFrom, sbjct_start, str(match)))
        sbjct_end = max( match.mSbjctTo, sbjct_end )
        matches.append(match)

    processChunk( last_sbjct_id, sbjct_start, sbjct_end, matches )

    for p in methods:
        p.finish()

    if options.loglevel >= 1:
        options.stdlog.write("# clusters: ninput=%i, noutput=%i, ninput_lines=%i, ncomponents=%i\n" %\
                                 (ninput, noutput, ninput_lines, ncomponents ) )
        options.stdlog.write("# skipped: queries=%i, components=%i, pre_reads=%i, post_reads=%i\n" %\
                                 (nskipped, nskipped_components,
                                  n_pre_filtered_reads, n_post_filtered_reads))

    E.Stop()
