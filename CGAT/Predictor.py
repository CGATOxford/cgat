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
Predictor.py - running gene predictions with exonerate
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, optparse, tempfile, time, subprocess

USAGE="""python %s [OPTIONS] peptide genome

Version: $Id: Predictor.py 2781 2009-09-10 11:33:14Z andreas $

Wrapper for running gene predictor on sequence chunks.

Format of gff file:

The last (9th) column containing the query identifier and an optional id. For example:

chr1    pseudogene.org  pseudogene (Processed)  31282205        31282680        .       +       .       Query Q9Y3D8; Id 74
chr1    pseudogene.org  pseudogene (Duplicated) 33225155        33225823        .       +       .       Query P15880; Id 77
chr1    pseudogene.org  pseudogene (Duplicated) 33226224        33226423        .       +       .       Query P15880; Id 77

or just

chr1    pseudogene.org  pseudogene (Processed)  31282205        31282680        .       +       .       Q9Y3D8
chr1    pseudogene.org  pseudogene (Duplicated) 33225155        33225823        .       +       .       P15880
chr1    pseudogene.org  pseudogene (Duplicated) 33226224        33226423        .       +       .       P15880

The file has to be sorted appropriately so that all the regions belonging to a group are adjacent.
"""

import Experiment as E
import Genomics
import IndexedFasta
import Prediction
import PredictionParser
try: import alignlib_lite
except ImportError: pass
import Cluster
import GTF
import threading

class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class ParsingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class Predictor:

    mLogLevel = 0
    
    def __init__(self):

        self.mDryRun = False

        self.mCluster = None

    def __call__( self,
                  filename_peptides,
                  filename_genome,
                  options = None):
        """run a single prediction using prediction program.

        return lines
        """

        options += " " + self.mInputOptions

        if not options:
            options += " " + self.mDefaultOptions

        statement = self.buildStatement( filename_peptides, filename_genome, options )

        if self.mLogLevel >= 4:
            print "# statement: %s" % statement
            sys.stdout.flush()

        if self.mDryRun:
            print statement
            return Prediction.Predictions()

        if self.mCluster:
            out, err = self.mCluster.runJob( statement )
        else:
            s = subprocess.Popen( statement,
                                  shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE,
                                  close_fds = True)                              

            (out, err) = s.communicate()

            if s.returncode != 0:
                print err
                print out
                raise IOError, "Error in running statement %s" % statement
        
        if self.mLogLevel >= 3:
            print"# received %i chars from %s" % (len(out), self.mExecutable)
            print out            

        results = self.mParser.parse( map( lambda x: x+ "\n", out.split("\n") ) )
        return results

    def setCluster( self, cluster ):
        """set cluster object."""
        self.mCluster = cluster

    def setOutputOptions(self, options):
        """set output options."""
        self.mOutputOptions = options 

    def setDryRun( self, flag = True):
        """setup dry run.

        Gene prediction is not run, but files are set up as if it were.
        """
        self.mDryRun = flag

    def setLogLevel( self, level ):
        """set loglevel."""
        self.mLogLevel = level


    def getParser( self ):
        """returns the parser object."""
        return self.mParser

    def predictFromSequences( self,
                              peptide_sequence,
                              genomic_sequence,
                              options = None ):
        """predict from two explicitely given sequences."""

        outfile, filename_peptides = tempfile.mkstemp()
        os.write( outfile, ">query\n%s" % peptide_sequence )
        os.close(outfile)

        outfile, filename_genome = tempfile.mkstemp()
        os.write( outfile, ">target\n%s" % genomic_sequence )
        os.close(outfile)

        self.mParser.addGenomicSequence( "target", genomic_sequence )
        result = self.__call__( filename_peptides, filename_genome, options )
        self.mParser.deleteGenomicSequence( "target" )
        os.remove( filename_peptides )
        os.remove( filename_genome )
        
        return result

    def predictFromRange( self,
                          peptide_sequence,
                          contig, strand, start, end,
                          fasta,
                          options = None ):
        """predict from a peptide sequence within a region."""

        
        genomic_sequence = fasta.getSequence( contig, strand, start, end )

        result = self.predictFromSequences( peptide_sequence, genomic_sequence, options )

        lcontig = fasta.getLength( contig )

        for r in result:

            # shift aligned region for fragment
            # this is strand dependent.
            r.shiftGenomicRegion( start, lcontig - start )
            ## set sbjct token to original
            r.mSbjctToken = contig
            
            ## sort out strand
            ## if match is on negative strand
            if Genomics.IsNegativeStrand( r.mSbjctStrand ):
                ## switch coordinates
                if Genomics.IsNegativeStrand( strand ):
                    strand = "+"
                else:
                    strand = "-"
                    
            r.mSbjctStrand = strand

        return result

class PredictorGenewise(Predictor):
    """run genewise as predictor.

    genewise only runs a single peptide sequence against a single
    DNA sequences. genewisedb runs collections of peptide sequences
    against collections of DNA sequences. Note that this might take a while
    so it is prudent to spread it over the cluster.
    """
    
    def __init__( self, parser = None ):
        Predictor.__init__(self)
        if parser:
            self.mParser = parser
        else:
            self.mParser = PredictionParser.PredictionParserGenewise()
            
        self.mExecutable = "genewise"
        self.mDefaultOptions="-pseudo -init endbias"
        self.mOutputOptions="-quiet -sum -gff -trans -pep -alb "
        self.mInputOptions = ""

    def buildStatement( self, filename_peptides, filename_genome, options ):
        return string.join( map(str, (
            self.mExecutable,
            options,
            self.mOutputOptions,
            filename_peptides,
            filename_genome,
            )), " " )

class PredictorGenewiseHMM(Predictor):
    
    def __init__( self, parser = None ):
        Predictor.__init__(self)
        if parser:
            self.mParser = parser
        else:
            self.mParser = PredictionParser.PredictionParserGenewiseHMM()
            
        self.mExecutable = "genewise -hmmer"
        self.mDefaultOptions="-pseudo -alg 623L "
        self.mInputOptions = ""
        self.mOutputOptions="-quiet -sum -gff -trans -pep -alb "

    def buildStatement( self, filename_peptides, filename_genome, options ):

        ## allow parser to extract query length from the
        ## submitted hmm.
        self.mParser.extractQueryLength( filename_peptides )
        
        return string.join( map(str, (
            self.mExecutable,
            options,
            self.mOutputOptions,
            filename_peptides,
            filename_genome,
            )), " " )
                        
class PredictorExonerate ( Predictor ):

    def __init__( self, parser = None ):
        
        Predictor.__init__(self)
        
        self.mExecutable = "exonerate"

        self.mInputOptions = "-Q protein -T dna --softmasktarget TRUE"
        self.mDefaultOptions = "-m p2g --forcegtag TRUE"
        self.mOutputOptions = '--forwardcoordinates FALSE --showalignment FALSE --showvulgar FALSE --showsugar FALSE --showcigar FALSE ' + \
                              ' --ryo "diy\t%S\t%ql\t%r\t%pi\t%ps\t%V\n" --showtargetgff FALSE --showquerygff FALSE'
        
        ## only predictions on this strand
        self.mRestrictToStrand = None

        if parser:
            self.mParser = parser
        else:
            self.mParser = PredictionParser.PredictionParserExonerate()

    def buildStatement( self, filename_peptides, filename_genome, options ):
        return string.join( map(str, (
            self.mExecutable,
            options,
            self.mOutputOptions,
            "--query",
            filename_peptides,
            "--target",
            filename_genome,
            )), " " )

##--------------------------------------------------------------------------
class PredictionParser :
    """base class of prediction parsers."""
    mLogLevel = 0
    
    def __init__( self,
                  peptides = None,
                  peptide = None,
                  genomes = None,
                  genome = None,
                  filename_genome = None,
                  filename_peptides = None,
                  indexed_genome = None,
                  ):
        """pass sequence data."""
        self.mPeptides = peptides
        self.mPeptide = peptide
        self.mGenomes = genomes
        self.mGenome = genome        
        self.mIndexedGenome = indexed_genome
        self.mFilenamePeptides = filename_peptides
        self.mFilenameGenome = filename_genome

        ## stop codons to ignore at boundaries
        self.mBorderStopCodons = 0

    def parse( self, lines ):
        """parse lines.
        """
        return Prediction.Predictions()

    def getPeptideSequence( self, query_token ):
        """returns peptide sequence for query_token."""


        if self.mPeptides:
            return self.mPeptides[query_token]
        elif self.mFilenamePeptides:
            self.mPeptides = Genomics.ReadPeptideSequences( open( self.mFilenamePeptides, "r") )
            return self.mPeptides[query_token]
        elif self.mPeptide:
            return self.mPeptide

    def addGenomicSequence( self, sbjct_token, sequence, strand = "+" ):
        """add a sequence so that the parser nkows about it.
        """
        if strand == "-":
            raise "not implemented, please add positive strand"

        if self.mGenomes:
            if sbjct_token in self.mGenomes:
                raise KeyError, "key %s already in genomes" % (sbjct_token)
            self.mGenomes[sbjct_token] = sequence
        elif self.mFilenameGenome:
            raise "not implemented"
        elif self.mIndexedGenome:
            raise "not implemented"
        else:
            self.mGenome = sequence

    def deleteGenomicSequence( self, sbjct_token ):
        """remove a genomic sequence from parser."""
        
        if self.mGenomes:
            if sbjct_token in self.mGenomes:
                del self.mGenomes[sbjct_token]
        elif self.mFilenameGenome:
            raise "not implemented"
        elif self.mIndexedGenome:
            raise "not implemented"

    def getGenomicSequence( self, sbjct_token, sbjct_strand = None, sbjct_from = None, sbjct_to = None ):
        """returns a genomic sequence for given coordinates."""

        if self.mGenomes:
            genomic_sequence = self.mGenomes[sbjct_token]
        elif self.mFilenameGenome:
            self.mGenomes = Genomics.ReadGenomicSequences( open(self.mFilenameGenome, "r"), do_reverse = False )
            genomic_sequence = self.mGenomes[sbjct_token]
        elif self.mIndexedGenome:
            return self.mIndexedGenome.getSequence( sbjct_token, sbjct_strand, sbjct_from, sbjct_to )
        else:
            genomic_sequence = self.mGenome

        if sbjct_strand in ("-", "0", "-1", -1):
            genomic_sequence = string.translate( genomic_sequence, string.maketrans("ACGTacgt", "TGCAtgca") )[::-1]

        if sbjct_from != None and sbjct_to != None:
            return genomic_sequence[sbjct_from:sbjct_to]
        else:
            return genomic_sequence

    def getGenomicSequenceLength( self, sbjct_token ):
        """returns a genomic sequence for given coordinates."""

        if self.mGenomes:
            return len(self.mGenomes[sbjct_token])
        elif self.mFilenameGenome:
            self.mGenomes = Genomics.ReadGenomicSequences( open(self.mFilenameGenome, "r") )
            return len(self.mGenomes[sbjct_token])            
        elif self.mIndexedGenome:
            return self.mIndexedGenome.getLength( sbjct_token )
        else:
            return len(self.mGenome)
    
##--------------------------------------------------------------------------
class PredictionParserGenewise( PredictionParser ):
    """Parse output from genewise.
    """
    
    def __init__(self, **kwargs):
        
        PredictionParser.__init__( self, **kwargs )

        self.mParsers = [ self.parseSummary, self.parseTranslation,  self.parsePeptide, self.parseGFF, self.parseAlb ]

        self.mQueryLength = None
        
    def parse( self, lines ):
        """call the various parsers for a chunk of data.
        """
        if len(lines) == 0: return 0

        matches = Prediction.Predictions()
        self.mMatches = matches
        self.mTranslations = []
        self.mParsedPeptides = []
        self.mAlignments = []
        self.mSummaries = []
        self.strands = []

        # for line in lines:
        # print line,
        
        ###############################################################################
        ## call the various parsers for chunks separated by "//"
        chunks = filter( lambda x: lines[x][:2] == "//", range( len(lines )) )

        if len(chunks) % len(self.mParsers) != 0:

            for x in range(len(chunks)-1):
                print "#########################################################"
                print "".join(lines[chunks[x]+1:chunks[x+1]])
            
            raise ParsingError( "lengths of sections do not not fit: found %i, expected %i" % (len(chunks), len(self.mParsers)))

        chunks = [-1] + chunks
        
        chunk = 0
        while chunk < len(chunks) - 1:
            
            for parser in self.mParsers:
                parser(lines[chunks[chunk]+1:chunks[chunk+1]])
                chunk += 1

        ###############################################################################
        # build the Prediction entries
        for x in range(len(self.mSummaries)):
            
            entry = Prediction.Prediction()

            summary = self.mSummaries[x]

            (entry.mQueryToken, entry.mSbjctToken, entry.score,
             entry.mQueryFrom, entry.mQueryTo, 
             entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
             entry.mNIntrons) = \
                            (summary.mQueryToken, summary.mSbjctToken, summary.score,
                             summary.mQueryFrom, summary.mQueryTo, 
                             summary.mSbjctGenomeFrom, summary.mSbjctGenomeTo,
                             summary.mIntrons)

            entry.mMapPeptide2Genome, entry.mMapPeptide2Translation, nmatches, nindels = self.mAlignments[x]
            
            ## in case of pseudogenes, translations are empty, but peptides are there,
            ## so using those.
            entry.mTranslation = self.mParsedPeptides[x][1]

            entry.mSbjctStrand = self.strands[x]
            peptide_sequence = self.getPeptideSequence( entry.mQueryToken )

            if peptide_sequence:
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                alignlib_lite.py_rescoreAlignment( entry.mMapPeptide2Translation, row_seq, col_seq )

                entry.mQueryLength = len(peptide_sequence)            
                entry.mPercentIdentity = alignlib_lite.py_calculatePercentIdentity( entry.mMapPeptide2Translation, row_seq, col_seq ) * 100
                entry.mPercentSimilarity = alignlib_lite.py_calculatePercentSimilarity( entry.mMapPeptide2Translation ) * 100
                entry.mQueryCoverage = ( entry.mMapPeptide2Translation.getRowTo() - \
                                         entry.mMapPeptide2Translation.getRowFrom() + 1 ) * 100 /\
                                         entry.mQueryLength
            elif self.mQueryLength:
                ## for hmms, query lenght can be extraced from the HMM
                entry.mQueryLength = self.mQueryLength
                entry.mQueryCoverage = ( entry.mMapPeptide2Translation.getRowTo() - \
                                         entry.mMapPeptide2Translation.getRowFrom() + 1 ) * 100 /\
                                         entry.mQueryLength
                ## percent identity is given by the number of matches versus alignment length
                entry.mPercentIdentity = float(nmatches) / (nmatches + nindels)
                entry.mPercentSimilarity = entry.mPercentIdentity

            entry.mAlignmentString = string.join( map( \
                                      lambda x: string.join(map(str, x), " "),
                                      entry.mMapPeptide2Genome), " ")

            if len(entry.mMapPeptide2Genome) == 0:
                print "### PARSING ERROR: empty alignment"
                print str(entry)
                print string.join(lines, "")
                sys.exit(1)

            # get regions matching
            # use summary information, as it is correct for genomic data.
            # when -u is set (in contrast to -alb format)
            # genewise starts counting at 1, thus subtract 1 from first position
            # Note: When position starts with a gap, genewise says alignments starts
            # from position 1 while the first aligned residue is position 2
            if entry.mMapPeptide2Genome[0][0] == "G": entry.mQueryFrom += entry.mMapPeptide2Genome[0][1]
            if entry.mMapPeptide2Genome[-1][0] == "G": entry.mQueryTo -= entry.mMapPeptide2Genome[0][1]

            # a sanity check
            if entry.mQueryFrom != entry.mMapPeptide2Translation.getRowFrom():
                print "## PARSING ERROR: wrong start at %s vs %s: %i <> %i" %\
                      (entry.mQueryToken, entry.mSbjctToken,
                       entry.mQueryFrom,
                       entry.mMapPeptide2Translation.getRowFrom() )
                print str(entry)
                print string.join(lines,"")
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                print str( alignlib_lite.py_AlignmentFormatExplicit( entry.mMapPeptide2Translation, row_seq, col_seq ))
                sys.exit(1)
            if entry.mQueryTo != entry.mMapPeptide2Translation.getRowTo() and \
                   entry.mMapPeptide2Genome[-1][0] != "G":
                print "## PARSING ERROR: wrong end at %s vs %s: %i <> %i" %\
                  (entry.mQueryToken, entry.mSbjctToken,
                   entry.mQueryTo,
                   entry.mMapPeptide2Translation.getRowTo() )
                print str(entry)                
                print string.join(lines,"")
                row_seq = alignlib_lite.py_makeSequence( peptide_sequence )
                col_seq = alignlib_lite.py_makeSequence( entry.mTranslation )
                print str( alignlib_lite.py_AlignmentFormatExplicit( entry.mMapPeptide2Translation, row_seq, col_seq ))
                sys.exit(1)

            ## fix coordinates
            
            ## sic: on neg strand genewise has high,low coordinates
            if entry.mSbjctStrand == "-":
                lgenome = self.getGenomicSequenceLength( entry.mSbjctToken )                
                entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo = lgenome - entry.mSbjctGenomeFrom, lgenome - entry.mSbjctGenomeTo
                entry.mSbjctGenomeTo += 1
            else:
                entry.mSbjctGenomeFrom -= 1

            genomic_sequence = self.getGenomicSequence( entry.mSbjctToken, entry.mSbjctStrand,
                                                        entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo)

            if genomic_sequence:
                (entry.mNIntrons, entry.mNFrameShifts, entry.mNGaps, entry.mNSplits, entry.mNStopCodons, disruptions) = \
                                  Genomics.CountGeneFeatures( 0,
                                                              entry.mMapPeptide2Genome,
                                                              genomic_sequence,                                                          
                                                              self.mBorderStopCodons )
                
            
            matches.append( entry )
            

        matches.sort( lambda x,y : cmp( -x.score, -y.score))
        
        return matches
    
    def parseGFF( self, lines ):
        """extract strands from the gff section."""
        for line in lines:
            data = re.split( "\s+", line[:-1])
            if data[2] == "match":
                self.strands.append( data[6] )

    def parseSummary( self, lines ):
        """parse summary entry of genewise."""

        class SummaryResult:
            pass

        ## first line is header
        for line in lines[1:]:
            data = re.split( "\s+", line[:-1])
            entry = SummaryResult()
            entry.score = float(data[0])
            entry.mQueryToken = data[1]
            entry.mQueryFrom = int(data[2])
            entry.mQueryTo = int(data[3])
            entry.mSbjctToken = data[4]
            entry.mSbjctGenomeFrom = int(data[5])
            entry.mSbjctGenomeTo = int(data[6])
            entry.mIndels = int(data[7])
            if data[8] == "":
                entry.mIntrons = 0
            else:
                entry.mIntrons = int(data[8])
                
            self.mSummaries.append( entry )

    def parseTranslation( self, lines ):
        """read a translation entry from genewise.        
        """
        ## There might be some problems with multiple entries.
        ##
        self.mTranslations += self.__parseFasta( lines )

    def parsePeptide( self, lines ):
        """read a peptide entry from genewise.
        """
        self.mParsedPeptides += self.__parseFasta( lines )
        
    def __parseFasta( self, lines ):
        """read a fasta entry from genewise.

        stores tuples of description,sequence in result.
        """

        results = []
        fragments = []
        description = ""
        for line in lines:
            if line[0] in (">", "#"):
                if fragments:
                    results.append( (description, "".join(fragments) ) )
                    
                if line[0] != "#":
                    description = line[1:-1]
                else:
                    description = "pseudogene"
                    
                fragments = []
                continue
            elif line[:6] == "Making":
                continue
            fragments.append(line[:-1])

        if fragments:
            results.append( ( description, "".join(fragments) ) )

        return results
        
    def parseAlb( self, lines ):
        """Parse genewise logical block alignment and return aligment.

        ( alignment of peptide to genome,
          alignment of peptide to translation )

        """

        map_genewise2code = \
                          { ( "MATCH_STATE", "CODON") : "M",
                            ( "INTRON_STATE", "5SS_PHASE_0"): "5",
                            ( "INTRON_STATE", "5SS_PHASE_1"): "5",                        
                            ( "INTRON_STATE", "5SS_PHASE_2"): "5",
                            ( "INTRON_STATE", "CENTRAL_INTRON"): "I",
                            ( "MATCH_STATE", "3SS_PHASE_0"): "3",
                            ( "MATCH_STATE", "3SS_PHASE_1"): "3",
                            ( "MATCH_STATE", "3SS_PHASE_2"): "3",                        
                            ( "INSERT_STATE", "3SS_PHASE_0") : "3",
                            ( "INSERT_STATE", "3SS_PHASE_1") : "3",
                            ( "INSERT_STATE", "3SS_PHASE_2") : "3",
                            ( "INTRON_MATCH_1", "5SS_PHASE_0"): "5",
                            ( "INTRON_MATCH_1", "5SS_PHASE_1"): "5",
                            ( "INTRON_MATCH_1", "5SS_PHASE_2"): "5",                            
                            ( "INTRON_MATCH_1", "CENTRAL_INTRON"): "I",
                            ( "INTRON_MATCH_1", "PYRIMIDINE_TRACT"): "I",
                            ( "DELETE_STATE", "INSERT") : "G",
                            ( "INSERT_STATE", "CODON") : "G",
                            ( "INSERT_STATE", "SEQUENCE_DELETION") : "F", # a codon with less than 3 nucleotides at an insert
                            ( "MATCH_STATE", "SEQUENCE_INSERTION") : "F", # a codon with more than 3 nucleotides
                            ( "MATCH_STATE", "SEQUENCE_DELETION")  : "F", # a codon with less than 3 nucleotides
                            ( "LOOP_STATE", "RANDOM_SEQUENCE") : None,
                            ( "GENOMIC_RND_STATE", "RANDOM_SEQUENCE") : None,
                            ( "GENOMIC_RND_STATE", "CODON" ): None,
                            ( "END", "END") : None,
                            }


        ## note: all the numbers including the ranges can be negative (-1 for first position) !
        rx = re.compile( "^(\S+)\s+\[(\S+):(\S+)\s+\"(\S+)\"\s+(\S+)\],\[(\S+):(\S+)\s+\"(\S+)\"\s+(\S+)\]" )
        
        query_last_state = None
        query_first_from = 0

        sbjct_last_state = None
        sbjct_genomic_residue_first_from = 0
        sbjct_peptide_residue_first_from = 0

        states = []
        total_bin_score = 0.0

        nintrons = 0
        nframeshifts = 0
        nphases = [0,0,0]
        ngaps = 0

        ## parsing errors
        nerrors = 0

        for line in lines:

            try:
                (bin_score,
                 query_from, query_to, query_state, query_score,
                 sbjct_from, sbjct_to, sbjct_state, sbjct_score ) =  rx.match( line[:-1] ).groups()
            except AttributeError:
                raise ParsingError, line

            query_from, query_to, sbjct_from, sbjct_to = map(int, (query_from, query_to, sbjct_from, sbjct_to))

            total_bin_score += string.atof( bin_score )

            if query_last_state != query_state or sbjct_last_state != sbjct_state:

                if query_last_state:
                    states.append( ( query_last_state, query_first_from, query_from,
                                     sbjct_last_state, sbjct_first_from, sbjct_from ) )

                
                query_first_from = query_from
                sbjct_first_from = sbjct_from
                query_last_state = query_state
                sbjct_last_state = sbjct_state

        states.append( ( query_last_state, query_first_from, query_to,
                         sbjct_last_state, sbjct_first_from, sbjct_to) )

        ## translate genewise alignment to exonerate alignment codes
        n = 0
        
        ## residue positions on the query/sbjct for peptide sequences
        ## note: genewise would start the first state at -1:0. Thus, add
        ## 2 to get first aligned position in peptide
        query_peptide_from = states[0][1] + 2
        sbjct_peptide_from = 1

        map_peptide2genome = []
        map_query2sbjct = alignlib_lite.py_makeAlignmentVector()

        ## count number of matches/indels as proxy for pid
        ## See above for a loop, where counting could be done
        ## by scores.
        nmatches, nindels = 0, 0
        
        for x in range( len(states)):

            query_state, query_from, query_to, sbjct_state, sbjct_from, sbjct_to = states[x]

            query_len = query_to - query_from
            sbjct_len = sbjct_to - sbjct_from
            
            ## build alignment between peptide sequence and genomic sequence
            try:
                code = map_genewise2code[ (query_state, sbjct_state)]
            except KeyError:
                print "## ERROR: unknown genewise code: %s %s" % (query_state, sbjct_state)
                nerrors += 1
                code = "?"

            if query_state == "LOOP_STATE" and x > 0:
                if map_peptide2genome:
                    self.mAlignments.append( (map_peptide2genome, map_query2sbjct, nmatches, nindels) )
                map_query2sbjct = alignlib_lite.py_makeAlignmentVector()
                map_peptide2genome = []
                query_peptide_from = states[x+1][1] + 2
                sbjct_peptide_from = 1
                nmatches = 0
                nindels = 0
                continue

            ## add code for split codons
            ## move matches to second split codon, if
            ## a split codon is on a gap: add a Gap
            if sbjct_state == "5SS_PHASE_1":
                map_peptide2genome.append( ("S", 0, 1) )
                sbjct_len -= 1
            elif sbjct_state == "5SS_PHASE_2":
                map_peptide2genome.append( ("S", 0, 2) )
                sbjct_len -= 2
            elif sbjct_state == "3SS_PHASE_1":
                query_len = 0
                sbjct_len -= 2
            elif sbjct_state == "3SS_PHASE_2":
                query_len = 0                
                sbjct_len -= 1

            if code:
                if query_from != query_peptide_from - 2:
                    print "Mismatch", query_state, query_from, query_to, sbjct_state,
                    sbjct_from, sbjct_to, query_peptide_from
                
                map_peptide2genome.append( (code, query_len, sbjct_len) )

            if sbjct_state == "3SS_PHASE_1":
                map_peptide2genome.append( ("S", 1, 2) )
                sbjct_len = 2
                if query_state == "MATCH_STATE":
                    query_len = 1
                else:
                    query_len = 0
                    
            elif sbjct_state == "3SS_PHASE_2":
                map_peptide2genome.append( ("S", 1, 1) )                
                sbjct_len = 1
                if query_state == "MATCH_STATE":
                    query_len = 1
                else:
                    query_len = 0
                
            ## build alignment between peptide sequences
            query_increment = 0
            sbjct_increment = 0

            if query_state == "MATCH_STATE":
                 query_increment = query_len
                 nmatches += query_len
                 
            elif query_state == "DELETE_STATE":
                query_increment = query_len
                nindels += query_len
                
            elif query_state == "INSERT_STATE":
                query_increment = query_len
                nindels += query_len
                
            elif query_state == "UNKNOWN_LABEL":
                query_peptide_from += query_len
                ## there are three X's in the sbjct
                sbjct_peptide_from += 3             
                continue
            elif query_state == "LOOP_STATE":
                query_increment += query_len

            if sbjct_state == "CODON":
                sbjct_increment = sbjct_len 
            elif sbjct_state in ("3SS_PHASE_1", "3SS_PHASE_2"):
                sbjct_increment = 3
            elif sbjct_state == "INSERT":
                sbjct_increment = sbjct_len
            elif sbjct_state == "RANDOM_SEQUENCE":
                sbjct_increment = 0

            ## for a frameshift: genewise outputs two characters in the
            ## translated sequence? I do not know. 3 makes more sense.
            elif sbjct_state == "SEQUENCE_DELETION":
                sbjct_increment = 3

            ## for a frameshift: if there is a sequence insertion, no
            ## character is output.
            elif sbjct_state == "SEQUENCE_INSERTION":
                if sbjct_len > 3:
                    sbjct_increment = 3
                else:
                    sbjct_increment = 0

            if code and query_increment and sbjct_increment:
                alignlib_lite.py_addDiagonal2Alignment( map_query2sbjct, 
                                               query_peptide_from, query_peptide_from + query_increment - 1,
                                               sbjct_peptide_from - query_peptide_from
                                               )
                
            query_peptide_from += query_increment
            sbjct_peptide_from += sbjct_increment / 3
            
        if map_peptide2genome:
            self.mAlignments.append( (map_peptide2genome, map_query2sbjct, nmatches, nindels) )
    
##--------------------------------------------------------------------------
class PredictionParserGenewiseHMM( PredictionParserGenewise ):
    """Parse output from genewise.

    This subclass reads the profile from the file to extract
    the length of a query.
    """
    
    def __init__(self, **kwargs):
        
        PredictionParserGenewise.__init__( self, **kwargs )
        
    def extractQueryLength( self, filename_peptides ):
        """extract query length from HMM."""

        infile = open(filename_peptides, "r")
        for line in infile:
            if line[:4] == "LENG":
                self.mQueryLength = int(line[5:-1])
                break
        else:
            self.mQueryLength = None
        infile.close()
        
##--------------------------------------------------------------------------
class PredictionParserExonerate (PredictionParser):
    """Parse output from exonerate.
    """
    
    def __init__(self, **kwargs):
        PredictionParser.__init__( self, **kwargs )

        self.mRestrictToStrand = None

    def SetRestrictToStrand( self, v ):
        self.mRestrictToStrand = v
            
    def parse( self, lines ):

        matches = Prediction.Predictions()
        
        nmatches = 0
        for line in lines:
            
            if not re.match("diy", line): continue

            nmatches += 1
            data = string.split( line[:-1], "\t")

            entry = Prediction.Prediction()
            entry.mRank = nmatches 

            (sugar, entry.mQueryLength, rank, percent_identity, percent_similarity, vulgar ) = data[1:]

            (entry.mQueryToken, entry.mQueryFrom, entry.mQueryTo, entry.mQueryStrand, 
             entry.mSbjctToken, entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
             entry.mSbjctStrand, entry.score) = string.split(sugar, " " )

            if self.mRestrictToStrand and entry.mSbjctStrand != self.mRestrictToStrand:
                continue

            try:
                (entry.mQueryFrom, entry.mQueryTo,
                 entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                 entry.score, entry.mQueryLength) = map( \
                    int , (entry.mQueryFrom, entry.mQueryTo,
                           entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                           entry.score, entry.mQueryLength))
                
            except ValueError:
                raise InputError, "error while parsing ints: %s" % " ".join( (entry.mQueryFrom, entry.mQueryTo,
                                                                            entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo,
                                                                            entry.score, entry.mQueryLength) )
            
            try:            
                entry.mPercentIdentity, entry.mPercentSimilarity = map( float,
                                                                        (percent_identity, percent_similarity ))
            except ValueError:
                raise InputError, "error while parsing floats: %s, %s" % (percent_identity, percent_similarity)
            
            entry.mAlignmentString = vulgar
            try:
                entry.mMapPeptide2Genome = Genomics.String2Alignment( vulgar )
            except ValueError:
                raise InputError, "error while parsing alignment %s" % vulgar
            
            entry.mQueryCoverage = (entry.mQueryTo - entry.mQueryFrom) * 100 / entry.mQueryLength

            genomic_sequence = self.getGenomicSequence( entry.mSbjctToken, entry.mSbjctStrand,
                                                        entry.mSbjctGenomeFrom, entry.mSbjctGenomeTo )

            entry.mMapPeptide2Translation, entry.mTranslation = Genomics.Alignment2PeptideAlignment( \
                entry.mMapPeptide2Genome, entry.mQueryFrom, 0, genomic_sequence )

            ## add 1 for peptide alignments
            entry.mQueryFrom += 1

            if genomic_sequence:
                (entry.mNIntrons, entry.mNFrameShifts, entry.mNGaps, entry.mNSplits, entry.mNStopCodons, disruptions) = \
                                  Genomics.CountGeneFeatures( 0,
                                                              entry.mMapPeptide2Genome,
                                                              genomic_sequence,
                                                              self.mBorderStopCodons )                                                              

            matches.append(entry)
            
        return matches

##--------------------------------------------------------------------------
class PredictionParserEcho (PredictionParser):
    """simply output predictions - do not parse."""

    
    def __init__(self, **kwargs):
        PredictionParser.__init__( self, **kwargs )

    def parse( self, lines ):
        """call the various parsers for a chunk of data.
        """

        for line in lines:
            print line[:-1]

        return Prediction.Predictions()        

##--------------------------------------------------------------------------
class TandemPredictor:
    """run a tandem prediction

    If a match is found, look for further related matches upstream or
    downstream.

    This is a hack that works for one additional query downstream or
    upstream. To be fully workable, sort out reverse strand predictions
    and overlapping predictions. Also, make sure that if several results
    are obtained, they are not spanning each other.
    """

    def __init__(self, predictor,
                 upstream_queries = [],
                 downstream_queries = []):
        """run predictor on queries."""
        self.mPredictor = predictor

        self.mDownStreamQueries = downstream_queries
        self.mUpStreamQueries = upstream_queries        

        ## maximal search region for additional elements
        self.mMaxExtension = 5000
        ## minimum distance between consecutive elements
        self.mMinDistance  = 0

        ## minimum length region to search
        self.mMinGenomeSize = 100

        self.mLogLevel = 0

    def setLogLevel( self, level ):
        self.mLogLevel = level

    def runTandem( self,
                   filename_peptides,
                   filename_genome,
                   options,
                   sbjct_token, sbjct_strand,
                   start, end ):

        """run a tandem prediction with region given by start, end."""

        parser = self.mPredictor.getParser()
        lgenome = parser.getGenomicSequenceLength( sbjct_token )

        start = max( 0, start )
        end = min( end, lgenome )

        if end - start < self.mMinGenomeSize:
            return None

        if self.mLogLevel >= 2:
            sys.stdout.write("# running tandem prediction on fragment %s:%s:%i:%i\n" % (sbjct_token,
                                                                                        sbjct_strand,
                                                                                        start,
                                                                                        end))

        ## patch - tempdir is loosing temporary directory assignment
        ## why is that?
        tempfile.tempdir = "."
        new_outfile_genome, new_filename_genome = tempfile.mkstemp()

        ## Note: genewise fails with long contig names, thus
        ## truncate the path.
        key = os.path.basename(new_filename_genome)
        genomic_sequence =  parser.getGenomicSequence( sbjct_token,
                                                       sbjct_strand,
                                                       start, end )
        os.write(new_outfile_genome, ">%s\n%s\n" % (key, genomic_sequence) )

        
        os.close( new_outfile_genome )
        parser.addGenomicSequence( key, genomic_sequence )
                                   
        result = self.mPredictor( filename_peptides, new_filename_genome, options )
        result.shiftGenomicRegion( start, lgenome - end )
        
        for r in result:
            r.mSbjctToken = sbjct_token
            
        parser.deleteGenomicSequence( key )

        if self.mLogLevel < 2:
            os.remove( new_filename_genome )
        
        if len(result) == 0:
            return None
        else:
            # return only the best match
            return result.getBestMatch()
    
    def __call__( self,
                  filename_peptides,
                  filename_genome,
                  options = None):
        """call runRecursive."""

        ## do predictions to get anchors
        results = self.mPredictor( filename_peptides, filename_genome, options )

        if not results: return results

        ## get parser object
        parser = self.mPredictor.getParser()

        ## place to store new results
        new_results = Prediction.Predictions()

        ## sort results by location
        results.sort( lambda x,y: cmp( (x.mSbjctToken,x.mSbjctStrand,x.mSbjctGenomeFrom),
                                       (y.mSbjctToken,y.mSbjctStrand,y.mSbjctGenomeFrom) ) )


        last_sbjct_token, last_sbjct_strand = None, None
        
        for result in results:

            if self.mLogLevel >= 2:
                sys.stdout.write("# extending prediction: %s\n" % (str(result) ))

            sbjct_token, sbjct_strand = result.mSbjctToken, result.mSbjctStrand
            
            if last_sbjct_token != sbjct_token or last_sbjct_strand != sbjct_strand:
                
                lgenome = parser.getGenomicSequenceLength( sbjct_token )
                start = 0
                end   = lgenome

            last_sbjct_token = sbjct_token
            last_sbjct_strand = sbjct_strand
            
            upstream_results = Prediction.Predictions()

            start = max( start, result.mSbjctGenomeFrom - self.mMinDistance + self.mMaxExtension )

            for filename_query in self.mUpStreamQueries:
                r = self.runTandem( filename_query,
                                    filename_genome,
                                    options,
                                    sbjct_token, sbjct_strand,
                                    start,
                                    max( start, result.mSbjctGenomeFrom - self.mMinDistance) )

                if r:
                    if self.mLogLevel >= 2:
                        sys.stdout.write("# upstream extension successful: %s\n" % (str(r) ))
                                             
                    upstream_results.add( r )
                else:
                    if self.mLogLevel >= 2:
                        sys.stdout.write("# upstream extension unsuccessful.\n" )
                                             
            downstream_results = Prediction.Predictions()
            for filename_query in self.mDownStreamQueries:
                r = self.runTandem( filename_query,
                                    filename_genome,
                                    options,
                                    sbjct_token, sbjct_strand,
                                    min( lgenome, result.mSbjctGenomeTo + self.mMinDistance),
                                    min( lgenome, result.mSbjctGenomeTo + self.mMinDistance + self.mMaxExtension ) )

                if r:
                    if self.mLogLevel >= 2:
                        sys.stdout.write("# downstream extension successful: %s\n" % (str(r) ))
                                             
                    downstream_results.add( r )
                else:
                    if self.mLogLevel >= 2:
                        sys.stdout.write("# downstream extension unsuccessful.\n" )
                                             
            ## now sort out the predictions by combining anchors
            ## with upstream and downstream predictions
            if upstream_results:
                upstream_results.sort( lambda x,y: cmp( (x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom), (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeTo) ) )
                new_result = upstream_results[0]

                for r in upstream_results[1:]:
                    new_result.Add( r, combine_queries = True )
                    
                new_result.Add( result, combine_queries = True )
            else:
                new_result = result

            if downstream_results:
                downstream_results.sort( lambda x,y: cmp( (x.mSbjctToken, x.mSbjctStrand, x.mSbjctGenomeFrom), (y.mSbjctToken, y.mSbjctStrand, y.mSbjctGenomeTo) ) )

                for r in downstream_results:
                    new_result.Add( r, combine_queries = True )

            new_results.add( new_result )

            start = new_result.mSbjctGenomeTo

        return new_results
    
    
##--------------------------------------------------------------------------
class RecursivePredictor:
    """run recursive prediction.

    Note: this will always run the full peptide set against the matched fragments,
    so only use this if you have only a single input sequence/hmm and the matches
    returned are not repetetive.
    """
    
    def __init__(self, predictor ):

        self.mPredictor = predictor
        
        # analyze at least 500 kb
        self.mMinGenomeSize = 500

        self.mLogLevel = 0

        self.mLevel = 0

    def setLogLevel( self, level ):
        self.mLogLevel = level
        
    def runRecursive( self,
                      filename_peptides,
                      filename_genome,
                      options ):
        
        ## do predictions
        results = self.mPredictor( filename_peptides, filename_genome, options )

        self.mLevel += 1

        infile = open(filename_genome, "r" )
        genome = Genomics.ReadGenomicSequences( infile, do_reverse = False )
        infile.close()
        
        ## collect new fragments to run
        fragments = {}

        for result in results:
            lgenome = len(genome[result.mSbjctToken])
            if result.mSbjctGenomeFrom >= self.mMinGenomeSize:
                key = "%s:%s:%i:%i" % (result.mSbjctToken, result.mSbjctStrand, 0, result.mSbjctGenomeFrom)
                fragments[key] = ( result.mSbjctToken, result.mSbjctStrand, 0, result.mSbjctGenomeFrom, lgenome)
            if lgenome - result.mSbjctGenomeTo >= self.mMinGenomeSize:
                key = "%s:%s:%i:%i" % (result.mSbjctToken, result.mSbjctStrand, result.mSbjctGenomeTo, lgenome)                
                fragments[key] = ( result.mSbjctToken, result.mSbjctStrand, result.mSbjctGenomeTo, lgenome, lgenome )

        if not fragments: return results

        new_genome_files = []
            
        first = True
        single = len(genome) == 1
        for key, data in fragments.items():

            sbjct_token, sbjct_strand, xfrom, xto, lgenome = data
            ## extract sequence to be read
            if Genomics.IsNegativeStrand( sbjct_strand ):
                xfrom, xto = lgenome - xto, lgenome - xfrom

            ## if single sequence genomes, keep them that way
            ## otherwise: combine them
            if single or first:
                if not first: os.close(new_outfile_genome)
                
                new_outfile_genome, new_filename_genome = tempfile.mkstemp()                
                new_genome_files.append( new_filename_genome )
                first = False
                                    
            os.write(new_outfile_genome, ">%s\n%s\n" % (key, genome[sbjct_token][xfrom:xto]) )

            if self.mLogLevel >= 2:
                sys.stdout.write( "# adding fragment to %s: %s for %s:%s-%s.\n" % (new_filename_genome, key, sbjct_token, xfrom, xto) )
            
        os.close(new_outfile_genome)

        ## run recursive predictions
        new_results = Prediction.Predictions()
        for new_filename_genome in new_genome_files:
            new_results.combine( self.runRecursive( filename_peptides,
                                                    new_filename_genome,
                                                    options ) )
            

        ## change coordinates of predictions so that they correspond to the original
        ## submitted contigs.
        for result in new_results:
            sbjct_token, sbjct_strand, sbjct_from, sbjct_to, lgenome = fragments[r.mSbjctToken]
            # shift aligned region for fragment
            # this is strand dependent.
            r.shiftGenomicRegion( sbjct_from, lgenome - sbjct_to )
            ## set sbjct token to original
            r.mSbjctToken = sbjct_token
            ## sort out strand
            ## if match is on negative strand
            if Genomics.IsNegativeStrand(r.mSbjctStrand):
                ## switch coordinates
                sbjct_from, sbjct_to = lgenome - sbjct_to, lgenome-sbjct_from
                if Genomics.IsNegativeStrand(sbjct_strand):
                    sbjct_strand = "+"
                else:
                    sbjct_strand = "-"
            
            r.mSbjctStrand = sbjct_strand

        ## clean up
        os.remove( new_filename_genome )

        results.combine( new_results )
        return results
    
    def __call__( self,
                  filename_peptides,
                  filename_genome,
                  options = None):
        """call runRecursive."""

        return self.runRecursive( filename_peptides,
                                  filename_genome,
                                  options )
                

def getPredictor( options, parser_options ):
    """setup predictor."""

    ## sort out parsers:
    output_options = None
    if options.pretty:
        parser_options = {}
        prediction_parser = PredictionParserEcho()
        if options.method == "genewise":
            output_options = "-pretty"
        elif options.method == "genewise-hmm":
            output_options = "-pretty"
        elif options.method == "exonerate":
            output_options = "--showalignment"
    else:
        if options.method == "genewise":
            prediction_parser = PredictionParserGenewise( **parser_options )
        elif options.method == "genewise-hmm":
            if 'filename_peptides' in parser_options:
                del parser_options['filename_peptides']
            prediction_parser = PredictionParserGenewiseHMM( **parser_options )
        elif options.method == "exonerate":
            prediction_parser = PredictionParserExonerate( **parser_options )                                    
        
    if options.method == "genewise":
        predictor = PredictorGenewise( parser = prediction_parser )
    elif options.method == "genewise-hmm":
        predictor = PredictorGenewiseHMM( parser = prediction_parser )
    elif options.method == "exonerate":
        predictor = PredictorExonerate( parser = prediction_parser )
    else:
        raise "Unknown method %s" % (options.method)

    if output_options:
        predictor.setOutputOptions( output_options )

    if options.dry_run:
        predictor.setDryRun( options.dry_run )

    predictor.setLogLevel( options.loglevel )

    if options.use_cluster:
        cluster = Cluster.ClusterSGE( options )
        predictor.setCluster( cluster )

        tempfile.tempdir = "."
        
    if options.recursive:
        predictor = RecursivePredictor( predictor )
        predictor.setLogLevel( options.loglevel )
        
    if options.downstream_queries or options.upstream_queries:
        predictor = TandemPredictor( predictor,
                                     upstream_queries = options.upstream_queries,
                                     downstream_queries = options.downstream_queries )
        predictor.setLogLevel( options.loglevel )
        
    return predictor

def runLocally( filename_peptides, filename_genome, options, parser_options ):
    """run jobs locally."""
    
    result = Prediction.Predictions()
    
    if options.genome_chunk_size:
        infile = open(filename_genome, "r")
        filename_genomes = IndexedFasta.splitFasta( infile,
                                                 options.genome_chunk_size,
                                                 dir=os.curdir )
        infile.close()
        
        if options.loglevel >= 1:
            options.stdlog.write("# starting %i jobs.\n" % len(filename_genomes) )

        x = 0
        for filename_genome in filename_genomes:
            
            x += 1
            
            if options.loglevel >= 1:
                options.stdlog.write("# running job %i on %s.\n" % (x, filename_genome ) )
            
            predictor = getPredictor( options, parser_options )

            r = predictor( filename_peptides,
                           filename_genome,
                           options = options.options )

            if options.loglevel >= 2:
                options.stdlog.write("# %s\n" % (str(r)))
            
            result.combine( r )
            
    else:
        predictor = getPredictor( options, parser_options )
        
        result = predictor( filename_peptides,
                            filename_genome,
                            options = options.options )
    
    return result

def runOnCluster( filename_peptides, filename_genome, options, parser_options ):
    """run jobs on cluster.

    Warning: does not work perfectly - there are some problems
    with the parsers which surface when called from inside a makefile,
    but (strangely) not when called from the shell directly.
    """
    all_results = Prediction.Predictions()

    def save_result( request, result ):
        if result:
            all_results.combine( result )            

    def report_error( request, exc_info ):
        options.stdlog.write( "# exception occured in request #%s: %s\n" % \
                              (request.requestID, exc_info[1]) )
        options.stdlog.flush()

    def predict( predictor, filename_peptides, filename_genome, options ):
        result = predictor( filename_peptides, filename_genome, options = options)
        return result
    
    if options.genome_chunk_size:
        infile = open(filename_genome, "r")
        filename_genomes = IndexedFasta.splitFasta( infile,
                                                    options.genome_chunk_size,
                                                    dir=os.curdir )
        infile.close()
        
        if options.loglevel >= 1:
            options.stdlog.write("# starting %i parallel jobs with a pool size of %i.\n" %\
                                 ( len(filename_genomes), options.cluster_num_jobs ) )
            options.stdlog.flush()
            
        ## run jobs in parallel

        import threadpool
        pool = threadpool.ThreadPool( options.cluster_num_jobs )

        for filename_genome in filename_genomes:
            predictor = getPredictor( options, parser_options )
            request = threadpool.WorkRequest( predict,
                                              (predictor, filename_peptides, filename_genome, options.options),
                                              (),
                                              None,
                                              save_result,
                                              report_error )
            
            pool.putRequest( request )

        pool.wait()
        
    return all_results

def predictFromGFF( options ):
    """run predictions from a gff file that associates queries with genomic regions.
    """
    
    if not options.genome_file or not options.filename_peptides or options.filename_peptides == "-":
        raise "the --from-pairs option requires both a genome and an indexed set of queries."""

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    p, e = os.path.splitext(options.filename_peptides)
    if e in (".fasta" , ".fa" ):
        peptides = Genomics.ReadPeptideSequences( open(options.filename_peptides, "r") )
    else:
        peptides = IndexedFasta.IndexedFasta( options.filename_peptides )
        
    ninput, noutput, nskipped = 0, 0, 0

    if options.from_gff == "-":
        infile = sys.stdin
    else:
        infile = open( options.from_gff, "r" )
        
    predictor = getPredictor( options, parser_options )

    options.stdout.write( "%s\n" % Prediction.Prediction().getHeader() )

    prediction_id = 0
    for gffs in GTF.joined_iterator( GTF.iterator( infile ) ): 

        ninput += 1

        gffs.sort( lambda x,y: cmp(x.start, y.start) )

        g = gffs[0]
        contig, strand = g.name, g.strand
        if "Query" in g:
            query = g["Query"]
        else:
            query = g.info

        if "Id" in g:
            id = g["Id"]
        else:
            prediction_id += 1
            id = prediction_id

        if contig not in fasta:
            if options.loglevel >= 1:
                options.stdlog.write( "# contig %s not in genome - query %s skipped\n" % (contig, query))
                options.stdlog.flush()
            nskipped += 1
            continue

        if query not in peptides:
            if options.loglevel >= 1:
                options.stdlog.write( "# query %s not in peptides - query %s skipped\n" % (query, query))
                options.stdlog.flush()
            nskipped += 1
            continue

        start, end = gffs[0].start, gffs[-1].end

        lcontig = fasta.getLength( contig )

        if end > lcontig:
            raise "region for query %s out of bounds: %s:%i > %i" % (query, contig, end, lcontig)

        if options.add_flank:
            start = max( 0, start - options.add_flank )
            end = min( lcontig, end + options.add_flank )

        if Genomics.IsNegativeStrand( strand):
            start, end = lcontig - end, lcontig - start
            
        if options.loglevel >= 2:
            options.stdlog.write( "# predicting query %s in region %s:%s:%i:%i\n" %\
                                      (query, contig, strand, start, end))
            options.stdlog.flush()

        result = predictor.predictFromRange( 
            peptides[query],
            contig,
            strand,
            start, 
            end,
            fasta,
            options = options.options )
        
        noutput += 1

        for r in result:
            r.mPredictionId = id
            r.mQueryToken = query
            options.stdout.write( "%s\n" % str(r))
            
        options.stdout.flush()

    if options.from_gff != "-": infile.close()

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

def predictFromPairs( options ):
    """run predictions from two fasta files, one with peptide sequences and the
    other with genomic sequences. 

    Sequences are associated via the identifier.

    The cds file should be indexd.
    """
    
    if not options.genome_file or not options.filename_peptides or options.filename_peptides == "-":
        raise "the --from-pairs option requires both a genome and an indexed set of queries."""


    try:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    except KeyError:
        fasta = Genomics.ReadGenomicSequences( open(options.genome_file, "r"), do_reverse = False )

    try:
        peptides = IndexedFasta.IndexedFasta( options.filename_peptides )
        keys = peptides.getContigs()
    except KeyError:
        peptides = Genomics.ReadPeptideSequences( open(options.filename_peptides, "r") )
        keys = peptides.keys()

    ninput, noutput, nskipped = 0, 0, 0

    predictor = getPredictor( options, parser_options )

    options.stdout.write( "%s\n" % Prediction.Prediction().getHeader() )

    for query in keys:

        ninput += 1
        peptide_sequence = peptides[query]

        if query not in fasta:
            E.warn( "query %s not cds database: skipped" % (query,))
            nskipped += 1
            continue

        genomic_sequence = fasta[query]

        E.debug( "predicting query %s against sequence of %i" % (query, len(genomic_sequence)))

        result = predictor.predictFromSequences( 
            peptide_sequence,
            genomic_sequence, 
            options = options.options )
        
        noutput += 1

        for r in result:
            r.mPredictionId = query
            r.mQueryToken = query
            options.stdout.write( "%s\n" % str(r))
            
        options.stdout.flush()

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: Predictor.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-p", "--filename-peptides", dest="filename_peptides", type ="string",
                       help="filename with peptide sequences. Use - for stdin.")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("genewise", "exonerate", "genewise-hmm"),
                      help="""gene prediction method to use. If genewise-hmm is chosen, the peptide sequence file
                      is assumed to contain a profile."""  )
    
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome pattern."  )
    
    parser.add_option("-r", "--range", dest="range", type="string",
                      help="range to scan."  )

    parser.add_option( "--from-gff", dest="from_gff", type="string",
                       help="""read a list of query identifiers and genomic locations from a tab-separated file (- for stdin) and run a 
predictor for each. 
""" )

    parser.add_option( "--from-pairs", dest="from_pairs", action="store_true",
                       help="""predict pairs of peptide/cds sequences from the two fasta files. The identifiers must match""")

    parser.add_option("-R", "--filename-ranges", dest="filename_ranges", type="string",
                      help="filename with genomic ranges to scan a set of queries against."  )

    parser.add_option("-o", "--options", dest="options", type="string",
                       help="options to pass on to predictor.")

    parser.add_option("-i", "--pretty", dest="pretty", action="store_true",
                      help="pretty print results directly.")
    
    parser.add_option( "--dry-run", dest="dry_run", action="store_true",
                      help="prepare files, but do not execute prediction program.")

    parser.add_option( "--forward-coordinates", dest="forward_coordinates", action="store_true",
                      help="return result in forward coordinates.")

    parser.add_option( "--both-strands", dest="both_strands", action="store_true",
                      help="analyze both strands.")

    parser.add_option( "--use-fasta-coordinates", dest="use_fasta_coordinates", action="store_true",
                       help="""use fasta coordinates. The format of the description line should thus be
                       >id contig:strand:from:to:lcontig in 0-based reverse coordinates.
                       """ )

    parser.add_option( "--recursive", dest="recursive", action="store_true",
                      help="analyze genomic stretches recursively.")

    parser.add_option( "--downstream-queries", dest="downstream_queries", type="string", action="append",
                      help="downstream queries, triggered for each match.")
    
    parser.add_option( "--upstream-queries", dest="upstream_queries", type="string", action="append",
                      help="upstream queries, triggered for each match.")

    parser.add_option( "--genome-chunk-size", dest="genome_chunk_size", type="int",
                      help="chunk size of genomic files.")
    
    parser.add_option( "--input-coordinate-format", dest="input_coordinate_format", type="choice",
                       choices=( "zero-both", "zero-forward" ),
                       help="coordinate format." )

    parser.add_option( "--add-flank", dest="add_flank", type="int",
                      help="add a flanking region to each region to be queried.")
    
    parser.set_defaults(
        method = "exonerate",
        options   = "",
        genome_file = "genome",
        filename_peptides = "-",
        filename_ranges = None,
        filename_profile = None,
        range = None,
        pretty = False,
        dry_run = False,
        forward_coordinates = False,
        both_strands = False,
        use_fasta_coordinates = False,
        recursive = False,
        genome_chunk_size = 0,
        downstream_queries = [],
        upstream_queries = [] ,
        input_coordinate_format= "zero-both",
        add_flank = 0,
        from_pairs = False,
        from_gff = None,
        )

    (options, args) = E.Start( parser,
                               add_pipe_options = True,
                               add_cluster_options = True )

    parser_options = {}

    ##########################################################################
    ##########################################################################
    ##########################################################################    
    ## setup input files
    ##########################################################################
    ## translation table between internal contig names
    ## and original locations:
    map_key2contig = {}    

    if options.from_gff:

        predictFromGFF( options )
        E.Stop()
        sys.exit(0)

    elif options.from_pairs:
        
        predictFromPairs( options )
        E.Stop()
        sys.exit(0)
        
    if options.filename_ranges:

        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        
        if options.filename_ranges == "-":
            infile = sys.stdin
        else:
            infile = open(options.filename_ranges)
            
        genomic_sequences = {}
        
        outfile, filename_genome = tempfile.mkstemp()
       
        for line in infile:
            if line[0] == "#": continue

            sbjct_token, sbjct_strand, sbjct_from, sbjct_to = line[:-1].split(":")

            sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)

            key = "%s:%s:%s:%s" % (sbjct_token, sbjct_strand, sbjct_from, sbjct_to)
            
            genomic_sequence = fasta.getSequence( sbjct_token, sbjct_strand,
                                                  sbjct_from, sbjct_to )
            
            genomic_sequences[key] = genomic_sequence

            map_key2contig[key] = (sbjct_token, sbjct_strand, sbjct_from, sbjct_to )
        
            os.write(outfile, ">%s\n%s\n" % (key, genomic_sequence))

        os.close(outfile)

        if options.filename_ranges == "-":
            infile.close()

        parser_options['genomes'] = genomic_sequences

    elif options.range:
        
        fasta = IndexedFasta.IndexedFasta( options.genome_file )

        sbjct_token, sbjct_strand, sbjct_from, sbjct_to = Genomics.String2Location( options.range )

        ## convert to my coordinate format
        converter = IndexedFasta.getConverter( options.input_coordinate_format )        
        sbjct_from, sbjct_to = converter( sbjct_from, sbjct_to,
                                          Genomics.IsPositiveStrand(sbjct_strand),
                                          fasta.getLength( sbjct_token ) )

        genomic_sequence = fasta.getSequence( sbjct_token, sbjct_strand,
                                              sbjct_from, sbjct_to )

        key = "%s:%s:%s:%s" % (sbjct_token, sbjct_strand, sbjct_from, sbjct_to)
        
        map_key2contig[key] = (sbjct_token, sbjct_strand, sbjct_from, sbjct_to )
        
        genomic_sequences = { key : genomic_sequence }
        
        outfile, filename_genome = tempfile.mkstemp()
        os.write(outfile, ">%s\n%s\n" % (key, genomic_sequence))
        os.close(outfile)
        parser_options['genomes'] = genomic_sequences
        
    else:
        filename_genome = options.genome_file
        parser_options['filename_genome'] = filename_genome
        
    peptide_sequences = {}
    if options.filename_peptides == "-":
        ## save peptides to temporary file
        outfile, filename_peptides = tempfile.mkstemp()
        for line in sys.stdin:
            os.write(outfile, line)
        os.close(outfile)
        parser_options['filename_peptides'] = filename_peptides
    else:
        filename_peptides = options.filename_peptides
        parser_options['filename_peptides'] = options.filename_peptides

    ##########################################################################
    ##########################################################################
    ##########################################################################    
    ## setup cluster run by splitting peptide and or genomic files
    ##########################################################################
    if options.use_cluster:
        result = runOnCluster( filename_peptides, filename_genome, options, parser_options )
    else:
        result = runLocally( filename_peptides, filename_genome, options, parser_options )        

    E.info( "obtained %i results" % len(result))

    ##########################################################################
    ##########################################################################
    ##########################################################################    
    ## post-processing
    ##########################################################################
    if options.range:
        if not options.dry_run:
            os.remove(filename_genome)
        for r in result:
            sbjct_token, sbjct_strand, sbjct_from, sbjct_to = map_key2contig[r.mSbjctToken]
            lgenome = fasta.getLength( sbjct_token )
            # shift aligned region for fragment
            # this is strand dependent.
            r.shiftGenomicRegion( sbjct_from, lgenome - sbjct_to )
            ## set sbjct token to original
            r.mSbjctToken = sbjct_token
            ## sort out strand
            ## if match is on negative strand
            if r.mSbjctStrand in ("-", "0", -1, "-1"):
                ## switch coordinates
                sbjct_from, sbjct_to = lgenome - sbjct_to, lgenome-sbjct_from
                if sbjct_strand in ("-", "0", -1, "-1"):
                    sbjct_strand = "+"
                else:
                    sbjct_strand = "-"
            
            r.mSbjctStrand = sbjct_strand

    elif options.use_fasta_coordinates:
        ## get coordinates from fasta file.
        ## the contig length should be supplied in the file
        ## for this to work.
        coordinates = {}
        infile = open(filename_genome,"r")
        for line in infile:
            if line[0] == ">":
                tokens = re.split("\s+", line[1:-1] )
                if len(tokens) == 1:
                    coord = tokens[0]
                    id = coord
                else:
                    id, coord = tokens[:2]

                try:
                    contig, strand, xfrom, xto, lcontig = coord.split(":")
                except ValueError:
                    E.warn( "skipped: can not parse coordinates for %s = %s" % (id, coord))
                    continue

                coordinates[id] = ( contig, strand, int(xfrom), int(xto), int(lcontig) )
                
        infile.close()
        
        for r in result:
            sbjct_token, sbjct_strand, sbjct_from, sbjct_to, lcontig = coordinates[r.mSbjctToken]
            
            # shift aligned region for fragment
            # this is strand dependent.
            r.shiftGenomicRegion( sbjct_from, lcontig - sbjct_to )
            ## set sbjct token to original
            r.mSbjctToken = sbjct_token
            
            ## sort out strand
            ## if match is on negative strand
            if r.mSbjctStrand in ("-", "0", -1, "-1"):
                ## switch coordinates
                if sbjct_strand in ("-", "0", -1, "-1"):
                    sbjct_strand = "+"
                else:
                    sbjct_strand = "-"
            
            r.mSbjctStrand = sbjct_strand

            ## convert to forward coordinates, if so desired.
            if options.forward_coordinates:
                r.mSbjctGenomeFrom, r.mSbjctGenomeTo = lcontig - r.mSbjctGenomeTo, lcontig - r.mSbjctGenomeFrom

        ## do not convert again (lcontig from fasta object does not reflect
        ## true contig size).
                
        options.forward_coordinates = False

    if options.forward_coordinates:
        for r in result:
            lgenome = fasta.getLength( sbjct_token )
            if r.mSbjctStrand in ("-", "0", -1, "-1"):
                r.mSbjctGenomeFrom, r.mSbjctGenomeTo = lgenome - r.mSbjctGenomeTo, lgenome - r.mSbjctGenomeFrom

    options.stdout.write( "%s\n" % Prediction.Prediction().getHeader() )
    x = 0
    for r in result:
        x += 1
        r.mPredictionId = x
        options.stdout.write( "%s\n" % str(r))

    ## clean up
    if options.filename_peptides == "-" and not options.dry_run:
        os.remove( filename_peptides)

    E.Stop()

