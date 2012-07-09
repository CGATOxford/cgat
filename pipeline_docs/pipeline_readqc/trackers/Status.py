import os, sys, re, types, itertools, math, numpy

from ReadqcReport import *

class FastqcStatus( Status ):
    '''status information for mapping stage.'''
    
    pattern = "(.*)_fastqc_status"

    def testBasicStatistics( self, track ):
        '''basic statistics.
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Basic Statistics'" ), 0

    def testPerBaseSequenceQuality( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per base sequence quality'" ), 0

    def testPerSequenceQualityScores( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per sequence quality scores'" ), 0
    
    def testPerBaseSequenceContent( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per base sequence content'" ), 0

    def testPerBaseGCContent( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per base GC content'" ), 0

    def testPerSequenceGCContent( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per sequence GC content'" ), 0

    def testPerBaseNContent( self, track ):
        '''per base sequence quality
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Per base N content'" ), 0

    def testSequenceLengthDistribution( self, track ):
        '''Sequence length distribution
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Sequence Length Distribution'" ), 0

    def testSequenceDuplicationLevels( self, track ):
        '''Sequence duplication levels
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Sequence Duplication Levels'" ), 0

    def testOverrepresentedSequences( self, track ):
        '''Overrepresented sequences
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Overrepresented sequences'" ), 0

    def testKmerContent( self, track ):
        '''Kmer content
        '''
        return self.getValue( "SELECT status FROM %(track)s_fastqc_status WHERE name = 'Kmer Content'" ), 0


        


