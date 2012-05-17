import os, sys, re, types, itertools, math, numpy

from RnaseqTranscriptsReport import *

class TranscriptStatus( Status ):
    '''status information on transcriptome building.'''
    
    pattern = "(.*)_transcript_counts"

    # minimum percentage of bases covered by reads for a gene to
    # to be considered as to be fully covered
    min_coverage = 90

    # minimum number of reads for genes to be considered expressed
    min_reads = 5

    # reference gene set to use
    reference = "refcoding"

    def testGeneCoverage( self, track ):
        '''test coverage of known protein coding transcript models with reads.

        PASS: >= 20% of expressed genes >90% covered
        WARN: <= 20% of expressed genes >90% covered
        FAIL: <= 10% of expressed genes >90% covered

        Only genes with at least 5 reads mapping to them are used in order to only take into
        account the number of genes that are expressen.

        Read directionality is ignored in this test.
        '''
        
        values = self.getValues( """SELECT max(c.coverage_anysense_pcovered) FROM 
                                            %(track)s_transcript_counts as c,
                                            %(reference)s_transcript2gene as i
                                         WHERE c.coverage_anysense_nval > %(min_reads)i
                                   AND i.transcript_id = c.transcript_id 
                                   GROUP BY i.gene_id""" )
        
        value = float(len( [x for x in values if x >= self.min_coverage] )) / len(values)
        
        if value >= 0.2: status= "PASS"
        elif value >= 0.1: status= "WARNING"
        else: status= "FAIL"

        return status, "%5.2f%%" % (100.0 * value)
        
    def testFragmentation( self, track ):
        '''test transcript construction.

        A large number of fragments are an indicator that transcript
        building failed.

        PASS: >= 25% of transcripts are complete
        WARN: >= 10% of transcripts are complete
        FAIL: < 10% of transicrpts are complete

        This test uses transcripts classification. The test statistic
        is the ratio of ``complete`` transcripts and ``complete+fragment``
        transcripts.
        '''

        complete = self.getValue( """SELECT COUNT(*) FROM %(track)s_class WHERE class = 'complete'""" )
        fragments = self.getValue( """SELECT COUNT(*) FROM %(track)s_class WHERE class like '%%fragment%%' """ )

        value = float(complete) / (complete + fragments )

        if value >= 0.25: status= "PASS"
        elif value >= 0.1: status= "WARNING"
        else: status= "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testDirectionality( self, track ):
        '''test proportion of antisense reads within known transcript models.

        PASS: <= 1% of reads are antisense
        WARN: <= 5% of reads are antisense
        FAIL: >5% of reads are antisense

        If the proportion is close to 50% (>40%), this test is ignored.

        Note that some level of antisense expression is expected.
        '''

        values = self.get( """SELECT SUM(c.coverage_sense_pcovered), 
                                     SUM(c.coverage_antisense_pcovered)
                              FROM
                                  %(track)s_transcript_counts as c
                              WHERE c.coverage_anysense_nval > %(min_reads)i""" )
        
        sense, antisense = values[0]
        value = float(antisense) / (sense + antisense)

        if value >= 0.4: status = "NA"
        elif value <= 0.05: status= "PASS"
        elif value <= 0.10: status= "WARNING"
        else: status= "FAIL"

        return status, "%5.2f%%" % (100.0 * value)





