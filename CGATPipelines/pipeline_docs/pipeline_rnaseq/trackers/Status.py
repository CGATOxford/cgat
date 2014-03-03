import os
import sys
import re
import types
import itertools
import math
import numpy

from RnaseqReport import *


class MappingStatus(Status):

    '''status information for mapping stage.'''

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM view_mapping")
        return tuple([x[0] for x in d])

    slices = ("Mapping", "RepetetiveRNA", "SplicedAlignments")

    def testMapping(self, track):
        '''proportion of reads mapped.

        PASS : >=60% reads mapped
        WARN : >=40% reads mapped
        FAIL : < 40% reads mapped

        '''

        value = self.getValue(
            "SELECT reads_mapped/CAST( reads_in AS FLOAT) from view_mapping WHERE track = '%(track)s'")
        if value is None:
            value, status = 0, "NA"
        elif value >= 0.6:
            status = "PASS"
        elif value >= 0.4:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testRepetetiveRNA(self, track):
        '''proportion of reads mapping to repetetive RNA.

        PASS : < 10% reads mapping to repetetive RNA
        WARN : < 50% reads mapping to repetetive RNA
        FAIL : >=50% reads mapping to repetetive RNA

        '''

        value = self.getValue(
            "SELECT 1.0 - reads_norna/CAST(reads_mapped AS FLOAT) from view_mapping WHERE track = '%(track)s'")
        if value < 0.05:
            status = "PASS"
        elif value < 0.1:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testSplicedAlignments(self, track):
        '''proportion of spliced alignments.

        PASS: >= 15% of alignments spliced
        WARN: >=  5% of alignments spliced
        FAIL: <   5% of alignments spliced

        '''

        # track = re.sub("-", "_", track)
        value = self.getValue(
            "SELECT spliced/CAST(input AS FLOAT) from exon_validation WHERE track = '%(track)s.accepted'")
        if value >= 0.15:
            status = "PASS"
        elif value >= 0.05:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)


class TranscriptStatus(Status):

    '''status information on transcriptome building.'''

    pattern = "(.*)_transcript_counts"

    # minimum percentage of bases covered by reads for a gene to
    # to be considered as to be fully covered
    min_coverage = 90

    # minimum number of reads for genes to be considered expressed
    min_reads = 5

    # reference gene set to use
    reference = "refcoding"

    def testGeneCoverage(self, track):
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

        value = float(
            len([x for x in values if x >= self.min_coverage])) / len(values)

        if value >= 0.2:
            status = "PASS"
        elif value >= 0.1:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testFragmentation(self, track):
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

        complete = self.getValue(
            """SELECT COUNT(*) FROM %(track)s_class WHERE class = 'complete'""" )
        fragments = self.getValue(
            """SELECT COUNT(*) FROM %(track)s_class WHERE class like '%%fragment%%' """ )

        value = float(complete) / (complete + fragments)

        if value >= 0.25:
            status = "PASS"
        elif value >= 0.1:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testDirectionality(self, track):
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

        if value >= 0.4:
            status = "NA"
        elif value <= 0.05:
            status = "PASS"
        elif value <= 0.10:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)


class ExpressionStatus(Status):

    '''status information - estimating transcription levels
    for protein coding transcripts.
    '''

    pattern = "(.*)_transcript_counts"

    tablename = "refcoding_cuffdiff_gene_levels"

    # minimum expression level for transcripts to be
    # considered expressed
    min_fpkm = 1.0

    tracks = [x.asTable() for x in EXPERIMENTS]

    def testErrorBars(self, track):
        '''test error bars of expression level measurements.

        PASS: median relative error <= 50%
        WARN: median relative error <= 100%
        FAIL: median relative error > 100%

        The relative error is the size of the confidence interval
        divided by the absolute expression level divided by two.

        This test fails if more than half the isoforms can not be
        estimated with an accuracy of at least half the expression 
        value.
        '''

        statement = '''SELECT (%(track)s_conf_hi - %(track)s_conf_lo ) / %(track)s_fpkm / 2
                       FROM %(tablename)s WHERE %(track)s_fpkm > %(min_fpkm)f'''

        values = self.getValues(statement)
        value = numpy.median(values)

        if value <= 0.5:
            status = "PASS"
        elif value <= 1.0:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)


class ExpressionStatusNoncoding(ExpressionStatus):
    tablename = "refnoncoding_cuffdiff_gene_levels"


class DifferentialExpressionStatus(Status):
    pattern = "(.*)_cuffdiff_gene_diff"

    def testTests(self, track):
        '''test if tests for differential expression are successful.

        Unsuccessful test have status NO_CALL or FAIL, meaning that the
        test failed or the expression levels of at least one gene was very
        low. If that is the case, some statistical tests become very unreliable.

        PASS: <= 10% of tests failed.
        WARN: <= 50% of tests failed
        FAIL: > 50% of tests failed

        '''

        failed = self.getValue(
            '''SELECT COUNT(*) FROM %(track)s_cuffdiff_gene_diff WHERE status = 'FAIL' OR status = 'NO_TEST' ''')
        tested = self.getValue(
            '''SELECT COUNT(*) FROM %(track)s_cuffdiff_gene_diff''')

        value = float(failed) / tested

        if value <= 0.1:
            status = "PASS"
        elif value <= 0.5:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)
