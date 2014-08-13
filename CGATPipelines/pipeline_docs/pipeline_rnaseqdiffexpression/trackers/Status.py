import itertools
import numpy

from RnaseqDiffExpressionReport import *


class ExpressionStatus(Status):

    '''status information - estimating transcription levels
    for protein coding transcripts.
    '''

    pattern = "(.*)_gene_diff"

    # minimum expression level for transcripts to be
    # considered expressed
    min_fpkm = 1.0

    def _testErrorBars(self, track, part):

        statement = '''SELECT (%(part)s_std) / %(part)s_mean / 2
        FROM %(track)s_gene_diff
        WHERE %(part)s_mean > %(min_fpkm)f'''

        values = self.getValues(statement)
        value = numpy.median(values)

        if value <= 0.5:
            status = "PASS"
        elif value <= 1.0:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testErrorBarsTreatment(self, track):
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
        return self._testErrorBars(track, "treatment")

    def testErrorBarsControl(self, track):
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
        return self._testErrorBars(track, "control")


class ExpressionStatusNoncoding(ExpressionStatus):
    tablename = "refnoncoding_cuffdiff_gene_levels"


class DifferentialExpressionStatus(Status):
    pattern = "(.*)_gene_diff"

    def testTests(self, track):
        '''test if tests for differential expression are successful.

        Unsuccessful test have status NO_CALL or FAIL, meaning that
        the test failed or the expression levels of at least one gene
        was very low. If that is the case, some statistical tests
        become very unreliable.

        PASS: <= 10% of tests failed.
        WARN: <= 50% of tests failed
        FAIL: > 50% of tests failed

        '''

        failed = self.getValue(
            '''SELECT COUNT(*) FROM %(track)s_gene_diff WHERE status = 'FAIL' OR status = 'NO_TEST' ''')
        tested = self.getValue( '''SELECT COUNT(*) FROM %(track)s_gene_diff''')

        value = float(failed) / tested

        if value <= 0.1:
            status = "PASS"
        elif value <= 0.5:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)
