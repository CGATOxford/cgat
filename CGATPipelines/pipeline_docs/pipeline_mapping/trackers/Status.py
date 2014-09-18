import os
import sys
import re
import types
import itertools
import math
import numpy

from MappingReport import *


class MappingStatus(Status):

    '''status information for mapping stage.'''

    @property
    def tracks(self):
        d = self.get("SELECT DISTINCT track FROM view_mapping")
        return tuple([x[0] for x in d])

    slices = ("Mapping", "PairMapping", "RepetetiveRNA", "SplicedAlignments")

    def testMapping(self, track):
        '''proportion of reads mapped.

        PASS : >=80% reads mapped
        WARN : >=40% reads mapped
        FAIL : < 40% reads mapped

        '''
        value = self.getValue(
            """SELECT reads_mapped/CAST(reads_total AS FLOAT)
            from view_mapping WHERE track = '%(track)s'""")

        if value >= 0.8:
            status = "PASS"
        elif value >= 0.4:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testPairMapping(self, track):
        '''proportion of pairs mapped uniquely and in proper pairs.

        PASS : >=80% pairs mapped
        WARN : >=40% pairs mapped
        FAIL : < 40% pairs mapped
        '''

        value = self.getValue("""SELECT
        pairs_proper_unique_alignments / CAST(pairs_total AS FLOAT)
        FROM view_mapping
        WHERE track = '%(track)s'
        AND pairs_total > 0""")

        if value is None:
            return "NA", 0

        if value >= 0.8:
            status = "PASS"
        elif value >= 0.4:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testRepetetiveRNA(self, track):
        '''proportion of reads mapping to repetetive RNA.

        PASS : < 10% alignments mapping to repetetive RNA
        WARN : < 50% alignments mapping to repetetive RNA
        FAIL : >=50% alignments mapping to repetetive RNA

        '''

        value = self.getValue(
            "SELECT 1.0 - alignments_no_rna/CAST(alignments_mapped AS FLOAT) from view_mapping WHERE track = '%(track)s'")
        if value < 0.05:
            status = "PASS"
        elif value < 0.1:
            status = "WARNING"
        else:
            status = "FAIL"

        return "NA", "NOT TESTED"
        return status, "%5.2f%%" % (100.0 * value)

    def testSplicedAlignments(self, track):
        '''proportion of spliced alignments.

        PASS: >= 15% of alignments spliced
        WARN: >=  5% of alignments spliced
        FAIL: <   5% of alignments spliced

        '''

        if not self.hasTable('exon_validation'):
            return "NA", 0

        value = self.getValue(
            "SELECT spliced/CAST(input AS FLOAT) from exon_validation WHERE track = '%(track)s'")
        if value >= 0.15:
            status = "PASS"
        elif value >= 0.05:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)
