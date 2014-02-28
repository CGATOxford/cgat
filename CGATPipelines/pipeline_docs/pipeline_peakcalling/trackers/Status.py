import os
import sys
import re
import types
import itertools
import math
import numpy

from PeakcallingReport import *


class PeakCallingStatus(Status):

    '''status information for mapping stage.'''

    def testCalling(self, track):
        '''number of peaks called. The number of peaks expected in a sample
        will of course vary wildly.

        PASS : >1000 peaks called
        WARN : >100 peaks called
        FAIL : <100 peaks called

        '''

        suffix = re.sub("^[^_]*_", "", self.pattern[:-1])

        value = self.getValue(
            """SELECT COUNT(*) FROM %(track)s_%(suffix)s""" )
        if value >= 1000:
            status = "PASS"
        elif value >= 100:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%i" % value


class PeakCallingStatusMACS(PeakCallingStatus):
    pattern = ("(.*)_macs_regions$")


class PeakCallingStatusMACS2(PeakCallingStatus):
    pattern = ("(.*)_macs2_regions$")


class PeakCallingStatusSPP(PeakCallingStatus):
    pattern = ("(.*)_spp_regions$")


class PeakCallingStatusSICER(PeakCallingStatus):
    pattern = ("(.*)_sicer_regions$")


class PeakCallingStatusZinba(PeakCallingStatus):
    pattern = ("(.*)_zinba_regions$")


class PeakCallingStatusPeakRanger(PeakCallingStatus):
    pattern = ("(.*)_peakranger_peaks$")


class PeakCallingStatusCCAT(PeakCallingStatus):
    pattern = ("(.*)_ccat_peaks$")


class EncodeQualityMetrics(Status):

    '''
    See http://code.google.com/p/phantompeakqualtools/ and :pmid:`22955991`
    '''

    tablename = "spp_quality"

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM %(tablename)s")

    def testNSC(self, track):
        '''Normalized strand correlation coefficient (NSC)

        PASS : >=1.1
        WARN : >=1.05
        FAIL : < 1.05

        '''

        value = self.getValue(
            """SELECT nsc FROM %(tablename)s WHERE track = '%(track)s'""" )
        if value >= 1.1:
            status = "PASS"
        elif value >= 1.05:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%f" % value

    def testRSC(self, track):
        '''Relative strand cross-correlation coefficient (RSC)

        PASS : >=1.0
        WARN : >= 0.8
        FAIL : <0.8
        '''

        value = self.getValue(
            """SELECT rsc FROM %(tablename)s WHERE track = '%(track)s'""" )
        if value >= 1.0:
            status = "PASS"
        elif value >= 0.8:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%f" % value

    def testMappedReads(self, track):
        '''ENCODE recommends 10Mio uniquely mapped reads for mammalian genomes
        for point-binding intervals and at least 20Mio mapped reads for broad peaks.

        PASS : >= 20000000
        WARN : >= 10000000
        FAIL : <10000000
        '''
        value = self.getValue(
            """SELECT mapped_reads FROM %(tablename)s WHERE track = '%(track)s'""" )
        if value >= 2e7:
            status = "PASS"
        elif value >= 1e7:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%i" % value
