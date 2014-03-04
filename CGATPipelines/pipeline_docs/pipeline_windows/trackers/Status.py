import os
import sys
import re
import types
import itertools
import math
import numpy

from MedipReport import *


class MappingStatus(Status):

    '''status information for mapping stage.'''

    def getTracks(self):
        d = self.get(
            "SELECT DISTINCT track FROM bam_stats WHERE track LIKE '%%.genome'")
        return tuple([re.sub(".genome", "", x[0]) for x in d])

    def testMapping(self, track):
        '''proportion of reads mapped.

        PASS : >=60% reads mapped
        WARN : >=40% reads mapped
        FAIL : < 40% reads mapped

        '''
        value = self.getValue( """SELECT reads_mapped/CAST( reads_total AS FLOAT) 
                                         FROM bam_stats 
                                         WHERE track = '%(track)s.genome'""" )
        if value >= 0.6:
            status = "PASS"
        elif value >= 0.4:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)

    def testPairing(self, track):
        '''proportion of reads mapped.

        PASS : >=95% reads in proper pairs
        WARN : >=80% reads in proper pairs
        FAIL : < 80% reads in proper pairs

        '''

        value = self.getValue( """SELECT proper_pair/CAST( mapped AS FLOAT) 
                                         FROM bam_stats 
                                         WHERE track = '%(track)s.genome'""" )
        if value >= 0.95:
            status = "PASS"
        elif value >= 0.8:
            status = "WARNING"
        else:
            status = "FAIL"

        return status, "%5.2f%%" % (100.0 * value)
