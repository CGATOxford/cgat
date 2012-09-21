import os, sys, re, types, itertools, math, numpy

from PeakcallingReport import *

class PeakCallingStatus( Status ):
    '''status information for mapping stage.'''

    def testCalling( self, track ):
        '''number of peaks called. The number of peaks expected in a sample
        will of course vary wildly.
        
        PASS : >1000 peaks called
        WARN : >100 peaks called
        FAIL : <100 peaks called

        '''

        suffix = re.sub( "^[^_]*_", "", self.pattern )
        value = self.getValue( """SELECT COUNT(*) FROM %(track)s_%(suffix)s""" )
        if value >= 1000: status= "PASS"
        elif value >= 100: status= "WARNING"
        else: status= "FAIL"

        return status, "%i" % value

class PeakCallingStatusMACS( PeakCallingStatus ):
    pattern = ("(.*)_macs_regions" )

class PeakCallingStatusSPP( PeakCallingStatus ):
    pattern = ("(.*)_spp_regions" )

class PeakCallingStatusSICER( PeakCallingStatus ):
    pattern = ("(.*)_sicer_regions" )

class PeakCallingStatusZinba( PeakCallingStatus ):
    pattern = ("(.*)_zinba_regions" )

class PeakCallingStatusPeakRanger( PeakCallingStatus ):
    pattern = ("(.*)_peakranger_regions" )

