import os, sys, re, types
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats

from SphinxReport.Tracker import *

from RnaseqReport import *
import math

##################################################################################
##################################################################################
##################################################################################
## Trackers that access reference statistics
##################################################################################

class ReferenceData(TrackerSQL):
    """Base class f or Trackers accessing reference table."""
    pattern = "(.*)_transcript_counts$" 
    reference = "refcoding"

class TranscriptCoverage(ReferenceData):
    """Coverage of reference transcripts."""
    mXLabel = "overlap / %"
    def __call__(self, track, slice = None ):
        data = self.getValues( """SELECT coverage_pcovered FROM %(track)s_transcript_counts WHERE coverage_nval > 0""" )
        return odict( (("covered", data ) ,) )

class GeneCoverage(ReferenceData):
    '''Coverage of reference genes - max transcript coverage per gene.'''
    mXLabel = "number of transcripts"
    def __call__(self, track, slice = None ):
        data = self.getValues( """SELECT max(c.coverage_pcovered) FROM 
                                            %(track)s_transcript_counts as c,
                                            %(reference)s_transcript2gene as i
                                         WHERE c.coverage_nval > 0
                                   AND i.transcript_id = c.transcript_id 
                                   GROUP BY i.gene_id """ )
        return odict( (("covered", data ) ,) )

class CoverageVsLengthByReadDepth(ReferenceData):
    """plot the absolute coverage of a known gene versus its length.
    Dots are colored by read depth.
    """

    mXLabel = "log(length)"

    def __call__(self, track, slice = None):
        reference = self.reference
        statement = """SELECT AVG(exons_sum) AS ref_length,
                                MIN(c.coverage_pcovered) AS coverage, 
                                AVG(c.coverage_mean) AS read_depth
                        FROM %(track)s_transcript_counts AS c,
                                %(reference)s_transcript2gene as i
                        WHERE i.transcript_id = c.transcript_id AND 
                                c.coverage_nval > 0
                        GROUP BY i.gene_id"""

        data = [ (math.log(x[0]), x[1], math.log(x[2]) ) for x in self.get( statement % locals() ) ]
        r = odict( zip( ("log(length)", "log(coverage)", "log(read_depth)"), zip(*data)))
        return odict( zip( ("log(length)", "log(coverage)", "log(read_depth)"), zip(*data)))

##=================================================================
## Coverage
##=================================================================
class MeanVsMaxReadDepth( ReferenceData ):
    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "maximum read depth"

    def __call__(self, track, slice = None ):
        reference = self.reference
        statement = "SELECT coverage_mean, coverage_max, exons_sum FROM %(track)s_transcript_counts" % locals()
        data = [ (x[0], x[1], math.log( x[2]) ) for x in self.get( statement) if x[2] > 0 ]
        return odict( zip( ("mean coverage", "max coverage", "length" ), zip(*data) ) )

class MeanVsMedianReadDepth( ReferenceData ):
    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "median read depth"

    def __call__(self, track, slice = None ):
        reference = self.reference
        statement = "SELECT coverage_mean, coverage_median, exons_sum FROM %(track)s_transcript_counts" % locals()
        data = [ (x[0], x[1], math.log( x[2]) ) for x in self.get( statement) if x[2] > 0 ]
        return odict( zip( ("mean coverage", "median coverage", "length" ), zip(*data) ) )



