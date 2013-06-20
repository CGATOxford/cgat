import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from ChipseqReport import *

class MacsSummary( DefaultTracker ):
    pattern = "(macs_summary)"

    def getTracks( self, subset = None ):
        return self.getValues( "SELECT track FROM macs_summary ORDER BY track" )
    
    def __call__(self, track, slice = None ):

        resultsdir = os.path.abspath( os.path.join( EXPORTDIR, "MACS" ) )
        
        fields = (
            "called_positive", "called_negative",
            "scan_window", "shift",
            "tag_treatment_total", "tag_treatment_filtered", 
            "tag_control_total", "tag_control_filtered", 
            "ncandidates_positive", "ncandidates_negative", 
            "min_tags",
            "paired_peaks", )

        f = ",".join(fields)
        data = self.getFirstRow( '''SELECT %(f)s FROM macs_summary WHERE track="%(track)s"''' % locals())
        result = odict( zip( fields,data) )

        if os.path.exists( resultsdir ):
            result["peakshape"] = "`pdf <%(resultsdir)s/%(track)s_model.pdf>`_" % locals()

        return result

class MacsDiagnostics(ChipseqTracker):
    """Closest distance of transcript models to gene models in the reference set."""

    pattern = "(.*)_macsdiag"

    def __call__(self, track, slice = None ):

        data = self.get( "SELECT fc,npeaks,p20,p30,p40,p50,p60,p70,p80,p90 FROM %(track)s_macsdiag" % locals() )

        result = odict()
        for fc,npeaks,p20,p30,p40,p50,p60,p70,p80,p90 in data:
            result[fc] = odict()
            result[fc]["npeaks"] = npeaks
            result[fc]["proportion of reads"] = range( 20, 100, 10)
            result[fc]["proportion of peaks"] = map( float, (p20,p30,p40,p50,p60,p70,p80,p90) ) 
            
        return result

class MacsFiltering(ChipseqTracker, SingleTableTrackerColumns ):
    column = "fdr"
    table = "macs_fdr"
