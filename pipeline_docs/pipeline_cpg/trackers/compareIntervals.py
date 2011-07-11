import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class UniqueIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_unique_intervals" % locals() )
        return odict( zip( ("Unique intervals", "mean_interval_length" ), data) )

##################################################################################
class UniqueIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT (stop-start) FROM %(track)s_unique_intervals" % locals() )
        return { "length" : data }

##################################################################################
class SharedIntervals( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_shared_intervals" % locals() )
        return odict( zip( ("Shared intervals", "mean_interval_length" ), data) )

##################################################################################
class sharedIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT (stop-start) FROM %(track)s_shared_intervals" % locals() )
        return { "length" : data }

##################################################################################
class UniqueIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_unique_intervals u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "peakval" : data }

##################################################################################
class UniqueIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_unique_intervals" u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "avgval" : data }

##################################################################################
class UniqueIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT fold FROM %(track)s_unique_intervals" u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "Fold Change" : data }


##################################################################################
class SharedIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_shared_intervals u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "peakval" : data }

##################################################################################
class SharedIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_shared_intervals" u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "avgval" : data }

##################################################################################
class SharedIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_shared_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT fold FROM %(track)s_shared_intervals" u, %(track)_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start 
                                  AND u.stop=i.end''' % locals() )
        return { "Fold Change" : data }

##################################################################################
class UniqueIntervalCoverage( cpgTracker ):
    """ """

    mPattern = "_unique_coverage$"

    def __call__(self, track, slice = None):
        fields = self.getValues("SELECT distinct track FROM bam_stats")
        sep = ""
        select_statement = ""
        out_fields = []
        for f in fields:
            if f.find("input") < 0 and f.find("2") < 0:
                select_statement += sep + f.replace("-","_")
                out_fields.append(f)
                sep = ","
        print out_fields
        statement = '''SELECT %(select_statement)s FROM %(track)s_unique_coverage''' % locals()
        data = self.get( statement )
        return result



