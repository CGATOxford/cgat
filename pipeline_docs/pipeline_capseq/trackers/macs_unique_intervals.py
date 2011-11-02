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

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getFirstRow( "SELECT COUNT(*) as number, round(AVG(stop-start),0) as length FROM %(track)s_merged_unique_intervals" % locals() )
        return odict( zip( ("Unique intervals", "mean_interval_length" ), data) )

##################################################################################
class UniqueIntervalLengths( cpgTracker ):
    """Distribution of interval length. """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( "SELECT (stop-start) FROM %(track)s_merged_unique_intervals" % locals() )
        return { "length" : data }

##################################################################################
class UniqueIntervalPeakValues( cpgTracker ):
    """Distribution of maximum interval coverage (the number of reads at peak). """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.peakval FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "peakval" : data }

##################################################################################
class UniqueIntervalAverageValues( cpgTracker ):
    """Distribution of average coverage (the average number of reads within the interval) """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT avgval FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return { "avgval" : data }

##################################################################################
class UniqueIntervalFoldChange( cpgTracker ):
    """Distribution of fold change """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT i.fold FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start''' % locals() )
        return odict( [("Fold Change", data)] )

##################################################################################
class UniqueIntervalTSS( cpgTracker ):
    """Distribution of distance to closest TSS """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getValues( '''SELECT closest_dist FROM %(track)s_merged_unique_intervals u, 
                                  %(track)s_macs_merged_intervals i, %(track)s_merged_tss t
                                  WHERE u.contig=i.contig
                                  AND u.start=i.start
                                  AND t.gene_id=i.interval_id''' % locals() )
        return { "distance" : data }

##################################################################################
class UniqueIntervalCpGDensity( cpgTracker ):
    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pCpG FROM %(track)s_merged_unique_intervals u, 
                               %(track)s_macs_merged_intervals i,%(track)s_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class UniqueIntervalCpGObsExp1( cpgTracker ):
    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp1 FROM %(track)s_merged_unique_intervals u, 
                               %(track)s_macs_merged_intervals i,%(track)s_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class UniqueIntervalCpGObsExp2( cpgTracker ):
    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT CpG_ObsExp2 FROM %(track)s_merged_unique_intervals u, 
                               %(track)s_macs_merged_intervals i,%(track)s_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class UniqueIntervalCpGNumber( cpgTracker ):
    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT nCpG FROM %(track)s_merged_unique_intervals u, 
                               %(track)s_macs_merged_intervals i,%(track)s_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
class UniqueIntervalGCContent( cpgTracker ):
    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( '''SELECT pGC FROM %(track)s_merged_unique_intervals u, 
                               %(track)s_macs_merged_intervals i,%(track)s_merged_composition c
                               WHERE u.contig=i.contig
                               AND u.start=i.start
                               AND c.gene_id=i.interval_id''' % locals() )
        return data

##################################################################################
##################################################################################
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
        result = odict(zip((out_fields),zip(*data))) 
        return result

##################################################################################
##################################################################################
##################################################################################
class UniqueIntervalLengthVsAverageValue( cpgTracker ):
    """Length vs average value. """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT length, avgval FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("length", "avgval"), zip(*data) ) )

##################################################################################
class UniqueIntervalLengthVsPeakValue( cpgTracker ):
    """Length vs peak value """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT length, peakval FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("length", "peakval"), zip(*data) ) )

##################################################################################
class UniqueIntervalLengthVsFoldChange( cpgTracker ):
    """Length vs fold change"""

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT length, fold FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("length", "foldchange"), zip(*data) ) )

##################################################################################
class UniqueIntervalAvgValVsPeakVal( cpgTracker ):
    """average value vs peak value """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT avgval, peakval FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("avgval", "peakval"), zip(*data) ) )

##################################################################################
class UniqueIntervalAvgValVsFoldChange( cpgTracker ):
    """average value vs fold change """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT avgval, fold FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("avgval", "foldchange"), zip(*data) ) )

##################################################################################
class UniqueIntervalPeakValVsFoldChange( cpgTracker ):
    """Peak value vs fold change """

    mPattern = "_merged_unique_intervals$"

    def __call__(self, track, slice = None):
        data = self.get( '''SELECT peakval, fold FROM %(track)s_merged_unique_intervals u, %(track)s_macs_merged_intervals i
                            WHERE u.contig=i.contig
                            AND u.start=i.start''' % locals() )
        return odict( zip( ("peakval", "foldchange"), zip(*data) ) )

