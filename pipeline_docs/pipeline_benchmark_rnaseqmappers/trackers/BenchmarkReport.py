import os, sys, re, types, itertools, glob, collections

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

import numpy

# get from config file
UCSC_DATABASE="hg19"

###################################################################
###################################################################
## parameterization

EXPORTDIR=P['benchmark_rnaseqmappers_exportdir']
DATADIR=P['benchmark_rnaseqmappers_datadir']
DATABASE=P['benchmark_rnaseqmappers_backend']

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "%s/*.bam" % DATADIR), "%s/(\S+).bam" % DATADIR) 

###########################################################################
class BenchmarkTracker( TrackerSQL ):
    '''Define convenience tracks for plots'''
    def __init__(self, *args, **kwargs ):
        TrackerSQL.__init__(self, *args, backend = DATABASE, **kwargs )
    

class CoverageTrackerExons( BenchmarkTracker ):
    pattern = "(.*)_exon_coverage"    
    section = "exon"
    
class CoverageTrackerRegions( BenchmarkTracker ):
    pattern = "(.*)_region_coverage"
    section = "region"

###########################################################################
class CoverageProportions( BenchmarkTracker ):
    def __call__(self, track ):
        statement = '''SELECT CAST(antisense_nreads AS FLOAT)/ (antisense_nreads + sense_nreads ) AS proportion
                              FROM %(track)s_%(section)s_coverage 
                              ORDER BY gene_id'''
        return odict( (("proportion", self.getValues(statement)),) )

###########################################################################
class CoverageCounts( BenchmarkTracker ):
    def __call__(self, track ):
        statement = '''SELECT sense_nreads 
                              FROM %(track)s_%(section)s_coverage 
                              ORDER BY gene_id'''
        return odict( (("nreads", self.getValues(statement)),) )

###########################################################################
class CoverageTotals( BenchmarkTracker ):
    '''percentages are those of mapped reads.

    ratio is percent antisense compared to sense.
    '''
    
    def __call__(self, track ):
        r = odict()
        r["total"] = self.getValue( '''SELECT reads_total FROM bam_stats WHERE track = '%(track)s' ''' )
        mapped = self.getValue( '''SELECT reads_mapped FROM bam_stats WHERE track = '%(track)s' ''' )
        r["mapped"] = mapped
        statement = '''SELECT SUM(sense_nreads) + SUM(antisense_nreads) AS anysense,
                              CAST((SUM(sense_nreads) + SUM(antisense_nreads)) AS FLOAT) / %(mapped)i AS anysense_percent,
                              SUM(antisense_nreads) AS antisense,
                              CAST(SUM(antisense_nreads) AS FLOAT) / %(mapped)i AS antisense_percent,
                              SUM(sense_nreads) AS sense,
                              CAST(SUM(sense_nreads) AS FLOAT) / %(mapped)i AS sense_percent,
                              CAST(SUM(antisense_nreads) AS FLOAT)/ (SUM(antisense_nreads) + SUM(sense_nreads) ) AS ratio
                              FROM %(track)s_%(section)s_coverage '''
        r.update( self.getRow(statement) )
        return r

###########################################################################
class CoverageCorrelation( BenchmarkTracker ):
    limit = 100000
    tracks = ["all"]

    def __call__(self, track ):
        
        tracks = self.getTracks()
        fields = ",".join( [ "x%i.sense_nreads AS %s" % (i,x) for i,x in enumerate( tracks ) ] )
        tables = ", ".join( [ "%s_%s_coverage AS x%i" %( x,self.section,i) for i,x in enumerate( tracks ) ] )
        where = " AND ".join( ["x0.gene_id = x%i.gene_id" % i for i in range(len(tracks))] )

        statement = '''SELECT %(fields)s FROM %(tables)s 
                              WHERE %(where)s LIMIT %(limit)i'''

        data = self.getAll(statement)
        return self.getAll( statement )

class CoverageProportionsExons( CoverageProportions, CoverageTrackerExons ): pass
class CoverageProportionsRegions( CoverageProportions, CoverageTrackerRegions ): pass
class CoverageTotalsExons( CoverageTotals, CoverageTrackerExons ): pass
class CoverageTotalsRegions( CoverageTotals, CoverageTrackerRegions ): pass
class CoverageCountsExons( CoverageCounts, CoverageTrackerExons ): pass
class CoverageCountsRegions( CoverageCounts, CoverageTrackerRegions ): pass
class CoverageCorrelationExons( CoverageCorrelation, CoverageTrackerExons ): pass

###########################################################################
###########################################################################
###########################################################################

class ExonValidationSummary( BenchmarkTracker, SingleTableTrackerRows ):
    table = "exon_validation"

class TranscriptomeValidationSummary( BenchmarkTracker, SingleTableTrackerRows ):
    table = "transcriptome_validation"

###########################################################################
###########################################################################
###########################################################################
class ReadCorrespondenceTracker( BenchmarkTracker ):
    pattern = "(.*)_overrun"
    table = "read_correspondence"
    limit = 100000

class QCFailedBasesInMappedReads( ReadCorrespondenceTracker ):

    def __call__(self, track):
        statement = '''SELECT nfailed FROM %(table)s WHERE %(track)s_nh > 0 LIMIT %(limit)i'''
        return odict( ((track, self.getValues( statement )),) )

def asMatrix( rownames, colnames, data ):
    
    nrows, ncols = len(rownames), len(colnames)
    matrix = numpy.zeros( ( nrows, ncols ) )
    for x, y in data: 
        try: matrix[x,y] += 1
        except IndexError: pass
    results = odict()
    for x in range(nrows):
        r = odict()
        for y in range(ncols): 
            r[str(colnames[y])] = matrix[x,y]
        results[str(rownames[x])] = r
        
    return results

class QCFailedVersusMatched( ReadCorrespondenceTracker ):
    '''returns a matrix: QCfailed bases vs found'''
    max_failed = 60
    max_nmatched = 20
    tracks = ["all"]

    def __call__(self, track):
        statement = '''SELECT nfailed, nmatched FROM %(table)s LIMIT %(limit)i'''

        return asMatrix( range( self.max_failed),
                         range( self.max_nmatched),
                         self.get(statement) )


class QCFailedVersusHits( ReadCorrespondenceTracker ):
    '''returns a matrix: the number hits versus QCfailed bases.'''
    max_failed = 60
    max_nh = 10

    def __call__(self, track):
        statement = '''SELECT nfailed, %(track)s_nh FROM %(table)s LIMIT %(limit)i'''

        return asMatrix( range( self.max_failed),
                         range( self.max_nh),
                         self.get(statement) )

class QCFailedVersusLocation( ReadCorrespondenceTracker ):
    '''returns a matrix: the number locations versus QCfailed bases.'''
    max_failed = 60
    max_locations = 20
    def __call__(self, track):
        statement = '''SELECT nfailed, nlocations FROM %(table)s WHERE %(track)s_nh > 0 LIMIT %(limit)i'''

        return asMatrix( range( self.max_failed),
                         range( self.max_locations),
                         self.get(statement) )

class MappedLocationsVersusHits( ReadCorrespondenceTracker ):

    max_nmatched = 20
    max_nh = 10
    def __call__(self, track):
        statement = '''SELECT nmatched, %(track)s_nh FROM %(table)s WHERE nmatched > 0 LIMIT %(limit)i'''
        
        return asMatrix( range( self.max_nmatched),
                         range( self.max_nh),
                         self.get(statement) )

class MatrixFound( ReadCorrespondenceTracker ):

    def getSlices( self ):
        return self.getTracks()

    def __call__(self, track, slice ):
        return self.getValue( '''SELECT COUNT(*) FROM %(table)s WHERE %(track)s_nh > 0 AND %(slice)s_nh > 0''' )

class ReadQualitiesVersusHits( ReadCorrespondenceTracker ):
    
    slices = map(str, range( 0, 10 ))
    
    def __call__(self, track, slice ):
        return self.getValues( '''SELECT nfailed FROM %(table)s WHERE %(track)s_nh > 0 LIMIT %(limit)i''' )
