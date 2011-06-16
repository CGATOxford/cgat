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

class CoverageProportionsExons( CoverageProportions, CoverageTrackerExons ): pass
class CoverageProportionsRegions( CoverageProportions, CoverageTrackerRegions ): pass
class CoverageTotalsExons( CoverageTotals, CoverageTrackerExons ): pass
class CoverageTotalsRegions( CoverageTotals, CoverageTrackerRegions ): pass
class CoverageCountsExons( CoverageCounts, CoverageTrackerExons ): pass
class CoverageCountsRegions( CoverageCounts, CoverageTrackerRegions ): pass

###########################################################################
###########################################################################
###########################################################################

class ExonValidationSummary( BenchmarkTracker, SingleTableTrackerRows ):
    table = "exon_validation"

###########################################################################
###########################################################################
###########################################################################
class ReadCorrespondenceTracker( BenchmarkTracker ):
    tracks = ['data_bfast',
              'data_bowtie',
              'data_bwa',
              'data_shrimp',
              'data_tophat',
              'data_trim3x5_bfast',
              'data_trim3x5_bowtie',
              'data_trim3x5_bwa',
              'data_trim3x5_shrimp',
              'data_trim3x5_tophat',
              'strictdata_bowtie',
              'strictdata_tophat' ]

    table = "read_correspondence"
    limit = 100000

class QCFailedBasesInMappedReads( ReadCorrespondenceTracker ):

    def __call__(self, track):
        statement = '''SELECT nfailed FROM %(table)s WHERE %(track)s_nh > 0 LIMIT %(limit)i'''
        return odict( ((track, self.getValues( statement )),) )

class QCFailedVersusMappedLocations( ReadCorrespondenceTracker ):
    '''returns a matrix: the number hits versus QCfailed bases.'''
    max_failed = 60
    max_nh = 10
    def __call__(self, track):
        statement = '''SELECT nfailed, %(track)s_nh FROM %(table)s LIMIT %(limit)i'''
        
        matrix = numpy.zeros( (self.max_failed, self.max_nh ))
        for nfailed, nh in self.get(statement): 
            try:
                matrix[nfailed,nh] += 1
            except IndexError: pass
        results = odict()
        for x in range(0,self.max_failed):
            r = odict()
            for y in range(0,self.max_nh): r[str(y)] = matrix[x,y]
            results[str(x)] = r
            
        return results
