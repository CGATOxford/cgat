import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram
import ChipseqReport


from SphinxReport.Tracker import *

class TrackerROI( TrackerSQL ):

    def getTracks(self, subset):
        return self.getValues( "SELECT DISTINCT(class) FROM regions_of_interest" )

##################################################################################
##################################################################################
##################################################################################
class AllRegionsOfInterest( TrackerROI ):
    '''return all regions of interest.'''
    
    def __call__(self, track, slice = None ):

        columns = self.getColumns( "regions_of_interest" )
        extra_columns = sorted( [ x for x in columns if x not in ("class", "contig", "start", "end", "roi_id", "pos" ) ] )
        extra =",".join( extra_columns)
        data = self.get( '''SELECT roi_id, contig, start, end, %(extra)s FROM
        regions_of_interest WHERE class = '%(track)s' ''' % locals())

        n = odict()
        i = 0
        for d in data:
            i += 1
            roi_id, contig, start, end = d[:4]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[roi_id] = odict( zip ( ("pos",) + tuple(extra_columns), (pos,)+d[3:] ) )

        return n

##################################################################################
##################################################################################
##################################################################################
class RegionSizes( TrackerROI ):
    '''return all regions of interest.'''
    
    def __call__(self, track, slice = None ):

        return odict( ( ("size",
                         self.getValues( "SELECT end-start FROM regions_of_interest WHERE class='%(track)s'" % locals())
                         ), )
                      )

##################################################################################
##################################################################################
##################################################################################
class ROIOverlapCounts( ChipseqReport.DefaultTracker ):
    mPattern = "_intervals$"
    
    def __call__( self, track, slice = None ):
        '''return a table linking regions of interest,
        their associated marker SNP and ChIP-seq intervals.

        Note that only associations are shown where both SNPS and
        intervals are present.
        '''

        columns = ("class", "nroi", "nsnps", "nintervals")
        statement = '''
        SELECT class,
               COUNT(DISTINCT roi.roi_id),
               COUNT(DISTINCT(snp.snp)),
               COUNT(DISTINCT ovl.interval_id)
        FROM regions_of_interest as roi,
             snps as snp,
             %(track)s_roi as ovl
        WHERE
             snp.roi_id = roi.roi_id AND
             ovl.roi_id = roi.roi_id 
        GROUP BY class
        '''

        data = self.get( statement % locals() )

        return odict( zip( columns,
                           zip(*data ) ) )
    

##################################################################################
##################################################################################
##################################################################################
class ROIOverlap( ChipseqReport.DefaultTracker ):
    mPattern = "_intervals$"
    
    def __call__( self, track, slice = None ):
        '''return a table linking regions of interest,
        their associated marker SNP and ChIP-seq intervals.
        '''

        columns = ("class", "contig", "roi_start", "roi_end", "snp", "pos", "iv_start", "iv_end", "distance" )
        statement = '''
        SELECT class, roi.contig, roi.start, roi.end, snp.snp, snp.pos, i.start, i.end,
        min( abs(i.start - snp.pos), abs(i.end - snp.pos) )
        FROM regions_of_interest as roi,
             %(track)s_roi as ovl,
             %(track)s_intervals as i
             LEFT JOIN snps as snp ON snp.roi_id = roi.roi_id 
        WHERE
             ovl.roi_id = roi.roi_id AND
             ovl.interval_id = i.interval_id
        ORDER BY class, roi.contig
        '''

        data = self.get( statement % locals() )

        return odict( zip( columns,
                           zip(*data ) ) )

##################################################################################
##################################################################################
##################################################################################
class ROIOverlapWithGenes( ChipseqReport.DefaultTracker ):
    mPattern = "_intervals$"

    def __call__( self, track, slice = None ):
        '''return a table linking regions of interest,
        their associated marker SNP and ChIP-seq intervals.

        '''

        columns = ("class", "contig", "roi_start", "roi_end", "snp", "pos", 
                   "genenames",
                   "iv_start", "iv_end", "iv_distance",
                   "gene_id", "gene_name", 
                   "gene_strand", "gene_distance" )

        statement = '''
        SELECT class, roi.contig, roi.start, roi.end, 
               snp.snp, snp.pos, 
               (SELECT group_concat( gene_name, ',') FROM roi_genes WHERE roi_genes.roi_id = roi.roi_id) AS genenames,
               i.start, i.end,
               min( abs(i.start - snp.pos), abs(i.end - snp.pos) ) AS iv_d,
               closest_id, info.gene_name, closest_strand, closest_dist
        FROM regions_of_interest as roi,
             %(track)s_roi as ovl,
             %(track)s_intervals as i,
             %(track)s_tss as tss,
             gene_info as info
             LEFT JOIN snps as snp ON snp.roi_id = roi.roi_id 
        WHERE
             ovl.roi_id = roi.roi_id AND
             ovl.interval_id = i.interval_id AND
             ovl.interval_id = tss.gene_id AND
             tss.closest_id = info.gene_id
        ORDER BY class, roi.contig
        '''

        data = self.get( statement % locals() )

        return odict( zip( columns,
                           zip(*data ) ) )

##################################################################################
##################################################################################
##################################################################################
class GWASOverlapWithGenes( ChipseqReport.DefaultTracker ):
    mPattern = "_intervals$"
    tablename = "gwas_merged"

    def __call__( self, track, slice = None ):
        '''return a table linking regions of interest,
        their associated marker SNP and ChIP-seq intervals.
        '''

        columns = ("class", "contig", 
                   "roi_start", "roi_end", 
                   "snp", "pos", 
                   "genenames",
                   "iv_start", "iv_end", "iv_distance",
                   "gene_id", "gene_name", 
                   "gene_strand", "gene_distance" )

        statement = '''
        SELECT class, regions.contig, regions.start, regions.end, 
               regions.snp, regions.pos, 
               (SELECT group_concat( gene_name, ',') FROM %(tablename)s_genes WHERE 
                       %(tablename)s_genes.roi_id = regions.roi_id) AS genenames,
               i.start, i.end,
               min( abs(i.start - regions.pos), abs(i.end - regions.pos) ) AS iv_d,
               closest_id, info.gene_name, closest_strand, closest_dist
        FROM %(tablename)s as regions,
             %(track)s_%(tablename)s as ovl,
             %(track)s_intervals as i,
             %(track)s_tss as tss,
             gene_info as info
        WHERE
             ovl.roi_id = regions.roi_id AND
             ovl.interval_id = i.interval_id AND
             ovl.interval_id = tss.gene_id AND
             tss.closest_id = info.gene_id
        ORDER BY class, regions.contig
        '''
        
        data = self.get( statement % self.members( locals() ) )

        return odict( zip( columns,
                           zip(*data ) ) )

    
##################################################################################
##################################################################################
##################################################################################
class GWASIntervalList( ChipseqReport.DefaultTracker ):
    '''return a list with all gwas intervals.'''
    
    def getTracks(self, subset = None ):
        return ["gwas_merged"]
    
    def __call__( self, track, slice = None ):
        
        statement = '''
        SELECT class, contig, start, end, snp FROM
             %(track)s ORDER BY class, contig, start
        ''' % (locals())
        
        data = self.get( statement % self.members( locals() ) )
        return odict( zip( ("class", "contig", "start", "end", "snps"),
                           zip(*data)))


##################################################################################
##################################################################################
##################################################################################
class SelectionOverlapWithGenes( ChipseqReport.DefaultTracker ):
    mPattern = "_intervals$"
    tablename = "selection"

    def __call__( self, track, slice = None ):
        '''return a table linking regions of interest,
        their associated marker SNP and ChIP-seq intervals.
        '''

        columns = ("class", "contig", 
                   "roi_start", "roi_end", 
                   "iv_start", "iv_end", 
                   "gene_id", "gene_name", 
                   "gene_strand", "gene_distance" )

        statement = '''
        SELECT class, regions.contig, regions.start, regions.end, 
               i.start, i.end,
               closest_id, info.gene_name, closest_strand, closest_dist
        FROM %(tablename)s as regions,
             %(track)s_%(tablename)s as ovl,
             %(track)s_intervals as i,
             %(track)s_tss as tss,
             gene_info as info
        WHERE
             ovl.roi_id = regions.roi_id AND
             ovl.interval_id = i.interval_id AND
             ovl.interval_id = tss.gene_id AND
             tss.closest_id = info.gene_id
        ORDER BY class, regions.contig
        '''
        
        data = self.get( statement % self.members( locals() ) )

        return odict( zip( columns,
                           zip(*data ) ) )
