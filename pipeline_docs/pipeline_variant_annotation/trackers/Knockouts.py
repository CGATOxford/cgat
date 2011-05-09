import os, sys, re, types
from VariantsReport import *

from Effects import TrackerEffects

#####################################################
#####################################################
#####################################################
class TranscriptsTruncatedStopsMin( TrackerEffects ):
    '''return number of truncated codons for transcipts with variants.'''

    def __call__(self, track, slice = None ):

        p = self.getPrefix(slice)
        field_select = p + "stop_min"
        field_where = p + "nalleles"

        statement = '''SELECT COUNT(*)
        FROM %(track)s_effects WHERE %(field_where)s > 0 AND %(field_select)s >= %(min_truncated)i''' % \
            self.members(locals())
        
        return odict((("transcripts with stops", self.getValue( statement )),))

#####################################################
#####################################################
#####################################################
class TranscriptsTruncatedStopsMax( TrackerEffects ):
    '''return number of truncated codons for transcipts with variants.'''

    def __call__(self, track, slice = None ):

        p = self.getPrefix(slice)
        field_select = p + "stop_max"
        field_where = p + "nalleles"

        statement = '''SELECT COUNT(*)
        FROM %(track)s_effects WHERE %(field_where)s > 0 AND %(field_select)s >= %(min_truncated)i''' %\
            self.members(locals())
        
        return odict((("transcripts with stops", self.getValue( statement )),))

#####################################################
#####################################################
#####################################################
class GenesTruncatedStops( TrackerEffects ):
    '''integrate stats for transcripts by genes.'''

    def __call__(self, track, slice = None ):

        p = self.getPrefix(slice)
        statement = '''
                 SELECT MIN(%(field)s) AS m
                 FROM
                 %(track)s_effects AS e,
                 annotations.transcript_info AS i
                 WHERE i.transcript_id = e.transcript_id
                 GROUP BY i.gene_id
                 HAVING m >= %(min_truncated)i
                 '''

        r = odict()

        field = p + "stop_min"
        r["genes with min stops"] = len(self.getValues( statement % self.members(locals()) ))
        field = p + "stop_max"
        r["genes with max stops"] = len(self.getValues( statement % self.members(locals()) ))
        return r

#####################################################
#####################################################
#####################################################
class GeneListTruncatedStopsMin( TrackerEffects ):
    '''output a genelist of genes with truncated stops.
    '''

    def __call__(self, track, slice = None ):

        headers = ("gene_id", "gene_name", "min(cds_len)", "ntranscripts", "truncated")
        
        field = self.getPrefix(slice) + "stop_min"        
        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            MIN(e.cds_len),
            COUNT(DISTINCT i.transcript_id),
            MIN(%(field)s) AS m
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        GROUP BY i.gene_id
        HAVING m >= %(min_truncated)i 
        ''' % self.members(locals())

        return odict( zip( headers,
                           zip(*self.get( statement ))) )

#####################################################
#####################################################
#####################################################
class GenesNMDMin( TrackerEffects ):
    '''return genes in which all transcripts are knocked out by :term:`nmd`.'''

    def __call__(self, track, slice = None ):

        raise NotImplementedError
        slice = "cds"
        p = self.getPrefix(slice)
        field_select = p + "stop_min"
        field_where = p + "nalleles"

        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            MIN(e.cds_len),
            COUNT(DISTINCT i.transcript_id),
            MIN(%(field)s) AS m
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        GROUP BY i.gene_id
        ''' % self.members(locals())
        
        return odict((("transcripts with nmd", self.getValue( statement )),))

#####################################################
#####################################################
#####################################################
class TranscriptListTruncatedStopsMin( TrackerEffects ):
    '''output a genelist of genes with truncated stops.'''

    def __call__(self, track, slice = None ):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "truncated")
        
        field = self.getPrefix(slice) + "stop_min"        
        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            %(field)s AS m
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id AND
              m > %(min_truncated)i 
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict( zip( headers,
                           zip(*self.get( statement ))) )

#####################################################
#####################################################
#####################################################
class TranscriptListTruncatedStopsMax( TrackerEffects ):
    '''output a genelist of genes with truncated stops.'''

    def __call__(self, track, slice = None ):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "truncated")
        
        field = self.getPrefix(slice) + "stop_max"        
        statement = '''
        SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            %(field)s AS m
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id AND
              m > %(min_truncated)i 
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict( zip( headers,
                           zip(*self.get( statement ))) )



#####################################################
#####################################################
#####################################################
class TranscriptsNMDMin( TrackerEffects ):
    '''return number of transcripts that are likely to be affected by :term:`nmd`.'''

    def __call__(self, track, slice = None ):

        slice = "cds"
        p = self.getPrefix(slice)
        field_select = p + "stop_min"
        field_where = p + "nalleles"

        statement = '''
        SELECT COUNT(*)
        FROM %(track)s_effects WHERE 
        %(field_where)s > 0 AND 
        %(field_select)s > 0 AND 
        cds_len - %(field_select)s * 3 < last_exon_start''' % \
            self.members(locals())
        
        return odict((("transcripts with nmd", self.getValue( statement )),))

#####################################################
#####################################################
#####################################################
class TranscriptsNMDMax( TrackerEffects ):
    '''return number of transcripts that are likely to be affected by :term:`nmd`.'''

    def __call__(self, track, slice = None ):

        p = self.getPrefix(slice)
        field_select = p + "stop_max"
        field_where = p + "nalleles"

        statement = '''
        SELECT COUNT(*)
        FROM %(track)s_effects WHERE 
        %(field_where)s > 0 AND 
        %(field_select)s > 0 AND 
        cds_len - %(field_select)s * 3 < last_exon_start''' % \
            self.members(locals())
        
        return odict((("transcripts with nmd", self.getValue( statement )),))

#####################################################
#####################################################
#####################################################
class TranscriptListNMDMin( TrackerEffects ):
    '''output a genelist of genes with truncated stops.'''

    def __call__(self, track, slice = None ):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "truncated", "last_exon_start" )
        
        field = self.getPrefix(slice) + "stop_min"        
        statement = '''
        SELECT
            DISTINCT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            %(field)s AS m,
            e.last_exon_start
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id AND
              m > 0 AND 
              cds_len - m * 3 < last_exon_start
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict( zip( headers,
                           zip(*self.get( statement ))) )

#####################################################
#####################################################
#####################################################
class TranscriptListNMDMax( TrackerEffects ):
    '''output a genelist of genes with truncated stops.'''

    def __call__(self, track, slice = None ):

        headers = ("gene_id", "gene_name", "transcript_id", "cds_len", "truncated", "last_exon_start")
        
        field = self.getPrefix(slice) + "stop_max"        
        statement = '''
        SELECT
            DISTINCT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.cds_len,
            %(field)s AS m,
            e.last_exon_start
        FROM
            %(track)s_effects AS e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id AND
              m > 0 AND 
              cds_len - m * 3 < last_exon_start
        ORDER BY i.gene_id
        ''' % self.members(locals())

        return odict( zip( headers,
                           zip(*self.get( statement ))) )

