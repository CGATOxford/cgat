import os, sys, re, types
from VariantsReport import *

#####################################################
#####################################################
#####################################################
class TrackerEffects( VariantsTracker ):
    
    # minimum number of truncated codons
    min_truncated = 5

    mPattern = "_effects$"

    def getPrefix( self, slice ):
        if slice == None or slice == "all": prefix = ""
        else: prefix = "%s_" % slice
        return prefix
    
    def getSlices(self, subset = None ):
        if subset == None:
            return []            
        elif "separate" in subset:
            return ("all", "splice", "cds" )
        return subset

#####################################################
#####################################################
#####################################################
class TranscriptsNumAlleles( TrackerEffects ):
    '''return number of variants.'''
    def __call__(self, track, slice = None ):

        field = self.getPrefix(slice) + "nalleles"
        
        statement = '''SELECT %(field)s, COUNT(*)
        FROM %(track)s_effects GROUP BY %(field)s''' % locals()
        
        return odict( self.get( statement ) )

#####################################################
#####################################################
#####################################################
class VariantsPerTranscipt( TrackerEffects ):
    '''number of variants in transcripts.'''
    
    mPattern = "_effects_cds$"
    
    def __call__(self, track, slice = None ):
        field = self.getPrefix(slice) + "nvariant_sites"
        return odict((("variant_sites", self.getValues("SELECT %(field)s FROM %(track)s_effects" % locals() )),))

#####################################################
#####################################################
#####################################################
class TranscriptsGenotypeResolvable( TrackerEffects ):
    '''return genotype counts.

    ``None`` are transrcipts that have no variant.
    '''
    def __call__(self, track, slice = None ):

        field = self.getPrefix(slice) + "genotype"
        
        statement = '''SELECT %(field)s, COUNT(*)
        FROM %(track)s_effects GROUP BY %(field)s''' % locals()

        data = odict(self.get( statement ))
        data["``-``"] = data[None]
        del data[None]

        return data

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
class TranscriptsPerGene( VariantsTracker ):
    '''number of transcripts per gene.'''
    tracks = ("annotations.transcript_info",)
    
    def __call__(self, track, slice = None ):
        statement = '''SELECT COUNT(*) FROM %(track)s GROUP BY gene_id''' % locals()
        return odict( ( ("ntranscripts", self.getValues( statement ) ), ) )

#####################################################
#####################################################
#####################################################
class SplicingCounts( VariantsTracker):
    '''number of transcripts per gene.'''
    mPattern = "_effects_splicing$"
    
    def __call__(self, track, slice = None ):
        columns = ("nintrons",
                   "ncanonical", 
                   "nunchanged", 
                   "ndisrupted",
                   "nnonsynonymous",
                   "nnovel", 
                   "nnunknown",       
                   "nsynonymous",
                   "nframeshifts",
                   "nunchanged_frames",
                   "ncorrected_frames",       
                   "nuncorrected_frames" )

        result = odict()
        for column in columns:
            result[column] = self.getValue('''SELECT SUM(%(column)s) FROM %(track)s_effects_splicing''' % locals() )
        return result

#####################################################
#####################################################
#####################################################
class FrameShiftCorrection( VariantsTracker ):
    '''number of frameshift introns and how many have been corrected.'''
    mPattern = "_effects_splicing$"

    def __call__(self, track, slice = None ):
        result = odict( zip(
            ( "nframeshifts", "nunchanged", "ncorrected", "nuncorrected"),
            self.getFirstRow(
            '''SELECT SUM(nframeshifts), SUM(nunchanged_frames), SUM(ncorrected_frames), SUM(nuncorrected_frames)
            FROM %(track)s_effects_splicing''' % locals() )))
        return result

#####################################################
#####################################################
#####################################################
class FrameShiftCorrectedTranscripts( VariantsTracker ):
    '''return the transcripts that have corrected frameshifts.'''
    mPattern = "_effects_splicing$"

    def __call__(self, track, slice = None ):
        return odict( self.get(
                '''SELECT transcript_id, ncorrected_frames
            FROM %(track)s_effects_splicing WHERE ncorrected_frames > 0''' % locals()) )

#####################################################
#####################################################
#####################################################
class CDSCountsTranscripts( VariantsTracker):
    '''number of transcripts with certain features.'''
    mPattern = "_effects_translation$"
    
    def __call__(self, track, slice = None ):
        columns = ("ncorrected_frames", 
                   "ndeleted_bases",  
                   "ndeleted_codons", 
                   "nincomplete_codons",      
                   "ninserted_bases", 
                   "ninserted_codons",   
                   "nnonsynonymous_codons",
                   "noframe_codons",
                   "nstop_codons",
                   "nstops",  
                   "nsynonymous_codons",
                   "nunaffected_codons", 
                   "nwrong_frames")

        result = odict()
        for column in columns:
            result[column] = self.getValue('''SELECT COUNT(*) FROM %(track)s_effects_translation
                                    WHERE %(column)s > 0''' % locals() )
        return result

#####################################################
#####################################################
#####################################################
class CDSCountsVariants( VariantsTracker ):
    '''number of transcripts with certain features.'''
    mPattern = "_effects_translation$"
    
    def __call__(self, track, slice = None ):
        columns = ("ncorrected_frames", 
                   "ndeleted_bases",  
                   "ndeleted_codons", 
                   "nincomplete_codons",      
                   "ninserted_bases", 
                   "ninserted_codons",   
                   "nnonsynonymous_codons",
                   "noframe_codons",
                   "nstop_codons",
                   "nstops",  
                   "nsynonymous_codons",
                   "nunaffected_codons", 
                   "nwrong_frames")

        select = ",".join(["SUM(%s)" % x for x in columns] )
        statement = '''SELECT %(select)s FROM %(track)s_effects_translation''' % locals() 
        
        return odict( zip( columns, self.getFirstRow( statement) ) )

#####################################################
#####################################################
#####################################################
class VariantsCounts( VariantsTracker ):
    '''number of variants in transcripts.'''
    
    mPattern = "_effects_cds$"
    
    def __call__(self, track, slice = None ):
        result = odict()
        result["cds"]      = self.getValue("SELECT COUNT(*) FROM %(track)s_effects_cds" % locals() )
        result["splicing"] = self.getValue("SELECT COUNT(*) FROM %(track)s_effects_splicing" % locals() )
        return result

#####################################################
#####################################################
#####################################################
class TrackerVariants( VariantsTracker ):
    '''type of cds variants.'''
    
    mPattern = "_effects_cds$"
    column = None
    aggregate = None

    def process(self, data):
        return data
    
    def __call__(self, track, slice = None ):
        if not self.column: raise NotImplementedError

        if self.aggregate == None or self.aggregate == "transcript":
            statement = 'SELECT %(column)s, COUNT(*) FROM %(track)s_effects_cds GROUP BY %(column)s'
        elif self.aggregate == "position":
            statement = 'SELECT %(column)s, COUNT(DISTINCT snp_position) FROM %(track)s_effects_cds GROUP BY %(column)s'

        data = odict( self.get( statement % \
                                       dict( self.members(), **locals())))
                    

        return self.process(data) 

class VariantsCDSEffectCodes( TrackerVariants ):
    column = "code"

class VariantsCDSEffectCodesPerStrain( VariantsCDSEffectCodes ):
    aggregate = "position"
    
class VariantsCDSEffectCodesPerPosition( VariantsTracker ):
    '''count effects per SNP.'''
    mPattern = "^effects_cds$"
    mAsTable = True

    def __call__(self, track, slice = None ):
        
        statement = """SELECT 
                  COUNT(*) AS 'all',
                  SUM(CASE WHEN X > 0 AND N = 0 AND S = 0 THEN 1 ELSE 0 END) AS 'X', 
                  SUM(CASE WHEN N > 0 AND X = 0 AND S = 0 THEN 1 ELSE 0 END) AS 'N',
                  SUM(CASE WHEN S > 0 AND X = 0 AND N = 0 THEN 1 ELSE 0 END) AS 'S' 
                  FROM %(track)s"""

        result = self.getFirstRow( statement % locals() )
        
        r = odict( zip( ("all", "X", "N", "S"), result ) )
        r["ambiguous"] = result[0] - sum(result[1:])
        del r["all"]
        return r

class VariantsCDSVariantCodes( TrackerVariants ):
    '''type of cds variants according variant code.'''
    column = "variant_code"
    
    def process(self,data):
        '''quote `+` and `-`.'''
        result = odict()
        for key, counts in data.iteritems():
            result["``%s``" % key ] = counts
        return result
        
class VariantsCDSVariantCodesPerPosition( VariantsCDSVariantCodes ):
    aggregate = "position"

class VariantsCDSVariantTypes( TrackerVariants ):
    '''type of cds variants according variant code.

    The counts are symmetrized, i.e., W,D = D,W.

    Note that the counts are not normalized by transcripts.
    If there is a variant affecting multiple transcripts,
    it will be counted multiple times.
    '''
    column = "variant_type"
    
    def process(self,data):
        '''symmetrize counts.'''
        result = collections.defaultdict(int)
        for key, counts in data.iteritems():
            try:
                if "," in key: key = ",".join( sorted(key.split(",")))
            except TypeError:
                pass
            result[key] += counts
        return result

##################################################################################
##################################################################################
##################################################################################
## Substitutions
##################################################################################
class SubstitutionMatrixNucleotides( VariantsTracker ):
    '''Matrix of single nucleotide substitution polymorphisms.

    Counts both homozygous ([ACGT]) and heterozygous substitutions ([^ACGT])
    '''
    pattern = "(.*)_annotations$"
    slices = ("A","C","G","T")

    def __call__(self, track, slice=None):
        snpcodes = "','".join(SNP_CODES)
        data = self.get( "SELECT genotype, COUNT(*) FROM %(track)s_annotations WHERE reference_base='%(slice)s' AND variant_type IN ('%(snpcodes)s') GROUP BY genotype" % locals() )
        d = odict( [ x for x in data if x[0] in "ACGT"] + [(x[0].lower(),x[1]) for x in data if x[0] not in "ACGT"] )
        return d
        

