import os, sys, re, types, math, itertools
from VariantsReport import *

class PolyphenTracker( VariantsTracker ):

    mPattern = "^polyphen_[^_]+$"
    mAsTables = True
    tablename_map = "polyphen_map"

    def getTracks( self, subset = None ):
        return [ x for x in VariantsTracker.getTracks( self, subset ) if x != "polyphen_map" ]

class PolyphenCounts( PolyphenTracker ):
    '''returns overview of polyphen results.'''

    def getSlices( self, subset = None ):
        return ["snps", "proteins", "genes" ]

    def __call__(self, track, slice = None):

        data = odict()
        join = "WHERE 1"
        if slice == "snps":
            s = "p.snp_id"
        elif slice == "proteins":
            s = "DISTINCT p.protein_id"
        elif slice == "genes":
            s = "DISTINCT gene_id"
            join = ", annotations.transcript_info AS t WHERE t.protein_id = p.protein_id"
        else:
            s = "p.snp_id"

        data["total"] = self.getValue( '''SELECT COUNT(%(s)s) FROM %(track)s AS p %(join)s
                                       ''' % self.members(locals()))

        data.update( odict( self.get( '''SELECT prediction, COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            GROUP BY prediction''' % self.members(locals()) ) ) )

        data.update( odict( self.get( '''SELECT pph2_class, COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            GROUP BY pph2_class''' % self.members(locals()) ) ) )

        data.update( odict( self.get( '''SELECT based_on, COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            GROUP BY based_on''' % self.members(locals()) ) ) )

        data["with structure"] = self.getValue( '''SELECT COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            AND Nstruct > 0''' % self.members(locals()) )

        data["with pfam"] = self.getValue( '''SELECT COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            AND PfamHit IS NOT NULL''' % self.members(locals()) )

        data["with site"] = self.getValue( '''SELECT COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            AND site IS NOT NULL''' % self.members(locals()) )

        data["with region"] = self.getValue( '''SELECT COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            AND region IS NOT NULL''' % self.members(locals()) )

        data["with PHAT"] = self.getValue( '''SELECT COUNT(%(s)s) 
                                            FROM %(track)s AS p %(join)s
                                            AND PHAT IS NOT NULL''' % self.members(locals()) )

        return data

class PolyphenResults( StrainTracker ):
    '''returns overview of polyphen results.'''

    tablename = "polyphen"
    tablename_map = "polyphen_map"

    def getSlices( self, subset = None):
        return [ "polyphen_HumVar", "polyphen_HumDiv" ]

    def __call__(self, track, slice = None):

        s = "DISTINCT locus_id"

        if track == "all": where = "1"
        else: where = "map.track = '%(track)s'" % locals()

        data = odict( self.get( '''
        SELECT prediction, COUNT(%(s)s) 
        FROM %(slice)s AS t, %(tablename_map)s as map
        WHERE map.snp_id = t.snp_id AND %(where)s GROUP BY prediction ''' % self.members(locals() ) ) )

        return data

class PolyphenDeleteriousGenesPerStrain( StrainTracker ):
    '''returns number of genes with deleterious SNPs per strain.'''

    tablename = "polyphen"
    tablename_map = "polyphen_map"

    def getSlices( self, subset = None):
        return [ "polyphen_HumVar", "polyphen_HumDiv" ]

    def __call__(self, track, slice = None):

        s = "DISTINCT i.gene_id"

        if track == "all": where = "1"
        else: where = "map.track = '%(track)s'" % locals()

        data = odict( self.get( '''
        SELECT prediction, COUNT(%(s)s) 
        FROM %(slice)s AS t, 
             %(tablename_map)s as map,
             annotations.transcript_info AS i
        WHERE map.snp_id = t.snp_id AND 
              i.transcript_id = map.transcript_id AND
              prediction != 'benign' AND
              %(where)s 
              GROUP BY prediction ''' % self.members(locals() ) ) )

        return data

class PolyphenDistribution( PolyphenTracker ):
    '''return distributions for various columns.'''

    def getSlices( self, subset = None ):
        return ("pph2_prob", 
                "pph2_TPR", 
                "pph2_FDR", 
                "pph2_FDR", 
                "dscore",
                "Nobs",
                "NStruct",
                "dVol",
                "dProp" )
    
    def __call__(self, track, slice = None ):
        
        return odict( ( 
                (slice, self.getValues( "SELECT %(slice)s FROM %(track)s WHERE %(slice)s IS NOT NULL" % locals () )), ) )

class PolyphenPDeleteriousPerGene( PolyphenTracker ):
    '''return the proportion of deleterious SNPs per gene
    '''
    tablename_map = "polyphen_map"

    def __call__(self, track, slice = None):

        data = self.getValues( '''
        SELECT 
               COUNT(DISTINCT case t.prediction when 'possiblydamaging' then map.locus_id when 'probablydamaging' then map.locus_id else NULL end) /
                   CAST(COUNT(DISTINCT locus_id) AS FLOAT)
               FROM %(track)s as t, 
                    %(tablename_map)s as map, 
                    protein_stats as s,
                    annotations.transcript_info as i 
        WHERE map.snp_id = t.snp_id AND 
              i.transcript_id = map.transcript_id AND
              s.protein_id = map.protein_id
        GROUP BY i.gene_id
        ''' % self.members(locals() ) )

        return odict( ( ("pdeleterious", data), ) )

class PolyphenSnpsPerGeneAndCategory( PolyphenTracker ):
    '''return the proportion of deleterious SNPs per gene
    '''
    tablename_map = "polyphen_map"

    def getSlices( self, subset = None ):
        return ( "unknown", "benign", "possiblydamaging", "probablydamaging")

    def __call__(self, track, slice = None):

        data = self.getValues( '''
        SELECT 
               COUNT(DISTINCT locus_id)
               FROM %(track)s as t, 
                    %(tablename_map)s as map, 
                    protein_stats as s,
                    annotations.transcript_info as i 
        WHERE map.snp_id = t.snp_id AND 
              i.transcript_id = map.transcript_id AND
              s.protein_id = map.protein_id AND
              prediction = '%(slice)s'
        GROUP BY i.gene_id
        ''' % self.members(locals() ) )

        return odict( ((slice, data),) )

class PolyphenEffectAndLength( PolyphenTracker ):
    '''plot number of snps in gene versus gene length. 
    
    Color dot by the proportion of deleterious SNPs found.
    '''
    tablename_map = "polyphen_map"

    def __call__(self, track, slice = None):

        data = self.get( '''
        SELECT 
               MAX(s.length),
               COUNT(DISTINCT locus_id) as nsnps, 
               COUNT(DISTINCT case t.prediction when 'possiblydamaging' then map.locus_id when 'probablydamaging' then map.locus_id else NULL end) /
                   CAST(COUNT(DISTINCT locus_id) AS FLOAT) as pdeleterious
               FROM %(track)s as t, 
                    %(tablename_map)s as map, 
                    protein_stats as s,
                    annotations.transcript_info as i 
        WHERE map.snp_id = t.snp_id AND 
              i.transcript_id = map.transcript_id AND
              s.protein_id = map.protein_id
        GROUP BY i.gene_id
        ''' % self.members(locals() ) )

        data = odict( zip( ("length", "nsnps", "pdeleterious" ), zip(*data) ) )
        data["length"] = [ math.log(x) for x in data["length"] ]
        data["nsnps"] = [ math.log(x) for x in data["nsnps"] ]
        
        return data

class PolyphenEnrichment( PolyphenTracker ):
    '''plot number of snps in gene versus gene length. 
    
    Color dot by the proportion of deleterious SNPs found.
    '''
    tablename_map = "polyphen_map"

    def __call__(self, track, slice = None):

        data = self.get( '''SELECT code, COUNT(*) 
        FROM %(track)s_genestats GROUP BY code
        ''' % self.members(locals() ) )

        return odict(data)

class PolyphenEnrichmentLength( PolyphenTracker ):
    '''plot number of snps in gene versus gene length. 
    
    Color dot by the proportion of deleterious SNPs found.
    '''
    tablename_map = "polyphen_map"

    def getSlices( self, subset = None ):
        return [ '---', ] + [ "".join( x ) for x in itertools.product( ("0","1"), ("0","1"), ("0","1") )  ]

    def __call__(self, track, slice = None):

        data = self.getValues( '''SELECT length
        FROM %(track)s_genestats WHERE code = '%(slice)s'
        ''' % self.members(locals() ) )

        return odict( (("length",data),))

class PolyphenPhase( PolyphenTracker ):
    '''return the number of SNPs in various categories 
       according to the codon position.'''

    def getSlices( self, subset = None ):
        return [ "benign", "possiblydamaging", "probablydamaging", "unknown" ]

    def __call__(self, track, slice = None ):

        data = self.get( '''SELECT phase, COUNT(distinct s.snp_id)
                                  FROM %(track)s AS s, %(tablename_map)s AS m
                                  WHERE m.snp_id = s.snp_id 
                                  AND s.prediction = '%(slice)s'
                                  GROUP BY phase
        ''' % self.members(locals() ) )
        
        return odict(data)


class PolyphenSpeciesDistribution( PolyphenTracker ):
    '''return the number of SNPs in various categories 
       according to the codon position.'''

    def getSlices( self, subset = None ):
        return [ "benign", "possiblydamaging", "probablydamaging", "unknown" ]

    def __call__(self, track, slice = None ):

        data = self.getValues( '''SELECT count(DISTINCT track) 
                               FROM %(track)s AS s, %(tablename_map)s as m 
                               WHERE s.snp_id = m.snp_id AND s.prediction = '%(slice)s' 
                               GROUP BY m.locus_id''' % (self.members(locals())))
        
        return odict((("counts", data),))
