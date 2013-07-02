import os, sys, re, types, math, itertools
from VariantsReport import *

class PolyphenHitsHumanDiv( VariantsTracker ):
    '''output a genelist of genes with deleterious variants from Polyphen'''

    mPattern = "^annotations$"

    def __call__(self, track, slice = None ):

        headers = ("Track", "gene_id", "gene_name", "transcript_id", "chromosome", "position", "reference_base", "genotype", "variant_type", "consensus_quality", "rms_mapping_quality", "Ref Amino acid", "Var amino acid","prediction", "pph2_class")
               
        statement = '''
        SELECT distinct
            a.track,
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            a.chromosome, 
            a.position,
            a.reference_base,
            a.genotype,
            a.variant_type,
            a.consensus_quality,
            a.rms_mapping_quality,
            phd.aa1,
            phd.aa2,
            phd.prediction,
            phd.pph2_class
        FROM
            annotations a, polyphen_map pm, polyphen_humdiv phd,
            annotations.transcript_info AS i
        WHERE i.transcript_id = pm.transcript_id
        AND pm.protein_id=phd.protein_id
        AND pm.snp_id=phd.snp_id
        AND a.chromosome=pm.contig
        AND a.position=pm.pos
        ORDER By a.track, i.gene_id
        '''

        data = self.getAll( statement )
        return data

class PolyphenHitsHumanVar( VariantsTracker ):
    '''output a genelist of genes with deleterious variants from Polyphen'''

    mPattern = "^annotations$"

    def __call__(self, track, slice = None ):

        headers = ("Track", "gene_id", "gene_name", "transcript_id", "chromosome", "position", "reference_base", "genotype", "variant_type", "consensus_quality", "rms_mapping_quality", "Ref Amino acid", "Var amino acid","prediction", "pph2_class")
        
        #field = self.getPrefix(slice) + "stop_min"        
        statement = '''
        SELECT distinct
            a.track,
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            a.chromosome, 
            a.position,
            a.reference_base,
            a.genotype,
            a.variant_type,
            a.consensus_quality,
            a.rms_mapping_quality,
            phd.aa1,
            phd.aa2,
            phd.prediction,
            phd.pph2_class
        FROM
            annotations a, polyphen_map pm, polyphen_humvar phd,
            annotations.transcript_info AS i
        WHERE i.transcript_id = pm.transcript_id
        AND pm.protein_id=phd.protein_id
        AND pm.snp_id=phd.snp_id
        AND a.chromosome=pm.contig
        AND a.position=pm.pos
        ORDER By a.track, i.gene_id
        ''' % self.members(locals())

        return odict( zip( headers, zip(*self.get( statement ))) )


